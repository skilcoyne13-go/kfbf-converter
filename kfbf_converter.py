#!/usr/bin/env python3
"""
kfbf_converter.py - Convert KFBio fluorescent whole-slide images (.kfbf) into
open computational-pathology formats.

Output modes
------------
ome    Multi-channel pyramidal OME-TIFF          <stem>.ome.tif
rgb    Composite 8-bit RGB pyramidal BigTIFF     <stem>.rgb.tif
split  One grayscale pyramidal BigTIFF/channel   <stem>_<channel>.tif
all    All of the above

All writers stream tiles level-by-level straight into ``tifffile.TiffWriter``;
full-resolution planes are never assembled in memory.

Usage
-----
    python kfbf_converter.py slide.kfbf --mode all -o out/
    python kfbf_converter.py slides_folder --mode ome -o out/
    python kfbf_converter.py slide.kfbf --info        # dump parsed structure

As a module
-----------
    from kfbf_converter import KFBioSlide, convert_slide
    with KFBioSlide("slide.kfbf") as slide:
        tile = slide.get_tile(level=0, col=3, row=5, channel=0)
    convert_slide("slide.kfbf", "out", mode="all")

Dependencies: numpy, tifffile, imagecodecs, pillow, tqdm
"""

from __future__ import annotations

import argparse
import io
import logging
import math
import mmap
import os
import re
import struct
import sys
import threading
from collections import OrderedDict, deque
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable, Iterator, Sequence

import numpy as np
import tifffile
from tqdm import tqdm

try:
    import imagecodecs
except ImportError:  # pragma: no cover
    imagecodecs = None

try:
    from PIL import Image
except ImportError:  # pragma: no cover
    Image = None

__version__ = "1.0.0"

log = logging.getLogger("kfbf_converter")

# --------------------------------------------------------------------------
# Format constants
# --------------------------------------------------------------------------

MAGIC = b"\xf1\x01\xee\xee"
SIGNATURE = b"KFBF"
TLV_MARKER = b"\xff\x01\xee\xee"  # 0xEEEE01FF little-endian
INDEX_MARKER = b"\xf1\x04\xee\xee"  # 0xEEEE04F1
RECORD_TAIL = 0xEEEE04FF
RECORD_DELIM = 0xEEEE04F1
RECORD_SIZE = 64
SRC_TILE = 256

TAG_CHANNEL_COUNT = 0x4B
TAG_CHANNEL_FIRST = 0x4C

DEFAULT_MPP = 0.2528
FALLBACK_CHANNEL_NAMES = ["DAPI", "FITC", "TRITC", "Cy5"]

# Synthetic downsample levels are generated only when the file itself does
# not contain a pyramid reaching below this size.
MIN_PYRAMID_DIM = 1024

# Upper bound on uncompressed tiles tifffile buffers for parallel compression.
WRITE_BUFFER_BYTES = 64 * 2**20
try:
    import inspect as _inspect

    _WRITE_HAS_BUFFERSIZE = "buffersize" in _inspect.signature(
        tifffile.TiffWriter.write).parameters
except (TypeError, ValueError):  # pragma: no cover
    _WRITE_HAS_BUFFERSIZE = False

_RECORD_DTYPE = np.dtype(
    {
        "names": ["col", "row", "gcols", "grows", "scale", "z",
                  "ptr_off", "ptr_len", "tail", "delim"],
        "formats": ["<u4", "<u4", "<u4", "<u4", "<f4", "<u4",
                    "<u8", "<u8", "<u4", "<u4"],
        "offsets": [0x00, 0x04, 0x08, 0x0C, 0x10, 0x14,
                    0x20, 0x28, 0x38, 0x3C],
        "itemsize": RECORD_SIZE,
    }
)

# (keywords, RGB) - checked in order, first match wins.
_COLOR_RULES: list[tuple[tuple[str, ...], tuple[int, int, int]]] = [
    (("CY7", "AF750", "AF790", "750", "770", "780", "NIR", "NEAR"), (0, 255, 255)),
    (("CY5", "AF647", "647", "670", "FAR", "RED"), (255, 0, 0)),
    (("TRITC", "CY3", "AF555", "AF568", "555", "568", "570", "ORANGE",
      "RHOD", "TXR", "PE"), (255, 165, 0)),
    (("FITC", "AF488", "488", "GFP", "GREEN", "520", "FAM"), (0, 255, 0)),
    (("DAPI", "HOECHST", "405", "460", "BLUE", "UV"), (0, 0, 255)),
]
_INDEX_COLORS = [(0, 0, 255), (0, 255, 0), (255, 165, 0), (255, 0, 0),
                 (0, 255, 255), (255, 0, 255), (255, 255, 0), (255, 255, 255)]


class KFBFError(Exception):
    """Raised for unreadable or malformed .kfbf files."""


def channel_color(name: str, index: int = 0) -> tuple[int, int, int]:
    """Default pseudo-colour for a fluorophore name (falls back to index)."""
    upper = re.sub(r"[^A-Z0-9]", "", name.upper())
    for keys, rgb in _COLOR_RULES:
        if any(k in upper for k in keys):
            return rgb
    return _INDEX_COLORS[index % len(_INDEX_COLORS)]


def ome_color(rgb: Sequence[int]) -> int:
    """RGB -> OME-XML Color (signed 32-bit RGBA)."""
    r, g, b = (int(v) & 0xFF for v in rgb)
    value = (r << 24) | (g << 16) | (b << 8) | 0xFF
    return value - (1 << 32) if value >= (1 << 31) else value


def decode_jpeg(data: bytes) -> np.ndarray:
    """Decode a JPEG payload to a 2-D uint8 array."""
    arr = None
    if imagecodecs is not None:
        try:
            arr = imagecodecs.jpeg_decode(data)
        except Exception:
            arr = None
    if arr is None:
        if Image is None:
            raise KFBFError("neither imagecodecs nor Pillow could decode JPEG")
        with Image.open(io.BytesIO(data)) as im:
            arr = np.asarray(im.convert("L") if im.mode not in ("L", "I;16") else im)
    if arr.ndim == 3:
        if arr.shape[2] == 1 or (np.array_equal(arr[..., 0], arr[..., 1])
                                 and np.array_equal(arr[..., 0], arr[..., -1])):
            arr = arr[..., 0]
        else:  # unexpected colour JPEG: luminance
            arr = (arr[..., :3].astype(np.float32) @ np.array(
                [0.299, 0.587, 0.114], np.float32)).round()
    if arr.dtype != np.uint8:
        arr = np.clip(arr, 0, 255).astype(np.uint8)
    return arr


# --------------------------------------------------------------------------
# Reader
# --------------------------------------------------------------------------

@dataclass
class Level:
    """One pyramid level (index 0 = highest resolution)."""

    index: int
    scale: float          # objective magnification of this level
    width: int
    height: int
    grid_cols: int
    grid_rows: int
    downsample: float     # relative to level 0
    synthetic: bool = False  # generated by 2x2 averaging, not in file
    coord_mode: str = ""     # how record x/y were interpreted
    raw_grid: tuple = ()     # raw record fields +0x08/+0x0C (diagnostic)


class KFBioSlide:
    """Reader for KFBio fluorescent ``.kfbf`` slides.

    Parameters
    ----------
    path : path to the .kfbf file
    mpp_fallback : pixel size (µm/px) used when the header value is invalid
    channel_names : optional override of the channel labels
    tile_cache : number of decoded 256x256 tiles kept in an LRU cache
    """

    def __init__(self, path: str | os.PathLike, mpp_fallback: float = DEFAULT_MPP,
                 channel_names: Sequence[str] | None = None, tile_cache: int = 1024):
        self.path = Path(path)
        self._fh = open(self.path, "rb")
        try:
            self.file_size = os.fstat(self._fh.fileno()).st_size
            if self.file_size < 0x60:
                raise KFBFError(f"file too small ({self.file_size} bytes)")
            self._mm = mmap.mmap(self._fh.fileno(), 0, access=mmap.ACCESS_READ)
        except Exception:
            self._fh.close()
            raise
        self._cache: OrderedDict = OrderedDict()
        self._cache_size = max(0, int(tile_cache))
        self._lock = threading.Lock()
        self.decode_failures = 0
        self.warnings: list[str] = []
        try:
            self._parse_header(mpp_fallback)
            self._parse_metadata()
            self._parse_index()
            self._resolve_channels(channel_names)
            self._build_levels()
        except Exception:
            self.close()
            raise

    # ---- lifecycle -------------------------------------------------------
    def close(self) -> None:
        mm = getattr(self, "_mm", None)
        if mm is not None and not mm.closed:
            mm.close()
        if not self._fh.closed:
            self._fh.close()

    def __enter__(self) -> "KFBioSlide":
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    def _warn(self, msg: str) -> None:
        self.warnings.append(msg)
        log.warning("%s: %s", self.path.name, msg)

    # ---- parsing ---------------------------------------------------------
    def _parse_header(self, mpp_fallback: float) -> None:
        mm = self._mm
        if mm[0:4] != MAGIC or mm[4:8] != SIGNATURE:
            raise KFBFError(f"not a KFBF file (bad magic/signature)")
        (self.version,) = struct.unpack_from("<f", mm, 0x0C)
        self.n_records, self.height, self.width, self.magnification = \
            struct.unpack_from("<4I", mm, 0x10)
        self.codec = bytes(mm[0x20:0x28]).split(b"\0", 1)[0].decode("ascii", "replace")
        (self.index_offset,) = struct.unpack_from("<Q", mm, 0x44)
        (mpp,) = struct.unpack_from("<f", mm, 0x4C)
        if math.isfinite(mpp) and 0.01 <= mpp <= 100.0:
            # strip float32 noise (0.2527999877... -> 0.2528)
            self.mpp, self.mpp_from_header = float(f"{mpp:.6g}"), True
        else:
            self.mpp, self.mpp_from_header = float(mpp_fallback), False
            self._warn(f"invalid MPP in header ({mpp!r}); using {mpp_fallback}")
        if self.width == 0 or self.height == 0:
            raise KFBFError(f"zero image dimensions in header")
        if self.codec.upper() not in ("JPEG", "JPG", ""):
            self._warn(f"unexpected codec {self.codec!r}; attempting JPEG decode anyway")
        if not (0x60 <= self.index_offset < self.file_size):
            raise KFBFError(f"tile index offset 0x{self.index_offset:X} "
                            f"outside file")

    def _parse_metadata(self) -> None:
        """Parse the Tag|Length|Value block that follows 0x5C."""
        mm = self._mm
        self.tlv: dict[int, list[bytes]] = {}
        pos = 0x5C
        if mm[pos:pos + 4] == TLV_MARKER:
            pos += 4
        else:
            self._warn("metadata TLV marker not found at 0x5C")
            return
        end = self.index_offset if self.index_offset > pos else self.file_size
        while pos + 8 <= end:
            tag, length = struct.unpack_from("<II", mm, pos)
            if (tag == 0 and length == 0) or (tag >> 16) == 0xEEEE:
                break  # terminator or next block marker
            if length > end - pos - 8:
                break
            self.tlv.setdefault(tag, []).append(bytes(mm[pos + 8:pos + 8 + length]))
            pos += 8 + length

    @staticmethod
    def _decode_label(raw: bytes) -> str | None:
        for part in raw.split(b"\0"):
            text = part.decode("latin-1").strip()
            if (text and len(text) <= 64 and any(ch.isalpha() for ch in text)
                    and all(ch.isprintable() for ch in text)):
                return text
        return None

    def _parse_index(self) -> None:
        mm = self._mm
        base = self.index_offset
        # Records normally start right after the 4-byte marker; pick whichever
        # alignment yields more valid tail markers.
        best = None
        for start in ((base + 4, base) if mm[base:base + 4] == INDEX_MARKER else (base, base + 4)):
            avail = (self.file_size - start) // RECORD_SIZE
            n = min(self.n_records, avail)
            if n <= 0:
                continue
            recs = np.frombuffer(mm, dtype=_RECORD_DTYPE, count=n, offset=start).copy()
            score = int(np.count_nonzero(recs["tail"] == RECORD_TAIL))
            if best is None or score > best[0]:
                best = (score, recs)
        if best is None:
            raise KFBFError(f"tile index table is empty/truncated")
        score, recs = best
        if len(recs) < self.n_records:
            self._warn(f"tile index truncated: {len(recs)} of {self.n_records} records")
        if score < 0.5 * len(recs):
            self._warn(f"only {score}/{len(recs)} tile records carry the 0xEEEE04FF "
                       f"tail marker; layout may differ from spec")
        ok = (np.isfinite(recs["scale"]) & (recs["scale"] > 0)
              & (recs["ptr_off"] < self.file_size) & (recs["ptr_len"] < self.file_size))
        if not ok.all():
            self._warn(f"ignoring {int((~ok).sum())} invalid tile records")
            recs = recs[ok]
        if len(recs) == 0:
            raise KFBFError(f"no valid tile records")
        self._records = recs

    def _resolve_channels(self, override: Sequence[str] | None) -> None:
        c = None
        if TAG_CHANNEL_COUNT in self.tlv and len(self.tlv[TAG_CHANNEL_COUNT][0]) >= 4:
            (c,) = struct.unpack_from("<I", self.tlv[TAG_CHANNEL_COUNT][0])
            if not 1 <= c <= 64:
                self._warn(f"implausible channel count {c} in metadata")
                c = None
        if c is None:
            # Infer from distance between the offset and length pointer arrays
            gaps = (self._records["ptr_len"][:64].astype(np.int64)
                    - self._records["ptr_off"][:64].astype(np.int64)) // 8
            gaps = gaps[(gaps >= 1) & (gaps <= 16)]
            c = int(np.median(gaps)) if len(gaps) else (len(override) if override else 3)
            self._warn(f"channel count tag missing; inferred C={c}")
        if override:
            if len(override) != c:
                raise KFBFError(f"{len(override)} channel names given but slide has {c} channels")
            names = list(override)
        else:
            names = []
            for i in range(c):
                raw = self.tlv.get(TAG_CHANNEL_FIRST + i)
                label = self._decode_label(raw[0]) if raw else None
                if label is None:
                    label = (FALLBACK_CHANNEL_NAMES[i] if i < len(FALLBACK_CHANNEL_NAMES)
                             else f"Channel{i}")
                names.append(label)
        self.n_channels = c
        self.channel_names = names
        self.channel_colors = [channel_color(n, i) for i, n in enumerate(names)]
        self._detect_length_width()

    def _detect_length_width(self) -> None:
        """Byte lengths may be stored as uint32 or uint64 per channel."""
        c, mm, size = self.n_channels, self._mm, self.file_size
        votes64 = votes32 = 0
        for rec in self._records[:64]:
            po, pl = int(rec["ptr_off"]), int(rec["ptr_len"])
            if po + 8 * c > size or pl + 4 * c > size:
                continue
            offs = struct.unpack_from(f"<{c}Q", mm, po)
            if pl + 8 * c <= size:
                lens = struct.unpack_from(f"<{c}Q", mm, pl)
                if all(0 < ln <= size and o + ln <= size for o, ln in zip(offs, lens)):
                    votes64 += 1
                    continue
            lens = struct.unpack_from(f"<{c}I", mm, pl)
            if all(0 < ln and o + ln <= size for o, ln in zip(offs, lens)):
                votes32 += 1
        self._len_fmt = f"<{c}Q" if votes64 >= votes32 and votes64 > 0 else f"<{c}I"

    def _build_levels(self) -> None:
        recs = self._records
        scales = np.unique(recs["scale"])[::-1]  # high -> low magnification
        s0 = float(scales[0])
        if self.magnification and abs(s0 - self.magnification) > 0.5:
            log.debug("top level scale %.3f differs from header magnification %d",
                      s0, self.magnification)

        # Z planes
        zvals = np.unique(recs["z"])
        if len(zvals) > 1 and zvals.max() < 4096:
            self.z_values = [int(v) for v in zvals]
        else:
            if len(zvals) > 1:
                self._warn("ignoring implausible Z-plane indices")
            self.z_values = [0]
        zmap = {v: i for i, v in enumerate(self.z_values)}
        multi_z = len(self.z_values) > 1

        # Record +0x00/+0x04 hold the tile's top-left *pixel* coordinate within
        # its level (real KFBio files), so col = x // 256, row = y // 256.
        # Tile-index and level-0-pixel layouts are detected too, per level.
        levels: list[Level] = []
        self._lookup: dict[tuple[int, int, int, int], int] = {}
        unaligned = collisions = out_of_range = 0
        rec_z = recs["z"].astype(np.int64)
        for s in scales:
            s = float(s)
            idx = np.nonzero(recs["scale"] == np.float32(s))[0]
            xs = recs["col"][idx].astype(np.int64)
            ys = recs["row"][idx].astype(np.int64)
            ds = s0 / s
            log2 = math.log2(ds)
            if abs(log2 - round(log2)) < 0.15:  # snap float32 noise to 2^k
                ds = float(2 ** round(log2))
            w = max(1, math.ceil(self.width / ds - 1e-6))
            h = max(1, math.ceil(self.height / ds - 1e-6))
            ncols, nrows = math.ceil(w / SRC_TILE), math.ceil(h / SRC_TILE)
            if levels and (w >= levels[-1].width and h >= levels[-1].height):
                self._warn(f"skipping non-decreasing level at scale {s:g}")
                continue
            xmax, ymax = int(xs.max()), int(ys.max())
            if xmax < ncols and ymax < nrows and (ncols > 1 or nrows > 1) \
                    and (xs % SRC_TILE).any():
                mode, div = "index", 1
            elif xmax < w and ymax < h:
                mode, div = "pixel", SRC_TILE
            elif xmax < self.width and ymax < self.height and ds > 1:
                mode, div = "base-pixel", int(round(SRC_TILE * ds))
            else:
                mode, div = "pixel", SRC_TILE
            if div > 1:
                unaligned += int(np.count_nonzero((xs % div) | (ys % div)))
            li = len(levels)
            raw = (int(recs["gcols"][idx[0]]), int(recs["grows"][idx[0]]))
            levels.append(Level(li, s, w, h, ncols, nrows, ds, coord_mode=mode,
                                raw_grid=raw))
            for i, x, y in zip(idx, xs // div, ys // div):
                if x >= ncols or y >= nrows:
                    out_of_range += 1
                    continue
                key = (li, zmap.get(int(rec_z[i]), 0) if multi_z else 0, int(x), int(y))
                collisions += key in self._lookup
                self._lookup[key] = int(i)
        if unaligned:
            self._warn(f"{unaligned} tile records not aligned to the {SRC_TILE}px grid")
        if out_of_range:
            self._warn(f"{out_of_range} tile records outside their level bounds ignored")
        if collisions:
            self._warn(f"{collisions} duplicate tile records (last one wins)")

        self.n_source_levels = len(levels)
        # Synthesize a pyramid if the file contains only high-resolution levels
        while max(levels[-1].width, levels[-1].height) > MIN_PYRAMID_DIM:
            p = levels[-1]
            levels.append(Level(len(levels), p.scale / 2, max(1, math.ceil(p.width / 2)),
                                max(1, math.ceil(p.height / 2)), 0, 0,
                                p.downsample * 2, synthetic=True))
        self.levels = levels

    # ---- properties ------------------------------------------------------
    @property
    def n_levels(self) -> int:
        return len(self.levels)

    @property
    def n_z(self) -> int:
        return len(self.z_values)

    @property
    def dimensions(self) -> tuple[int, int]:
        """(width, height) of level 0."""
        return self.levels[0].width, self.levels[0].height

    def level_mpp(self, level: int) -> float:
        return self.mpp * self.levels[level].downsample

    def summary(self, min_level_size: int = 256, max_levels: int | None = None) -> str:
        lines = [
            f"File          : {self.path}",
            f"Version       : {self.version:.3g}   codec: {self.codec}",
            f"Dimensions    : {self.width} x {self.height} px  @ {self.magnification}x",
            f"MPP           : {self.mpp:.5g} um/px"
            + ("" if self.mpp_from_header else " (fallback)"),
            f"Channels ({self.n_channels}) : "
            + ", ".join(f"{n} #{r:02X}{g:02X}{b:02X}" for n, (r, g, b)
                        in zip(self.channel_names, self.channel_colors)),
            f"Z planes      : {self.n_z} {self.z_values if self.n_z > 1 else ''}",
            f"Tile records  : {len(self._records)} (header: {self.n_records}), "
            f"length fmt {self._len_fmt[-1]}",
            f"Metadata tags : " + ", ".join(f"0x{t:X}" for t in sorted(self.tlv)),
            "Levels:",
        ]
        exported = {lv.index for lv in select_export_levels(self.levels, min_level_size,
                                                             max_levels)}
        lines[-1] = f"Levels ('*' = exported, min size {min_level_size}px):"
        for lv in self.levels:
            lines.append(
                f" {'*' if lv.index in exported else ' '}[{lv.index:2d}] "
                f"{lv.width:>7d} x {lv.height:<7d} scale {lv.scale:<9.4g} "
                f"ds {lv.downsample:<8g}"
                + (" (synthetic)" if lv.synthetic else
                   f" grid {lv.grid_cols}x{lv.grid_rows} coords={lv.coord_mode}"
                   f" raw+08/+0C={lv.raw_grid[0]},{lv.raw_grid[1]}"))
        return "\n".join(lines)

    # ---- pixel access ----------------------------------------------------
    def _payload(self, rec_index: int, channel: int) -> bytes | None:
        rec = self._records[rec_index]
        po, pl = int(rec["ptr_off"]), int(rec["ptr_len"])
        c = self.n_channels
        lsize = struct.calcsize(self._len_fmt)
        if po + 8 * c > self.file_size or pl + lsize > self.file_size:
            return None
        off = struct.unpack_from("<Q", self._mm, po + 8 * channel)[0]
        ln = struct.unpack_from(self._len_fmt, self._mm, pl)[channel]
        if ln == 0 or off + ln > self.file_size:
            return None
        return self._mm[off:off + ln]

    def get_tile(self, level: int, col: int, row: int, channel: int, z: int = 0) -> np.ndarray:
        """Return one decoded 256x256 uint8 source tile.

        Missing / out-of-range tiles return zeros. ``z`` is the Z-plane
        *index* (0..n_z-1).
        """
        if not 0 <= channel < self.n_channels:
            raise IndexError(f"channel {channel} out of range (0..{self.n_channels - 1})")
        if not 0 <= level < self.n_source_levels:
            raise IndexError(f"level {level} is not a source level "
                             f"(0..{self.n_source_levels - 1}); use read_region()")
        key = (level, z, col, row, channel)
        with self._lock:
            hit = self._cache.get(key)
            if hit is not None:
                self._cache.move_to_end(key)
                return hit
        ri = self._lookup.get((level, z, col, row))
        tile = np.zeros((SRC_TILE, SRC_TILE), np.uint8)
        if ri is not None:
            data = self._payload(ri, channel)
            if data is not None:
                try:
                    arr = decode_jpeg(data)
                    h, w = min(arr.shape[0], SRC_TILE), min(arr.shape[1], SRC_TILE)
                    tile[:h, :w] = arr[:h, :w]
                except Exception as exc:
                    with self._lock:
                        self.decode_failures += 1
                        first = self.decode_failures == 1
                    if first:
                        log.warning("%s: tile decode failed (L%d c%d r%d ch%d): %s",
                                    self.path.name, level, col, row, channel, exc)
        tile.flags.writeable = False
        if self._cache_size:
            with self._lock:
                self._cache[key] = tile
                if len(self._cache) > self._cache_size:
                    self._cache.popitem(last=False)
        return tile

    def read_region(self, level: int, x: int, y: int, width: int, height: int,
                    channel: int, z: int = 0) -> np.ndarray:
        """Read a (height, width) uint8 region of one channel at ``level``.

        Pixels outside the level bounds are zero.
        """
        lv = self.levels[level]
        out = np.zeros((height, width), np.uint8)
        x0, y0 = max(x, 0), max(y, 0)
        x1, y1 = min(x + width, lv.width), min(y + height, lv.height)
        if x1 <= x0 or y1 <= y0:
            return out
        if lv.synthetic:
            pw, ph = 2 * (x1 - x0), 2 * (y1 - y0)
            parent = self.levels[level - 1]
            src = self.read_region(level - 1, 2 * x0, 2 * y0, pw, ph, channel, z)
            vw, vh = min(pw, parent.width - 2 * x0), min(ph, parent.height - 2 * y0)
            if vw < pw or vh < ph:  # replicate edge instead of averaging with zeros
                src = np.pad(src[:vh, :vw], ((0, ph - vh), (0, pw - vw)), mode="edge")
            small = src.reshape(ph // 2, 2, pw // 2, 2).mean(axis=(1, 3), dtype=np.float32)
            out[y0 - y:y1 - y, x0 - x:x1 - x] = (small + 0.5).astype(np.uint8)
            return out
        for row in range(y0 // SRC_TILE, (y1 - 1) // SRC_TILE + 1):
            ty = row * SRC_TILE
            sy0, sy1 = max(y0, ty), min(y1, ty + SRC_TILE)
            for col in range(x0 // SRC_TILE, (x1 - 1) // SRC_TILE + 1):
                tx = col * SRC_TILE
                sx0, sx1 = max(x0, tx), min(x1, tx + SRC_TILE)
                tile = self.get_tile(level, col, row, channel, z)
                out[sy0 - y:sy1 - y, sx0 - x:sx1 - x] = tile[sy0 - ty:sy1 - ty, sx0 - tx:sx1 - tx]
        return out

    def read_region_channels(self, level: int, x: int, y: int, width: int, height: int,
                             z: int = 0) -> np.ndarray:
        """Read all channels -> (C, height, width) uint8."""
        return np.stack([self.read_region(level, x, y, width, height, c, z)
                         for c in range(self.n_channels)])

    def channel_histograms(self, level: int | None = None, max_pixels: int = 4_000_000
                           ) -> np.ndarray:
        """Exact per-channel 256-bin histograms over a low-resolution level
        (all Z planes). Returns (C, 256) int64."""
        if level is None:
            level = next((lv.index for lv in self.levels
                          if lv.width * lv.height <= max_pixels), self.n_levels - 1)
        lv = self.levels[level]
        hist = np.zeros((self.n_channels, 256), np.int64)
        step = 1024
        for z in range(self.n_z):
            for c in range(self.n_channels):
                for y in range(0, lv.height, step):
                    for x in range(0, lv.width, step):
                        w, h = min(step, lv.width - x), min(step, lv.height - y)
                        hist[c] += np.bincount(
                            self.read_region(level, x, y, w, h, c, z).ravel(), minlength=256)
        return hist


# --------------------------------------------------------------------------
# Writers
# --------------------------------------------------------------------------

_COMPRESSION = {"deflate": "zlib", "zlib": "zlib", "lzw": "lzw", "jpeg": "jpeg", "none": None}


def select_export_levels(levels: Sequence[Level], min_level_size: int = 256,
                         max_levels: int | None = None) -> list[Level]:
    """Practical pyramid: level 0 plus every level whose larger side is at least
    ``min_level_size`` px (KFBio stores levels down to 1x1, which Bio-Formats /
    QuPath handle poorly), optionally capped at ``max_levels``."""
    out = [lv for lv in levels if lv.index == 0 or max(lv.width, lv.height) >= min_level_size]
    return out[:max_levels] if max_levels else out


def _ordered_map(fn: Callable, items: Iterable, workers: int) -> Iterator:
    """Like ThreadPoolExecutor.map but with a bounded in-flight window."""
    if workers <= 1:
        for it in items:
            yield fn(it)
        return
    window = workers * 4
    with ThreadPoolExecutor(workers) as ex:
        pending: deque = deque()
        for it in items:
            pending.append(ex.submit(fn, it))
            if len(pending) >= window:
                yield pending.popleft().result()
        while pending:
            yield pending.popleft().result()


def percentile_stretch(hist: np.ndarray, low_pct: float, high_pct: float,
                       ignore_zero: bool = True) -> list[tuple[int, int]]:
    """Return per-channel (lo, hi) intensities from 256-bin histograms."""
    limits = []
    for h in hist:
        h = h.astype(np.float64).copy()
        if ignore_zero and h[1:].sum() > 0:
            h[0] = 0  # zero = background / missing tiles
        total = h.sum()
        if total == 0:
            limits.append((0, 255))
            continue
        cdf = np.cumsum(h) / total
        lo = int(np.searchsorted(cdf, low_pct / 100.0))
        hi = int(np.searchsorted(cdf, high_pct / 100.0))
        lo, hi = min(lo, 254), min(max(hi, lo + 1), 255)
        limits.append((lo, hi))
    return limits


def _composite_luts(limits, colors) -> np.ndarray:
    """(C, 256, 3) float32 LUTs: stretched intensity * channel colour."""
    v = np.arange(256, dtype=np.float32)
    luts = []
    for (lo, hi), rgb in zip(limits, colors):
        s = np.clip((v - lo) / float(hi - lo), 0.0, 1.0)
        luts.append(s[:, None] * np.asarray(rgb, np.float32)[None, :])
    return np.stack(luts)


class _PyramidWriter:
    """Shared machinery: writes one pyramidal BigTIFF from a tile function."""

    def __init__(self, slide: KFBioSlide, tile_size: int, compression: str,
                 workers: int, jpeg_quality: int = 90, show_progress: bool = True,
                 levels: Sequence[Level] | None = None):
        if tile_size % 16 or tile_size < 16:
            raise ValueError("tile size must be a positive multiple of 16")
        if compression not in _COMPRESSION:
            raise ValueError(f"unknown compression {compression!r}")
        self.slide = slide
        self.levels = list(levels) if levels is not None else select_export_levels(slide.levels)
        self.ts = tile_size
        self.compression = _COMPRESSION[compression]
        self.workers = max(1, workers)
        self.jpeg_quality = jpeg_quality
        self.show_progress = show_progress

    def _compression_kwargs(self, rgb: bool) -> dict:
        kw: dict = {"compression": self.compression}
        if self.compression == "jpeg":
            kw["compressionargs"] = {"level": self.jpeg_quality}
        elif self.compression in ("zlib", "lzw"):
            kw["predictor"] = True
            if self.compression == "zlib":
                kw["compressionargs"] = {"level": 6}
        return kw

    def write(self, out_path: Path, *, planes: Sequence, lead_shape: tuple,
              sample_shape: tuple,
              fetch: Callable[[int, object, int, int], np.ndarray],
              photometric: str, metadata: dict | None, ome: bool, desc: str) -> Path:
        """fetch(level, plane, tx, ty) -> tile array (ts, ts, *sample_shape)."""
        slide, ts = self.slide, self.ts
        levels = self.levels
        lead = tuple(d for d in lead_shape if d > 1)
        assert math.prod(lead) == len(planes), (lead, len(planes))
        total = sum(len(planes) * math.ceil(lv.width / ts) * math.ceil(lv.height / ts)
                    for lv in levels)
        tmp = out_path.with_name(out_path.name + ".part")
        tile_shape = (ts, ts) + sample_shape
        bar = tqdm(total=total, desc=desc, unit="tile", leave=False,
                   disable=not self.show_progress, dynamic_ncols=True)
        try:
            with tifffile.TiffWriter(tmp, bigtiff=True, ome=ome) as tif:
                for lv in levels:
                    ncols, nrows = math.ceil(lv.width / ts), math.ceil(lv.height / ts)
                    jobs = ((p, tx, ty) for p in planes
                            for ty in range(nrows) for tx in range(ncols))

                    def make(job, _lv=lv.index):
                        p, tx, ty = job
                        tile = fetch(_lv, p, tx, ty)
                        if tile.shape != tile_shape:  # pad edge tiles
                            padded = np.zeros(tile_shape, np.uint8)
                            padded[:tile.shape[0], :tile.shape[1]] = tile
                            tile = padded
                        return tile

                    def tiles(jobs=jobs, make=make):
                        for t in _ordered_map(make, jobs, self.workers):
                            bar.update(1)
                            yield t

                    res = 1e4 / slide.level_mpp(lv.index)  # pixels per cm
                    kw = dict(
                        shape=lead + (lv.height, lv.width) + sample_shape,
                        dtype=np.uint8, tile=(ts, ts), photometric=photometric,
                        resolution=(res, res), resolutionunit="CENTIMETER",
                        **self._compression_kwargs(photometric == "rgb"),
                    )
                    if photometric == "rgb":
                        kw["planarconfig"] = "contig"
                    if lv is levels[0]:
                        kw["subifds"] = len(levels) - 1
                        kw["metadata"] = metadata
                        kw["software"] = f"kfbf_converter {__version__}"
                    else:
                        kw["subfiletype"] = 1
                        kw["metadata"] = None
                    # tifffile reads ahead buffersize bytes of tiles to compress
                    # in parallel (default 512 MB); bound it to keep RAM flat.
                    # Older tifffile lacks `buffersize`: compress serially there,
                    # which streams tile-by-tile without read-ahead.
                    if _WRITE_HAS_BUFFERSIZE:
                        kw.update(maxworkers=self.workers, buffersize=WRITE_BUFFER_BYTES)
                    else:
                        kw.update(maxworkers=1)
                    tif.write(tiles(), **kw)
            os.replace(tmp, out_path)
        except BaseException:
            try:
                tmp.unlink()
            except OSError:
                pass
            raise
        finally:
            bar.close()
        return out_path


def _region(slide: KFBioSlide, level: int, tx: int, ty: int, ts: int) -> tuple[int, int, int, int]:
    lv = slide.levels[level]
    x, y = tx * ts, ty * ts
    return x, y, min(ts, lv.width - x), min(ts, lv.height - y)


def write_ome(slide: KFBioSlide, out_path: Path, writer: _PyramidWriter) -> Path:
    """Multi-channel pyramidal OME-TIFF: (C,Y,X) or (Z,C,Y,X) per level."""
    planes = [(z, c) for z in range(slide.n_z) for c in range(slide.n_channels)]
    axes = ("ZCYX" if slide.n_z > 1 else "CYX") if slide.n_channels > 1 else \
        ("ZYX" if slide.n_z > 1 else "YX")
    metadata = {
        "axes": axes,
        "Name": slide.path.stem,
        "PhysicalSizeX": slide.mpp, "PhysicalSizeXUnit": "µm",
        "PhysicalSizeY": slide.mpp, "PhysicalSizeYUnit": "µm",
        "Channel": {
            "Name": list(slide.channel_names),
            "Color": [ome_color(c) for c in slide.channel_colors],
        },
    }
    if slide.n_z > 1:
        # Z spacing is not stored in the header; z indices are written as-is
        metadata["Plane"] = {"TheZ": [z for z, _ in planes], "TheC": [c for _, c in planes],
                             "TheT": [0] * len(planes)}

    def fetch(level, plane, tx, ty):
        z, c = plane
        x, y, w, h = _region(slide, level, tx, ty, writer.ts)
        return slide.read_region(level, x, y, w, h, c, z)

    return writer.write(out_path, planes=planes, lead_shape=(slide.n_z, slide.n_channels),
                        sample_shape=(), fetch=fetch,
                        photometric="minisblack", metadata=metadata, ome=True,
                        desc=f"OME  {slide.path.stem}")


def write_rgb(slide: KFBioSlide, out_path: Path, writer: _PyramidWriter,
              low_pct: float = 0.1, high_pct: float = 99.9) -> Path:
    """Composite 8-bit RGB pyramidal BigTIFF with per-channel percentile stretch.

    Stretch limits are computed once per channel from an exact histogram of a
    low-resolution pyramid level, so every tile uses identical limits (no
    tile-seam brightness jumps).
    """
    limits = percentile_stretch(slide.channel_histograms(), low_pct, high_pct)
    for name, (lo, hi) in zip(slide.channel_names, limits):
        log.info("%s: stretch %s -> [%d, %d]", slide.path.stem, name, lo, hi)
    luts = _composite_luts(limits, slide.channel_colors)

    def fetch(level, z, tx, ty):
        x, y, w, h = _region(slide, level, tx, ty, writer.ts)
        acc = np.zeros((h, w, 3), np.float32)
        for c in range(slide.n_channels):
            acc += luts[c][slide.read_region(level, x, y, w, h, c, z)]
        np.clip(acc, 0, 255, out=acc)
        return (acc + 0.5).astype(np.uint8)

    metadata = None  # Z planes (if any) are written as consecutive pages
    return writer.write(out_path, planes=list(range(slide.n_z)), lead_shape=(slide.n_z,),
                        sample_shape=(3,),
                        fetch=fetch, photometric="rgb", metadata=metadata, ome=False,
                        desc=f"RGB  {slide.path.stem}")


def safe_channel_names(names: Sequence[str]) -> list[str]:
    out, seen = [], set()
    for i, n in enumerate(names):
        s = re.sub(r"[^A-Za-z0-9._-]+", "_", n).strip("._") or f"ch{i}"
        base, k = s, 2
        while s.lower() in seen:
            s, k = f"{base}_{k}", k + 1
        seen.add(s.lower())
        out.append(s)
    return out


def write_split(slide: KFBioSlide, out_path: Path, channel: int,
                writer: _PyramidWriter) -> Path:
    """Single-channel grayscale pyramidal BigTIFF."""

    def fetch(level, z, tx, ty):
        x, y, w, h = _region(slide, level, tx, ty, writer.ts)
        return slide.read_region(level, x, y, w, h, channel, z)

    metadata = None  # Z planes (if any) are written as consecutive pages
    return writer.write(out_path, planes=list(range(slide.n_z)), lead_shape=(slide.n_z,),
                        sample_shape=(),
                        fetch=fetch, photometric="minisblack", metadata=metadata,
                        ome=False, desc=f"{slide.channel_names[channel]:<5.5s} {slide.path.stem}")


# --------------------------------------------------------------------------
# High-level conversion
# --------------------------------------------------------------------------

MODES = ("ome", "rgb", "split", "all")


def _slide_stem(path: Path) -> str:
    return path.name[: -len(path.suffix)] if path.suffix else path.name


def planned_outputs(slide: KFBioSlide, output_dir: Path, mode: str) -> list[tuple[str, Path, int]]:
    stem = _slide_stem(slide.path)
    jobs: list[tuple[str, Path, int]] = []
    if mode in ("ome", "all"):
        jobs.append(("ome", output_dir / f"{stem}.ome.tif", -1))
    if mode in ("rgb", "all"):
        jobs.append(("rgb", output_dir / f"{stem}.rgb.tif", -1))
    if mode in ("split", "all"):
        for i, n in enumerate(safe_channel_names(slide.channel_names)):
            jobs.append(("split", output_dir / f"{stem}_{n}.tif", i))
    return jobs


def convert_slide(path: str | os.PathLike, output_dir: str | os.PathLike | None = None,
                  mode: str = "ome", compression: str = "deflate", tile_size: int = 512,
                  mpp: float = DEFAULT_MPP, overwrite: bool = False, workers: int | None = None,
                  channel_names: Sequence[str] | None = None, low_pct: float = 0.1,
                  high_pct: float = 99.9, jpeg_quality: int = 90,
                  show_progress: bool = True, min_level_size: int = 256,
                  max_levels: int | None = None) -> list[Path]:
    """Convert one .kfbf slide. Returns the list of files written."""
    if mode not in MODES:
        raise ValueError(f"mode must be one of {MODES}")
    path = Path(path)
    out_dir = Path(output_dir) if output_dir else path.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    workers = workers or min(8, os.cpu_count() or 1)
    written: list[Path] = []
    with KFBioSlide(path, mpp_fallback=mpp, channel_names=channel_names) as slide:
        export = select_export_levels(slide.levels, min_level_size, max_levels)
        log.info("%s: %dx%d, %d ch %s, %d Z, exporting %d of %d levels (%s)", path.name,
                 slide.width, slide.height, slide.n_channels, slide.channel_names,
                 slide.n_z, len(export), slide.n_levels,
                 ", ".join(f"{lv.width}x{lv.height}" for lv in export))
        writer = _PyramidWriter(slide, tile_size, compression, workers, jpeg_quality,
                                show_progress, levels=export)
        for kind, out, ch in planned_outputs(slide, out_dir, mode):
            if out.exists() and not overwrite:
                log.info("skip existing %s (use --overwrite)", out)
                continue
            if kind == "ome":
                write_ome(slide, out, writer)
            elif kind == "rgb":
                write_rgb(slide, out, writer, low_pct, high_pct)
            else:
                write_split(slide, out, ch, writer)
            log.info("wrote %s", out)
            written.append(out)
        if slide.decode_failures:
            log.warning("%s: %d tile(s) failed to decode and were written as zeros",
                        path.name, slide.decode_failures)
    return written


def find_slides(input_path: Path, recursive: bool = False) -> list[Path]:
    if input_path.is_file():
        return [input_path]
    if input_path.is_dir():
        it = input_path.rglob("*") if recursive else input_path.iterdir()
        return sorted(p for p in it if p.is_file() and p.suffix.lower() == ".kfbf")
    raise FileNotFoundError(f"input not found: {input_path}")


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

class _TqdmHandler(logging.Handler):
    def emit(self, record: logging.LogRecord) -> None:
        try:
            tqdm.write(self.format(record), file=sys.stderr)
        except Exception:  # pragma: no cover
            self.handleError(record)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="kfbf_converter",
        description="Convert KFBio fluorescent slides (.kfbf) to pyramidal OME-TIFF, "
                    "composite RGB BigTIFF, and/or per-channel grayscale BigTIFFs.")
    p.add_argument("input", type=Path, nargs="?", default=None,
                   help=".kfbf file or directory containing .kfbf files "
                        "(omit to pick a format interactively for the current folder)")
    p.add_argument("-m", "--mode", choices=MODES, default=None,
                   help="output type (default: ome)")
    p.add_argument("-o", "--output-dir", type=Path, default=None,
                   help="output directory (default: alongside the input)")
    p.add_argument("--compression", choices=["deflate", "lzw", "jpeg"], default="deflate",
                   help="TIFF tile compression (default: deflate)")
    p.add_argument("--jpeg-quality", type=int, default=90, help="quality for --compression jpeg")
    p.add_argument("--tile-size", type=int, default=512, help="output tile size (default: 512)")
    p.add_argument("--mpp", type=float, default=DEFAULT_MPP,
                   help=f"fallback µm/px if missing from header (default: {DEFAULT_MPP})")
    p.add_argument("--overwrite", action="store_true", help="overwrite existing outputs")
    p.add_argument("--recursive", action="store_true", help="search sub-directories for .kfbf")
    p.add_argument("--workers", type=int, default=None,
                   help="tile decode threads (default: min(8, CPUs))")
    p.add_argument("--channel-names", nargs="+", default=None, metavar="NAME",
                   help="override channel names (must match channel count)")
    p.add_argument("--low-pct", type=float, default=0.1, help="RGB stretch low percentile")
    p.add_argument("--high-pct", type=float, default=99.9, help="RGB stretch high percentile")
    p.add_argument("--min-level-size", type=int, default=256,
                   help="drop pyramid levels whose larger side is below this (default: 256)")
    p.add_argument("--max-levels", type=int, default=None,
                   help="cap the number of exported pyramid levels (incl. full resolution)")
    p.add_argument("--info", action="store_true", help="print parsed slide structure and exit")
    p.add_argument("-q", "--quiet", action="store_true", help="hide progress bars")
    p.add_argument("-v", "--verbose", action="store_true", help="debug logging")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    return p


_MENU = [("ome", "OME-TIFF (Default / QuPath)"),
         ("rgb", "Composite RGB BigTIFF"),
         ("split", "Split Grayscale Channels"),
         ("all", "All of the above")]


def prompt_mode(count: int) -> str | None:
    """Interactive output-format menu. Returns a mode, or None if cancelled."""
    print(f"Found {count} slide(s).")
    print("Select output format:")
    for i, (_, label) in enumerate(_MENU, 1):
        print(f"  [{i}] {label}")
    attempted = False
    while True:
        try:
            choice = input("Enter choice [1-4] (Press Enter for 1): ")
        except EOFError:
            print()
            if attempted:  # input ran out after an invalid answer: don't guess
                print("No valid choice entered; cancelled.")
                return None
            return _MENU[0][0]  # non-interactive stdin with no answer: default
        except KeyboardInterrupt:
            print("\nCancelled.")
            return None
        # strip a byte-order mark (e.g. text piped from PowerShell)
        for bom in (chr(0xFEFF), chr(0xEF) + chr(0xBB) + chr(0xBF)):
            choice = choice.replace(bom, "")
        choice = choice.strip()
        attempted = True
        if choice == "":
            return _MENU[0][0]
        if choice in {"1", "2", "3", "4"}:
            return _MENU[int(choice) - 1][0]
        print("Invalid choice, please enter a number from 1 to 4.")


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    handler = _TqdmHandler()
    handler.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
    log.addHandler(handler)
    log.setLevel(logging.DEBUG if args.verbose else logging.INFO)
    log.propagate = False

    if args.tile_size < 16 or args.tile_size % 16:
        log.error("--tile-size must be a positive multiple of 16")
        return 2
    if not 0 <= args.low_pct < args.high_pct <= 100:
        log.error("require 0 <= --low-pct < --high-pct <= 100")
        return 2

    interactive = args.input is None
    if interactive:
        args.input = Path.cwd()
    try:
        slides = find_slides(args.input, args.recursive)
    except FileNotFoundError as exc:
        log.error("%s", exc)
        return 2
    if interactive and not slides:
        print(f"No .kfbf files found in the current folder ({args.input}).")
        return 0
    if interactive and args.mode is None and not args.info:
        args.mode = prompt_mode(len(slides))
        if args.mode is None:
            return 130
    args.mode = args.mode or "ome"
    if not slides:
        log.error("no .kfbf files found in %s", args.input)
        return 1

    if args.info:
        rc = 0
        for s in slides:
            try:
                with KFBioSlide(s, args.mpp, args.channel_names) as slide:
                    print(slide.summary(args.min_level_size, args.max_levels), end="\n\n")
            except Exception as exc:
                log.error("%s: %s", s.name, exc)
                rc = 1
        return rc

    failures: list[tuple[Path, str]] = []
    n_written = 0
    for s in tqdm(slides, desc="Slides", unit="slide", disable=args.quiet or len(slides) < 2,
                  dynamic_ncols=True):
        out_dir = args.output_dir or s.parent
        try:
            n_written += len(convert_slide(
                s, out_dir, mode=args.mode, compression=args.compression,
                tile_size=args.tile_size, mpp=args.mpp, overwrite=args.overwrite,
                workers=args.workers, channel_names=args.channel_names,
                low_pct=args.low_pct, high_pct=args.high_pct,
                jpeg_quality=args.jpeg_quality, show_progress=not args.quiet,
                min_level_size=args.min_level_size, max_levels=args.max_levels))
        except KeyboardInterrupt:
            log.error("interrupted")
            return 130
        except Exception as exc:
            log.error("%s: conversion failed: %s", s.name, exc,
                      exc_info=args.verbose)
            failures.append((s, str(exc)))

    log.info("done: %d slide(s), %d file(s) written, %d failure(s)",
             len(slides), n_written, len(failures))
    for s, msg in failures:
        log.info("  FAILED %s: %s", s, msg)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
