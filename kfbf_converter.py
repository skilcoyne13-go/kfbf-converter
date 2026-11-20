#!/usr/bin/env python3
"""
kfbf_converter.py - Convert whole-slide images (KFBio .kfbf, Aperio .svs,
PerkinElmer/Akoya .qptiff and plain .tif/.tiff) into open, QuPath-ready
pyramidal formats.

Created by Scott Kilcoyne, DigitCells.
Copyright (c) 2026 DigitCells

Inputs and what each becomes
----------------------------
.kfbf       KFBio fluorescence slide      -> OME-TIFF / RGB / split (see modes)
.svs        Aperio brightfield slide      -> RGB pyramidal BigTIFF   <stem>.tif
.qptiff     PerkinElmer/Akoya QPTIFF      -> multi-channel pyramidal OME-TIFF
                                             <stem>.ome.tif (an RGB QPTIFF
                                             gets the RGB BigTIFF instead)
.tif/.tiff  plain RGB or grayscale image  -> pyramidal BigTIFF <stem>.tif, or
                                             <stem>.bigtiff.tif if that name
                                             is taken (e.g. by the input)

Given a folder, the tool counts each type it finds and asks what to convert
(--types / --yes skip the question). Files this tool wrote are never taken
as inputs.

Output modes (fluorescence .kfbf)
---------------------------------
ome    Multi-channel pyramidal OME-TIFF          <stem>.ome.tif
rgb    Composite 8-bit RGB pyramidal BigTIFF     <stem>.rgb.tif
split  One grayscale pyramidal BigTIFF/channel   <stem>_<channel>.tif
all    All of the above

Brightfield .svs slides always produce one RGB pyramidal BigTIFF,
<stem>.tif. By default the stored JPEG tiles are copied unchanged (no extra
quality loss); --compression / --tile-size re-encode them instead.

All writers stream tiles level-by-level straight into ``tifffile.TiffWriter``;
full-resolution planes are never assembled in memory. Each output is checked
before writing ("Data is too big for the file format chosen."), written to
<name>.part, verified by reading it back, then renamed.

Usage
-----
    python kfbf_converter.py slide.kfbf --mode all -o out/
    python kfbf_converter.py slide.svs -o out/
    python kfbf_converter.py slide.qptiff -o out/
    python kfbf_converter.py slides_folder --types qptiff tif --yes -o out/
    python kfbf_converter.py slide.kfbf --info        # dump parsed structure

As a module
-----------
    from kfbf_converter import KFBioSlide, SVSSlide, open_slide, convert_slide
    with open_slide("slide.qptiff") as slide:       # -> QPTIFFSlide
        dapi = slide.read_region(level=3, x=0, y=0, width=512, height=512, channel=0)
    with KFBioSlide("slide.kfbf") as slide:
        tile = slide.get_tile(level=0, col=3, row=5, channel=0)
    with open_slide("slide.svs") as slide:          # -> SVSSlide
        rgb = slide.read_region(level=2, x=0, y=0, width=512, height=512)
    convert_slide("slide.kfbf", "out", mode="all")

Dependencies: numpy, tifffile, imagecodecs, pillow, tqdm
"""

from __future__ import annotations

import argparse
import errno
import io
import logging
import math
import mmap
import os
import re
import shutil
import struct
import sys
import threading
import xml.etree.ElementTree as ElementTree
from collections import Counter, OrderedDict, deque
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

__version__ = "1.1.0"
__author__ = "Scott Kilcoyne"
__company__ = "DigitCells"

CREDIT = f"{__author__}, {__company__}"           # TIFF Artist tag, --version
MENU_BANNER = f"KFBio Slide Converter - {__company__}"

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


class SlideFormatError(KFBFError):
    """Raised for unreadable, malformed or unsupported slide files."""


class DataTooBigError(KFBFError):
    """An output cannot be written: file-format, file-system or free-space
    limit. ``lines()`` gives the three user-facing lines."""

    HEADLINE = "Data is too big for the file format chosen."

    def __init__(self, reason: str, fix: str):
        self.reason, self.fix = reason, fix
        super().__init__(f"{self.HEADLINE} {reason} {fix}")

    def lines(self) -> list[str]:
        return [self.HEADLINE, f"  {self.reason}", f"  Fix: {self.fix}"]


class VerificationError(KFBFError):
    """A written output could not be read back correctly (it was deleted)."""

    PREFIX = "Output failed verification and was deleted: "

    def __init__(self, reason: str):
        self.reason = reason
        super().__init__(self.PREFIX + reason)

    def lines(self) -> list[str]:
        return [str(self)]


def jpeg_sof_info(data: bytes) -> dict | None:
    """Parse JPEG markers up to SOS: size, components, Adobe/JFIF presence.

    Returns ``{"width", "height", "components": [(id, h, v), ...],
    "adobe_transform": int | None, "jfif": bool}`` or None if no SOF found.
    """
    info: dict = {"adobe_transform": None, "jfif": False}
    i, n = 2, len(data)
    if data[:2] != b"\xff\xd8":
        return None
    while i + 4 <= n and data[i] == 0xFF:
        m = data[i + 1]
        if m == 0xFF:  # fill byte
            i += 1
            continue
        if m in (0xD8, 0x01) or 0xD0 <= m <= 0xD7:
            i += 2
            continue
        seg_len = int.from_bytes(data[i + 2:i + 4], "big")
        seg = data[i + 4:i + 2 + seg_len]
        if m == 0xE0 and seg[:5] == b"JFIF\0":
            info["jfif"] = True
        elif m == 0xEE and seg[:5] == b"Adobe" and len(seg) >= 12:
            info["adobe_transform"] = seg[11]
        elif m in (0xC0, 0xC1, 0xC2, 0xC3) and len(seg) >= 6:
            nc = seg[5]
            info["height"] = int.from_bytes(seg[1:3], "big")
            info["width"] = int.from_bytes(seg[3:5], "big")
            info["components"] = [(seg[6 + 3 * k], seg[7 + 3 * k] >> 4, seg[7 + 3 * k] & 15)
                                  for k in range(nc) if 8 + 3 * k <= len(seg)]
            return info
        elif m == 0xDA:
            break
        i += 2 + seg_len
    return None


def _merge_jpeg_tables(tables: bytes | None, data: bytes) -> bytes:
    """Splice TIFF JPEGTables (abbreviated table spec) into an abbreviated tile."""
    if not tables:
        return data
    return bytes(tables[:-2]) + bytes(data[2:])  # drop tables' EOI and tile's SOI


def decode_jpeg_rgb(data: bytes, tables: bytes | None = None,
                    rgb_stored: bool = False) -> np.ndarray:
    """Decode a colour JPEG payload to an (h, w, 3) uint8 RGB array.

    Unlike :func:`decode_jpeg`, colour is kept. ``rgb_stored`` means the
    container says the samples were stored as RGB (TIFF photometric RGB),
    i.e. *no* YCbCr transform, even when the stream has no Adobe marker, as
    in some Aperio files. Otherwise the JPEG's own markers decide (JFIF /
    component IDs 1-2-3 = YCbCr, Adobe transform 0 or IDs 'R','G','B' = RGB).
    imagecodecs is used first, Pillow as a fallback.
    """
    arr = None
    if imagecodecs is not None:
        kw = {"colorspace": "RGB", "outcolorspace": "RGB"} if rgb_stored else {}
        try:
            arr = imagecodecs.jpeg_decode(data, tables=tables or None, **kw)
        except Exception:
            arr = None
    if arr is None:
        if Image is None:
            raise KFBFError("neither imagecodecs nor Pillow could decode JPEG")
        with Image.open(io.BytesIO(_merge_jpeg_tables(tables, data))) as im:
            if rgb_stored and im.mode == "RGB":
                # Without an Adobe marker libjpeg would treat the stored RGB
                # samples as YCbCr and convert them; asking for "YCbCr"
                # output returns the stored component values unchanged.
                im.draft("YCbCr", None)
            im.load()
            if im.mode != "RGB" and not (rgb_stored and im.mode == "YCbCr"):
                im = im.convert("RGB")
            arr = np.asarray(im)
    if arr.ndim == 2:
        arr = np.repeat(arr[..., None], 3, axis=2)
    elif arr.shape[2] == 1:
        arr = np.repeat(arr, 3, axis=2)
    elif arr.shape[2] > 3:
        arr = arr[..., :3]
    if arr.dtype != np.uint8:
        arr = np.clip(arr, 0, 255).astype(np.uint8)
    return arr


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


MIN_MAGNIFICATION, MAX_MAGNIFICATION = 0.5, 250.0


def format_magnification(value: float) -> str:
    """40.0 -> '40', 2.5 -> '2.5'."""
    return f"{value:g}"


class _SlideBase:
    """Shared reader machinery: read-only mmap, lifecycle, LRU tile cache,
    synthetic pyramid levels and tiled region reads.

    Subclasses set ``levels``, ``n_source_levels``, ``z_values``, ``mpp`` and
    implement their own ``get_tile`` / ``read_region``.
    """

    kind = ""              # "fluorescence" | "brightfield"
    is_rgb = False
    fill_value = 0         # value of missing tiles / out-of-bounds pixels
    tile_size = SRC_TILE   # source tile size
    dtype = np.dtype(np.uint8)  # pixel type (QPTIFF / TIFF may be 16-bit)
    input_type = ""        # "KFBF" | "SVS" | "QPTIFF" | "TIFF": what the user gave us
    min_file_size = 8

    def __init__(self, path: str | os.PathLike, tile_cache: int = 1024):
        self.path = Path(path)
        self._fh = open(self.path, "rb")
        try:
            self.file_size = os.fstat(self._fh.fileno()).st_size
            if self.file_size < self.min_file_size:
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

    # ---- lifecycle -------------------------------------------------------
    def close(self) -> None:
        mm = getattr(self, "_mm", None)
        if mm is not None and not mm.closed:
            mm.close()
        if not self._fh.closed:
            self._fh.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    def _warn(self, msg: str) -> None:
        self.warnings.append(msg)
        log.warning("%s: %s", self.path.name, msg)

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

    def _resolve_magnification(self, from_file, source: str,
                               override: float | None) -> None:
        """Objective magnification: ``override`` (--magnification) wins, else
        the file's value if plausible, else 0 (= unknown)."""
        if override is not None:
            self.magnification, self.magnification_source = float(override), "--magnification"
            return
        value = _float_or_none(from_file)
        if value is not None and MIN_MAGNIFICATION <= value <= MAX_MAGNIFICATION:
            self.magnification, self.magnification_source = value, source
        else:
            self.magnification, self.magnification_source = 0.0, "missing"
            if value:  # present but implausible
                self._warn(f"implausible magnification {value:g} in file; ignored "
                           f"(set --magnification)")

    def magnification_text(self) -> str:
        """e.g. '40x (header)' or 'unknown (set --magnification)'."""
        if not self.magnification:
            return "unknown (set --magnification)"
        return f"{format_magnification(self.magnification)}x ({self.magnification_source})"

    # ---- helpers ---------------------------------------------------------
    def _add_synthetic_levels(self, levels: list[Level]) -> list[Level]:
        """Append 2x-downsampled levels until the pyramid is small enough."""
        while max(levels[-1].width, levels[-1].height) > MIN_PYRAMID_DIM:
            p = levels[-1]
            levels.append(Level(len(levels), p.scale / 2, max(1, math.ceil(p.width / 2)),
                                max(1, math.ceil(p.height / 2)), 0, 0,
                                p.downsample * 2, synthetic=True))
        return levels

    def _cache_get(self, key):
        with self._lock:
            hit = self._cache.get(key)
            if hit is not None:
                self._cache.move_to_end(key)
            return hit

    def _cache_put(self, key, tile: np.ndarray) -> np.ndarray:
        tile.flags.writeable = False
        if self._cache_size:
            with self._lock:
                self._cache[key] = tile
                if len(self._cache) > self._cache_size:
                    self._cache.popitem(last=False)
        return tile

    def _note_decode_failure(self, what: str, exc: Exception) -> None:
        with self._lock:
            self.decode_failures += 1
            first = self.decode_failures == 1
        if first:
            log.warning("%s: tile decode failed (%s): %s", self.path.name, what, exc)

    def _read_plane(self, level: int, x: int, y: int, width: int, height: int,
                    tile_fn: Callable[[int, int, int], np.ndarray],
                    sample_shape: tuple = ()) -> np.ndarray:
        """Assemble a (height, width, *sample_shape) region from source tiles
        (``tile_fn(level, col, row)``) or, for synthetic levels, by 2x2
        averaging the parent level with edge replication."""
        out = np.full((height, width) + sample_shape, self.fill_value, self.dtype)
        lv = self.levels[level]
        x0, y0 = max(x, 0), max(y, 0)
        x1, y1 = min(x + width, lv.width), min(y + height, lv.height)
        if x1 <= x0 or y1 <= y0:
            return out
        if lv.synthetic:
            pw, ph = 2 * (x1 - x0), 2 * (y1 - y0)
            parent = self.levels[level - 1]
            src = self._read_plane(level - 1, 2 * x0, 2 * y0, pw, ph, tile_fn, sample_shape)
            vw, vh = min(pw, parent.width - 2 * x0), min(ph, parent.height - 2 * y0)
            if vw < pw or vh < ph:  # replicate edge instead of averaging with fill
                pad = ((0, ph - vh), (0, pw - vw)) + ((0, 0),) * len(sample_shape)
                src = np.pad(src[:vh, :vw], pad, mode="edge")
            small = src.reshape((ph // 2, 2, pw // 2, 2) + sample_shape).mean(
                axis=(1, 3), dtype=np.float32)
            out[y0 - y:y1 - y, x0 - x:x1 - x] = (small + 0.5).astype(self.dtype)
            return out
        ts = self.tile_size
        for row in range(y0 // ts, (y1 - 1) // ts + 1):
            ty = row * ts
            sy0, sy1 = max(y0, ty), min(y1, ty + ts)
            for col in range(x0 // ts, (x1 - 1) // ts + 1):
                tx = col * ts
                sx0, sx1 = max(x0, tx), min(x1, tx + ts)
                tile = tile_fn(level, col, row)
                out[sy0 - y:sy1 - y, sx0 - x:sx1 - x] = tile[sy0 - ty:sy1 - ty, sx0 - tx:sx1 - tx]
        return out


class KFBioSlide(_SlideBase):
    """Reader for KFBio fluorescent ``.kfbf`` slides.

    Parameters
    ----------
    path : path to the .kfbf file
    mpp_fallback : pixel size (µm/px) used when the header value is invalid
    channel_names : optional override of the channel labels
    tile_cache : number of decoded 256x256 tiles kept in an LRU cache
    """

    kind = "fluorescence"
    input_type = "KFBF"
    min_file_size = 0x60

    def __init__(self, path: str | os.PathLike, mpp_fallback: float = DEFAULT_MPP,
                 channel_names: Sequence[str] | None = None, tile_cache: int = 1024,
                 magnification: float | None = None):
        super().__init__(path, tile_cache)
        try:
            self._parse_header(mpp_fallback)
            self._resolve_magnification(self.header_magnification, "header", magnification)
            self._parse_metadata()
            self._parse_index()
            self._resolve_channels(channel_names)
            self._build_levels()
        except Exception:
            self.close()
            raise

    # ---- parsing ---------------------------------------------------------
    def _parse_header(self, mpp_fallback: float) -> None:
        mm = self._mm
        if mm[0:4] != MAGIC or mm[4:8] != SIGNATURE:
            raise KFBFError(f"not a KFBF file (bad magic/signature)")
        (self.version,) = struct.unpack_from("<f", mm, 0x0C)
        self.n_records, self.height, self.width, self.header_magnification = \
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
        if self.header_magnification and abs(s0 - self.header_magnification) > 0.5:
            log.debug("top level scale %.3f differs from header magnification %d",
                      s0, self.header_magnification)

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
        self.levels = self._add_synthetic_levels(levels)

    def summary(self, min_level_size: int = 256, max_levels: int | None = None) -> str:
        lines = [
            f"File          : {self.path}",
            f"Version       : {self.version:.3g}   codec: {self.codec}",
            f"Dimensions    : {self.width} x {self.height} px  @ "
            f"{format_magnification(self.magnification)}x",
            f"Magnification : {self.magnification_text()}",
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
        hit = self._cache_get(key)
        if hit is not None:
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
                    self._note_decode_failure(
                        f"L{level} c{col} r{row} ch{channel}", exc)
        return self._cache_put(key, tile)

    def read_region(self, level: int, x: int, y: int, width: int, height: int,
                    channel: int, z: int = 0) -> np.ndarray:
        """Read a (height, width) uint8 region of one channel at ``level``.

        Pixels outside the level bounds are zero.
        """
        return self._read_plane(level, x, y, width, height,
                                lambda lv, col, row: self.get_tile(lv, col, row, channel, z))

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


TIFF_MAGICS = (b"II*\0", b"MM\0*", b"II+\0", b"MM\0+")
_TIFF_JPEG, _PHOTOMETRIC_RGB, _PHOTOMETRIC_YCBCR = 7, 2, 6


def _parse_aperio_description(desc: str) -> dict[str, str]:
    """``'Aperio ... |AppMag = 20|MPP = 0.4859|...'`` -> {'AppMag': '20', ...}."""
    props: dict[str, str] = {}
    for part in desc.split("|")[1:]:
        key, sep, value = part.partition("=")
        if sep and key.strip():
            props[key.strip()] = value.strip()
    return props


def _float_or_none(text) -> float | None:
    try:
        v = float(text)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def _exact_step(stored: int, parent: int, ratio: int, tolerance: float = 0.02) -> int:
    """Size of a level that is exactly ``parent / ratio`` (floor or ceil).

    A larger stored size is padding and is always cropped (as before). A
    smaller one is extended with white only if it is within ``tolerance``
    (or 2 px) of the exact size; otherwise the stored size is kept."""
    lo, hi = parent // ratio, -(-parent // ratio)
    if stored > hi:
        return hi
    if stored >= lo or lo - stored > max(2, tolerance * lo):
        return stored
    return max(lo, 1)


@dataclass
class _SVSLevelInfo:
    """Where one SVS pyramid level's tiles live and how they are encoded."""

    page: object            # tifffile.TiffPage
    offsets: np.ndarray     # int64, row-major tile order
    counts: np.ndarray      # int64, 0 = missing tile
    compression: int
    photometric: int
    jpegtables: bytes | None
    rgb_stored: bool        # JPEG samples stored as RGB (no YCbCr transform)
    stride: int             # tiles per row in the stored page
    stored_size: tuple = ()  # (width, height) of the stored page


class SVSSlide(_SlideBase):
    """Reader for brightfield Aperio ``.svs`` slides (as written by KFBio).

    Same interface as :class:`KFBioSlide`, but tiles and regions are RGB:
    ``get_tile`` returns ``(tile_size, tile_size, 3)`` and ``read_region``
    ``(height, width, 3)`` uint8. Missing tiles and pixels outside the image
    are white (255), the colour of blank glass. Only the pyramid is read; the
    thumbnail, label and macro images are ignored.

    Parameters
    ----------
    path : path to the .svs file
    mpp_fallback : pixel size (µm/px) used when the file does not state one
    tile_cache : number of decoded RGB tiles kept in an LRU cache (default 256)
    """

    kind = "brightfield"
    input_type = "SVS"
    is_rgb = True
    fill_value = 255
    n_channels = 3
    channel_names = ["Red", "Green", "Blue"]
    channel_colors = [(255, 0, 0), (0, 255, 0), (0, 0, 255)]

    def __init__(self, path: str | os.PathLike, mpp_fallback: float = DEFAULT_MPP,
                 tile_cache: int = 256, magnification: float | None = None):
        # RGB tiles are 3x the size of fluorescence tiles: 256 x 192 KB = 48 MB
        super().__init__(path, tile_cache)
        self._tif = None
        try:
            if bytes(self._mm[:4]) not in TIFF_MAGICS:
                raise SlideFormatError("not a TIFF/SVS file (bad magic)")
            self._tif = tifffile.TiffFile(self.path)
            self._parse(mpp_fallback, magnification)
        except Exception:
            self.close()
            raise

    def close(self) -> None:
        tif = getattr(self, "_tif", None)
        if tif is not None:
            tif.close()
            self._tif = None
        super().close()

    # ---- parsing ---------------------------------------------------------
    def _parse(self, mpp_fallback: float, magnification: float | None = None) -> None:
        tif = self._tif
        self.is_bigtiff = tif.is_bigtiff
        self.is_aperio = tif.is_svs
        series = tif.series[0]
        base = series.levels[0].keyframe
        self.description = base.description or ""
        self.properties = (_parse_aperio_description(self.description)
                           if self.is_aperio else {})
        self.width, self.height = int(base.imagewidth), int(base.imagelength)
        self.codec = getattr(base.compression, "name", str(base.compression))
        self._resolve_magnification(self.properties.get("AppMag"), "Aperio AppMag",
                                    magnification)

        # Pixel size: Aperio "MPP =", else TIFF resolution tags, else fallback
        mpp = _float_or_none(self.properties.get("MPP"))
        self.mpp_source = "description"
        if mpp is None or not 0.01 <= mpp <= 100:
            mpp, self.mpp_source = self._mpp_from_resolution(base), "resolution tags"
        if mpp is None or not 0.01 <= mpp <= 100:
            mpp, self.mpp_source = float(mpp_fallback), "fallback"
            self._warn(f"no pixel size in file; using {mpp_fallback} um/px (set --mpp)")
        self.mpp = float(f"{mpp:.6g}")
        self.mpp_from_header = self.mpp_source != "fallback"

        if base.dtype != np.uint8 or base.samplesperpixel not in (3, 4):
            raise SlideFormatError(
                f"unsupported SVS pixel format ({base.samplesperpixel} samples of "
                f"{base.dtype}); only 8-bit RGB brightfield is supported")
        if not base.is_tiled:
            raise SlideFormatError("SVS level 0 is not tiled")
        if base.tilewidth != base.tilelength:
            raise SlideFormatError(f"non-square tiles ({base.tilewidth}x{base.tilelength})")
        if int(base.planarconfig) != 1:
            raise SlideFormatError("only contiguous (chunky) RGB SVS files are supported")
        self.tile_size = int(base.tilewidth)

        levels: list[Level] = []
        self._info: list[_SVSLevelInfo] = []
        bad = 0
        for lvs in series.levels:
            pg = lvs.keyframe
            if len(lvs.pages) > 1 and not levels:
                self._warn(f"{len(lvs.pages)} planes per level; only the first is read")
            w, h = int(pg.imagewidth), int(pg.imagelength)
            if not pg.is_tiled or pg.tilewidth != self.tile_size \
                    or pg.tilelength != self.tile_size or int(pg.planarconfig) != 1:
                self._warn(f"level {w}x{h} has a different tile layout; "
                           f"lower levels will be computed instead")
                break
            if levels and (w >= levels[-1].width or h >= levels[-1].height):
                self._warn(f"skipping non-decreasing level {w}x{h}")
                continue
            stored = (w, h)
            ds = (self.width / w + self.height / h) / 2
            log2 = math.log2(ds)
            if abs(log2 - round(log2)) < 0.15:  # Aperio pads level sizes
                ds = float(2 ** round(log2))
                if levels:
                    # KFBio level sizes are not exact steps: some are padded
                    # with white (1628x2305 stored as 1632x2320), some are a
                    # pixel or two short (2303 for half of 4609). Declare an
                    # exact step so readers see one pyramid: crop the padding,
                    # or extend short levels with white. Tiles are unchanged.
                    r = ds / levels[-1].downsample
                    if r == int(r):
                        w = _exact_step(w, levels[-1].width, int(r))
                        h = _exact_step(h, levels[-1].height, int(r))
            stride = math.ceil(stored[0] / self.tile_size)
            ncols, nrows = math.ceil(w / self.tile_size), math.ceil(h / self.tile_size)
            offsets = np.asarray(pg.dataoffsets, np.int64)
            counts = np.asarray(pg.databytecounts, np.int64)
            need = stride * math.ceil(stored[1] / self.tile_size)
            if len(offsets) < need or len(counts) != len(offsets):
                raise SlideFormatError(f"level {stored[0]}x{stored[1]}: {len(offsets)} "
                                       f"tiles, expected {need}")
            outside = (offsets + counts > self.file_size) | (offsets < 0)
            if outside.any():
                bad += int(outside.sum())
                counts = np.where(outside, 0, counts)
            comp, phot = int(pg.compression), int(pg.photometric)
            levels.append(Level(len(levels), (self.magnification or 1) / ds, w, h,
                                ncols, nrows, ds, coord_mode="tiff"))
            self._info.append(_SVSLevelInfo(
                pg, offsets[:need], counts[:need], comp, phot, pg.jpegtables or None,
                comp == _TIFF_JPEG and phot == _PHOTOMETRIC_RGB, stride, stored))
        if bad:
            self._warn(f"{bad} tile(s) point outside the file (truncated?); shown as white")
        self.z_values = [0]
        self.n_source_levels = len(levels)
        self.levels = self._add_synthetic_levels(levels)
        self.associated = [s.name for s in tif.series[1:]]

    @staticmethod
    def _mpp_from_resolution(page) -> float | None:
        try:
            xres = page.tags["XResolution"].value
            unit = int(page.tags["ResolutionUnit"].value)
        except (KeyError, TypeError, ValueError):
            return None
        ppu = xres[0] / xres[1] if isinstance(xres, tuple) else float(xres)
        if not ppu or unit not in (2, 3):
            return None
        return (2.54e4 if unit == 2 else 1e4) / ppu  # inch / cm -> µm per pixel

    # ---- info ------------------------------------------------------------
    def missing_tiles(self) -> int:
        return int(sum(np.count_nonzero(i.counts == 0) for i in self._info))

    def stored_bytes(self, level: int) -> tuple[int, int]:
        """(encoded bytes, missing tiles) of the tiles a level export covers."""
        info, lv = self._info[level], self.levels[level]
        grid = info.counts.reshape(-1, info.stride)[:lv.grid_rows, :lv.grid_cols]
        beyond = lv.grid_rows * lv.grid_cols - grid.size  # declared past the stored grid
        return int(grid.sum()), int(np.count_nonzero(grid == 0)) + beyond

    def stored_size(self, level: int) -> tuple[int, int]:
        """(width, height) of the level as stored, including any padding."""
        return self._info[level].stored_size

    def _size_note(self, lv: Level) -> str:
        sw, sh = self.stored_size(lv.index)
        if (sw, sh) == (lv.width, lv.height):
            return ""
        how = ("padding cropped" if sw >= lv.width and sh >= lv.height else
               "extended with white" if sw <= lv.width and sh <= lv.height else "adjusted")
        return f" (stored {sw}x{sh}, {how})"

    def summary(self, min_level_size: int = 256, max_levels: int | None = None) -> str:
        info0 = self._info[0]
        sub = self.jpeg_subsampling(0) if info0.compression == _TIFF_JPEG else None
        scanner = ((self.description.split("|", 1)[0].strip().splitlines() or ["-"])[0]
                   if self.is_aperio else "-")
        lines = [
            f"File          : {self.path}",
            f"Type          : brightfield RGB ({'Aperio SVS' if self.is_aperio else 'TIFF'}, "
            f"{'BigTIFF' if self.is_bigtiff else 'classic TIFF'})",
            f"Scanner       : {scanner}",
            f"Dimensions    : {self.width} x {self.height} px  @ "
            f"{format_magnification(self.magnification)}x",
            f"Magnification : {self.magnification_text()}",
            f"MPP           : {self.mpp:.5g} um/px"
            + ("" if self.mpp_source == "description" else f" ({self.mpp_source})"),
            f"Tiles         : {self.tile_size}x{self.tile_size} {self.codec}, photometric "
            f"{'RGB' if info0.photometric == _PHOTOMETRIC_RGB else 'YCbCr'}"
            + (f", subsampling {sub[0]}x{sub[1]}" if sub else "")
            + f", {self.missing_tiles()} missing",
            f"Associated    : {', '.join(self.associated) or '-'} (not exported)",
        ]
        exported = {lv.index for lv in select_export_levels(self.levels, min_level_size,
                                                             max_levels)}
        lines.append(f"Levels ('*' = exported, min size {min_level_size}px):")
        for lv in self.levels:
            lines.append(
                f" {'*' if lv.index in exported else ' '}[{lv.index:2d}] "
                f"{lv.width:>7d} x {lv.height:<7d} ds {lv.downsample:<8g}"
                + (" (synthetic)" if lv.synthetic else f" grid {lv.grid_cols}x{lv.grid_rows}"
                   + self._size_note(lv)))
        return "\n".join(lines)

    # ---- pixel access ----------------------------------------------------
    def _decode(self, info: _SVSLevelInfo, data: bytes, index: int) -> np.ndarray:
        if info.compression == _TIFF_JPEG:
            return decode_jpeg_rgb(data, info.jpegtables, info.rgb_stored)
        seg = info.page.decode(data, index, jpegtables=info.jpegtables)[0]
        return np.asarray(seg).reshape(self.tile_size, self.tile_size, -1)

    def raw_tile(self, level: int, col: int, row: int) -> bytes | None:
        """The stored (still encoded) bytes of one source tile, or None if missing."""
        info, lv = self._info[level], self.levels[level]
        stored_rows = len(info.offsets) // info.stride
        if not (0 <= col < min(lv.grid_cols, info.stride)
                and 0 <= row < min(lv.grid_rows, stored_rows)):
            return None  # outside the image, or white padding added past the stored grid
        i = row * info.stride + col
        off, cnt = int(info.offsets[i]), int(info.counts[i])
        return self._mm[off:off + cnt] if cnt > 0 else None

    def get_tile(self, level: int, col: int, row: int, z: int = 0) -> np.ndarray:
        """Return one decoded (tile_size, tile_size, 3) uint8 RGB source tile.

        Missing / out-of-range tiles are white (255).
        """
        if not 0 <= level < self.n_source_levels:
            raise IndexError(f"level {level} is not a source level "
                             f"(0..{self.n_source_levels - 1}); use read_region()")
        key = (level, z, col, row)
        hit = self._cache_get(key)
        if hit is not None:
            return hit
        ts = self.tile_size
        tile = np.full((ts, ts, 3), 255, np.uint8)
        data = self.raw_tile(level, col, row)
        if data is not None:
            try:
                arr = self._decode(self._info[level], data,
                                   row * self._info[level].stride + col)
                h, w = min(arr.shape[0], ts), min(arr.shape[1], ts)
                tile[:h, :w] = arr[:h, :w, :3]
                # Pixels past the stored image are encoder padding, not
                # image: make them white (matters where a short level was
                # extended to the exact pyramid size).
                sw, sh = self._info[level].stored_size
                vx, vy = sw - col * ts, sh - row * ts
                if vx < ts:
                    tile[:, max(vx, 0):] = 255
                if vy < ts:
                    tile[max(vy, 0):, :] = 255
            except Exception as exc:
                self._note_decode_failure(f"L{level} c{col} r{row}", exc)
        return self._cache_put(key, tile)

    def read_region(self, level: int, x: int, y: int, width: int, height: int,
                    z: int = 0) -> np.ndarray:
        """Read a (height, width, 3) uint8 RGB region at ``level``.

        Pixels outside the level bounds are white.
        """
        return self._read_plane(level, x, y, width, height,
                                lambda lv, col, row: self.get_tile(lv, col, row, z), (3,))

    # ---- tile copy (lossless passthrough) ---------------------------------
    def jpeg_subsampling(self, level: int) -> tuple[int, int] | None:
        """(h, v) luma sampling factors read from the level's JPEG SOF header."""
        lv, found = self.levels[level], None
        n = lv.grid_cols * lv.grid_rows
        for i in sorted({0, n // 2, n - 1}):
            data = self.raw_tile(level, i % lv.grid_cols, i // lv.grid_cols)
            sof = jpeg_sof_info(bytes(data[:4096])) if data is not None else None
            if sof is None or len(sof.get("components", ())) != 3:
                continue
            (_, h, v), *chroma = sof["components"]
            if any((ch, cv) != (1, 1) for _, ch, cv in chroma):
                return None
            if found not in (None, (h, v)):
                return None
            found = (h, v)
        return found

    def copy_tiles_plan(self, levels: Sequence[Level]) -> tuple[dict | None, str]:
        """Can ``levels`` be written by copying the stored JPEG tiles as-is?

        Returns ``({level_index: tifffile kwargs}, "")`` or ``(None, reason)``.
        """
        plan: dict[int, dict] = {}
        for lv in levels:
            if lv.synthetic:
                return None, f"level {lv.index} is computed, not stored"
            info = self._info[lv.index]
            if info.compression != _TIFF_JPEG:
                return None, f"source tiles are {self.codec}, not JPEG"
            if info.photometric != _PHOTOMETRIC_YCBCR:
                return None, "source JPEG tiles are not YCbCr"
            sub = self.jpeg_subsampling(lv.index)
            if sub is None:
                return None, f"could not determine JPEG subsampling of level {lv.index}"
            plan[lv.index] = {"subsampling": sub, "jpegtables": info.jpegtables}
        return plan, ""

    def blank_jpeg(self, subsampling: tuple[int, int]) -> bytes:
        """A white JPEG tile matching the source's tile size and subsampling."""
        key = ("blank", subsampling)
        with self._lock:
            hit = self._cache.get(key)
        if hit is None:
            hit = imagecodecs.jpeg_encode(
                np.full((self.tile_size, self.tile_size, 3), 255, np.uint8),
                level=90, subsampling=subsampling)
            with self._lock:
                self._cache[key] = hit  # tiny; LRU eviction is harmless
        return hit


class _PageSource:
    """One stored image plane of one pyramid level: a TIFF page's tiles or
    strips ("segments"), where they are in the file and how to decode them."""

    def __init__(self, page, byteorder: str):
        self.page = page
        self.width, self.height = int(page.imagewidth), int(page.imagelength)
        self.spp = int(page.samplesperpixel)
        self.tiled = bool(page.is_tiled)
        if self.tiled:
            self.seg_w, self.seg_h = int(page.tilewidth), int(page.tilelength)
        else:
            self.seg_w = self.width
            self.seg_h = max(1, min(int(page.rowsperstrip or self.height), self.height))
        self.ncols = math.ceil(self.width / self.seg_w)
        self.nrows = math.ceil(self.height / self.seg_h)
        self.offsets = np.asarray(page.dataoffsets, np.int64)
        self.counts = np.asarray(page.databytecounts, np.int64)
        self.compression = int(page.compression)
        self.photometric = int(page.photometric)
        self.predictor = int(getattr(page, "predictor", 1) or 1)
        self.jpegtables = page.jpegtables or None
        self.dtype = np.dtype(page.dtype)
        self.file_dtype = self.dtype.newbyteorder(byteorder)
        self.decode = page.decode  # may seek the file: fetch once, up front

    @property
    def n_segments(self) -> int:
        return self.ncols * self.nrows


class _TiffPagesSlide(_SlideBase):
    """Shared reader for slides stored as plain TIFF pages (QPTIFF, TIFF):
    each pyramid level is one page per plane (channel), read segment by
    segment (tiles or strips) through a small byte-bounded cache.

    Stored levels are used while each is an exact 2x (or 4x, ...) step of the
    one above; from the first that is not, lower levels are computed by 2x2
    averaging instead, so readers always see one consistent pyramid.
    """

    SEGMENT_CACHE_BYTES = 96 * 2 ** 20
    COMPUTED_TILE_CACHE_BYTES = 64 * 2 ** 20
    VIRTUAL_TILE = 512   # tile grid used for stripped (untiled) images

    def __init__(self, path: str | os.PathLike):
        super().__init__(path, tile_cache=0)  # sized in _set_pyramid
        self._tif = None
        self._segs: OrderedDict = OrderedDict()
        self._seg_bytes = 0
        self._src: list[list[_PageSource]] = []
        self._blank: dict = {}

    def close(self) -> None:
        tif = getattr(self, "_tif", None)
        if tif is not None:
            tif.close()
            self._tif = None
        super().close()

    def _damaged(self, why: str) -> SlideFormatError:
        log.debug("%s: %s", self.path.name, why)
        return SlideFormatError(f"{self.path.name}: this {self.input_type} is damaged or "
                                f"incomplete and cannot be converted")

    def _open_tiff(self) -> None:
        if bytes(self._mm[:4]) not in TIFF_MAGICS:
            raise SlideFormatError(f"{self.path.name}: not a TIFF file (bad magic)")
        try:
            self._tif = tifffile.TiffFile(self.path)
            n_pages = len(self._tif.pages)
        except Exception as exc:
            raise self._damaged(f"tifffile cannot open it: {exc}") from exc
        if not n_pages:
            raise self._damaged("no image pages")
        if not self._tif.is_bigtiff and self.file_size > CLASSIC_TIFF_MAX_BYTES:
            # classic TIFF offsets are 32-bit: data past 4 GB is unreachable
            raise self._damaged(f"classic TIFF of {_gb(self.file_size)} (over the 4 GB limit)")

    # ---- pyramid ---------------------------------------------------------
    def _source(self, page) -> _PageSource:
        src = _PageSource(page, self._tif.byteorder)
        if src.spp > 1 and int(page.planarconfig) != 1:
            raise SlideFormatError(f"{self.path.name}: separate colour planes "
                                   f"(planar configuration 2) are not supported")
        if len(src.offsets) < src.n_segments or len(src.counts) != len(src.offsets):
            raise self._damaged(f"page {page.index}: {len(src.offsets)} segments, "
                                f"expected {src.n_segments}")
        ends = src.offsets + src.counts
        if ((src.counts > 0) & ((ends > self.file_size) | (src.offsets < 0))).any():
            raise self._damaged(f"page {page.index}: data past the end of the file")
        nonzero = np.flatnonzero(src.counts[:src.n_segments])
        if src.compression != 1 and len(nonzero):  # try one segment now
            i = int(nonzero[len(nonzero) // 2])
            try:
                self._decode_segment(src, i)
            except Exception as exc:
                raise self._damaged(f"page {page.index} segment {i}: {exc}") from exc
        return src

    def _set_pyramid(self, stacks: Sequence[Sequence], tile_size: int | None = None) -> None:
        """``stacks`` = [[page per plane] per stored level], full resolution first."""
        first = stacks[0][0]
        if tile_size is None:
            tile_size = (int(first.tilewidth) if first.is_tiled
                         and first.tilewidth == first.tilelength else self.VIRTUAL_TILE)
        self.tile_size = tile_size
        levels: list[Level] = []
        for pages in stacks:
            srcs = [self._source(p) for p in pages]
            w, h = srcs[0].width, srcs[0].height
            if any((s.width, s.height) != (w, h) for s in srcs):
                break
            ds = 1.0
            if levels:
                prev = levels[-1]
                r = round(prev.width / w) if w else 0
                ok = r >= 2 and r & (r - 1) == 0 \
                    and w in (prev.width // r, -(-prev.width // r)) \
                    and h in (prev.height // r, -(-prev.height // r))
                if not ok:
                    log.info("%s: stored level %dx%d is not an exact 2x step of %dx%d; "
                             "computing the lower levels from %dx%d instead", self.path.name,
                             w, h, prev.width, prev.height, prev.width, prev.height)
                    break
                ds = prev.downsample * r
            levels.append(Level(len(levels), (self.magnification or 1) / ds, w, h,
                                math.ceil(w / tile_size), math.ceil(h / tile_size), ds,
                                coord_mode="tiff"))
            self._src.append(srcs)
        self.n_source_levels = len(levels)
        self.levels = self._add_synthetic_levels(levels)
        tile_bytes = tile_size * tile_size * self.dtype.itemsize * (3 if self.is_rgb else 1)
        self._cache_size = max(16, self.COMPUTED_TILE_CACHE_BYTES // tile_bytes)

    # ---- segment access ----------------------------------------------------
    def _decode_segment(self, src: _PageSource, index: int) -> np.ndarray:
        """(rows, cols, samples) array of one stored tile / strip."""
        off, cnt = int(src.offsets[index]), int(src.counts[index])
        data = self._mm[off:off + cnt]
        if src.compression == _TIFF_JPEG and src.dtype == np.uint8 and src.spp == 1 \
                and src.jpegtables is None:
            a = decode_jpeg(bytes(data))
            return a.reshape(a.shape[0], a.shape[1], 1)
        seg = src.decode(bytes(data), index, jpegtables=src.jpegtables)[0]
        a = np.asarray(seg)
        return a.reshape(a.shape[-3], a.shape[-2], a.shape[-1])

    def _segment(self, src: _PageSource, index: int) -> np.ndarray | None:
        if index >= len(src.counts) or src.counts[index] <= 0:
            return None
        if src.compression == 1 and not src.tiled:  # uncompressed strip: view the file
            row0 = (index // src.ncols) * src.seg_h
            rows = min(src.seg_h, src.height - row0)
            n = rows * src.width * src.spp
            if n * src.dtype.itemsize > src.counts[index]:
                return None
            a = np.frombuffer(self._mm, src.file_dtype, n, int(src.offsets[index]))
            return a.reshape(rows, src.width, src.spp)
        key = (id(src), index)
        with self._lock:
            hit = self._segs.get(key)
            if hit is not None:
                self._segs.move_to_end(key)
                return hit
        try:
            a = self._decode_segment(src, index)
        except Exception as exc:
            self._note_decode_failure(f"page {src.page.index} segment {index}", exc)
            return None
        a.flags.writeable = False
        with self._lock:
            self._segs[key] = a
            self._seg_bytes += a.nbytes
            while self._seg_bytes > self.SEGMENT_CACHE_BYTES and len(self._segs) > 1:
                _, old = self._segs.popitem(last=False)
                self._seg_bytes -= old.nbytes
        return a

    @property
    def _sample_shape(self) -> tuple:
        return (3,) if self.is_rgb else ()

    def get_tile(self, level: int, col: int, row: int, channel: int = 0) -> np.ndarray:
        """One (tile_size, tile_size[, 3]) tile of a stored level; missing
        tiles and pixels outside the image are ``fill_value``."""
        if not 0 <= level < self.n_source_levels:
            raise IndexError(f"level {level} is not a source level "
                             f"(0..{self.n_source_levels - 1}); use read_region()")
        src = self._src[level][channel]
        ts = self.tile_size
        tile = np.full((ts, ts) + self._sample_shape, self.fill_value, self.dtype)
        x0, y0 = col * ts, row * ts
        x1, y1 = min(x0 + ts, src.width), min(y0 + ts, src.height)
        if x1 <= x0 or y1 <= y0:
            return tile
        for sr in range(y0 // src.seg_h, (y1 - 1) // src.seg_h + 1):
            for sc in range(x0 // src.seg_w, (x1 - 1) // src.seg_w + 1):
                seg = self._segment(src, sr * src.ncols + sc)
                if seg is None:
                    continue
                gx, gy = sc * src.seg_w, sr * src.seg_h          # segment origin
                ax0, ay0 = max(x0, gx), max(y0, gy)
                ax1 = min(x1, gx + src.seg_w, gx + seg.shape[1])
                ay1 = min(y1, gy + src.seg_h, gy + seg.shape[0])
                if ax1 <= ax0 or ay1 <= ay0:
                    continue
                part = seg[ay0 - gy:ay1 - gy, ax0 - gx:ax1 - gx]
                part = part[..., :3] if self.is_rgb else part[..., 0]
                tile[ay0 - y0:ay1 - y0, ax0 - x0:ax1 - x0] = part
        return tile

    def read_region(self, level: int, x: int, y: int, width: int, height: int,
                    channel: int = 0, z: int = 0) -> np.ndarray:
        """(height, width) region of one channel, or (height, width, 3) for
        RGB slides (``channel`` is then ignored)."""
        plane = 0 if self.is_rgb else channel
        if not 0 <= plane < len(self._src[0]):
            raise IndexError(f"channel {channel} out of range (0..{len(self._src[0]) - 1})")
        return self._assemble(level, x, y, width, height, plane)

    def _assemble(self, level: int, x: int, y: int, width: int, height: int,
                  plane: int) -> np.ndarray:
        out = np.full((height, width) + self._sample_shape, self.fill_value, self.dtype)
        lv, ts = self.levels[level], self.tile_size
        x0, y0 = max(x, 0), max(y, 0)
        x1, y1 = min(x + width, lv.width), min(y + height, lv.height)
        for row in range(y0 // ts, (y1 - 1) // ts + 1) if y1 > y0 and x1 > x0 else ():
            ty = row * ts
            sy0, sy1 = max(y0, ty), min(y1, ty + ts)
            for col in range(x0 // ts, (x1 - 1) // ts + 1):
                tx = col * ts
                sx0, sx1 = max(x0, tx), min(x1, tx + ts)
                tile = self._level_tile(level, col, row, plane)
                out[sy0 - y:sy1 - y, sx0 - x:sx1 - x] = tile[sy0 - ty:sy1 - ty, sx0 - tx:sx1 - tx]
        return out

    def _level_tile(self, level: int, col: int, row: int, plane: int) -> np.ndarray:
        """One tile of any level. Computed levels are built tile by tile from
        the 2x2 tiles above them (2x2 average, edge pixels replicated, same
        result as the other readers) and cached, so writing a computed level
        reads each tile of the level above only once."""
        if not self.levels[level].synthetic:
            return self.get_tile(level, col, row, plane)
        key = (level, plane, col, row)
        hit = self._cache_get(key)
        if hit is not None:
            return hit
        ts, ss = self.tile_size, self._sample_shape
        parent = self.levels[level - 1]
        px, py = 2 * col * ts, 2 * row * ts
        vw, vh = min(2 * ts, parent.width - px), min(2 * ts, parent.height - py)
        tile = np.full((ts, ts) + ss, self.fill_value, self.dtype)
        if vw > 0 and vh > 0:
            src = self._assemble(level - 1, px, py, vw, vh, plane)
            pw, ph = vw + vw % 2, vh + vh % 2
            if (pw, ph) != (vw, vh):
                src = np.pad(src, ((0, ph - vh), (0, pw - vw)) + ((0, 0),) * len(ss),
                             mode="edge")
            small = src.reshape((ph // 2, 2, pw // 2, 2) + ss).mean(axis=(1, 3),
                                                                   dtype=np.float32)
            tile[:ph // 2, :pw // 2] = (small + 0.5).astype(self.dtype)
        return self._cache_put(key, tile)

    # ---- tile copy (passthrough of the stored, compressed tiles) ------------
    # Lossless codecs whose tiles are copied as-is (with the source predictor).
    # Other codecs (uncompressed, JPEG 2000, WebP, ZSTD, ...) are re-encoded so
    # that every TIFF reader can open the output.
    _COPY_LOSSLESS = {5: "lzw", 8: "zlib", 32946: "zlib"}

    def _jpeg_sampling(self, src: _PageSource) -> tuple[int, int] | None:
        """(h, v) luma sampling of a 3-component JPEG level, read from the SOF
        of a few tiles; None if unreadable, inconsistent or chroma-upsampled."""
        found = None
        n = src.n_segments
        for i in sorted({0, n // 2, n - 1}):
            if src.counts[i] <= 0:
                continue
            off = int(src.offsets[i])
            sof = jpeg_sof_info(bytes(self._mm[off:off + min(int(src.counts[i]), 4096)]))
            if sof is None or len(sof.get("components", ())) != 3:
                return None
            (_, h, v), *chroma = sof["components"]
            if any((ch, cv) != (1, 1) for _, ch, cv in chroma) or found not in (None, (h, v)):
                return None
            found = (h, v)
        return found if found in ((1, 1), (2, 1), (2, 2), (4, 1)) else None

    def _level_copy_kwargs(self, srcs: Sequence[_PageSource]) -> tuple[dict | None, str]:
        """tifffile kwargs that reproduce one stored level's encoding, so its
        tiles can be written unchanged; or (None, why not)."""
        s0, ts = srcs[0], self.tile_size
        if any(not s.tiled or (s.seg_w, s.seg_h) != (ts, ts) for s in srcs):
            return None, "it is not tiled like the output"
        if len({(s.compression, s.photometric, s.predictor, s.jpegtables) for s in srcs}) != 1:
            return None, "its channels are stored differently"
        if s0.compression == _TIFF_JPEG:
            if s0.dtype != np.uint8:
                return None, "its JPEG tiles are not 8-bit"
            kw: dict = {"compression": "jpeg"}
            if self.is_rgb:
                if s0.spp != 3 or s0.photometric not in (_PHOTOMETRIC_RGB, _PHOTOMETRIC_YCBCR):
                    return None, "its JPEG tiles are not plain RGB/YCbCr"
                sub = self._jpeg_sampling(s0)
                if sub is None:
                    return None, "its JPEG subsampling could not be read"
                kw["subsampling"] = sub
                if s0.photometric == _PHOTOMETRIC_RGB:  # keep "RGB" (no YCbCr transform)
                    kw["compressionargs"] = {"outcolorspace": "RGB"}
            elif s0.spp != 1:
                return None, "its JPEG tiles are not single-channel"
            if s0.jpegtables:
                kw["jpegtables"] = s0.jpegtables
            return kw, ""
        codec = self._COPY_LOSSLESS.get(s0.compression)
        if codec is None:
            name = "uncompressed" if s0.compression == 1 else \
                s0.page.compression.name.replace("_", " ").lower()
            return None, f"its tiles are {name}"
        if any(np.count_nonzero(s.counts[:s.n_segments] <= 0) for s in srcs):
            return None, "it has missing tiles"
        if s0.predictor not in (1, 2):
            return None, "it uses an unusual predictor"
        if s0.dtype.itemsize > 1 and self._tif.byteorder != "<":
            return None, "it is big-endian"
        if self.is_rgb and (s0.spp != 3 or s0.photometric != _PHOTOMETRIC_RGB):
            return None, "it is not plain RGB"
        if not self.is_rgb and (s0.spp != 1 or s0.photometric != 1):
            return None, "it is not plain grayscale"
        return {"compression": codec, "predictor": s0.predictor == 2}, ""

    def copy_tiles_plan(self, levels: Sequence[Level]) -> tuple[dict | None, str]:
        """Which stored ``levels`` can be written by copying their compressed
        tiles unchanged? Returns ``({level_index: tifffile kwargs}, notes)``:
        levels left out (computed, or stored in a way that can't be copied)
        are decoded and encoded again, and ``notes`` says why for stored
        ones. ``(None, reason)`` if no level can be copied."""
        plan: dict[int, dict] = {}
        notes: list[str] = []
        for lv in levels:
            if lv.synthetic:
                continue
            kw, why = self._level_copy_kwargs(self._src[lv.index])
            if kw is None:
                notes.append(f"level {lv.index} is re-encoded because {why}")
            else:
                plan[lv.index] = kw
        if not plan:
            return None, notes[0].split(" because ", 1)[1] if notes else "no stored level"
        return plan, "; ".join(notes)

    def rest_compression(self, plan: dict) -> str:
        """Compression for the levels a tile copy re-encodes: JPEG next to
        copied JPEG tiles (8-bit only), otherwise deflate."""
        if self.dtype == np.uint8 and any(kw["compression"] == "jpeg" for kw in plan.values()):
            return "jpeg"
        return DEFAULT_COMPRESSION

    def raw_tile(self, level: int, col: int, row: int, channel: int = 0) -> bytes | None:
        """The stored (still encoded) bytes of one tile, or None if missing."""
        src = self._src[level][channel]
        if not (0 <= col < src.ncols and 0 <= row < src.nrows):
            return None
        i = row * src.ncols + col
        off, cnt = int(src.offsets[i]), int(src.counts[i])
        return self._mm[off:off + cnt] if cnt > 0 else None

    def copied_tile(self, level: int, col: int, row: int, channel: int, kw: dict) -> bytes:
        """A stored tile's bytes for a tile copy; a missing tile becomes a
        blank (fill-colour) JPEG encoded like the level (``kw`` from the
        plan). Lossless levels are only copied when no tile is missing."""
        data = self.raw_tile(level, col, row, channel)
        return bytes(data) if data is not None else self.blank_tile(kw)

    def blank_tile(self, kw: dict) -> bytes:
        sub = kw.get("subsampling")
        rgb_coded = "compressionargs" in kw
        key = ("blank", sub, rgb_coded)
        with self._lock:
            hit = self._blank.get(key)
        if hit is None:
            shape = (self.tile_size, self.tile_size) + ((3,) if self.is_rgb else ())
            fill = np.full(shape, self.fill_value, np.uint8)
            args = dict(level=90)
            if self.is_rgb:
                args["subsampling"] = sub
                if rgb_coded:
                    args.update(colorspace="RGB", outcolorspace="RGB")
            hit = imagecodecs.jpeg_encode(fill, **args)
            with self._lock:
                self._blank[key] = hit
        return hit

    def blank_jpeg(self) -> bytes:
        """A blank (fill-colour) JPEG tile for missing grayscale tiles."""
        return self.blank_tile({})

    def copy_bytes(self, levels: Sequence[Level], jpeg_quality: int = 90,
                   compression: str | None = "jpeg") -> tuple[int, bool]:
        """(bytes, exact) of a tile copy of ``levels``: stored tile bytes for
        copied levels, an estimate (``compression``) for re-encoded ones."""
        total, exact, ts = 0, True, self.tile_size
        planes = len(self._src[0])
        samples = 3 if self.is_rgb else 1
        plan = self.copy_tiles_plan(levels)[0] or {}
        for lv in levels:
            if lv.index not in plan:
                cells = math.ceil(lv.width / ts) * math.ceil(lv.height / ts)
                total += int(cells * ts * ts * planes * samples * self.dtype.itemsize
                             * compression_factor(compression, jpeg_quality))
                exact = False
                continue
            for s in self._src[lv.index]:
                counts = s.counts[:s.n_segments]
                total += int(counts.sum()) + int(np.count_nonzero(counts == 0)) * 1024
        return total, exact

    def _pixel_text(self) -> str:
        return f"{self.dtype.itemsize * 8}-bit"

    def _levels_lines(self, min_level_size: int, max_levels: int | None) -> list[str]:
        exported = {lv.index for lv in select_export_levels(self.levels, min_level_size,
                                                             max_levels)}
        lines = [f"Levels ('*' = exported, min size {min_level_size}px):"]
        for lv in self.levels:
            lines.append(f" {'*' if lv.index in exported else ' '}[{lv.index:2d}] "
                         f"{lv.width:>7d} x {lv.height:<7d} ds {lv.downsample:<8g}"
                         + (" (computed)" if lv.synthetic else ""))
        return lines


_QPI_TEXT_TAGS = ("Biomarker", "Name", "Color", "ImageType", "Objective",
                  "ObjectiveMagnification", "Magnification", "PixelSizeMicrons",
                  "AcquisitionSoftware", "SlideID")


def parse_qptiff_description(desc: str | bytes | None) -> dict[str, str]:
    """Top-level fields of a PerkinElmer/Akoya QPTIFF page description
    (``<PerkinElmer-QPI-ImageDescription>``): ImageType, Name, Biomarker,
    Color, Objective, ... Nested elements (e.g. the filter ``Name`` inside
    ``Responsivity``) are ignored. Returns {} if it isn't QPTIFF XML."""
    if not desc:
        return {}
    if isinstance(desc, bytes):
        desc = desc.decode("utf-8", "replace")
    text = re.sub(r"^\s*<\?xml[^>]*\?>", "", desc.strip("\0 \r\n\t"))
    try:
        root = ElementTree.fromstring(text)
        out = {}
        for child in root:
            if len(child) == 0 and child.text is not None:
                out[child.tag] = child.text.strip()
        return out
    except ElementTree.ParseError:  # malformed XML: first occurrence of each tag
        out = {}
        for tag in _QPI_TEXT_TAGS:
            m = re.search(rf"<{tag}>([^<]*)</{tag}>", text)
            if m:
                out[tag] = m.group(1).strip()
        return out


def _parse_rgb_text(text: str | None) -> tuple[int, int, int] | None:
    """'0,128,255' -> (0, 128, 255); None if not three 0-255 numbers."""
    nums = re.findall(r"\d+", text or "")
    if len(nums) != 3 or any(int(n) > 255 for n in nums):
        return None
    return tuple(int(n) for n in nums)


def _mpp_from_resolution(page) -> float | None:
    return SVSSlide._mpp_from_resolution(page)


class QPTIFFSlide(_TiffPagesSlide):
    """Reader for PerkinElmer / Akoya ``.qptiff`` slides (Vectra Polaris,
    PhenoImager; also written by KFBio scanners), using tifffile only.

    Pages are classified by the ``<ImageType>`` in each page's XML
    description: ``FullResolution`` (level 0, one page per channel),
    ``ReducedResolution`` (lower levels), and ``Thumbnail`` / ``Overview`` /
    ``Label`` / ``Macro``, which are not part of the image and are skipped.
    Channel names come from ``<Biomarker>`` (or ``<Name>``), colours from
    ``<Color>``; the objective from ``<Objective>``; the pixel size from the
    XML or the TIFF resolution tags. 8- and 16-bit data keep their bit depth.

    Fluorescence QPTIFFs read like :class:`KFBioSlide` (``read_region(level,
    x, y, w, h, channel)``); brightfield (RGB) ones like :class:`SVSSlide`.
    """

    input_type = "QPTIFF"
    _LEVEL_TYPES = ("FullResolution", "ReducedResolution")

    def __init__(self, path: str | os.PathLike, mpp_fallback: float = DEFAULT_MPP,
                 channel_names: Sequence[str] | None = None, tile_cache: int = 0,
                 magnification: float | None = None):
        super().__init__(path)
        try:
            self._open_tiff()
            self._parse(mpp_fallback, channel_names, magnification)
        except Exception:
            self.close()
            raise

    def _parse(self, mpp_fallback: float, channel_names, magnification) -> None:
        tif = self._tif
        self.is_bigtiff = tif.is_bigtiff
        pages, infos = [], []
        try:
            for pg in tif.pages:
                pages.append(pg)
                infos.append(parse_qptiff_description(pg.description))
        except Exception as exc:
            raise self._damaged(f"cannot read page list: {exc}") from exc
        kinds = [i.get("ImageType", "") for i in infos]
        full = [p for p, k in zip(pages, kinds) if k == "FullResolution"]
        if not full:
            raise SlideFormatError(f"{self.path.name}: no full-resolution image found "
                                   f"(not a PerkinElmer/Akoya QPTIFF?)")
        self.associated = [k for k in kinds if k and k not in self._LEVEL_TYPES]
        base, info0 = full[0], infos[pages.index(full[0])]
        self.description = base.description or ""
        self.scanner = info0.get("AcquisitionSoftware", "")
        self.width, self.height = int(base.imagewidth), int(base.imagelength)
        self.dtype = np.dtype(base.dtype)
        self.is_rgb = int(base.samplesperpixel) in (3, 4) and \
            int(base.photometric) in (_PHOTOMETRIC_RGB, _PHOTOMETRIC_YCBCR)
        if self.is_rgb:
            if self.dtype != np.uint8:
                raise SlideFormatError(f"{self.path.name}: {self.dtype} RGB QPTIFF is not "
                                       f"supported; only 8-bit RGB")
            self.kind, self.fill_value = "brightfield", 255
            full = full[:1]
            self.channel_names = ["Red", "Green", "Blue"]
            self.channel_colors = [(255, 0, 0), (0, 255, 0), (0, 0, 255)]
        else:
            if self.dtype not in (np.uint8, np.uint16):
                raise SlideFormatError(
                    f"{self.path.name}: {self.dtype} pixels are not supported (8- or "
                    f"16-bit only; unmixed component data is 32-bit float)")
            if any(int(p.samplesperpixel) != 1 for p in full):
                raise SlideFormatError(f"{self.path.name}: unexpected samples per pixel "
                                       f"in a fluorescence QPTIFF")
            self.kind, self.fill_value = "fluorescence", 0
            names, colors = [], []
            for i, p in enumerate(full):
                info = infos[pages.index(p)]
                name = info.get("Biomarker") or info.get("Name") or f"Channel {i + 1}"
                names.append(name)
                colors.append(_parse_rgb_text(info.get("Color")) or channel_color(name, i))
            if channel_names:
                if len(channel_names) != len(names):
                    raise ValueError(f"--channel-names gives {len(channel_names)} names "
                                     f"for {len(names)} channels")
                names = list(channel_names)
            self.channel_names, self.channel_colors = names, colors
        self.n_channels = len(full) if not self.is_rgb else 3

        objective = info0.get("ObjectiveMagnification") or info0.get("Magnification") \
            or info0.get("Objective")
        m = re.search(r"\d+(?:\.\d+)?", objective or "")
        self._resolve_magnification(m.group(0) if m else None, "QPTIFF XML", magnification)

        mpp, self.mpp_source = _float_or_none(info0.get("PixelSizeMicrons")), "QPTIFF XML"
        if mpp is None or not 0.01 <= mpp <= 100:
            mpp, self.mpp_source = _mpp_from_resolution(base), "resolution tags"
        if mpp is None or not 0.01 <= mpp <= 100:
            mpp, self.mpp_source = float(mpp_fallback), "fallback"
            self._warn(f"no pixel size in file; using {mpp_fallback} um/px (set --mpp)")
        self.mpp = float(f"{mpp:.6g}")
        self.mpp_from_header = self.mpp_source != "fallback"

        # Pyramid: level 0, then ReducedResolution pages in groups of one
        # page per plane (same size, in file order).
        n = len(full)
        reduced = [p for p, k in zip(pages, kinds) if k == "ReducedResolution"
                   and int(p.samplesperpixel) == int(base.samplesperpixel)]
        stacks = [full]
        for i in range(0, len(reduced) - n + 1, n):
            group = reduced[i:i + n]
            if len({(int(p.imagewidth), int(p.imagelength)) for p in group}) != 1:
                break
            stacks.append(group)
        self.z_values = [0]
        self._set_pyramid(stacks)

    def summary(self, min_level_size: int = 256, max_levels: int | None = None) -> str:
        src = self._src[0][0]
        lines = [
            f"File          : {self.path}",
            f"Type          : QPTIFF {self.kind} "
            f"({'RGB' if self.is_rgb else f'{self.n_channels} channels'}, {self._pixel_text()}, "
            f"{'BigTIFF' if self.is_bigtiff else 'classic TIFF'})",
            f"Scanner       : {self.scanner or '-'}",
            f"Dimensions    : {self.width} x {self.height} px",
            f"Magnification : {self.magnification_text()}",
            f"MPP           : {self.mpp:.5g} um/px"
            + ("" if self.mpp_source == "QPTIFF XML" else f" ({self.mpp_source})"),
            f"Tiles         : {src.seg_w}x{src.seg_h} {src.page.compression.name}"
            + ("" if src.tiled else " (strips)"),
        ]
        if not self.is_rgb:
            lines.append(f"Channels ({self.n_channels}) : " + ", ".join(
                f"{n} (RGB {r},{g},{b})" for n, (r, g, b) in
                zip(self.channel_names, self.channel_colors)))
        lines.append(f"Associated    : {', '.join(self.associated) or '-'} (not exported)")
        lines += self._levels_lines(min_level_size, max_levels)
        return "\n".join(lines)


class TiffSlide(_TiffPagesSlide):
    """Reader for plain ``.tif`` / ``.tiff`` images: RGB (8-bit) or
    single-channel grayscale (8- or 16-bit), tiled or stripped, flat or
    already pyramidal. Pixel size comes from the TIFF resolution tags,
    magnification from an ``AppMag = `` in the description, if present.

    A TIFF that cannot be read -- truncated, corrupt, or a scanner TIFF over
    4 GB written as classic TIFF -- raises SlideFormatError("... this TIFF is
    damaged or incomplete and cannot be converted").
    """

    input_type = "TIFF"

    def __init__(self, path: str | os.PathLike, mpp_fallback: float = DEFAULT_MPP,
                 tile_cache: int = 0, magnification: float | None = None):
        super().__init__(path)
        try:
            self._open_tiff()
            self._parse(mpp_fallback, magnification)
        except Exception:
            self.close()
            raise

    def _parse(self, mpp_fallback: float, magnification: float | None) -> None:
        tif, name = self._tif, self.path.name
        self.is_bigtiff = tif.is_bigtiff
        try:
            series = tif.series[0]
            stacks = [list(lv.pages) for lv in series.levels]
            base = stacks[0][0]
            self.description = base.description or ""
        except Exception as exc:
            raise self._damaged(f"cannot read the image list: {exc}") from exc
        if base is None or len(stacks[0]) != 1:
            raise SlideFormatError(
                f"{name}: this TIFF holds {len(stacks[0])} images per level (several "
                f"channels or a stack); only RGB or single-channel grayscale TIFFs "
                f"can be converted")
        spp, phot = int(base.samplesperpixel), int(base.photometric)
        self.dtype = np.dtype(base.dtype)
        self.width, self.height = int(base.imagewidth), int(base.imagelength)
        if spp in (3, 4) and phot in (_PHOTOMETRIC_RGB, _PHOTOMETRIC_YCBCR):
            if self.dtype != np.uint8:
                raise SlideFormatError(f"{name}: {self.dtype} RGB TIFFs are not supported; "
                                       f"only 8-bit RGB")
            self.is_rgb, self.kind, self.fill_value, self.n_channels = \
                True, "brightfield", 255, 3
            self.channel_names = ["Red", "Green", "Blue"]
            self.channel_colors = [(255, 0, 0), (0, 255, 0), (0, 0, 255)]
        elif spp == 1 and phot == 1:  # minisblack
            if self.dtype not in (np.uint8, np.uint16):
                raise SlideFormatError(f"{name}: {self.dtype} pixels are not supported "
                                       f"(8- or 16-bit only)")
            self.is_rgb, self.kind, self.fill_value, self.n_channels = \
                False, "grayscale", 0, 1
            self.channel_names, self.channel_colors = ["Gray"], [(255, 255, 255)]
        else:
            what = {3: "palette (indexed colour)", 0: "inverted grayscale",
                    5: "CMYK"}.get(phot, f"{spp}-sample, photometric {phot}")
            raise SlideFormatError(f"{name}: {what} TIFFs are not supported; only RGB or "
                                   f"single-channel grayscale")
        m = re.search(r"AppMag\s*=\s*([\d.]+)", self.description)
        self._resolve_magnification(m.group(1) if m else None, "description", magnification)
        mpp, self.mpp_source = _mpp_from_resolution(base), "resolution tags"
        if mpp is None or not 0.01 <= mpp <= 100:
            mpp, self.mpp_source = float(mpp_fallback), "fallback"
            self._warn(f"no pixel size in file; using {mpp_fallback} um/px (set --mpp)")
        self.mpp = float(f"{mpp:.6g}")
        self.mpp_from_header = self.mpp_source != "fallback"
        self.is_pyramidal = len(stacks) > 1
        self.z_values = [0]
        self.associated = [s.name for s in tif.series[1:] if s.name]
        self._set_pyramid([pages[:1] for pages in stacks])

    def summary(self, min_level_size: int = 256, max_levels: int | None = None) -> str:
        src = self._src[0][0]
        lines = [
            f"File          : {self.path}",
            f"Type          : TIFF {'RGB' if self.is_rgb else 'grayscale'} "
            f"({self._pixel_text()}, {'BigTIFF' if self.is_bigtiff else 'classic TIFF'}, "
            f"{'pyramidal' if self.is_pyramidal else 'single resolution'})",
            f"Dimensions    : {self.width} x {self.height} px",
            f"Magnification : {self.magnification_text()}",
            f"MPP           : {self.mpp:.5g} um/px ({self.mpp_source})",
            f"Storage       : {'tiles' if src.tiled else 'strips'} {src.seg_w}x{src.seg_h} "
            f"{src.page.compression.name}",
            f"Associated    : {', '.join(self.associated) or '-'} (not exported)",
        ]
        lines += self._levels_lines(min_level_size, max_levels)
        return "\n".join(lines)


SLIDE_SUFFIXES = {".kfbf": "fluorescence", ".svs": "brightfield"}


def detect_slide_type(path: str | os.PathLike) -> str | None:
    """'fluorescence' / 'brightfield' from the file's magic bytes, else None."""
    with open(path, "rb") as fh:
        head = fh.read(8)
    if head[:4] == MAGIC and head[4:8] == SIGNATURE:
        return "fluorescence"
    if head[:4] in TIFF_MAGICS:
        return "brightfield"
    return None


def _peek_tiff(path: Path) -> dict | None:
    """First-page facts of a TIFF without reading its pyramid: written by
    us?, QPTIFF?, Aperio SVS? None if tifffile cannot open it."""
    try:
        with tifffile.TiffFile(path) as t:
            pg = t.pages[0]
            desc = pg.description or ""
            software = pg.software or ""
            return {"own": desc.startswith(APERIO_HEADER)
                    or software.startswith("kfbf_converter"),
                    "qptiff": bool(t.is_qpi), "svs": bool(t.is_svs)}
    except Exception:
        return None


def open_slide(path: str | os.PathLike, **kw) -> KFBioSlide | SVSSlide | QPTIFFSlide:
    """Open a slide with the right reader: magic bytes first, extension second.

    Keyword arguments (``mpp_fallback``, ``channel_names``, ``tile_cache``,
    ``magnification``) are passed to the reader; ``channel_names`` is ignored
    for brightfield. ``magnification`` overrides the objective magnification
    stored in the file. Any TIFF that tifffile identifies as QPTIFF opens as
    :class:`QPTIFFSlide`, whatever its extension.
    """
    path = Path(path)
    by_magic = detect_slide_type(path)
    by_ext = SLIDE_SUFFIXES.get(path.suffix.lower())
    if by_magic == "brightfield":
        peek = _peek_tiff(path)
        if peek is None:
            what = "QPTIFF" if path.suffix.lower() == ".qptiff" else "TIFF"
            raise SlideFormatError(f"{path.name}: this {what} is damaged or incomplete and "
                                   f"cannot be converted")
        if peek["own"]:  # our own output: never re-convert it
            raise SlideFormatError(f"{path.name}: already converted by kfbf_converter")
        if peek["qptiff"] or path.suffix.lower() == ".qptiff":
            return QPTIFFSlide(path, **{k: v for k, v in kw.items()
                                        if k in ("mpp_fallback", "channel_names",
                                                 "magnification")})
        is_tif = path.suffix.lower() in (".tif", ".tiff")
        if is_tif and not peek["svs"]:
            return TiffSlide(path, **{k: v for k, v in kw.items()
                                      if k in ("mpp_fallback", "magnification")})
        if by_ext != "brightfield" and not peek["svs"]:
            raise SlideFormatError(
                f"{path.name}: TIFF file is not an Aperio .svs slide")
    kind = by_magic or by_ext
    if kind is None:
        raise SlideFormatError(
            f"{path.name}: unrecognised slide file; expected a KFBio fluorescence "
            f"slide (.kfbf) or an Aperio brightfield slide (.svs)")
    if by_magic and by_ext and by_magic != by_ext:
        log.warning("%s: extension says %s but the file content is a %s slide; "
                    "reading it as %s", path.name, by_ext, by_magic, by_magic)
    if kind == "fluorescence":
        return KFBioSlide(path, **{k: v for k, v in kw.items()
                                   if k in ("mpp_fallback", "channel_names", "tile_cache",
                                            "magnification")})
    slide = SVSSlide(path, **{k: v for k, v in kw.items()
                              if k in ("mpp_fallback", "tile_cache", "magnification")})
    if path.suffix.lower() in (".tif", ".tiff"):
        slide.input_type = "TIFF"   # an Aperio-style .tif: read like .svs, named as TIFF
    return slide


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

    def __init__(self, slide: "_SlideBase", tile_size: int, compression: str,
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
        self.bigtiff = True  # always: classic TIFF stops at 4 GB

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
              photometric: str, metadata: dict | None, ome: bool, desc: str,
              tile_size: int | None = None,
              encoded: Callable[[Level], dict] | None = None,
              finalize: Callable[[Path], None] | None = None,
              thumbnail: Callable[[], np.ndarray] | None = None) -> Path:
        """fetch(level, plane, tx, ty) -> tile array (ts, ts, *sample_shape).

        With ``encoded``, ``encoded(level)`` gives the tifffile kwargs
        (``subsampling``, ``jpegtables``) of a level whose tiles are copied:
        fetch then returns already-compressed JPEG tile bytes for it. If
        ``encoded(level)`` is None that level is encoded normally (e.g. a
        computed level below copied ones). ``tile_size`` overrides the
        writer's tile size (it must equal the source tile size then).
        ``finalize(part_path)`` may edit the finished file (e.g. metadata)
        before it is verified and renamed.

        With ``thumbnail`` (single plane only) the file uses the Aperio
        layout that OpenSlide reads as a pyramid with objective power:
        level 0, then the thumbnail (``thumbnail()`` -> small array) as page
        1, then the other levels as top-level pages, each with an
        ``Aperio ...|AppMag = ..|MPP = ..`` description. Otherwise the
        lower levels are SubIFDs of level 0.
        """
        slide, ts = self.slide, tile_size or self.ts
        levels = self.levels
        lead = tuple(d for d in lead_shape if d > 1)
        assert math.prod(lead) == len(planes), (lead, len(planes))
        total = sum(len(planes) * math.ceil(lv.width / ts) * math.ceil(lv.height / ts)
                    for lv in levels)
        tmp = out_path.with_name(out_path.name + ".part")
        tile_shape = (ts, ts) + sample_shape
        bar = tqdm(total=total, desc=desc, unit="tile", leave=False,
                   disable=not self.show_progress, dynamic_ncols=True)
        level_shapes = [(lv.height, lv.width) + sample_shape for lv in levels]
        dtype = np.dtype(getattr(slide, "dtype", np.uint8))
        aperio = thumbnail is not None
        if aperio and len(planes) != 1:
            raise ValueError("the Aperio layout holds a single plane")
        enc0 = encoded(levels[0]) if encoded is not None else None
        comp0 = enc0.get("compression", "jpeg") if enc0 is not None else self.compression
        codec = "raw" if comp0 is None else (
            "JPEG" if comp0 == "jpeg" else _COMPRESSION_LABEL[comp0])
        try:
            with tifffile.TiffWriter(tmp, bigtiff=self.bigtiff, ome=ome) as tif:
                for lv in levels:
                    ncols, nrows = math.ceil(lv.width / ts), math.ceil(lv.height / ts)
                    jobs = ((p, tx, ty) for p in planes
                            for ty in range(nrows) for tx in range(ncols))
                    enc = encoded(lv) if encoded is not None else None

                    def make(job, _lv=lv.index, _raw=enc is not None):
                        p, tx, ty = job
                        tile = fetch(_lv, p, tx, ty)
                        if _raw:
                            return tile
                        if tile.shape != tile_shape:  # pad edge tiles
                            padded = np.full(tile_shape, slide.fill_value, dtype)
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
                        dtype=dtype, tile=(ts, ts), photometric=photometric,
                        resolution=(res, res), resolutionunit="CENTIMETER",
                        **({"compression": "jpeg", **enc} if enc is not None
                           else self._compression_kwargs(photometric == "rgb")),
                    )
                    if photometric == "rgb":
                        kw["planarconfig"] = "contig"
                    if aperio:
                        kw["metadata"] = None
                        kw["description"] = aperio_description(
                            slide, levels[0], lv, ts, codec, first=lv is levels[0])
                        # every level page is subfiletype 0, as in Aperio/KFBio
                        # files (1 = label, 9 = macro for SVS readers)
                        if lv is levels[0]:
                            kw["software"] = f"kfbf_converter {__version__}"
                            kw["extratags"] = [(315, "s", 0, CREDIT, True)]  # Artist
                    elif lv is levels[0]:
                        kw["subifds"] = len(levels) - 1
                        kw["metadata"] = metadata
                        kw["software"] = f"kfbf_converter {__version__}"
                        kw["extratags"] = [(315, "s", 0, CREDIT, True)]  # Artist
                    else:
                        kw["subfiletype"] = 1
                        kw["metadata"] = None
                    # tifffile reads ahead buffersize bytes of tiles to compress
                    # in parallel (default 512 MB); bound it to keep RAM flat.
                    # Older tifffile lacks `buffersize`: compress serially there,
                    # which streams tile-by-tile without read-ahead.
                    if enc is not None:
                        kw.update(maxworkers=1)  # nothing to compress
                    elif _WRITE_HAS_BUFFERSIZE:
                        kw.update(maxworkers=self.workers, buffersize=WRITE_BUFFER_BYTES)
                    else:
                        kw.update(maxworkers=1)
                    tif.write(tiles(), **kw)
                    if aperio and lv is levels[0]:  # Aperio: thumbnail is page 1
                        thumb = np.ascontiguousarray(thumbnail(), dtype=dtype)
                        tkw = dict(compression=self.compression or None)
                        if encoded is not None and dtype == np.uint8:
                            tkw = dict(compression="jpeg", compressionargs={"level": 90})
                        tif.write(thumb, photometric=photometric, metadata=None,
                                  description=f"Aperio Format thumbnail "
                                              f"{thumb.shape[1]}x{thumb.shape[0]}",
                                  **({"planarconfig": "contig"} if thumb.ndim == 3 else {}),
                                  **tkw)
            bar.close()
            if finalize is not None:
                finalize(tmp)
            # Never leave a file the user can't open: read it back first.
            verify_output(tmp, level_shapes, n_pages=len(planes), ome=ome,
                          layout="aperio" if aperio else "subifd")
            os.replace(tmp, out_path)
        except BaseException as exc:
            try:
                written = tmp.stat().st_size
            except OSError:
                written = 0
            try:
                tmp.unlink()
            except OSError:
                pass
            mapped = map_write_error(exc, out_path, written)
            if mapped is not None:
                raise mapped from exc
            raise
        finally:
            bar.close()
        return out_path


APERIO_HEADER = "Aperio Format written by kfbf_converter"
THUMBNAIL_MAX = 1024


def aperio_description(slide, level0: Level, lv: Level, tile: int, codec: str,
                       first: bool) -> str:
    """Aperio-style ImageDescription. OpenSlide reads ``AppMag`` as
    openslide.objective-power and ``MPP`` (level-0 µm/px) as openslide.mpp-x/y
    from the first page; every level repeats them, as Aperio/KFBio files do."""
    w0, h0 = level0.width, level0.height
    size = (f"{w0}x{h0} [0,0 {w0}x{h0}] ({tile}x{tile})" if first
            else f"{w0}x{h0} -> {lv.width}x{lv.height}")
    text = f"{APERIO_HEADER} {__version__}\r\n{size} {codec}"
    if slide.magnification:
        text += f"|AppMag = {format_magnification(slide.magnification)}"
    return text + f"|MPP = {slide.mpp:.6g}"


def thumbnail_level(slide) -> int:
    """Largest pyramid level (stored or synthetic) that fits THUMBNAIL_MAX."""
    for lv in slide.levels:
        if max(lv.width, lv.height) <= THUMBNAIL_MAX:
            return lv.index
    return slide.n_levels - 1


def _region(slide: "_SlideBase", level: int, tx: int, ty: int, ts: int) -> tuple[int, int, int, int]:
    lv = slide.levels[level]
    x, y = tx * ts, ty * ts
    return x, y, min(ts, lv.width - x), min(ts, lv.height - y)


OME_NS = "http://www.openmicroscopy.org/Schemas/OME/2016-06"
# OME 2016-06 schema order: OME = ... Experimenter*, ExperimenterGroup*, Instrument*, Image* ...;
# Image = ..., InstrumentRef?, ObjectiveSettings?, ImagingEnvironment?, StageLabel?, Pixels, ...
_OME_BEFORE_INSTRUMENT_REFS = ("<ImagingEnvironment", "<StageLabel", "<Pixels")


def ome_add_objective(xml: str, magnification: float) -> str:
    """Return OME-XML with ``Instrument > Objective NominalMagnification`` and
    every ``Image`` linked to it (``InstrumentRef`` + ``ObjectiveSettings``),
    inserted where the OME 2016-06 schema requires them."""
    if "<Instrument" in xml or "<ObjectiveSettings" in xml:
        raise ValueError("OME-XML already has an Instrument / ObjectiveSettings")
    first_image = re.search(r"<Image[\s>]", xml)
    if first_image is None:
        raise ValueError("OME-XML has no Image element")
    mag = f"{float(magnification):.1f}" if float(magnification).is_integer() \
        else f"{float(magnification):g}"
    instrument = (f'<Instrument ID="Instrument:0"><Objective ID="Objective:0:0" '
                  f'NominalMagnification="{mag}"/></Instrument>')
    refs = '<InstrumentRef ID="Instrument:0"/><ObjectiveSettings ID="Objective:0:0"/>'

    def link(m: re.Match) -> str:
        image = m.group(0)
        at = min((i for i in (image.find(t) for t in _OME_BEFORE_INSTRUMENT_REFS) if i > 0),
                 default=-1)
        if at < 0:
            raise ValueError("OME Image element has no Pixels")
        return image[:at] + refs + image[at:]

    i = first_image.start()
    xml = xml[:i] + instrument + re.sub(r"<Image[\s>].*?</Image>", link, xml[i:], flags=re.S)
    ElementTree.fromstring(xml.encode("utf-8"))  # must still be well-formed
    return xml


def set_ome_magnification(path: Path, magnification: float) -> None:
    """Rewrite the OME-XML ImageDescription of an OME-TIFF in place so it
    records the objective magnification (tifffile's writer can't)."""
    with tifffile.TiffFile(path) as tif:
        tag = tif.pages[0].tags["ImageDescription"]
        tif.filehandle.seek(tag.valueoffset)
        raw = tif.filehandle.read(tag.count)
    xml = raw.rstrip(b"\0").decode("utf-8")
    tifffile.tiffcomment(path, ome_add_objective(xml, magnification).encode("utf-8"))


def write_ome(slide: KFBioSlide | QPTIFFSlide, out_path: Path, writer: _PyramidWriter,
              copy_tiles: bool = False) -> Path:
    """Multi-channel pyramidal OME-TIFF: (C,Y,X) or (Z,C,Y,X) per level.

    With ``copy_tiles`` (QPTIFF with 8-bit JPEG tiles) the stored tiles of
    each stored level are copied unchanged; computed levels are encoded with
    the writer's compression. The writer's tile size must be the source's.
    """
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

    plan: dict = {}
    if copy_tiles:
        plan, why = slide.copy_tiles_plan(writer.levels)
        if plan is None:
            raise SlideFormatError(f"cannot copy source tiles: {why}")
        if writer.ts != slide.tile_size:
            raise ValueError("tile copy needs the source tile size")

    def fetch(level, plane, tx, ty):
        z, c = plane
        if level in plan:
            return slide.copied_tile(level, tx, ty, c, plan[level])
        x, y, w, h = _region(slide, level, tx, ty, writer.ts)
        return slide.read_region(level, x, y, w, h, c, z)

    def finalize(tmp: Path) -> None:
        if slide.magnification:
            set_ome_magnification(tmp, slide.magnification)

    return writer.write(out_path, planes=planes, lead_shape=(slide.n_z, slide.n_channels),
                        sample_shape=(), fetch=fetch,
                        photometric="minisblack", metadata=metadata, ome=True,
                        desc=f"OME  {slide.path.stem}", finalize=finalize,
                        encoded=(lambda lv: plan.get(lv.index)) if plan else None)


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

    def composite(level, z, x, y, w, h):
        acc = np.zeros((h, w, 3), np.float32)
        for c in range(slide.n_channels):
            acc += luts[c][slide.read_region(level, x, y, w, h, c, z)]
        np.clip(acc, 0, 255, out=acc)
        return (acc + 0.5).astype(np.uint8)

    def fetch(level, z, tx, ty):
        return composite(level, z, *_region(slide, level, tx, ty, writer.ts))

    def thumb():
        lv = slide.levels[thumbnail_level(slide)]
        return composite(lv.index, 0, 0, 0, lv.width, lv.height)

    metadata = None  # Z planes (if any) are written as consecutive pages
    return writer.write(out_path, planes=list(range(slide.n_z)), lead_shape=(slide.n_z,),
                        sample_shape=(3,),
                        fetch=fetch, photometric="rgb", metadata=metadata, ome=False,
                        desc=f"RGB  {slide.path.stem}",
                        thumbnail=thumb if slide.n_z == 1 else None)


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
                writer: _PyramidWriter, copy_tiles: bool = False) -> Path:
    """Single-channel grayscale pyramidal BigTIFF.

    With ``copy_tiles`` (a grayscale TIFF) the stored tiles of the levels in
    the slide's copy plan are written unchanged; other levels are encoded.
    """
    plan: dict = {}
    if copy_tiles:
        plan, why = slide.copy_tiles_plan(writer.levels)
        if plan is None:
            raise SlideFormatError(f"cannot copy source tiles: {why}")
        if writer.ts != slide.tile_size:
            raise ValueError("tile copy needs the source tile size")

    def fetch(level, z, tx, ty):
        if level in plan:
            return slide.copied_tile(level, tx, ty, channel, plan[level])
        x, y, w, h = _region(slide, level, tx, ty, writer.ts)
        return slide.read_region(level, x, y, w, h, channel, z)

    def thumb():
        lv = slide.levels[thumbnail_level(slide)]
        return slide.read_region(lv.index, 0, 0, lv.width, lv.height, channel)

    metadata = None  # Z planes (if any) are written as consecutive pages
    return writer.write(out_path, planes=list(range(slide.n_z)), lead_shape=(slide.n_z,),
                        sample_shape=(),
                        fetch=fetch, photometric="minisblack", metadata=metadata,
                        ome=False, desc=f"{slide.channel_names[channel]:<5.5s} {slide.path.stem}",
                        thumbnail=thumb if slide.n_z == 1 else None,
                        encoded=(lambda lv: plan.get(lv.index)) if plan else None)


def write_brightfield(slide: SVSSlide, out_path: Path, writer: _PyramidWriter,
                      copy_tiles: bool = False) -> Path:
    """RGB pyramidal BigTIFF for a brightfield slide (photometric RGB,
    contiguous samples, resolution in pixels/cm, Aperio page layout with
    AppMag/MPP so OpenSlide sees the pyramid and objective power).

    With ``copy_tiles`` the stored tiles are copied unchanged (no
    re-compression, so no extra quality loss). For .svs this requires
    JPEG/YCbCr tiles on every exported level, else SlideFormatError; for
    TIFF/QPTIFF pages, levels that can't be copied are encoded instead.
    """
    planes = list(range(slide.n_z))
    desc = f"BF   {slide.path.stem}"

    def thumb():
        lv = slide.levels[thumbnail_level(slide)]
        return slide.read_region(lv.index, 0, 0, lv.width, lv.height)

    if copy_tiles:
        plan, why = slide.copy_tiles_plan(writer.levels)
        if plan is None:
            raise SlideFormatError(f"cannot copy source tiles: {why}")
        ts = slide.tile_size

        def fetch_raw(level, z, tx, ty):
            if level not in plan:  # TIFF level that can't be copied: encode it
                x, y, w, h = _region(slide, level, tx, ty, ts)
                return slide.read_region(level, x, y, w, h)
            if isinstance(slide, SVSSlide):
                data = slide.raw_tile(level, tx, ty)
                return bytes(data) if data is not None else \
                    slide.blank_jpeg(plan[level]["subsampling"])
            return slide.copied_tile(level, tx, ty, 0, plan[level])

        return writer.write(out_path, planes=planes, lead_shape=(slide.n_z,),
                            sample_shape=(3,), fetch=fetch_raw, photometric="rgb",
                            metadata=None, ome=False, desc=desc,
                            tile_size=ts, encoded=lambda lv: plan.get(lv.index),
                            thumbnail=thumb)

    def fetch(level, z, tx, ty):
        x, y, w, h = _region(slide, level, tx, ty, writer.ts)
        return slide.read_region(level, x, y, w, h, z)

    return writer.write(out_path, planes=planes, lead_shape=(slide.n_z,), sample_shape=(3,),
                        fetch=fetch, photometric="rgb", metadata=None, ome=False, desc=desc,
                        thumbnail=thumb)


# --------------------------------------------------------------------------
# Size checks ("Data is too big for the file format chosen.")
# --------------------------------------------------------------------------

GiB = 1024 ** 3
CLASSIC_TIFF_MAX_BYTES = 4 * GiB - 1   # uint32 offsets
FAT_MAX_FILE_BYTES = 4 * GiB - 1       # FAT12/16/32 per-file limit
TIFF_MAX_DIM = 2 ** 32 - 1             # ImageWidth / ImageLength are uint32
JPEG_MAX_DIM = 65535                   # per JPEG stream, i.e. per tile
FAT_FILESYSTEMS = {"FAT", "FAT12", "FAT16", "FAT32", "VFAT", "MSDOS"}

# Expected compressed/uncompressed ratio, calibrated on real KFBio slides
# (worst seen: deflate 0.26, lzw 0.30, jpeg q90 0.073) with ~1.5x headroom
# because tissue-dense slides compress worse than the samples did.
COMPRESSION_FACTORS = {None: 1.0, "zlib": 0.40, "lzw": 0.45, "jpeg": 0.12}
_COMPRESSION_LABEL = {None: "none", "zlib": "deflate", "lzw": "lzw", "jpeg": "jpeg"}


def _gb(n: float) -> str:
    return f"{n / GiB:.2f} GB"


def compression_factor(compression: str | None, jpeg_quality: int = 90) -> float:
    if compression == "jpeg" and jpeg_quality > 90:  # grows steeply near q100
        return 0.12 + (min(jpeg_quality, 100) - 90) / 10 * 0.23
    return COMPRESSION_FACTORS.get(compression, 1.0)


@dataclass
class SizeEstimate:
    """Planned size of one output file (bytes)."""

    estimate: int        # expected file size
    uncompressed: int    # worst case: uncompressed pixel data + TIFF overhead
    exact: bool          # True for tile copy: estimate is the real size
    n_tiles: int         # tiles in the whole file
    tiles_per_ifd: int   # largest tile count of one page
    max_dim: int         # largest image width/height
    tile_size: int
    compression: str | None
    overhead: int


def _output_layout(slide, kind: str) -> tuple[int, int]:
    """(samples per pixel, planes) of one output kind."""
    if kind == "ome":
        return 1, slide.n_z * slide.n_channels
    if kind in ("rgb", "brightfield"):
        return 3, slide.n_z
    return 1, slide.n_z  # split


def estimate_output(slide, kind: str, writer: "_PyramidWriter",
                    copy_tiles: bool = False) -> SizeEstimate:
    """Estimate one output's size: uncompressed = sum over exported levels of
    tiles x tile_size^2 x samples x planes; expected = uncompressed x
    compression factor (exact source bytes for a tile copy) + TIFF overhead
    (tile offset/bytecount tables, IFDs, OME-XML)."""
    samples, planes = _output_layout(slide, kind)
    ts = slide.tile_size if copy_tiles else writer.ts
    grids = [math.ceil(lv.width / ts) * math.ceil(lv.height / ts) for lv in writer.levels]
    n_tiles = sum(grids) * planes
    entry = 16 if writer.bigtiff else 8   # one offset + one bytecount
    overhead = n_tiles * entry + len(writer.levels) * planes * 1024 + 256 * 1024
    if kind == "ome":
        overhead += planes * 512
    elif planes == 1:  # Aperio layout: thumbnail page (<= 1024 x 1024 x samples)
        thumb = THUMBNAIL_MAX * THUMBNAIL_MAX * samples
        overhead += int(thumb * compression_factor(
            "jpeg" if copy_tiles else writer.compression, writer.jpeg_quality))
    itemsize = np.dtype(getattr(slide, "dtype", np.uint8)).itemsize
    uncompressed = n_tiles * ts * ts * samples * itemsize + overhead
    if copy_tiles and hasattr(slide, "copy_bytes"):   # TIFF/QPTIFF: copied + encoded levels
        data, exact = slide.copy_bytes(writer.levels, writer.jpeg_quality, writer.compression)
        estimate = data + overhead
    elif copy_tiles:
        data = 0
        for lv in writer.levels:
            stored, missing = slide.stored_bytes(lv.index)
            data += stored + missing * 2048  # blank white tile is ~1-2 KB
        estimate, exact = data * planes + overhead, True
    else:
        estimate = int(uncompressed * compression_factor(writer.compression,
                                                         writer.jpeg_quality)) + overhead
        exact = False
    return SizeEstimate(estimate, uncompressed, exact, n_tiles, max(grids),
                        max(max(lv.width, lv.height) for lv in writer.levels), ts,
                        writer.compression, overhead)


def _volume_root(path: Path) -> str:
    """Drive / mount point that holds ``path`` (for messages)."""
    path = Path(os.path.abspath(path))
    if os.name == "nt":
        buf = _volume_path_name(path)
        if buf:
            return buf
        return path.anchor or str(path)
    p = path
    while not os.path.ismount(p) and p.parent != p:
        p = p.parent
    return str(p)


def _volume_path_name(path: Path) -> str | None:
    try:
        import ctypes
        k32 = ctypes.WinDLL("kernel32", use_last_error=True)
        buf = ctypes.create_unicode_buffer(1024)
        if k32.GetVolumePathNameW(ctypes.c_wchar_p(str(path)), buf, 1024):
            return buf.value
    except Exception:  # pragma: no cover - non-Windows or ctypes failure
        pass
    return None


def filesystem_type(path: str | os.PathLike) -> str | None:
    """File-system name of the volume holding ``path`` ('NTFS', 'exFAT',
    'FAT32', 'ext4', ...), or None if it cannot be determined."""
    path = Path(os.path.abspath(path))
    if os.name == "nt":
        root = _volume_path_name(path)
        if not root:
            return None
        try:
            import ctypes
            k32 = ctypes.WinDLL("kernel32", use_last_error=True)
            fs = ctypes.create_unicode_buffer(261)
            ok = k32.GetVolumeInformationW(ctypes.c_wchar_p(root), None, 0, None, None,
                                           None, fs, 261)
            return (fs.value or None) if ok else None
        except Exception:  # pragma: no cover
            return None
    if sys.platform.startswith("linux"):
        try:
            real = os.path.realpath(path)
            best, fstype = "", None
            with open("/proc/self/mounts", encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    parts = line.split()
                    if len(parts) < 3:
                        continue
                    mnt = parts[1].replace("\\040", " ")
                    if (real == mnt or real.startswith(mnt.rstrip("/") + "/")) \
                            and len(mnt) >= len(best):
                        best, fstype = mnt, parts[2]
            return fstype
        except OSError:
            return None
    return None  # macOS etc.: not detected; the free-space check still runs


def is_fat(fstype: str | None) -> bool:
    return bool(fstype) and fstype.upper() in FAT_FILESYSTEMS


def free_space_margin(estimate: int) -> int:
    """Space kept free on top of the estimate: 256 MiB + 5 %."""
    return 256 * 2 ** 20 + estimate // 20


def preflight_check(est: SizeEstimate, out_path: Path, bigtiff: bool = True,
                    jpeg_estimate: int | None = None) -> list[str]:
    """Refuse an output that cannot fit; return warnings for ones that might not.

    Raises DataTooBigError with the reason and a suggested fix.
    ``jpeg_estimate`` is the size with --compression jpeg, used to suggest it.
    """
    name = out_path.name
    warnings: list[str] = []

    def jpeg_hint(limit: int) -> str:
        """Suggest --compression jpeg only if it would bring the file under ``limit``."""
        if jpeg_estimate is not None and jpeg_estimate < est.estimate \
                and jpeg_estimate <= limit:
            return f", or use --compression jpeg (about {_gb(jpeg_estimate)})"
        return ""

    # --- TIFF format limits ---------------------------------------------
    if est.max_dim > TIFF_MAX_DIM:
        raise DataTooBigError(
            f"{name}: the image is {est.max_dim} px wide/high; TIFF stores "
            f"sizes as 32-bit numbers (max {TIFF_MAX_DIM} px).",
            "export fewer or smaller levels (--max-levels).")
    if not bigtiff:
        if est.estimate > CLASSIC_TIFF_MAX_BYTES:
            raise DataTooBigError(
                f"{name}: needs about {_gb(est.estimate)}, but classic TIFF files "
                f"are limited to {_gb(CLASSIC_TIFF_MAX_BYTES + 1)}.",
                "write BigTIFF" + jpeg_hint(CLASSIC_TIFF_MAX_BYTES) + ".")
        if est.tiles_per_ifd > 2 ** 32 - 1:
            raise DataTooBigError(
                f"{name}: {est.tiles_per_ifd} tiles in one image exceed the classic "
                f"TIFF limit of {2 ** 32 - 1}.", "write BigTIFF or use a larger --tile-size.")
        if est.uncompressed > CLASSIC_TIFF_MAX_BYTES:
            warnings.append(f"{name}: could exceed the 4 GB classic TIFF limit if it "
                            f"compresses poorly (up to {_gb(est.uncompressed)}).")
    if est.compression == "jpeg" and est.tile_size > JPEG_MAX_DIM:
        raise DataTooBigError(
            f"{name}: JPEG tiles are limited to {JPEG_MAX_DIM} px, tile size is "
            f"{est.tile_size} px.", "use a smaller --tile-size.")

    out_dir = out_path.parent
    volume = _volume_root(out_dir)

    # --- destination file system --------------------------------------------
    fstype = filesystem_type(out_dir)
    if is_fat(fstype):
        if est.estimate > FAT_MAX_FILE_BYTES:
            raise DataTooBigError(
                f"{name}: would be about {_gb(est.estimate)}, but the output drive "
                f"{volume} is {fstype}, which cannot store files larger than "
                f"{_gb(FAT_MAX_FILE_BYTES + 1)}.",
                "choose an NTFS or exFAT drive with -o" + jpeg_hint(FAT_MAX_FILE_BYTES) + ".")
        if not est.exact and est.uncompressed > FAT_MAX_FILE_BYTES:
            warnings.append(
                f"{name}: the output drive {volume} is {fstype} (max 4 GB per file). "
                f"Expected size is {_gb(est.estimate)}, but it could reach "
                f"{_gb(est.uncompressed)} if the image compresses poorly. If it fails, "
                f"choose an NTFS or exFAT drive with -o.")

    # --- free space -------------------------------------------------------
    try:
        free = shutil.disk_usage(out_dir).free
    except OSError:
        free = None
    if free is not None:
        margin = free_space_margin(est.estimate)
        if est.estimate + margin > free:
            need = est.estimate + margin - free
            raise DataTooBigError(
                f"{name}: needs about {_gb(est.estimate)} (plus a {_gb(margin)} safety "
                f"margin), but only {_gb(free)} is free on {volume}.",
                f"free up {_gb(need)} on {volume}, or choose another drive with -o"
                + jpeg_hint(free - free_space_margin(jpeg_estimate or 0)) + ".")
    return warnings


# OSError codes meaning "disk full" / "file too large for this file system"
_ENOSPC_WINERRORS = {112, 39}          # ERROR_DISK_FULL, ERROR_HANDLE_DISK_FULL
_EFBIG_WINERRORS = {223, 665}          # ERROR_FILE_TOO_LARGE, ..._SYSTEM_LIMITATION


def map_write_error(exc: BaseException, out_path: Path, written: int) -> DataTooBigError | None:
    """Turn a disk-full / file-too-large error raised while writing into a
    DataTooBigError with the real cause, or return None."""
    name, out_dir = out_path.name, out_path.parent
    volume = _volume_root(out_dir)
    if isinstance(exc, OSError):
        winerr = getattr(exc, "winerror", None)
        fstype = filesystem_type(out_dir)
        if exc.errno == errno.ENOSPC or winerr in _ENOSPC_WINERRORS:
            try:
                free = _gb(shutil.disk_usage(out_dir).free)
            except OSError:
                free = "?"
            return DataTooBigError(
                f"{name}: the drive {volume} ran out of space while writing "
                f"({_gb(written)} written, {free} free).",
                f"free up space on {volume}, or choose another drive with -o.")
        near_fat_limit = is_fat(fstype) and written >= FAT_MAX_FILE_BYTES - 256 * 2 ** 20
        if exc.errno == errno.EFBIG or winerr in _EFBIG_WINERRORS or near_fat_limit:
            limit = (f"{fstype}, which cannot store files larger than "
                     f"{_gb(FAT_MAX_FILE_BYTES + 1)}" if is_fat(fstype)
                     else f"{fstype or 'its file system'}, which refused a larger file")
            return DataTooBigError(
                f"{name}: reached {_gb(written)}, but the output drive {volume} is {limit}.",
                "choose an NTFS or exFAT drive with -o, or use --compression jpeg.")
    if isinstance(exc, (ValueError, struct.error)):
        msg = str(exc).lower()
        if ("too large" in msg or "bigtiff" in msg or "exceeds" in msg
                or ("out of range" in msg and isinstance(exc, struct.error))):
            return DataTooBigError(f"{name}: tifffile refused the data ({exc}).",
                                   "write BigTIFF, or use --compression jpeg.")
    return None


def verify_output(path: Path, level_shapes: Sequence[tuple], n_pages: int,
                  ome: bool = False, layout: str = "subifd") -> None:
    """Re-open a written file and check it is usable: page and SubIFD counts,
    the shape of every level, tile data inside the file, and one decoded
    tile per level. Raises VerificationError.

    ``layout="subifd"``: ``n_pages`` top-level pages (planes), each with the
    lower levels as SubIFDs. ``layout="aperio"``: one plane; level 0, a
    thumbnail, then the other levels as top-level pages."""
    try:
        size = path.stat().st_size
        with tifffile.TiffFile(path) as tif:
            pages = tif.pages
            if layout == "aperio":
                expected = len(level_shapes) + 1
                if len(pages) != expected:
                    raise VerificationError(f"expected {expected} image page(s) "
                                            f"(levels + thumbnail), found {len(pages)}")
                if not tif.is_svs:
                    raise VerificationError("Aperio description is missing")
                thumb = pages[1].asarray()
                if thumb.size == 0 or pages[1].is_tiled:
                    raise VerificationError("thumbnail page is missing or invalid")
                plane_levels = [[pages[0]] + [pages[i] for i in range(2, expected)]]
            else:
                if len(pages) != n_pages:
                    raise VerificationError(f"expected {n_pages} image page(s), "
                                            f"found {len(pages)}")
                plane_levels = []
                for pi in range(n_pages):
                    top = pages[pi]
                    subs = tifffile.TiffPages(top) if top.subifds else []
                    plane_levels.append([top] + list(subs))
            if ome and not tif.is_ome:
                raise VerificationError("OME-XML metadata is missing")
            fh = tif.filehandle
            for pi, levels in enumerate(plane_levels):
                if len(levels) != len(level_shapes):
                    raise VerificationError(f"expected {len(level_shapes)} pyramid level(s), "
                                            f"found {len(levels)}")
                for li, (pg, shape) in enumerate(zip(levels, level_shapes)):
                    if tuple(pg.shape) != tuple(shape):
                        raise VerificationError(f"level {li} has shape {tuple(pg.shape)}, "
                                                f"expected {tuple(shape)}")
                    offs = np.asarray(pg.dataoffsets, np.int64)
                    cnts = np.asarray(pg.databytecounts, np.int64)
                    if len(offs) == 0 or (cnts <= 0).any() or (offs + cnts > size).any():
                        raise VerificationError(f"level {li} tile data is missing or past "
                                                f"the end of the file (truncated?)")
                    if pi == 0:  # decode one tile per level
                        j = len(offs) // 2
                        decode, tables = pg.decode, pg.jpegtables  # may seek: get first
                        fh.seek(int(offs[j]))
                        data = fh.read(int(cnts[j]))
                        seg = decode(data, j, jpegtables=tables)[0]
                        if seg is None or np.asarray(seg).size == 0:
                            raise VerificationError(f"level {li} tile {j} decoded to nothing")
    except VerificationError:
        raise
    except Exception as exc:
        raise VerificationError(f"could not read it back ({type(exc).__name__}: {exc})") \
            from exc


# --------------------------------------------------------------------------
# High-level conversion
# --------------------------------------------------------------------------

MODES = ("ome", "rgb", "split", "all")


def _slide_stem(path: Path) -> str:
    return path.name[: -len(path.suffix)] if path.suffix else path.name


def planned_outputs(slide: "_SlideBase", output_dir: Path, mode: str | None
                    ) -> list[tuple[str, Path, int]]:
    """[(kind, output path, channel)] for a slide. Brightfield slides always
    give one RGB BigTIFF; ``mode`` applies to fluorescence slides only."""
    stem = _slide_stem(slide.path)
    if getattr(slide, "input_type", "") == "TIFF":
        out = tiff_output_path(slide.path, output_dir)
        return [("brightfield", out, -1) if slide.kind == "brightfield" else ("gray", out, 0)]
    if slide.kind == "brightfield":
        return [("brightfield", output_dir / f"{stem}.tif", -1)]
    if getattr(slide, "input_type", "") == "QPTIFF":   # fluorescence QPTIFF: OME-TIFF only
        return [("ome", output_dir / f"{stem}.ome.tif", -1)]
    mode = mode or "ome"
    jobs: list[tuple[str, Path, int]] = []
    if mode in ("ome", "all"):
        jobs.append(("ome", output_dir / f"{stem}.ome.tif", -1))
    if mode in ("rgb", "all"):
        jobs.append(("rgb", output_dir / f"{stem}.rgb.tif", -1))
    if mode in ("split", "all"):
        for i, n in enumerate(safe_channel_names(slide.channel_names)):
            jobs.append(("split", output_dir / f"{stem}_{n}.tif", i))
    return jobs


def _same_path(a: Path, b: Path) -> bool:
    try:
        if a.exists() and b.exists():
            return os.path.samefile(a, b)
    except OSError:
        pass
    return os.path.normcase(os.path.abspath(a)) == os.path.normcase(os.path.abspath(b))


def tiff_output_path(src: Path, output_dir: Path) -> Path:
    """Output name for a .tif/.tiff input: ``<stem>.tif``, or
    ``<stem>.bigtiff.tif`` when that name would clash -- with the input
    itself, with a file this tool didn't write, or with the output of a
    ``<stem>.svs`` / ``<stem>.qptiff`` next to the input."""
    stem = _slide_stem(src)
    plain, alt = output_dir / f"{stem}.tif", output_dir / f"{stem}.bigtiff.tif"
    if _same_path(plain, src):
        return alt
    if plain.exists() and not (_peek_tiff(plain) or {}).get("own"):
        return alt
    try:
        siblings = {p.name.lower() for p in src.parent.iterdir()}
    except OSError:
        siblings = set()
    if any(f"{stem}{ext}".lower() in siblings for ext in (".svs", ".qptiff")):
        return alt
    return plain


DEFAULT_COMPRESSION = "deflate"
DEFAULT_TILE_SIZE = 512


def report_output_error(exc: DataTooBigError | VerificationError) -> None:
    """Print a user-facing failure message (no traceback) above progress bars."""
    for line in exc.lines():
        tqdm.write(line, file=sys.stderr)


def _preflight_output(slide, kind: str, out: Path, writer: _PyramidWriter,
                      copy_tiles: bool) -> SizeEstimate:
    est = estimate_output(slide, kind, writer, copy_tiles)
    jpeg_est = None
    if not copy_tiles and writer.compression != "jpeg":
        jpeg_est =int((est.uncompressed - est.overhead) * compression_factor("jpeg")) \
            + est.overhead
    log.debug("%s: preflight estimate %s (uncompressed bound %s, %s, %d tiles)", out.name,
              _gb(est.estimate), _gb(est.uncompressed),
              "exact" if est.exact else _COMPRESSION_LABEL.get(writer.compression), est.n_tiles)
    for msg in preflight_check(est, out, writer.bigtiff, jpeg_est):
        log.warning("%s", msg)
    return est


def _setup_tile_copy(slide: "_TiffPagesSlide", export: Sequence[Level], path: Path,
                     writer_ts: int, writer_comp: str) -> tuple[bool, int, str]:
    """Plan a tile copy for a TIFF / QPTIFF: returns (copy_tiles, tile size,
    compression for the levels that are re-encoded) and logs the outcome."""
    plan, why = slide.copy_tiles_plan(export)
    if plan is None:
        log.info("%s: the stored tiles can't be copied (%s); re-encoding them",
                 path.name, why)
        return False, writer_ts, writer_comp
    comp = slide.rest_compression(plan)
    copied = sorted({kw["compression"] for kw in plan.values()})
    names = "/".join("JPEG" if c == "jpeg" else _COMPRESSION_LABEL[c] for c in copied)
    redo = [lv for lv in export if lv.index not in plan]
    log.info("%s: copying the stored %s tiles unchanged (no re-compression)%s", path.name,
             names, f"; {len(redo)} level(s) computed or re-encoded as "
                    f"{'JPEG' if comp == 'jpeg' else comp}" if redo else "")
    if why:
        log.info("%s: %s", path.name, why)
    return True, slide.tile_size, comp


def convert_slide(path: str | os.PathLike, output_dir: str | os.PathLike | None = None,
                  mode: str | None = "ome", compression: str | None = None,
                  tile_size: int | None = None,
                  mpp: float = DEFAULT_MPP, overwrite: bool = False, workers: int | None = None,
                  channel_names: Sequence[str] | None = None, low_pct: float = 0.1,
                  high_pct: float = 99.9, jpeg_quality: int = 90,
                  show_progress: bool = True, min_level_size: int = 256,
                  max_levels: int | None = None, preflight: bool = True,
                  errors: list | None = None,
                  magnification: float | None = None,
                  skipped: list | None = None) -> list[Path]:
    """Convert one slide: .kfbf (KFBio fluorescence), .svs (brightfield),
    .qptiff (PerkinElmer/Akoya) or a plain .tif/.tiff.

    Returns the list of files written. ``mode`` applies to .kfbf slides;
    .svs slides always produce ``<stem>.tif``, fluorescence .qptiff slides
    ``<stem>.ome.tif``, and .tif/.tiff images ``<stem>.tif`` (or
    ``<stem>.bigtiff.tif`` when that name is taken). Outputs that already
    exist are skipped (appended to ``skipped`` if it is a list) unless
    ``overwrite``.
    ``compression`` / ``tile_size`` default to deflate / 512 px; for
    brightfield, leaving both unset copies the stored JPEG tiles unchanged.

    Before each output is written its size is checked against the TIFF
    format, the destination file system and the free space
    (``preflight=False`` skips this). An output that is too big
    (DataTooBigError) or fails read-back verification (VerificationError)
    is deleted. If ``errors`` is a list, such per-output failures are
    printed and appended to it as ``(output_path, exception)`` and the
    remaining outputs are still written; otherwise the exception is raised.

    ``magnification`` overrides the objective magnification stored in the
    slide (written to the OME-XML ``Objective``).
    """
    if mode is not None and mode not in MODES:
        raise ValueError(f"mode must be one of {MODES}")
    path = Path(path)
    out_dir = Path(output_dir) if output_dir else path.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    workers = workers or min(8, os.cpu_count() or 1)
    written: list[Path] = []
    with open_slide(path, mpp_fallback=mpp, channel_names=channel_names,
                    magnification=magnification) as slide:
        export = select_export_levels(slide.levels, min_level_size, max_levels)
        levels_txt = ", ".join(f"{lv.width}x{lv.height}" for lv in export)
        itype = slide.input_type or "slide"
        copy_tiles = False
        writer_ts = tile_size or DEFAULT_TILE_SIZE
        writer_comp = compression or DEFAULT_COMPRESSION
        if slide.dtype.itemsize > 1 and writer_comp == "jpeg":
            raise SlideFormatError(f"{path.name}: this {itype} has {slide.dtype.itemsize * 8}"
                                   f"-bit pixels, which JPEG cannot store; use "
                                   f"--compression deflate or lzw")
        if itype == "TIFF" and getattr(slide, "is_pyramidal", False):
            log.info("%s: this TIFF is already pyramidal; writing it again with every "
                     "level and the metadata", path.name)
        if slide.kind == "grayscale":
            log.info("%s: TIFF grayscale %dx%d, %d-bit, %.4g um/px, exporting %d of %d "
                     "levels (%s)", path.name, slide.width, slide.height,
                     slide.dtype.itemsize * 8, slide.mpp, len(export), slide.n_levels,
                     levels_txt)
            if mode is not None:
                log.info("%s: --mode %s applies to KFBF slides only; writing the TIFF as "
                         "a grayscale BigTIFF", path.name, mode)
            if compression is None and tile_size is None:
                copy_tiles, writer_ts, writer_comp = _setup_tile_copy(
                    slide, export, path, writer_ts, writer_comp)
        elif slide.kind == "brightfield":
            log.info("%s: %sbrightfield %dx%d RGB, %.4g um/px, exporting %d of %d levels (%s)",
                     path.name, "" if itype == "SVS" else f"{itype} ", slide.width,
                     slide.height, slide.mpp, len(export), slide.n_levels, levels_txt)
            if itype == "QPTIFF":
                log.info("%s: this QPTIFF is brightfield (RGB); writing an RGB BigTIFF "
                         "instead of an OME-TIFF", path.name)
            if mode is not None:
                if itype == "SVS":
                    log.info("%s: --mode %s applies to fluorescence (.kfbf) slides only; "
                             "writing the brightfield RGB BigTIFF", path.name, mode)
                else:
                    log.info("%s: --mode %s applies to KFBF slides only; writing the "
                             "%s as an RGB BigTIFF", path.name, mode, itype)
            if channel_names:
                log.info("%s: --channel-names ignored for brightfield slides", path.name)
            if compression is None and tile_size is None and isinstance(slide, SVSSlide):
                plan, why = slide.copy_tiles_plan(export)
                copy_tiles = plan is not None
                if copy_tiles:
                    log.info("%s: copying the stored JPEG tiles unchanged "
                             "(no re-compression)", path.name)
                else:
                    log.info("%s: %s; re-encoding with %s", path.name, why,
                             DEFAULT_COMPRESSION)
            elif compression is None and tile_size is None:  # TIFF / RGB QPTIFF
                copy_tiles, writer_ts, writer_comp = _setup_tile_copy(
                    slide, export, path, writer_ts, writer_comp)
            if not copy_tiles and compression is None and itype == "TIFF":
                # Re-encoding a JPEG-sized slide with deflate makes it ~10x bigger.
                writer_comp = "jpeg"
                log.info("%s: writing JPEG tiles (quality %d) to keep the output close to "
                         "the input size; use --compression deflate for lossless",
                         path.name, jpeg_quality)
        elif itype == "QPTIFF":
            log.info("%s: QPTIFF fluorescence %dx%d, %d ch %s, %d-bit, %.4g um/px, "
                     "exporting %d of %d levels (%s)", path.name, slide.width, slide.height,
                     slide.n_channels, slide.channel_names, slide.dtype.itemsize * 8,
                     slide.mpp, len(export), slide.n_levels, levels_txt)
            if mode not in (None, "ome"):
                log.info("%s: --mode %s applies to KFBF slides only; a QPTIFF is written "
                         "as an OME-TIFF", path.name, mode)
            if compression is None and tile_size is None:
                copy_tiles, writer_ts, writer_comp = _setup_tile_copy(
                    slide, export, path, writer_ts, writer_comp)
        else:
            log.info("%s: %s fluorescence %dx%d, %d ch %s, %d Z, exporting %d of %d levels "
                     "(%s)", path.name, itype, slide.width, slide.height, slide.n_channels,
                     slide.channel_names, slide.n_z, len(export), slide.n_levels, levels_txt)
        writer = _PyramidWriter(slide, writer_ts, writer_comp, workers, jpeg_quality,
                                show_progress, levels=export)
        jobs = planned_outputs(slide, out_dir, mode)
        for _, out, _ in jobs:  # leftovers of an earlier crashed/killed run
            stale = out.with_name(out.name + ".part")
            if stale.exists():
                try:
                    stale.unlink()
                    log.info("removed stale partial file %s left by an earlier "
                             "interrupted run", stale)
                except OSError as exc:
                    log.warning("could not remove stale partial file %s: %s", stale, exc)
        for kind, out, ch in jobs:
            if out.resolve() == path.resolve():
                raise SlideFormatError(f"{path.name}: output {out.name} would overwrite "
                                       f"the input slide; choose another -o folder")
            if out.exists() and not overwrite:
                log.info("skip existing %s (use --overwrite)", out)
                if skipped is not None:
                    skipped.append(out)
                continue
            try:
                if preflight:
                    _preflight_output(slide, kind, out, writer, copy_tiles and
                                      kind in ("brightfield", "ome", "gray"))
                if kind == "brightfield":
                    write_brightfield(slide, out, writer, copy_tiles=copy_tiles)
                elif kind == "ome":
                    write_ome(slide, out, writer, copy_tiles=copy_tiles)
                elif kind == "rgb":
                    write_rgb(slide, out, writer, low_pct, high_pct)
                else:  # "split" channel of a .kfbf, or a grayscale TIFF
                    write_split(slide, out, ch, writer,
                                copy_tiles=copy_tiles and kind == "gray")
            except (DataTooBigError, VerificationError) as exc:
                if errors is None:
                    raise
                report_output_error(exc)
                log.debug("%s failed", out, exc_info=True)
                errors.append((out, exc))
                continue
            log.info("wrote %s", out)
            written.append(out)
        if slide.decode_failures:
            log.warning("%s: %d tile(s) failed to decode and were written as %s",
                        path.name, slide.decode_failures,
                        "white" if slide.kind == "brightfield" else "zeros")
    return written


def find_slides(input_path: Path, recursive: bool = False) -> list[Path]:
    """A file is returned as-is; a directory is searched for .kfbf and .svs
    slides (exact suffix match, case-insensitive)."""
    if input_path.is_file():
        return [input_path]
    if input_path.is_dir():
        it = input_path.rglob("*") if recursive else input_path.iterdir()
        return sorted(p for p in it if p.is_file() and p.suffix.lower() in SLIDE_SUFFIXES)
    raise FileNotFoundError(f"input not found: {input_path}")


def slide_kind_by_suffix(path: Path) -> str | None:
    """'fluorescence' / 'brightfield' from the file extension (for counting)."""
    return SLIDE_SUFFIXES.get(path.suffix.lower())


# Input types, in menu order
INPUT_TYPES = ("TIFF", "SVS", "QPTIFF", "KFBF")
INPUT_SUFFIXES = {".tif": "TIFF", ".tiff": "TIFF", ".svs": "SVS", ".qptiff": "QPTIFF",
                  ".kfbf": "KFBF"}
TYPE_FLAGS = {"tif": "TIFF", "tiff": "TIFF", "svs": "SVS", "qptiff": "QPTIFF", "kfbf": "KFBF"}
SUPPORTED_TEXT = ".tif/.tiff, .svs, .qptiff or .kfbf"
# names this tool writes (".part" = unfinished output)
OUTPUT_NAME_ENDINGS = (".ome.tif", ".rgb.tif", ".bigtiff.tif", ".part")


def own_output_reason(path: Path, sibling_names: set[str] | None = None) -> str | None:
    """Why ``path`` is one of this tool's outputs (so never an input), or None:
    an output name (.ome.tif, .rgb.tif, .bigtiff.tif, .part), ``<stem>.tif``
    next to ``<stem>.svs``, a split channel ``<stem>_<channel>.tif`` next to
    ``<stem>.kfbf``, or a Software tag starting with "kfbf_converter" (catches
    renamed outputs)."""
    name = path.name.lower()
    if name.endswith(OUTPUT_NAME_ENDINGS):
        return "output file name"
    if path.suffix.lower() not in (".tif", ".tiff"):
        return None
    if sibling_names is None:
        try:
            sibling_names = {p.name.lower() for p in path.parent.iterdir()}
        except OSError:
            sibling_names = set()
    stem = _slide_stem(path).lower()
    if path.suffix.lower() == ".tif" and f"{stem}.svs" in sibling_names:
        return f"made from {stem}.svs"
    for i, ch in enumerate(stem):
        if ch == "_" and f"{stem[:i]}.kfbf" in sibling_names:
            return f"split channel of {stem[:i]}.kfbf"
    if (_peek_tiff(path) or {}).get("own"):
        return "written by kfbf_converter (Software tag)"
    return None


def input_type_of(path: Path) -> str | None:
    """'TIFF' / 'SVS' / 'QPTIFF' / 'KFBF' for counting and messages: the
    extension, except that any TIFF tifffile identifies as QPTIFF is QPTIFF."""
    kind = INPUT_SUFFIXES.get(path.suffix.lower())
    if kind in (None, "TIFF"):
        try:
            magic = detect_slide_type(path)
        except OSError:
            magic = None
        if magic == "fluorescence":
            return "KFBF"
        if magic == "brightfield":
            peek = _peek_tiff(path) or {}
            if peek.get("qptiff"):
                return "QPTIFF"
            if kind is None:
                return "SVS" if peek.get("svs") else "TIFF"
    return kind


def discover_slides(input_path: Path, recursive: bool = False) -> list[tuple[Path, str]]:
    """[(path, input type)]. A file is returned as-is; a folder is searched
    for .tif/.tiff, .svs, .qptiff and .kfbf files (case-insensitive), skipping
    anything this tool wrote (see :func:`own_output_reason`)."""
    if input_path.is_file():
        return [(input_path, input_type_of(input_path) or "?")]
    if not input_path.is_dir():
        raise FileNotFoundError(f"input not found: {input_path}")
    it = input_path.rglob("*") if recursive else input_path.iterdir()
    siblings: dict[Path, set[str]] = {}
    found = []
    for p in sorted(it):
        if not p.is_file() or p.suffix.lower() not in INPUT_SUFFIXES:
            continue
        names = siblings.get(p.parent)
        if names is None:
            names = siblings[p.parent] = {q.name.lower() for q in p.parent.iterdir()}
        why = own_output_reason(p, names)
        if why:
            log.debug("ignoring %s: %s", p.name, why)
            continue
        found.append((p, input_type_of(p) or "?"))
    return found


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
        description="Convert slides to open pyramidal formats: KFBio fluorescence "
                    "(.kfbf) to OME-TIFF, composite RGB BigTIFF and/or per-channel "
                    "grayscale BigTIFFs; PerkinElmer/Akoya .qptiff to OME-TIFF; "
                    "brightfield .svs and plain .tif/.tiff images to a pyramidal BigTIFF.")
    p.add_argument("input", type=Path, nargs="?", default=None,
                   help="a slide file, or a folder of .tif/.tiff, .svs, .qptiff and .kfbf "
                        "files (omit to use the current folder and choose from a menu)")
    p.add_argument("-m", "--mode", choices=MODES, default=None,
                   help="output type for .kfbf slides (default: ome); the other types "
                        "always give one file each")
    p.add_argument("--types", nargs="+", choices=sorted(TYPE_FLAGS), default=None,
                   metavar="TYPE",
                   help="only convert these slide types: tif svs qptiff kfbf "
                        "(default: all found); skips the 'what to convert' menu")
    p.add_argument("-y", "--yes", action="store_true",
                   help="don't ask anything: convert every slide found with the defaults")
    p.add_argument("-o", "--output-dir", type=Path, default=None,
                   help="output directory (default: alongside the input)")
    p.add_argument("--compression", choices=["deflate", "lzw", "jpeg"], default=None,
                   help="TIFF tile compression (default: deflate; for .svs the original "
                        "JPEG tiles are copied unchanged unless this or --tile-size is set)")
    p.add_argument("--jpeg-quality", type=int, default=90, help="quality for --compression jpeg")
    p.add_argument("--tile-size", type=int, default=None,
                   help=f"output tile size (default: {DEFAULT_TILE_SIZE}; .svs tile copy "
                        f"keeps the source tile size)")
    p.add_argument("--mpp", type=float, default=DEFAULT_MPP,
                   help=f"fallback µm/px if missing from header (default: {DEFAULT_MPP})")
    p.add_argument("--magnification", type=float, default=None, metavar="X",
                   help="objective magnification, e.g. 40 (default: from the slide); "
                        "use if the slide's value is missing or wrong")
    p.add_argument("--overwrite", action="store_true", help="overwrite existing outputs")
    p.add_argument("--recursive", action="store_true",
                   help="also search sub-folders")
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
    p.add_argument("--no-preflight", action="store_true",
                   help="skip the size check before writing (use only if it wrongly "
                        "refuses an output; written files are still verified)")
    p.add_argument("--info", action="store_true", help="print parsed slide structure and exit")
    p.add_argument("-q", "--quiet", action="store_true", help="hide progress bars")
    p.add_argument("-v", "--verbose", action="store_true", help="debug logging")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__} - {CREDIT}")
    return p


_MENU = [("ome", "OME-TIFF (Default / QuPath)"),
         ("rgb", "Composite RGB BigTIFF"),
         ("split", "Split Grayscale Channels"),
         ("all", "All of the above")]


def _read_answer(prompt: str) -> str:
    """input() with byte-order marks (e.g. text piped from PowerShell) removed."""
    choice = input(prompt)
    for bom in (chr(0xFEFF), chr(0xEF) + chr(0xBB) + chr(0xBF)):
        choice = choice.replace(bom, "")
    return choice.strip()


def _found_line(count: int, n_brightfield: int) -> str:
    n_fluor = count - n_brightfield
    if not n_brightfield:
        return f"Found {count} slide(s)."
    if not n_fluor:
        return f"Found {n_brightfield} brightfield (.svs) slide(s)."
    return (f"Found {n_fluor} fluorescence (.kfbf) and {n_brightfield} "
            f"brightfield (.svs) slide(s).")


def stdin_is_console() -> bool:
    """True when a person can answer prompts (stdin is an interactive console)."""
    try:
        return bool(sys.stdin) and sys.stdin.isatty()
    except Exception:
        return False


def _menu_choice(n: int) -> int | None:
    """Ask for a number 1..n; Enter = 1. Returns the 0-based choice, or None
    if cancelled (Ctrl+C, or input ran out after an invalid answer). Input
    that ends before any answer (nobody at the keyboard) picks the default."""
    attempted = False
    while True:
        try:
            choice = _read_answer(f"Enter choice [1-{n}] (Press Enter for 1): ")
        except EOFError:
            print()
            if attempted:  # input ran out after an invalid answer: don't guess
                print("No valid choice entered; cancelled.")
                return None
            return 0  # non-interactive stdin with no answer: default
        except KeyboardInterrupt:
            print("\nCancelled.")
            return None
        attempted = True
        if choice == "":
            return 0
        if choice.isdigit() and 1 <= int(choice) <= n:
            return int(choice) - 1
        print(f"Invalid choice, please enter a number from 1 to {n}.")


def _confirm(found_line: str, what: str) -> bool:
    print(found_line)
    print(what)
    try:
        answer = _read_answer("Press Enter to start (or type q to cancel): ")
    except EOFError:
        print()
        return True  # non-interactive stdin: proceed
    except KeyboardInterrupt:
        print("\nCancelled.")
        return False
    if answer.lower() in ("q", "quit", "n", "no"):
        print("Cancelled.")
        return False
    return True


def confirm_brightfield(count: int) -> bool:
    """Brightfield-only folder: one output type, so just confirm."""
    return _confirm(_found_line(count, count),
                    "Each one will be converted to an RGB pyramidal BigTIFF (<name>.tif).")


_TYPE_RESULT = {
    "TIFF": "a pyramidal BigTIFF (<name>.tif, or <name>.bigtiff.tif if that name is taken)",
    "SVS": "an RGB pyramidal BigTIFF (<name>.tif)",
    "QPTIFF": "a multi-channel pyramidal OME-TIFF (<name>.ome.tif)",
}


def confirm_type(input_type: str, count: int) -> bool:
    """A folder with one slide type (not .kfbf): say what will happen, confirm."""
    if input_type == "SVS":
        return confirm_brightfield(count)
    return _confirm(f"Found {count} {input_type} slide(s).",
                    f"Each one will be converted to {_TYPE_RESULT.get(input_type, 'BigTIFF')}.")


_TYPE_MENU = {"TIFF": "BigTIFF", "SVS": "BigTIFF", "QPTIFF": "OME-TIFF",
              "KFBF": "choose the format next"}


def prompt_types(counts: dict) -> list[str] | None:
    """'What do you want to convert?' menu for a folder with several slide
    types. ``counts`` = {input type: number of files}; only types found are
    listed. Returns the chosen types, or None if cancelled."""
    found = [t for t in INPUT_TYPES if counts.get(t)]
    print("Found " + ", ".join(f"{counts[t]} {t}" for t in found) + " slide(s).")
    print("What do you want to convert?")
    print("  [1] Convert all (Default - press Enter)")
    for i, t in enumerate(found, 2):
        print(f"  [{i}] Only {t:<6} -> {_TYPE_MENU[t]}")
    choice = _menu_choice(len(found) + 1)
    if choice is None:
        return None
    return found if choice == 0 else [found[choice - 1]]


def prompt_mode(count: int, n_brightfield: int = 0, header: str | None = None) -> str | None:
    """Interactive output-format menu for .kfbf slides. Returns a mode, or
    None if cancelled.

    ``count`` is the total number of slides, ``n_brightfield`` how many of
    them are brightfield (.svs); the menu applies to the fluorescence ones.
    ``header`` replaces the "Found ..." line.
    """
    if header:
        print(header)
    else:
        print(_found_line(count, n_brightfield))
        if n_brightfield:
            print("Brightfield (.svs) slides are always converted to an RGB BigTIFF "
                  "(<name>.tif);\nthe format below applies to the fluorescence (.kfbf) slides.")
    print("Select output format:")
    for i, (_, label) in enumerate(_MENU, 1):
        print(f"  [{i}] {label}")
    choice = _menu_choice(len(_MENU))
    return None if choice is None else _MENU[choice][0]


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    handler = _TqdmHandler()
    handler.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
    for old in [h for h in log.handlers if isinstance(h, _TqdmHandler)]:
        log.removeHandler(old)  # main() may be called repeatedly in-process
    log.addHandler(handler)
    log.setLevel(logging.DEBUG if args.verbose else logging.INFO)
    log.propagate = False

    if args.tile_size is not None and (args.tile_size < 16 or args.tile_size % 16):
        log.error("--tile-size must be a positive multiple of 16")
        return 2
    if args.magnification is not None and not (
            MIN_MAGNIFICATION <= args.magnification <= MAX_MAGNIFICATION):
        log.error("--magnification must be between %g and %g",
                  MIN_MAGNIFICATION, MAX_MAGNIFICATION)
        return 2
    if not 0 <= args.low_pct < args.high_pct <= 100:
        log.error("require 0 <= --low-pct < --high-pct <= 100")
        return 2

    no_arg = args.input is None
    if no_arg:
        args.input = Path.cwd()
    try:
        found = discover_slides(args.input, args.recursive)
    except FileNotFoundError as exc:
        log.error("%s", exc)
        return 2
    if args.types:
        wanted = {TYPE_FLAGS[t] for t in args.types}
        found = [(p, t) for p, t in found if t in wanted]
    if not found:
        if no_arg:
            print(f"No {SUPPORTED_TEXT} files found in the current folder ({args.input}).")
            return 0
        log.error("no %s files found in %s", SUPPORTED_TEXT, args.input)
        return 1

    # Ask what to convert when a person is there: always for the current
    # folder (no argument), and for a folder passed on the command line or
    # dropped on run_converter.bat when stdin is an interactive console.
    ask = (not args.yes and not args.info and args.input.is_dir()
           and (no_arg or stdin_is_console()))
    if ask:
        types = [t for t in INPUT_TYPES if any(ft == t for _, ft in found)]
        menu_shown = False
        if (len(types) > 1 and not args.types) or (
                args.mode is None and ("KFBF" in types or len(types) == 1)):
            print(MENU_BANNER)  # only when a question follows
        if len(types) > 1 and not args.types:
            chosen = prompt_types(Counter(t for _, t in found))
            if chosen is None:
                return 130
            found = [(p, t) for p, t in found if t in chosen]
            types, menu_shown = [t for t in types if t in chosen], True
        if args.mode is None:
            n_kfbf = sum(t == "KFBF" for _, t in found)
            if "KFBF" in types:
                header = None if len(types) == 1 and not menu_shown else \
                    f"Output format for the {n_kfbf} KFBF slide(s):"
                args.mode = prompt_mode(n_kfbf, header=header)
                if args.mode is None:
                    return 130
            elif len(types) == 1 and not menu_shown:
                if not confirm_type(types[0], len(found)):
                    return 130

    if args.info:
        rc = 0
        for s, _ in found:
            try:
                with open_slide(s, mpp_fallback=args.mpp, channel_names=args.channel_names,
                                magnification=args.magnification) as slide:
                    print(slide.summary(args.min_level_size, args.max_levels), end="\n\n")
            except Exception as exc:
                log.error("%s", _plain_error(s, exc))
                rc = 1
        return rc

    failures: list[tuple[Path, str]] = []
    stats = {t: Counter() for t in INPUT_TYPES + ("?",)}
    n_written = 0
    for s, itype in tqdm(found, desc="Slides", unit="slide",
                         disable=args.quiet or len(found) < 2, dynamic_ncols=True):
        out_dir = args.output_dir or s.parent
        output_errors: list = []
        skipped: list = []
        written: list = []
        failed = False
        try:
            written = convert_slide(
                s, out_dir, mode=args.mode, compression=args.compression,
                tile_size=args.tile_size, mpp=args.mpp, overwrite=args.overwrite,
                workers=args.workers, channel_names=args.channel_names,
                low_pct=args.low_pct, high_pct=args.high_pct,
                jpeg_quality=args.jpeg_quality, show_progress=not args.quiet,
                min_level_size=args.min_level_size, max_levels=args.max_levels,
                preflight=not args.no_preflight, errors=output_errors,
                magnification=args.magnification, skipped=skipped)
        except KeyboardInterrupt:
            log.error("interrupted")
            return 130
        except (DataTooBigError, VerificationError) as exc:  # raised outside an output
            report_output_error(exc)
            log.debug("traceback", exc_info=True)
            failures.append((s, exc.lines()[0] if isinstance(exc, VerificationError)
                             else f"{exc.HEADLINE} {exc.reason}"))
            failed = True
        except SlideFormatError as exc:  # unreadable / unsupported / damaged input
            log.error("%s", _plain_error(s, exc))
            log.debug("traceback", exc_info=True)
            failures.append((s, _plain_error(s, exc, with_name=False)))
            failed = True
        except Exception as exc:
            log.error("%s: conversion failed: %s", s.name, exc,
                      exc_info=args.verbose)
            failures.append((s, str(exc)))
            failed = True
        for out, exc in output_errors:
            failures.append((out, str(exc) if isinstance(exc, VerificationError)
                             else f"{exc.HEADLINE} {exc.reason}"))
        n_written += len(written)
        status = ("failed" if failed or output_errors else
                  "converted" if written else "skipped")
        stats[itype if itype in stats else "?"][status] += 1

    log.info("done: %d slide(s), %d file(s) written, %d failure(s)",
             len(found), n_written, len(failures))
    total = sum(stats.values(), Counter())
    log.info("Converted %d, skipped %d (already done), failed %d",
             total["converted"], total["skipped"], total["failed"])
    if sum(1 for c in stats.values() if c) > 1:
        for t, c in stats.items():
            if c:
                log.info("  %-6s: converted %d, skipped %d (already done), failed %d",
                         t, c["converted"], c["skipped"], c["failed"])
    for s, msg in failures:
        log.info("  FAILED %s: %s", s, msg)
    return 1 if failures else 0


def _plain_error(path: Path, exc: Exception, with_name: bool = True) -> str:
    """'<name>: <message>' without repeating the name the message starts with."""
    msg = str(exc)
    prefix = f"{path.name}: "
    if msg.startswith(prefix):
        msg = msg[len(prefix):]
    return f"{path.name}: {msg}" if with_name else msg


if __name__ == "__main__":
    sys.exit(main())
