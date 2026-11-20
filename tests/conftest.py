"""Synthetic KFBio fixtures built from the confirmed file layouts.

No real slide data is used or committed: every fixture is generated into
pytest's tmp_path.

- ``.kfbf``: the reverse-engineered layout documented in README
  (header, TLV metadata, 64-byte tile index records, per-channel
  grayscale JPEG payloads with C x uint64 offset / length arrays).
- ``.svs``: an Aperio-style classic TIFF like the KFBio sample (level 0,
  thumbnail, lower levels as top-level pages, label, macro; 256 px YCbCr
  4:2:2 JPEG tiles; MPP in the ImageDescription; small levels padded).
"""

from __future__ import annotations

import math
import struct
import sys
from dataclasses import dataclass, field
from pathlib import Path

import imagecodecs
import numpy as np
import pytest
import tifffile

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import kfbf_converter as kc  # noqa: E402

TS = 256


def _pattern(h: int, w: int, seed: int) -> np.ndarray:
    """Smooth, JPEG-friendly but position-dependent uint8 image."""
    yy, xx = np.mgrid[0:h, 0:w].astype(np.float32)
    v = 128 + 60 * np.sin(xx / (17 + seed)) + 50 * np.cos(yy / (23 + 2 * seed)) + seed * 3
    return np.clip(v, 0, 255).astype(np.uint8)


# --------------------------------------------------------------------------
# .kfbf
# --------------------------------------------------------------------------

@dataclass
class KfbfInfo:
    path: Path
    width: int
    height: int
    channels: list[str]
    n_z: int
    mpp: float
    levels: list[tuple[int, int]]                       # (w, h) per stored level
    missing: set = field(default_factory=set)           # {(level, z, col, row)}
    tiles: dict = field(default_factory=dict)           # (level, z, col, row, ch) -> decoded
    images: dict = field(default_factory=dict)          # (level, z, ch) -> full plane (decoded)


def make_kfbf(path: Path, width: int = 700, height: int = 520,
              channels=("DAPI", "FITC", "TRITC"), n_z: int = 1, mpp: float = 0.2528,
              missing=((0, 0, 1, 1),), quality: int = 95,
              magnification: int = 40) -> KfbfInfo:
    """Write a synthetic .kfbf. ``missing`` = tiles (level, z, col, row) with no record."""
    c = len(channels)
    info = KfbfInfo(path, width, height, list(channels), n_z, mpp, [], set(missing))
    levels, w, h = [], width, height
    while True:
        levels.append((w, h))
        if w <= 1 and h <= 1:
            break
        w, h = max(1, math.ceil(w / 2)), max(1, math.ceil(h / 2))
    info.levels = levels

    # TLV metadata block at 0x5C
    tlv = bytearray(b"\xff\x01\xee\xee")
    tlv += struct.pack("<II", 0x4B, 4) + struct.pack("<I", c)
    for i, name in enumerate(channels):
        raw = name.encode("latin-1") + b"\0"
        tlv += struct.pack("<II", 0x4C + i, len(raw)) + raw
    tlv += struct.pack("<II", 0, 0)

    body = bytearray(0x5C) + tlv
    records = []
    for li, (lw, lh) in enumerate(levels):
        scale = 40.0 / (2 ** li)
        for z in range(n_z):
            planes = [_pattern(lh, lw, seed=ch + 4 * z + li) for ch in range(c)]
            dec_planes = [np.zeros((lh, lw), np.uint8) for _ in range(c)]
            for row in range(math.ceil(lh / TS)):
                for col in range(math.ceil(lw / TS)):
                    if (li, z, col, row) in info.missing:
                        continue
                    offs, lens = [], []
                    for ch in range(c):
                        t = np.zeros((TS, TS), np.uint8)
                        part = planes[ch][row * TS:(row + 1) * TS, col * TS:(col + 1) * TS]
                        t[:part.shape[0], :part.shape[1]] = part
                        enc = imagecodecs.jpeg_encode(t, level=quality)
                        dec = imagecodecs.jpeg_decode(enc)
                        info.tiles[(li, z, col, row, ch)] = dec
                        dec_planes[ch][row * TS:row * TS + part.shape[0],
                                       col * TS:col * TS + part.shape[1]] = \
                            dec[:part.shape[0], :part.shape[1]]
                        offs.append(len(body))
                        lens.append(len(enc))
                        body += enc
                    ptr_off = len(body)
                    body += struct.pack(f"<{c}Q", *offs)
                    ptr_len = len(body)
                    body += struct.pack(f"<{c}Q", *lens)
                    vw, vh = min(TS, lw - col * TS), min(TS, lh - row * TS)
                    records.append((col * TS, row * TS, vw, vh, scale, z, ptr_off, ptr_len))
            for ch in range(c):
                info.images[(li, z, ch)] = dec_planes[ch]

    index_offset = len(body)
    body += b"\xf1\x04\xee\xee"
    for x, y, vw, vh, scale, z, po, pl in records:
        rec = bytearray(64)
        struct.pack_into("<IIIIfI", rec, 0, x, y, vw, vh, scale, z)
        struct.pack_into("<QQ", rec, 0x20, po, pl)
        struct.pack_into("<II", rec, 0x38, 0xEEEE04FF, 0xEEEE04F1)
        body += rec

    struct.pack_into("<4s4s", body, 0, b"\xf1\x01\xee\xee", b"KFBF")
    struct.pack_into("<f", body, 0x0C, 2.1)
    struct.pack_into("<4I", body, 0x10, len(records), height, width, magnification)
    body[0x20:0x28] = b"JPEG\0\0\0\0"
    struct.pack_into("<Q", body, 0x44, index_offset)
    struct.pack_into("<f", body, 0x4C, mpp)
    path.write_bytes(bytes(body))
    return info


# --------------------------------------------------------------------------
# .svs
# --------------------------------------------------------------------------

@dataclass
class SvsInfo:
    path: Path
    width: int
    height: int
    mpp: float
    level_sizes: list[tuple[int, int]]   # declared (true) sizes
    stored_sizes: list[tuple[int, int]]  # page sizes in file
    missing: set = field(default_factory=set)  # {(level, col, row)}


def _rgb_pattern(h: int, w: int, seed: int = 0) -> np.ndarray:
    """H&E-like: white background with pink/purple 'tissue' blobs."""
    yy, xx = np.mgrid[0:h, 0:w].astype(np.float32)
    tissue = (np.sin(xx / 40 + seed) * np.cos(yy / 30) > 0.2)
    img = np.full((h, w, 3), 245, np.uint8)
    img[tissue] = (200, 120, 190)
    img[..., 2] = np.clip(img[..., 2].astype(np.int16) - (xx % 50 < 3) * 60, 0, 255)
    return img


def _zero_tile_bytecounts(path: Path, page_index: int, tile_indices) -> None:
    """Mark tiles as missing by setting their TileByteCounts entries to 0."""
    with tifffile.TiffFile(path) as t:
        tag = t.pages[page_index].tags["TileByteCounts"]
        dtype = {3: "<H", 4: "<I", 16: "<Q"}[int(tag.dtype)]
        base = tag.valueoffset
    size = struct.calcsize(dtype)
    with open(path, "r+b") as fh:
        for i in tile_indices:
            fh.seek(base + i * size)
            fh.write(struct.pack(dtype, 0))


def make_svs(path: Path, width: int = 1400, height: int = 1000, mpp: float = 0.5,
             n_levels: int = 3, pad_last: bool = True, missing=((0, 1, 1),),
             quality: int = 90, aperio: bool = True, shrink_last: int = 0) -> SvsInfo:
    """Write a synthetic KFBio-style Aperio .svs (classic TIFF).

    ``pad_last`` pads the last level with white (as in KFBio's Sample2);
    ``shrink_last`` instead stores it that many pixels too narrow (as in a
    KFBio H&E slide whose 4609 px level is followed by a 2303 px one).
    """
    img = _rgb_pattern(height, width)
    desc = ("Aperio Image Library Synthetic\r\n"
            f"{width}x{height} [0,0 {width}x{height}] (256x256) JPEG/RGB Q={quality}"
            f"|AppMag = 20|MPP = {mpp}") if aperio else "synthetic"
    level_sizes, stored_sizes, pages = [], [], []
    w, h = width, height
    for li in range(n_levels):
        if li:
            w, h = math.ceil(w / 2), math.ceil(h / 2)
        level = img[::2 ** li, ::2 ** li][:h, :w].copy() if li else img
        sw, sh = w, h
        if shrink_last and li == n_levels - 1 and li > 0:
            sw = w - shrink_last
            level = level[:, :sw].copy()
        elif pad_last and li == n_levels - 1 and li > 0:
            sw, sh = -(-w // 32) * 32, -(-h // 32) * 32   # pad like KFBio, with white
            padded = np.full((sh, sw, 3), 255, np.uint8)
            padded[:h, :w] = level
            level = padded
        level_sizes.append((w, h))
        stored_sizes.append((sw, sh))
        pages.append(level)
    kw = dict(compression="jpeg", compressionargs={"level": quality}, photometric="rgb",
              subsampling=(2, 1), metadata=None)
    with tifffile.TiffWriter(path, bigtiff=False) as t:
        t.write(pages[0], tile=(TS, TS), description=desc, **kw)
        t.write(img[::8, ::8].copy(), description="thumbnail", **kw)          # page 1
        for lv in pages[1:]:
            t.write(lv, tile=(TS, TS),
                    description=("Aperio Image Library Synthetic"
                                 f"|AppMag = 20|MPP = {mpp}") if aperio else "x", **kw)
        t.write(img[:64, :96].copy(), compression="lzw", photometric="rgb", subfiletype=1,
                description="label 96x64", metadata=None)
        t.write(img[:64, :192].copy(), subfiletype=9, description="macro 192x64", **kw)
    info = SvsInfo(path, width, height, mpp, level_sizes, stored_sizes, set(missing))
    page_of_level = [0] + [2 + i for i in range(n_levels - 1)]
    for li, col, row in missing:
        stride = math.ceil(stored_sizes[li][0] / TS)
        _zero_tile_bytecounts(path, page_of_level[li], [row * stride + col])
    return info


# --------------------------------------------------------------------------
# pytest fixtures
# --------------------------------------------------------------------------

@pytest.fixture
def kfbf(tmp_path) -> KfbfInfo:
    return make_kfbf(tmp_path / "fluo.kfbf")


@pytest.fixture
def svs(tmp_path) -> SvsInfo:
    return make_svs(tmp_path / "bright.svs")


@pytest.fixture
def mixed_dir(tmp_path) -> Path:
    d = tmp_path / "slides"
    d.mkdir()
    make_kfbf(d / "fluo.kfbf")
    make_svs(d / "bright.svs")
    return d


@pytest.fixture(autouse=True)
def _quiet_logger():
    """Tests may call main() many times; keep the module logger clean."""
    yield
    for h in list(kc.log.handlers):
        kc.log.removeHandler(h)
    kc.log.propagate = True   # main() turns it off; caplog needs it on
    kc.log.setLevel(0)


def tiff_levels(path: Path) -> dict:
    """Read structure of a written pyramid: level 0 + SubIFD levels, or for
    the Aperio layout level 0 + top-level pages 2.. (page 1 = thumbnail)."""
    with tifffile.TiffFile(path) as t:
        top = t.pages[0]
        aperio = not top.subifds and t.is_svs
        if aperio:
            subs = [t.pages[i] for i in range(2, len(t.pages))]
        else:
            subs = list(tifffile.TiffPages(top)) if top.subifds else []
        levels = []
        for pg in [top] + subs:
            xres = pg.tags["XResolution"].value
            levels.append(dict(
                shape=tuple(pg.shape), photometric=pg.photometric.name,
                compression=pg.compression.name, subfiletype=int(pg.subfiletype),
                planarconfig=int(pg.planarconfig), tile=(pg.tilelength, pg.tilewidth),
                px_per_cm=xres[0] / xres[1], unit=int(pg.tags["ResolutionUnit"].value),
                subsampling=(pg.tags["YCbCrSubSampling"].value
                             if "YCbCrSubSampling" in pg.tags else None)))
        return dict(n_pages=len(t.pages), is_ome=t.is_ome, ome_xml=t.ome_metadata,
                    aperio=aperio, description=top.description,
                    thumbnail_shape=tuple(t.pages[1].shape) if aperio else None,
                    levels=levels,
                    arrays=[lv.asarray() for lv in t.series[0].levels],
                    n_series_levels=len(t.series[0].levels))


# --------------------------------------------------------------------------
# .qptiff (PerkinElmer / Akoya) and plain .tif
# --------------------------------------------------------------------------

def qpi_description(image_type: str, name: str | None = None, color=None,
                    objective: str = "20x", biomarker: str | None = None) -> str:
    """Per-page XML like a Vectra Polaris / KFBio QPTIFF."""
    parts = ['<?xml version="1.0" encoding="utf-16"?>',
             "<PerkinElmer-QPI-ImageDescription>",
             "  <DescriptionVersion>2</DescriptionVersion>",
             "  <AcquisitionSoftware>Synthetic 1.0</AcquisitionSoftware>",
             f"  <ImageType>{image_type}</ImageType>",
             f"  <Objective>{objective}</Objective>"]
    if name is not None:
        parts.append(f"  <Name>{name}</Name>")
        if biomarker:
            parts.append(f"  <Biomarker>{biomarker}</Biomarker>")
        parts.append(f"  <Color>{','.join(map(str, color))}</Color>")
        parts += ["  <Responsivity>", "    <Filter>",           # nested Name: ignored
                  f"      <Name>filter-{name}</Name>", "    </Filter>", "  </Responsivity>"]
    parts.append("</PerkinElmer-QPI-ImageDescription>")
    return "\r\n".join(parts)


@dataclass
class QptiffInfo:
    path: Path
    width: int
    height: int
    channels: list[str]
    colors: list[tuple[int, int, int]]
    mpp: float
    level_sizes: list[tuple[int, int]]   # stored pyramid levels (w, h)
    planes: dict = field(default_factory=dict)  # (level, ch) -> stored array


QPTIFF_CHANNELS = (("DAPI", (0, 0, 255)), ("Opal 520", (0, 255, 0)),
                   ("Opal 690", (255, 0, 0)))


def make_qptiff(path: Path, width: int = 1300, height: int = 1000, channels=QPTIFF_CHANNELS,
                n_levels: int = 3, dtype=np.uint8, mpp: float = 0.4977,
                objective: str = "20x", bad_last: bool = False, rgb: bool = False,
                biomarkers: bool = False) -> QptiffInfo:
    """Write a synthetic QPTIFF: per channel a FullResolution page, an RGB
    Thumbnail page, ReducedResolution pages per level, then Overview and
    Label. 8-bit planes are JPEG tiles, 16-bit ones deflate. ``bad_last``
    stores the last level 40 rows short (not an exact 2x step)."""
    dtype = np.dtype(dtype)
    chans = [("RGB", (0, 0, 0))] if rgb else list(channels)
    info = QptiffInfo(path, width, height, [c for c, _ in chans], [c for _, c in chans],
                      mpp, [])
    res = 1e4 / mpp
    if rgb:
        comp = dict(compression="jpeg", compressionargs={"level": 90}, photometric="rgb")
    elif dtype == np.uint8:
        comp = dict(compression="jpeg", compressionargs={"level": 95}, photometric="minisblack")
    else:
        comp = dict(compression="zlib", photometric="minisblack")
    sizes, w, h = [], width, height
    for li in range(n_levels):
        sizes.append((w, h))
        w, h = math.ceil(w / 2), math.ceil(h / 2)
    if bad_last and n_levels > 1:
        sizes[-1] = (sizes[-1][0], sizes[-1][1] - 40)
    info.level_sizes = sizes
    with tifffile.TiffWriter(path, bigtiff=True) as t:
        def page(arr, desc, subfiletype=0, tiled=True, kw=comp, r=res):
            t.write(arr, tile=(TS, TS) if tiled else None, description=desc,
                    software="PerkinElmer-QPI", metadata=None, subfiletype=subfiletype,
                    resolution=(r, r), resolutionunit="CENTIMETER", **kw)

        for li, (lw, lh) in enumerate(sizes):
            for ci, (name, color) in enumerate(chans):
                if rgb:
                    arr = _rgb_pattern(lh, lw, seed=li)
                else:
                    arr = _pattern(lh, lw, seed=ci + 3 * li).astype(dtype)
                    if dtype == np.uint16:
                        arr = arr * np.uint16(200) + np.uint16(ci)
                info.planes[(li, ci)] = arr
                kind = "FullResolution" if li == 0 else "ReducedResolution"
                page(arr, qpi_description(kind, None if rgb else name, color, objective,
                                          f"CD{ci + 3}" if biomarkers else None),
                     subfiletype=0 if li == 0 else 1, r=res / 2 ** li)
            if li == 0:  # thumbnail sits after level 0, as in real files
                page(np.full((64, 80, 3), 200, np.uint8), qpi_description("Thumbnail"),
                     tiled=False, kw=dict(photometric="rgb"))
        page(np.full((40, 90, 3), 120, np.uint8), qpi_description("Overview"),
             tiled=False, kw=dict(photometric="rgb"))
        page(np.full((30, 50, 3), 60, np.uint8), qpi_description("Label"),
             tiled=False, kw=dict(photometric="rgb"))
    return info


def make_tiff(path: Path, width: int = 900, height: int = 700, rgb: bool = True,
              dtype=np.uint8, tiled: bool = False, pyramid: int = 1,
              mpp: float | None = 0.25, compression: str | None = "zlib") -> np.ndarray:
    """Write a plain (non-Aperio) TIFF; returns the full-resolution image.
    ``pyramid`` > 1 adds that many levels in total as SubIFDs."""
    img = _rgb_pattern(height, width) if rgb else _pattern(height, width, 1).astype(dtype)
    if not rgb and np.dtype(dtype) == np.uint16:
        img = img * np.uint16(250)
    kw = dict(photometric="rgb" if rgb else "minisblack", compression=compression,
              metadata=None)
    if mpp:
        kw.update(resolution=(1e4 / mpp, 1e4 / mpp), resolutionunit="CENTIMETER")
    with tifffile.TiffWriter(path, bigtiff=False) as t:
        t.write(img, tile=(TS, TS) if tiled else None, rowsperstrip=None if tiled else 16,
                subifds=pyramid - 1 if pyramid > 1 else None, **kw)
        for li in range(1, pyramid):
            t.write(img[::2 ** li, ::2 ** li], tile=(TS, TS) if tiled else None,
                    subfiletype=1, **kw)
    return img
