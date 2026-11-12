"""Readers: .kfbf fluorescence, .svs brightfield, format detection, RGB JPEG decode."""

import shutil

import imagecodecs
import numpy as np
import pytest
import tifffile

import kfbf_converter as kc
from conftest import TS, make_kfbf, make_svs


# ---- .kfbf -----------------------------------------------------------------

def test_kfbf_header_and_levels(kfbf):
    with kc.KFBioSlide(kfbf.path) as s:
        assert s.kind == "fluorescence" and not s.is_rgb
        assert s.dimensions == (kfbf.width, kfbf.height)
        assert s.channel_names == kfbf.channels
        assert s.n_channels == 3 and s.n_z == 1
        assert s.mpp == pytest.approx(kfbf.mpp)
        assert [(lv.width, lv.height) for lv in s.levels[:s.n_source_levels]] == kfbf.levels
        assert all(lv.coord_mode == "pixel" for lv in s.levels)


def test_kfbf_tiles_regions_and_missing(kfbf):
    with kc.KFBioSlide(kfbf.path) as s:
        for (li, z, col, row, ch), ref in list(kfbf.tiles.items())[:20]:
            np.testing.assert_array_equal(s.get_tile(li, col, row, ch, z), ref)
        miss = s.get_tile(0, 1, 1, 0)                     # fixture leaves (0,0,1,1) out
        assert miss.shape == (TS, TS) and not miss.any()  # fluorescence: black
        ref = kfbf.images[(0, 0, 1)].copy()
        ref[TS:2 * TS, TS:2 * TS] = 0
        np.testing.assert_array_equal(s.read_region(0, 0, 0, kfbf.width, kfbf.height, 1), ref)
        # region across tile borders, partly outside the image (outside = 0)
        r = s.read_region(1, 300, 200, 100, 100, 2)
        lv1 = s.levels[1]
        assert r.shape == (100, 100)
        np.testing.assert_array_equal(r[:lv1.height - 200, :50],
                                      kfbf.images[(1, 0, 2)][200:, 300:350])
        assert not r[:, lv1.width - 300:].any()
        assert s.read_region_channels(2, 0, 0, 64, 32).shape == (3, 32, 64)


def test_kfbf_synthetic_levels_2d(tmp_path):
    """A slide whose stored pyramid stops above 1024 px gets 2x2-averaged levels."""
    info = make_kfbf(tmp_path / "big.kfbf", width=2300, height=300, channels=("DAPI",),
                     missing=())
    with kc.KFBioSlide(info.path) as s:
        assert s.n_source_levels == len(info.levels)
        # keep only the stored full-resolution level, as if the file had no pyramid
        s.levels = s._add_synthetic_levels(s.levels[:1])
        s.n_source_levels = 1
        assert s.levels[1].synthetic and s.levels[1].width == 1150
        full = s.read_region(0, 0, 0, 2300, 300, 0).astype(np.float32)
        exp = (full.reshape(150, 2, 1150, 2).mean(axis=(1, 3)) + 0.5).astype(np.uint8)
        np.testing.assert_array_equal(s.read_region(1, 0, 0, 1150, 150, 0), exp)


# ---- .svs ------------------------------------------------------------------

def test_svs_structure(svs):
    with kc.open_slide(svs.path) as s:
        assert isinstance(s, kc.SVSSlide)
        assert s.kind == "brightfield" and s.is_rgb and s.n_channels == 3
        assert s.dimensions == (svs.width, svs.height)
        assert s.mpp == pytest.approx(svs.mpp) and s.mpp_source == "description"
        assert s.magnification == 20
        assert [(lv.width, lv.height) for lv in s.levels] == svs.level_sizes
        assert s.stored_size(2) == svs.stored_sizes[2] != svs.level_sizes[2]  # padding cropped
        assert [lv.downsample for lv in s.levels] == [1, 2, 4]
        assert s.associated == ["Thumbnail", "Label", "Macro"]
        assert s.jpeg_subsampling(0) == (2, 1)
        assert s.missing_tiles() == 1
        text = s.summary()
        assert "brightfield RGB" in text and "padding cropped" in text


def test_svs_tiles_match_tifffile_and_missing_is_white(svs, tmp_path):
    # reference pixels from an identical file without the missing tile
    # (tifffile itself cannot read zero-bytecount JPEG tiles)
    full = make_svs(tmp_path / "full.svs", missing=())
    with tifffile.TiffFile(full.path) as t:
        ref0 = t.pages[0].asarray()
        ref1 = t.pages[2].asarray()
    with kc.open_slide(svs.path) as s:
        tile = s.get_tile(0, 2, 0)
        assert tile.shape == (TS, TS, 3) and tile.dtype == np.uint8
        np.testing.assert_array_equal(tile, ref0[:TS, 2 * TS:3 * TS])
        miss = s.get_tile(0, 1, 1)                # fixture zeroes (0, col 1, row 1)
        assert (miss == 255).all()
        region = s.read_region(1, 100, 50, 400, 300)
        assert region.shape == (300, 400, 3)
        np.testing.assert_array_equal(region, ref1[50:350, 100:500])
        whole = s.read_region(0, 0, 0, svs.width, svs.height)
        expect = ref0.copy()
        expect[TS:2 * TS, TS:2 * TS] = 255
        np.testing.assert_array_equal(whole, expect)
        out = s.read_region(2, 300, 200, 100, 100)  # past the true (cropped) edge
        assert (out[50:, :] == 255).all() and (out[:, 50:] == 255).all()


def test_svs_short_level_extended_with_white(tmp_path):
    """Level stored 2 px narrower than half its parent (257 -> 255) is declared
    at 257 px, which needs an extra tile column that reads as white."""
    info = make_svs(tmp_path / "short.svs", width=1028, height=600, missing=(),
                    pad_last=False, shrink_last=2)
    assert info.stored_sizes[2] == (255, 150)
    with kc.open_slide(info.path) as s:
        lv = s.levels[2]
        assert (lv.width, lv.height, lv.grid_cols) == (257, 150, 2)
        assert "extended with white" in s.summary()
        assert s.raw_tile(2, 1, 0) is None             # beyond the stored grid
        assert s.stored_bytes(2)[1] == 1                # counted as one blank tile
        region = s.read_region(2, 0, 0, 257, 150)
        assert (region[:, 255:] == 255).all()
        with tifffile.TiffFile(info.path) as t:
            np.testing.assert_array_equal(region[:, :255], t.pages[3].asarray())
    out = tmp_path / "out"
    (p,) = kc.convert_slide(info.path, out, show_progress=False)   # tile copy
    with tifffile.TiffFile(p) as t:
        assert len(t.series) == 1
        assert [lv.shape for lv in t.series[0].levels] == [(600, 1028, 3), (300, 514, 3),
                                                           (150, 257, 3)]
        last = t.series[0].levels[2].asarray()
    assert np.abs(last[:, :255].astype(int) - region[:, :255]).max() == 0
    assert (last[:, 256:] >= 254).all()   # added tile column: blank white JPEG
    # column 255 is the scanner's own JPEG padding, copied unchanged (KFBio
    # repeats the edge pixel there; this synthetic file pads with black)
    (p2,) = kc.convert_slide(info.path, tmp_path / "re", compression="deflate",
                             show_progress=False)
    with tifffile.TiffFile(p2) as t:
        assert (t.series[0].levels[2].asarray()[:, 255:] == 255).all()   # re-encode: white


def test_svs_synthetic_levels_rgb(tmp_path):
    info = make_svs(tmp_path / "one.svs", width=1100, height=300, n_levels=1, missing=())
    with kc.open_slide(info.path) as s:
        assert s.n_source_levels == 1
        assert [(lv.width, lv.height, lv.synthetic) for lv in s.levels] == \
            [(1100, 300, False), (550, 150, True)]
        full = s.read_region(0, 0, 0, 1100, 300).astype(np.float32)
        exp = (full.reshape(150, 2, 550, 2, 3).mean(axis=(1, 3)) + 0.5).astype(np.uint8)
        got = s.read_region(1, 0, 0, 550, 150)
        assert got.shape == (150, 550, 3)
        np.testing.assert_array_equal(got, exp)
        with pytest.raises(IndexError):
            s.get_tile(1, 0, 0)


def test_svs_odd_size_synthetic_edge_replication(tmp_path):
    info = make_svs(tmp_path / "odd.svs", width=1101, height=301, n_levels=1, missing=())
    with kc.open_slide(info.path) as s:
        lv1 = s.levels[1]
        assert (lv1.width, lv1.height) == (551, 151)
        full = s.read_region(0, 0, 0, 1101, 301)
        pad = np.pad(full, ((0, 1), (0, 1), (0, 0)), mode="edge").astype(np.float32)
        exp = (pad.reshape(151, 2, 551, 2, 3).mean(axis=(1, 3)) + 0.5).astype(np.uint8)
        np.testing.assert_array_equal(s.read_region(1, 0, 0, 551, 151), exp)


# ---- detection -------------------------------------------------------------

def test_open_slide_by_magic_then_extension(tmp_path, kfbf, svs, caplog):
    with kc.open_slide(kfbf.path) as s:
        assert isinstance(s, kc.KFBioSlide)
    renamed = tmp_path / "really_svs.kfbf"       # magic wins over extension
    shutil.copy(svs.path, renamed)
    with kc.open_slide(renamed) as s:
        assert isinstance(s, kc.SVSSlide)
    assert "file content is a brightfield slide" in caplog.text
    junk = tmp_path / "junk.dat"
    junk.write_bytes(b"hello world, not a slide")
    with pytest.raises(kc.SlideFormatError, match="unrecognised slide file"):
        kc.open_slide(junk)
    bad_kfbf = tmp_path / "broken.kfbf"            # extension known, magic wrong
    bad_kfbf.write_bytes(b"\0" * 200)
    with pytest.raises(kc.KFBFError, match="bad magic"):
        kc.open_slide(bad_kfbf)


def test_open_slide_rejects_plain_tiff(tmp_path):
    plain = tmp_path / "already_converted.tif"
    tifffile.imwrite(plain, np.zeros((64, 64, 3), np.uint8), photometric="rgb")
    with pytest.raises(kc.SlideFormatError, match="not an Aperio .svs"):
        kc.open_slide(plain)


def test_find_slides_suffixes(tmp_path):
    for name in ("a.kfbf", "b.KFBF", "c.svs", "d.SVS", "e.kfb", "f.tif", "g.svs.part",
                 "h.kfbfx", "i.svsx"):
        (tmp_path / name).write_bytes(b"x")
    (tmp_path / "sub").mkdir()
    (tmp_path / "sub" / "j.svs").write_bytes(b"x")
    names = [p.name for p in kc.find_slides(tmp_path)]
    assert names == ["a.kfbf", "b.KFBF", "c.svs", "d.SVS"]
    names = sorted(p.name for p in kc.find_slides(tmp_path, recursive=True))
    assert names == ["a.kfbf", "b.KFBF", "c.svs", "d.SVS", "j.svs"]


# ---- RGB JPEG decode -------------------------------------------------------

def _strip_adobe_and_number_components(jpeg: bytes) -> bytes:
    """Make an 'Aperio RGB' JPEG: RGB samples, no APP14, component IDs 1,2,3."""
    out, i = bytearray(jpeg[:2]), 2
    while i < len(jpeg):
        m, ln = jpeg[i + 1], int.from_bytes(jpeg[i + 2:i + 4], "big")
        seg = bytearray(jpeg[i:i + 2 + ln])
        if m == 0xDA:  # SOS header lists component IDs too
            for k in range(seg[4]):
                seg[5 + 2 * k] = k + 1
            out += seg + jpeg[i + 2 + ln:]
            break
        if m == 0xEE:
            i += 2 + ln
            continue
        if m in (0xC0, 0xC1):
            for k in range(seg[9]):
                seg[10 + 3 * k] = k + 1
        out += seg
        i += 2 + ln
    return bytes(out)


def _maxerr(a, b):
    """Max error away from the colour edge at row 32 (chroma subsampling blurs it)."""
    keep = np.r_[0:24, 40:64]
    return int(np.abs(a.astype(int)[keep] - b.astype(int)[keep]).max())


def test_decode_jpeg_rgb_paths():
    rgb = np.zeros((64, 64, 3), np.uint8)
    rgb[:32] = (200, 30, 40)
    rgb[32:] = (20, 160, 220)
    ycc = imagecodecs.jpeg_encode(rgb, level=98)
    assert _maxerr(kc.decode_jpeg_rgb(ycc), rgb) < 8
    assert kc.jpeg_sof_info(ycc)["jfif"] is True

    raw_rgb = imagecodecs.jpeg_encode(rgb, level=98, colorspace="RGB", outcolorspace="RGB")
    aperio = _strip_adobe_and_number_components(raw_rgb)
    info = kc.jpeg_sof_info(aperio)
    assert info["adobe_transform"] is None and [c[0] for c in info["components"]] == [1, 2, 3]
    wrong = imagecodecs.jpeg_decode(aperio)          # treated as YCbCr: wrong colours
    assert _maxerr(wrong, rgb) > 40
    right = kc.decode_jpeg_rgb(aperio, rgb_stored=True)
    assert _maxerr(right, rgb) < 8

    gray = imagecodecs.jpeg_encode(np.ascontiguousarray(rgb[..., 0]), level=98)
    assert kc.decode_jpeg_rgb(gray).shape == (64, 64, 3)


def test_decode_jpeg_rgb_pillow_fallback(monkeypatch):
    rgb = np.full((32, 32, 3), (180, 60, 120), np.uint8)
    enc = imagecodecs.jpeg_encode(rgb, level=98)
    monkeypatch.setattr(kc, "imagecodecs", None)
    assert np.abs(kc.decode_jpeg_rgb(enc).astype(int) - rgb).max() < 8
    # Pillow path for 'Aperio RGB' (no Adobe marker) undoes the wrong YCbCr step
    raw = imagecodecs.jpeg_encode(np.full((32, 32, 3), (200, 30, 40), np.uint8), level=98,
                                  colorspace="RGB", outcolorspace="RGB")
    got = kc.decode_jpeg_rgb(_strip_adobe_and_number_components(raw), rgb_stored=True)
    assert np.abs(got.astype(int) - (200, 30, 40)).max() < 10
