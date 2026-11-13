"""BigTIFF outputs in the Aperio page layout, so OpenSlide (and QuPath through
it) sees the whole pyramid, the objective magnification and the pixel size."""

import shutil

import pytest
import tifffile

import kfbf_converter as kc
from conftest import make_kfbf, make_svs, tiff_levels

try:
    import openslide
except Exception:  # pragma: no cover - optional test dependency
    openslide = None

needs_openslide = pytest.mark.skipif(openslide is None, reason="openslide-python not installed")


@needs_openslide
@pytest.mark.parametrize("which", ["svs-copy", "svs-deflate", "kfbf-rgb", "kfbf-split"])
def test_openslide_reads_pyramid_magnification_and_mpp(which, tmp_path):
    out = tmp_path / "out"
    if which.startswith("svs"):
        info = make_svs(tmp_path / "b.svs")
        comp = None if which == "svs-copy" else "deflate"
        (p,) = kc.convert_slide(info.path, out, compression=comp, show_progress=False)
        mag, mpp, n_levels = "20", 0.5, 3
    else:
        info = make_kfbf(tmp_path / "f.kfbf")
        mode = which.split("-")[1]
        written = kc.convert_slide(info.path, out, mode=mode, show_progress=False)
        p = written[0]
        mag, mpp, n_levels = "40", 0.2528, 2
    assert openslide.OpenSlide.detect_format(str(p)) == "aperio"
    with openslide.OpenSlide(str(p)) as s:
        assert s.level_count == n_levels
        assert s.properties["openslide.objective-power"] == mag
        assert float(s.properties["openslide.mpp-x"]) == pytest.approx(mpp)
        assert [round(d) for d in s.level_downsamples] == [2 ** i for i in range(n_levels)]
        assert "thumbnail" in s.associated_images
        low = s.read_region((0, 0), n_levels - 1, s.level_dimensions[-1])
        assert low.size == s.level_dimensions[-1]


@needs_openslide
def test_openslide_without_magnification(tmp_path):
    info = make_kfbf(tmp_path / "nomag.kfbf", magnification=0)
    (p,) = kc.convert_slide(info.path, tmp_path / "out", mode="rgb", show_progress=False)
    assert "AppMag" not in tiff_levels(p)["description"]
    with openslide.OpenSlide(str(p)) as s:
        assert "openslide.objective-power" not in s.properties
        assert s.level_count == 2


def test_magnification_override_in_description(tmp_path):
    info = make_svs(tmp_path / "b.svs")
    (p,) = kc.convert_slide(info.path, tmp_path / "out", magnification=40,
                            show_progress=False)
    assert "|AppMag = 40|MPP = 0.5" in tiff_levels(p)["description"]


def test_zstack_bigtiff_keeps_subifd_layout(tmp_path):
    info = make_kfbf(tmp_path / "z.kfbf", n_z=2, missing=())
    (p,) = kc.convert_slide(info.path, tmp_path / "out", mode="rgb", show_progress=False)
    got = tiff_levels(p)
    assert not got["aperio"] and got["n_pages"] == 2      # one page per Z plane


def test_own_output_is_not_a_slide(tmp_path):
    info = make_svs(tmp_path / "b.svs")
    (p,) = kc.convert_slide(info.path, tmp_path / "out", show_progress=False)
    renamed = tmp_path / "again.svs"
    shutil.copy(p, renamed)
    for path in (p, renamed):
        with pytest.raises(kc.SlideFormatError, match="already converted by kfbf_converter"):
            kc.open_slide(path)


def test_output_never_overwrites_input(tmp_path, monkeypatch):
    info = make_svs(tmp_path / "b.svs")
    monkeypatch.setattr(kc, "planned_outputs",
                        lambda slide, out_dir, mode: [("brightfield", info.path, -1)])
    with pytest.raises(kc.SlideFormatError, match="would overwrite the input"):
        kc.convert_slide(info.path, tmp_path, show_progress=False)
    assert info.path.stat().st_size > 0


def test_verify_detects_missing_thumbnail(tmp_path):
    info = make_svs(tmp_path / "b.svs", missing=())
    (p,) = kc.convert_slide(info.path, tmp_path / "out", show_progress=False)
    shapes = [(1000, 1400, 3), (500, 700, 3), (250, 350, 3)]
    kc.verify_output(p, shapes, n_pages=1, layout="aperio")
    with pytest.raises(kc.VerificationError, match="levels \\+ thumbnail"):
        kc.verify_output(p, shapes + [(125, 175, 3)], n_pages=1, layout="aperio")
    plain = tmp_path / "plain.tif"                  # same pages, no Aperio description
    with tifffile.TiffFile(p) as t:
        arrays = [pg.asarray() for pg in t.pages]
    with tifffile.TiffWriter(plain) as w:
        for i, a in enumerate(arrays):              # page 1 (thumbnail) stays stripped
            w.write(a, tile=None if i == 1 else (256, 256), photometric="rgb", metadata=None)
    with pytest.raises(kc.VerificationError, match="Aperio description is missing"):
        kc.verify_output(plain, shapes, n_pages=1, layout="aperio")
