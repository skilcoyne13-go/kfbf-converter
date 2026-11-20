"""QPTIFF -> OME-TIFF, plain TIFF -> BigTIFF, the slide-type menu, never
re-reading our own outputs, damaged TIFFs and the end-of-run summary."""

import builtins
import dataclasses
import logging
import os
import re
import shutil
import xml.etree.ElementTree as ET

import numpy as np
import pytest
import tifffile

import kfbf_converter as kc
from conftest import make_kfbf, make_qptiff, make_svs, make_tiff, tiff_levels

OME = {"o": kc.OME_NS}


def _ome(path):
    with tifffile.TiffFile(path) as t:
        root = ET.fromstring(t.ome_metadata)
        series = t.series[0]
        return dict(root=root, n_pages=len(t.pages), axes=series.axes, dtype=series.dtype,
                    shapes=[tuple(lv.shape) for lv in series.levels],
                    compression=t.pages[0].compression.name)


def _channels(root):
    px = root.find(".//o:Pixels", OME)
    return [(c.get("Name"), int(c.get("Color"))) for c in px.findall("o:Channel", OME)]


# ---- QPTIFF reader -----------------------------------------------------------------

def test_qptiff_reader_metadata(tmp_path):
    q = make_qptiff(tmp_path / "q.qptiff")
    with kc.open_slide(q.path) as s:
        assert isinstance(s, kc.QPTIFFSlide) and s.input_type == "QPTIFF"
        assert s.kind == "fluorescence" and s.dtype == np.uint8
        assert s.channel_names == ["DAPI", "Opal 520", "Opal 690"]   # not the filter Names
        assert s.channel_colors == [(0, 0, 255), (0, 255, 0), (255, 0, 0)]
        assert s.mpp == pytest.approx(0.4977) and s.mpp_source == "resolution tags"
        assert s.magnification == 20 and s.magnification_source == "QPTIFF XML"
        assert [(lv.width, lv.height) for lv in s.levels] == q.level_sizes
        assert s.n_source_levels == 3
        assert s.associated == ["Thumbnail", "Overview", "Label"]   # skipped
        with tifffile.TiffFile(q.path) as t:      # level 1, channel 2 as stored (JPEG)
            ref = [pg for pg in t.pages if pg.imagewidth == 650][2].asarray()
        np.testing.assert_array_equal(s.read_region(1, 0, 0, 650, 500, 2), ref)


def test_qptiff_biomarker_wins_over_name(tmp_path):
    q = make_qptiff(tmp_path / "q.qptiff", biomarkers=True)
    with kc.open_slide(q.path) as s:
        assert s.channel_names == ["CD3", "CD4", "CD5"]


def test_qptiff_to_ome_copies_tiles_and_keeps_metadata(tmp_path):
    q = make_qptiff(tmp_path / "q.qptiff")
    (p,) = kc.convert_slide(q.path, tmp_path / "out", show_progress=False)
    assert p.name == "q.ome.tif"
    got = _ome(p)
    assert got["axes"] == "CYX" and got["dtype"] == np.uint8
    assert got["shapes"] == [(3, h, w) for w, h in q.level_sizes]
    assert got["n_pages"] == 3                         # no thumbnail / overview / label
    assert got["compression"] == "JPEG"
    assert _channels(got["root"]) == [("DAPI", kc.ome_color((0, 0, 255))),
                                      ("Opal 520", kc.ome_color((0, 255, 0))),
                                      ("Opal 690", kc.ome_color((255, 0, 0)))]
    px = got["root"].find(".//o:Pixels", OME)
    assert float(px.get("PhysicalSizeX")) == pytest.approx(0.4977)
    obj = got["root"].find(".//o:Instrument/o:Objective", OME)
    assert obj.get("NominalMagnification") == "20.0"
    # stored JPEG tiles are copied byte for byte
    with tifffile.TiffFile(q.path) as src, tifffile.TiffFile(p) as out:
        full = [pg for pg in src.pages if "FullResolution" in (pg.description or "")]
        for c in range(3):
            sp, b = full[c], out.pages[c]
            src.filehandle.seek(sp.dataoffsets[1])
            out.filehandle.seek(b.dataoffsets[1])
            assert src.filehandle.read(sp.databytecounts[1]) == \
                out.filehandle.read(b.databytecounts[1])


def test_qptiff_16bit_keeps_bit_depth(tmp_path):
    q = make_qptiff(tmp_path / "q16.qptiff", dtype=np.uint16, n_levels=2)
    (p,) = kc.convert_slide(q.path, tmp_path / "out", show_progress=False)
    with tifffile.TiffFile(p) as t:
        lv0 = t.series[0].levels[0].asarray()
    assert lv0.dtype == np.uint16 and lv0.max() > 255
    for c in range(3):
        np.testing.assert_array_equal(lv0[c], q.planes[(0, c)])
    with pytest.raises(kc.SlideFormatError, match="16-bit pixels, which JPEG cannot store"):
        kc.convert_slide(q.path, tmp_path / "out2", compression="jpeg", show_progress=False)


def test_qptiff_inexact_level_is_computed(tmp_path, caplog):
    q = make_qptiff(tmp_path / "q.qptiff", width=2600, height=2000, bad_last=True)
    caplog.set_level(logging.INFO, logger="kfbf_converter")
    with kc.open_slide(q.path) as s:
        assert s.n_source_levels == 2      # stored 650x460 is not half of 1300x1000
        assert [(lv.width, lv.height, lv.synthetic) for lv in s.levels] == \
            [(2600, 2000, False), (1300, 1000, False), (650, 500, True)]
        parent = s.read_region(1, 0, 0, 1300, 1000, 0).astype(np.float32)
        expect = (parent.reshape(500, 2, 650, 2).mean(axis=(1, 3)) + 0.5).astype(np.uint8)
        np.testing.assert_array_equal(s.read_region(2, 0, 0, 650, 500, 0), expect)
    assert "not an exact 2x step" in caplog.text
    (p,) = kc.convert_slide(q.path, tmp_path / "out", show_progress=False)
    assert _ome(p)["shapes"][-1] == (3, 500, 650)


def test_qptiff_identified_by_content_any_extension(tmp_path):
    q = make_qptiff(tmp_path / "scan.tif")
    assert kc.input_type_of(q.path) == "QPTIFF"
    with kc.open_slide(q.path) as s:
        assert isinstance(s, kc.QPTIFFSlide)
    (p,) = kc.convert_slide(q.path, tmp_path / "out", show_progress=False)
    assert p.name == "scan.ome.tif"


def test_qptiff_brightfield_gets_rgb_bigtiff(tmp_path, caplog):
    q = make_qptiff(tmp_path / "bf.qptiff", rgb=True, n_levels=2)
    caplog.set_level(logging.INFO, logger="kfbf_converter")
    (p,) = kc.convert_slide(q.path, tmp_path / "out", show_progress=False)
    assert p.name == "bf.tif"
    info = tiff_levels(p)
    assert info["aperio"] and info["levels"][0]["shape"] == (1000, 1300, 3)
    assert "|AppMag = 20|" in info["description"]
    assert "this QPTIFF is brightfield (RGB); writing an RGB BigTIFF" in caplog.text


def test_qptiff_messages_never_say_kfbf(tmp_path, capsys):
    q = make_qptiff(tmp_path / "q.qptiff")
    assert kc.main([str(q.path), "-o", str(tmp_path / "out"), "-q", "-v"]) == 0
    err = capsys.readouterr().err
    assert "QPTIFF fluorescence" in err
    assert "kfbf" not in err.replace("kfbf_converter", "").lower()


def _too_big_lines(err):
    """The reason and Fix lines printed right under the headline."""
    lines = err.splitlines()
    assert lines.count(kc.DataTooBigError.HEADLINE) == 1, err
    i = lines.index(kc.DataTooBigError.HEADLINE)
    return lines[i + 1], lines[i + 2]


def test_qptiff_preflight_too_big(tmp_path, monkeypatch, capsys):
    q = make_qptiff(tmp_path / "q.qptiff")
    out = tmp_path / "out"
    out.mkdir()
    monkeypatch.setattr(kc, "filesystem_type", lambda p: "NTFS")
    monkeypatch.setattr(kc.shutil, "disk_usage",
                        lambda p: shutil._ntuple_diskusage(10 ** 12, 10 ** 12, 1000))
    assert kc.main([str(q.path), "-o", str(out), "-q"]) == 1
    err = capsys.readouterr().err
    reason, fix = _too_big_lines(err)
    vol = re.escape(str(kc._volume_root(out)))
    gb = r"\d+\.\d\d GB"
    assert re.fullmatch(rf"  q\.ome\.tif: needs about {gb} \(plus a {gb} safety margin\), "
                        rf"but only 0\.00 GB is free on {vol}\.", reason), reason
    assert re.fullmatch(rf"  Fix: free up {gb} on {vol}, or choose another drive with -o\.",
                        fix), fix
    assert "kfbf" not in err.replace("kfbf_converter", "").lower()
    assert "Traceback" not in err
    assert os.listdir(out) == []                          # no output, no .part


def test_tiff_preflight_too_big(tmp_path, monkeypatch, capsys):
    src_dir = tmp_path / "in"
    src_dir.mkdir()
    make_tiff(src_dir / "photo.tif")
    before = (src_dir / "photo.tif").read_bytes()
    out = tmp_path / "out"
    out.mkdir()
    monkeypatch.setattr(kc, "filesystem_type", lambda p: "FAT32")
    monkeypatch.setattr(kc.shutil, "disk_usage",
                        lambda p: shutil._ntuple_diskusage(10 ** 15, 0, 10 ** 15))
    real = kc.estimate_output

    def fake(slide, kind, writer, copy_tiles=False):
        est = real(slide, kind, writer, copy_tiles)
        size = int(5.2 * kc.GiB)
        return dataclasses.replace(est, estimate=size, uncompressed=max(size, est.uncompressed),
                                   exact=False)
    monkeypatch.setattr(kc, "estimate_output", fake)
    assert kc.main([str(src_dir / "photo.tif"), "-o", str(out), "-q"]) == 1
    err = capsys.readouterr().err
    reason, fix = _too_big_lines(err)
    assert reason == (f"  photo.tif: would be about 5.20 GB, but the output drive "
                      f"{kc._volume_root(out)} is FAT32, which cannot store files larger "
                      f"than 4.00 GB.")
    assert fix.startswith("  Fix: choose an NTFS or exFAT drive with -o")
    assert "kfbf" not in err.replace("kfbf_converter", "").lower()
    assert "Traceback" not in err
    assert os.listdir(out) == []                          # no output, no .part
    assert (src_dir / "photo.tif").read_bytes() == before  # input untouched


# ---- plain TIFF -------------------------------------------------------------------

def test_tiff_rgb_to_bigtiff_never_overwrites_input(tmp_path):
    img = make_tiff(tmp_path / "photo.tif")              # stripped, deflate, 0.25 um/px
    before = (tmp_path / "photo.tif").read_bytes()
    (p,) = kc.convert_slide(tmp_path / "photo.tif", None, show_progress=False)
    assert p == tmp_path / "photo.bigtiff.tif"            # same name as the input -> .bigtiff
    assert (tmp_path / "photo.tif").read_bytes() == before
    info = tiff_levels(p)
    assert info["aperio"] and info["levels"][0]["photometric"] == "RGB"
    assert info["levels"][0]["px_per_cm"] == pytest.approx(1e4 / 0.25)
    np.testing.assert_array_equal(info["arrays"][0], img)


def test_tiff_extension_tiff_gets_plain_name(tmp_path):
    make_tiff(tmp_path / "photo.tiff")
    (p,) = kc.convert_slide(tmp_path / "photo.tiff", None, show_progress=False)
    assert p.name == "photo.tif"


@pytest.mark.parametrize("neighbour", ["svs", "foreign"])
def test_tiff_output_name_clash_uses_bigtiff(tmp_path, neighbour):
    make_tiff(tmp_path / "x.tiff")
    if neighbour == "svs":
        make_svs(tmp_path / "x.svs")                      # x.svs will write x.tif
    else:
        make_tiff(tmp_path / "x.tif")                     # someone else's x.tif
    (p,) = kc.convert_slide(tmp_path / "x.tiff", None, show_progress=False)
    assert p.name == "x.bigtiff.tif"


def test_tiff_grayscale_16bit(tmp_path):
    img = make_tiff(tmp_path / "g.tiff", rgb=False, dtype=np.uint16, tiled=True)
    (p,) = kc.convert_slide(tmp_path / "g.tiff", tmp_path / "out", show_progress=False)
    info = tiff_levels(p)
    assert info["aperio"] and info["arrays"][0].dtype == np.uint16
    np.testing.assert_array_equal(info["arrays"][0], img)


def test_tiff_already_pyramidal_is_rewritten(tmp_path, caplog):
    img = make_tiff(tmp_path / "pyr.tiff", width=1400, height=1200, tiled=True, pyramid=3)
    caplog.set_level(logging.INFO, logger="kfbf_converter")
    with kc.open_slide(tmp_path / "pyr.tiff") as s:
        assert s.is_pyramidal and s.n_source_levels == 3
    (p,) = kc.convert_slide(tmp_path / "pyr.tiff", tmp_path / "out", show_progress=False)
    assert "this TIFF is already pyramidal" in caplog.text
    info = tiff_levels(p)
    assert [lv["shape"] for lv in info["levels"]] == [(1200, 1400, 3), (600, 700, 3),
                                                      (300, 350, 3)]
    np.testing.assert_array_equal(info["arrays"][0], img)


def test_tiff_without_pixel_size_uses_fallback(tmp_path, caplog):
    make_tiff(tmp_path / "nores.tiff", mpp=None)
    with kc.open_slide(tmp_path / "nores.tiff") as s:
        assert s.mpp_source == "fallback" and s.mpp == kc.DEFAULT_MPP
    assert "no pixel size in file" in caplog.text


def _damage(tmp_path):
    make_tiff(tmp_path / "truncated.tif", tiled=True)
    data = (tmp_path / "truncated.tif").read_bytes()
    (tmp_path / "truncated.tif").write_bytes(data[: len(data) // 2])   # cut in half
    (tmp_path / "garbage.tiff").write_bytes(b"II*\0" + os.urandom(3000))  # header, then junk
    return ["truncated.tif", "garbage.tiff"]


def test_damaged_tiff_is_reported(tmp_path):
    for name in _damage(tmp_path):
        with pytest.raises(kc.SlideFormatError,
                           match=f"{name}: this TIFF is damaged or incomplete and "
                                 f"cannot be converted"):
            kc.open_slide(tmp_path / name)


def test_damaged_tiff_does_not_stop_the_batch(tmp_path, capsys):
    _damage(tmp_path)
    make_tiff(tmp_path / "good.tiff")
    make_svs(tmp_path / "slide.svs")
    out = tmp_path / "out"
    assert kc.main([str(tmp_path), "-o", str(out), "-q"]) == 1
    assert sorted(os.listdir(out)) == ["good.tif", "slide.tif"]
    err = capsys.readouterr().err
    assert "truncated.tif: this TIFF is damaged or incomplete and cannot be converted" in err
    assert "garbage.tiff: this TIFF is damaged or incomplete and cannot be converted" in err
    assert "Converted 2, skipped 0 (already done), failed 2" in err
    assert "TIFF  : converted 1, skipped 0 (already done), failed 2" in err
    assert "SVS   : converted 1, skipped 0 (already done), failed 0" in err


def test_summary_counts_skipped(tmp_path, capsys):
    make_tiff(tmp_path / "a.tiff")
    assert kc.main([str(tmp_path), "-q"]) == 0
    assert kc.main([str(tmp_path), "-q"]) == 0       # a.tif is ours now: ignored as input
    err = capsys.readouterr().err
    assert "Converted 1, skipped 0 (already done), failed 0" in err
    assert "Converted 0, skipped 1 (already done), failed 0" in err   # 2nd run


# ---- never take our own outputs as inputs -----------------------------------------------

def _own_output(path, src_kind="svs"):
    """Write a real output of this tool to ``path`` (then rename it there)."""
    tmp = path.parent / "_src"
    tmp.mkdir(exist_ok=True)
    src = make_svs(tmp / "o.svs").path
    (p,) = kc.convert_slide(src, tmp, show_progress=False)
    shutil.move(p, path)
    shutil.rmtree(tmp)


@pytest.mark.parametrize("name,setup", [
    ("s.ome.tif", None),
    ("s.rgb.tif", None),
    ("s.bigtiff.tif", None),
    ("s.tif.part", None),
    ("k_DAPI.tif", "k.kfbf"),          # split channel next to its .kfbf
    ("k_Opal_520.tif", "k.kfbf"),      # channel names may contain '_'
    ("b.tif", "b.svs"),                # brightfield output next to its .svs
])
def test_skip_rules_by_name(tmp_path, name, setup):
    make_tiff(tmp_path / name)
    if setup:
        (make_kfbf if setup.endswith(".kfbf") else make_svs)(tmp_path / setup)
    assert kc.own_output_reason(tmp_path / name)
    found = [p.name for p, _ in kc.discover_slides(tmp_path)]
    assert name not in found and (setup is None or setup in found)


def test_skip_rule_software_tag_catches_renamed_output(tmp_path):
    _own_output(tmp_path / "renamed.tiff")
    assert "Software tag" in kc.own_output_reason(tmp_path / "renamed.tiff")
    assert kc.discover_slides(tmp_path) == []


def test_user_tiffs_with_similar_names_are_inputs(tmp_path):
    make_tiff(tmp_path / "liver_HE.tif")       # no liver.kfbf: a real input
    make_tiff(tmp_path / "b.tiff")             # .tiff next to b.svs: not our name
    make_svs(tmp_path / "b.svs")
    found = sorted((p.name, t) for p, t in kc.discover_slides(tmp_path))
    assert found == [("b.svs", "SVS"), ("b.tiff", "TIFF"), ("liver_HE.tif", "TIFF")]


# ---- menu ------------------------------------------------------------------------

def _run(monkeypatch, argv, answers, console=False):
    """main() with scripted answers; ``console`` = a person at the keyboard."""
    monkeypatch.setattr(kc, "stdin_is_console", lambda: console)
    it = iter(answers)
    asked = []

    def fake_input(prompt=""):
        asked.append(prompt)
        print(prompt, end="")
        try:
            return next(it)
        except StopIteration:
            raise EOFError
    monkeypatch.setattr(builtins, "input", fake_input)
    return kc.main(argv), asked


@pytest.fixture
def three_types(tmp_path):
    d = tmp_path / "slides"
    d.mkdir()
    make_tiff(d / "t1.tiff")
    make_tiff(d / "t2.tiff", rgb=False)
    make_svs(d / "b.svs")
    make_qptiff(d / "q.qptiff")
    return d


def test_type_menu_lists_found_types_and_enter_converts_all(three_types, monkeypatch,
                                                           capsys):
    rc, asked = _run(monkeypatch, [str(three_types), "-q"], [""], console=True)
    assert rc == 0
    out = capsys.readouterr().out
    assert out.startswith(kc.MENU_BANNER + "\n"
                          "Found 2 TIFF, 1 SVS, 1 QPTIFF slide(s).\n"
                          "What do you want to convert?\n"
                          "  [1] Convert all (Default - press Enter)\n"
                          "  [2] Only TIFF   -> BigTIFF\n"
                          "  [3] Only SVS    -> BigTIFF\n"
                          "  [4] Only QPTIFF -> OME-TIFF\n"
                          "Enter choice [1-4] (Press Enter for 1): ")
    assert "KFBF" not in out and len(asked) == 1        # no .kfbf: no format menu
    assert sorted(p.name for p in three_types.glob("*.tif")) == \
        ["b.tif", "q.ome.tif", "t1.tif", "t2.tif"]


def test_type_menu_only_one_type(three_types, monkeypatch, capsys):
    rc, _ = _run(monkeypatch, [str(three_types), "-q"], ["4"], console=True)
    assert rc == 0
    assert sorted(p.name for p in three_types.glob("*.tif")) == ["q.ome.tif"]


def test_type_menu_invalid_then_cancel(three_types, monkeypatch, capsys):
    rc, _ = _run(monkeypatch, [str(three_types), "-q"], ["9"], console=True)
    out = capsys.readouterr().out
    assert rc == 130 and "Invalid choice, please enter a number from 1 to 4." in out
    assert list(three_types.glob("*.tif")) == []


def test_type_menu_with_kfbf_then_format_menu(tmp_path, monkeypatch, capsys):
    make_kfbf(tmp_path / "f.kfbf")
    make_qptiff(tmp_path / "q.qptiff")
    rc, _ = _run(monkeypatch, [str(tmp_path), "-q"], ["", "3"], console=True)
    out = capsys.readouterr().out
    assert rc == 0
    assert "  [2] Only QPTIFF -> OME-TIFF\n  [3] Only KFBF   -> choose the format next\n" in out
    assert "Output format for the 1 KFBF slide(s):\nSelect output format:" in out
    assert sorted(p.name for p in tmp_path.glob("*.tif")) == \
        ["f_DAPI.tif", "f_FITC.tif", "f_TRITC.tif", "q.ome.tif"]


def test_nobody_at_keyboard_converts_everything_without_asking(three_types, monkeypatch):
    rc, asked = _run(monkeypatch, [str(three_types), "-q"], [], console=False)
    assert rc == 0 and asked == []
    assert len(list(three_types.glob("*.tif"))) == 4


def test_yes_and_types_skip_the_menu(three_types, monkeypatch, capsys):
    rc, asked = _run(monkeypatch, [str(three_types), "-q", "--types", "svs", "tif"], [],
                     console=True)
    assert rc == 0 and asked == []
    assert sorted(p.name for p in three_types.glob("*.tif")) == ["b.tif", "t1.tif", "t2.tif"]
    rc, asked = _run(monkeypatch, [str(three_types), "-q", "--yes", "--overwrite"], [],
                     console=True)
    assert rc == 0 and asked == []
    assert "q.ome.tif" in [p.name for p in three_types.glob("*.tif")]


def test_single_type_folder_skips_menu_but_confirms(tmp_path, monkeypatch, capsys):
    make_qptiff(tmp_path / "q.qptiff")
    monkeypatch.chdir(tmp_path)
    rc, asked = _run(monkeypatch, ["-q"], [""])
    out = capsys.readouterr().out
    assert rc == 0 and len(asked) == 1
    assert out.startswith(kc.MENU_BANNER + "\n"
                          "Found 1 QPTIFF slide(s).\nEach one will be converted to a "
                          "multi-channel pyramidal OME-TIFF (<name>.ome.tif).\n")
    assert "What do you want to convert?" not in out
    assert (tmp_path / "q.ome.tif").exists()


def test_menu_ignores_bom_and_ctrl_c(three_types, monkeypatch, capsys):
    rc, _ = _run(monkeypatch, [str(three_types), "-q"], ["﻿3"], console=True)
    assert rc == 0 and [p.name for p in three_types.glob("*.tif")] == ["b.tif"]

    def interrupt(prompt=""):
        raise KeyboardInterrupt
    monkeypatch.setattr(builtins, "input", interrupt)
    assert kc.main([str(three_types), "-q"]) == 130
