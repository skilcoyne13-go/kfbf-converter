"""Writers, output modes, mixed-folder routing, CLI and interactive menus."""

import builtins
import os
import time
from pathlib import Path

import numpy as np
import pytest
import tifffile

import kfbf_converter as kc
from conftest import make_kfbf, make_svs, tiff_levels


def _out(tmp_path):
    d = tmp_path / "out"
    d.mkdir(exist_ok=True)
    return d


# ---- fluorescence ----------------------------------------------------------

@pytest.mark.parametrize("mode,names", [
    ("ome", ["fluo.ome.tif"]),
    ("rgb", ["fluo.rgb.tif"]),
    ("split", ["fluo_DAPI.tif", "fluo_FITC.tif", "fluo_TRITC.tif"]),
    ("all", ["fluo.ome.tif", "fluo.rgb.tif", "fluo_DAPI.tif", "fluo_FITC.tif",
             "fluo_TRITC.tif"]),
])
def test_fluorescence_modes(kfbf, tmp_path, mode, names):
    out = _out(tmp_path)
    written = kc.convert_slide(kfbf.path, out, mode=mode, show_progress=False)
    assert sorted(p.name for p in written) == sorted(names)
    assert sorted(os.listdir(out)) == sorted(names)        # no .part left over
    with kc.KFBioSlide(kfbf.path) as s:
        export = kc.select_export_levels(s.levels)
        assert [(lv.width, lv.height) for lv in export] == [(700, 520), (350, 260)]
        for p in written:
            info = tiff_levels(p)
            samples = (3,) if p.name.endswith(".rgb.tif") else ()
            assert len(info["levels"]) == len(export)
            for lv, got in zip(export, info["levels"]):
                assert got["shape"] == (lv.height, lv.width) + samples
                assert got["px_per_cm"] == pytest.approx(1e4 / (s.mpp * lv.downsample),
                                                         rel=1e-6)
                assert got["unit"] == 3   # centimetre
                assert got["tile"] == (512, 512)
            if p.name.endswith(".ome.tif"):
                assert [lv["subfiletype"] for lv in info["levels"]] == [0, 1]
                assert not info["aperio"]
            else:  # Aperio layout: levels are top-level pages, thumbnail is page 1
                assert info["aperio"] and info["n_pages"] == len(export) + 1
                assert [lv["subfiletype"] for lv in info["levels"]] == [0, 0]
                assert "|AppMag = 40|MPP = 0.2528" in info["description"]
            if p.name.endswith(".ome.tif"):
                assert info["is_ome"] and info["n_pages"] == 3
                xml = info["ome_xml"]
                for name in kfbf.channels:
                    assert f'Name="{name}"' in xml
                assert 'PhysicalSizeX="0.2528"' in xml and 'Color="65535"' in xml
                np.testing.assert_array_equal(info["arrays"][0][1],
                                              s.read_region(0, 0, 0, 700, 520, 1))
            elif "_" in p.stem:
                ch = kfbf.channels.index(p.stem.split("_")[1])
                np.testing.assert_array_equal(info["arrays"][0],
                                              s.read_region(0, 0, 0, 700, 520, ch))
            else:
                assert info["levels"][0]["photometric"] == "RGB"


def test_fluorescence_z_stack_ome(tmp_path):
    info = make_kfbf(tmp_path / "z.kfbf", n_z=2, missing=())
    out = _out(tmp_path)
    (p,) = kc.convert_slide(info.path, out, mode="ome", show_progress=False)
    got = tiff_levels(p)
    assert got["n_pages"] == 6 and got["arrays"][0].shape == (2, 3, 520, 700)
    assert 'SizeZ="2"' in got["ome_xml"]


def test_convert_slide_old_signature_still_works(kfbf, tmp_path):
    written = kc.convert_slide(kfbf.path, _out(tmp_path), mode="ome", compression="deflate",
                               tile_size=512, show_progress=False)
    assert [p.name for p in written] == ["fluo.ome.tif"]


# ---- brightfield -------------------------------------------------------------

def _check_bf_geometry(info, svs):
    assert info["aperio"] and info["n_pages"] == 4      # 3 levels + thumbnail
    assert info["n_series_levels"] == 3
    assert info["thumbnail_shape"] == (500, 700, 3)      # largest level <= 1024 px
    assert info["description"].startswith(kc.APERIO_HEADER)
    assert "|AppMag = 20|MPP = 0.5" in info["description"]
    for (w, h), got, ds in zip(svs.level_sizes, info["levels"], (1, 2, 4)):
        assert got["shape"] == (h, w, 3)
        assert got["planarconfig"] == 1    # contig
        assert got["px_per_cm"] == pytest.approx(1e4 / (svs.mpp * ds), rel=1e-6)
        assert got["unit"] == 3
    assert [lv["subfiletype"] for lv in info["levels"]] == [0, 0, 0]


def test_brightfield_default_copies_tiles(svs, tmp_path):
    out = _out(tmp_path)
    written = kc.convert_slide(svs.path, out, show_progress=False)
    assert [p.name for p in written] == ["bright.tif"]
    assert os.listdir(out) == ["bright.tif"]
    info = tiff_levels(written[0])
    _check_bf_geometry(info, svs)
    for lv in info["levels"]:
        assert (lv["photometric"], lv["compression"]) == ("YCBCR", "JPEG")
        assert lv["subsampling"] == (2, 1) and lv["tile"] == (256, 256)
    with kc.open_slide(svs.path) as s:
        for li, arr in enumerate(info["arrays"]):
            ref = s.read_region(li, 0, 0, arr.shape[1], arr.shape[0])
            diff = np.abs(arr.astype(int) - ref)
            assert diff.max() <= 1        # identical; blank white tile may be 254/255
        # stored bytes are copied unchanged
        with tifffile.TiffFile(written[0]) as t:
            pg = t.pages[0]
            t.filehandle.seek(pg.dataoffsets[0])
            assert t.filehandle.read(pg.databytecounts[0]) == bytes(s.raw_tile(0, 0, 0))


@pytest.mark.parametrize("compression", ["deflate", "lzw", "jpeg"])
def test_brightfield_reencode(svs, tmp_path, compression):
    out = _out(tmp_path)
    (p,) = kc.convert_slide(svs.path, out, compression=compression, show_progress=False)
    info = tiff_levels(p)
    _check_bf_geometry(info, svs)
    lv0 = info["levels"][0]
    assert lv0["tile"] == (512, 512)
    with kc.open_slide(svs.path) as s:
        ref = s.read_region(0, 0, 0, svs.width, svs.height)
    if compression == "jpeg":
        assert lv0["photometric"] == "YCBCR"
        assert np.abs(info["arrays"][0].astype(int) - ref).mean() < 3
    else:
        assert lv0["photometric"] == "RGB" and lv0["compression"] != "JPEG"
        np.testing.assert_array_equal(info["arrays"][0], ref)   # missing tile -> white
        assert (info["arrays"][0][256:512, 256:512] == 255).all()


def test_brightfield_tile_size_forces_reencode(svs, tmp_path):
    (p,) = kc.convert_slide(svs.path, _out(tmp_path), tile_size=256, show_progress=False)
    lv0 = tiff_levels(p)["levels"][0]
    assert lv0["tile"] == (256, 256) and lv0["photometric"] == "RGB"


# ---- routing / CLI ---------------------------------------------------------------

def test_mixed_folder_routing(mixed_dir, tmp_path, capsys):
    out = _out(tmp_path)
    assert kc.main([str(mixed_dir), "-o", str(out), "-q"]) == 0
    assert sorted(os.listdir(out)) == ["bright.tif", "fluo.ome.tif"]
    assert kc.main([str(mixed_dir), "-o", str(out), "-q", "--mode", "all"]) == 0
    assert sorted(os.listdir(out)) == ["bright.tif", "fluo.ome.tif", "fluo.rgb.tif",
                                       "fluo_DAPI.tif", "fluo_FITC.tif", "fluo_TRITC.tif"]
    err = capsys.readouterr().err
    assert err.count("--mode all applies to fluorescence (.kfbf) slides only") == 1
    assert "--mode ome applies" not in err   # default mode is not 'ignored'


def test_skip_existing_and_overwrite(mixed_dir, tmp_path, capsys):
    out = _out(tmp_path)
    assert kc.main([str(mixed_dir), "-o", str(out), "-q"]) == 0
    stamps = {p.name: p.stat().st_mtime_ns for p in out.iterdir()}
    capsys.readouterr()
    assert kc.main([str(mixed_dir), "-o", str(out), "-q"]) == 0
    err = capsys.readouterr().err
    assert err.count("skip existing") == 2 and "0 file(s) written" in err
    assert stamps == {p.name: p.stat().st_mtime_ns for p in out.iterdir()}
    time.sleep(0.01)
    assert kc.main([str(mixed_dir), "-o", str(out), "-q", "--overwrite"]) == 0
    assert "2 file(s) written" in capsys.readouterr().err
    assert all(p.stat().st_mtime_ns != stamps[p.name] for p in out.iterdir())


def test_info_both_types(mixed_dir, capsys):
    assert kc.main([str(mixed_dir), "--info"]) == 0
    text = capsys.readouterr().out
    assert "Channels (3) : DAPI" in text
    assert "brightfield RGB (Aperio SVS, classic TIFF)" in text
    assert "MPP           : 0.5 um/px" in text


def test_version_output(capsys):
    with pytest.raises(SystemExit) as ei:
        kc.main(["--version"])
    assert ei.value.code == 0
    assert capsys.readouterr().out.strip() == (
        f"kfbf_converter {kc.__version__} - {kc.__author__}, {kc.__company__}")


def test_bat_window_title_and_crlf():
    bat = (Path(kc.__file__).parent / "run_converter.bat").read_bytes()
    assert f"\r\ntitle KFBio Slide Converter - {kc.__company__}\r\n".encode() in bat
    assert bat.count(b"\n") == bat.count(b"\r\n")          # CRLF only


def test_artist_and_software_tags_on_every_output(kfbf, svs, tmp_path):
    from conftest import make_qptiff, make_tiff
    out = _out(tmp_path)
    q = make_qptiff(tmp_path / "q.qptiff")
    make_tiff(tmp_path / "photo.tif")
    written = (kc.convert_slide(kfbf.path, out, mode="all", show_progress=False)
               + kc.convert_slide(svs.path, out, show_progress=False)
               + kc.convert_slide(svs.path, out / "re", compression="deflate",
                                  show_progress=False)
               + kc.convert_slide(q.path, out, show_progress=False)
               + kc.convert_slide(tmp_path / "photo.tif", out, show_progress=False))
    assert len(written) == 9
    for p in written:
        with tifffile.TiffFile(p) as t:
            tags = t.pages[0].tags
            assert tags["Artist"].value == f"{kc.__author__}, {kc.__company__}", p.name
            assert tags["Software"].value.startswith("kfbf_converter "), p.name
        assert kc.own_output_reason(p) is not None, p.name     # still skipped as ours


def test_unknown_file_is_a_clear_error(tmp_path, capsys):
    bad = tmp_path / "x.dat"
    bad.write_bytes(b"nope" * 10)
    assert kc.main([str(bad), "-q"]) == 1
    assert "unrecognised slide file" in capsys.readouterr().err


MENU = ("Select output format:\n"
        "  [1] OME-TIFF (Default / QuPath)\n"
        "  [2] Composite RGB BigTIFF\n"
        "  [3] Split Grayscale Channels\n"
        "  [4] All of the above\n")


def _interactive(monkeypatch, cwd, answers):
    monkeypatch.chdir(cwd)
    it = iter(answers)

    def fake_input(prompt=""):
        print(prompt, end="")
        try:
            return next(it)
        except StopIteration:
            raise EOFError
    monkeypatch.setattr(builtins, "input", fake_input)
    return kc.main(["-q"])


def test_menu_fluorescence_only_unchanged(tmp_path, monkeypatch, capsys):
    make_kfbf(tmp_path / "a.kfbf")
    assert _interactive(monkeypatch, tmp_path, ["1"]) == 0
    out = capsys.readouterr().out
    assert out.startswith(kc.MENU_BANNER + "\nFound 1 slide(s).\n" + MENU)
    assert (tmp_path / "a.ome.tif").exists()


def test_menu_mixed_counts(mixed_dir, monkeypatch, capsys):
    # several types: first "what to convert" (Enter = all), then the .kfbf format
    assert _interactive(monkeypatch, mixed_dir, ["", "2"]) == 0
    out = capsys.readouterr().out
    assert out.startswith(kc.MENU_BANNER + "\n"
                          "Found 1 SVS, 1 KFBF slide(s).\n"
                          "What do you want to convert?\n"
                          "  [1] Convert all (Default - press Enter)\n"
                          "  [2] Only SVS    -> BigTIFF\n"
                          "  [3] Only KFBF   -> choose the format next\n")
    assert "Output format for the 1 KFBF slide(s):\n" + MENU in out
    assert sorted(p.name for p in mixed_dir.glob("*.tif")) == ["bright.tif", "fluo.rgb.tif"]


def test_menu_brightfield_only_confirm(tmp_path, monkeypatch, capsys):
    make_svs(tmp_path / "b.svs")
    assert _interactive(monkeypatch, tmp_path, ["q"]) == 130
    out = capsys.readouterr().out
    assert out.startswith(kc.MENU_BANNER + "\nFound 1 brightfield (.svs) slide(s).\n")
    assert "Select output format" not in out
    assert not (tmp_path / "b.tif").exists()
    assert _interactive(monkeypatch, tmp_path, [""]) == 0
    assert (tmp_path / "b.tif").exists()
