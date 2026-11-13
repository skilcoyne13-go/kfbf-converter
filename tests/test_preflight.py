"""'Data is too big for the file format chosen.': preflight, write errors,
read-back verification, stale .part cleanup."""

import collections
import dataclasses
import errno
import logging
import os
import sys
from pathlib import Path

import pytest

import kfbf_converter as kc

HEADLINE = "Data is too big for the file format chosen."
GiB = kc.GiB
Usage = collections.namedtuple("Usage", "total used free")


def _out(tmp_path):
    d = tmp_path / "out"
    d.mkdir(exist_ok=True)
    return d


def _lines(err):
    return err.splitlines()


def _message_after(lines, headline=HEADLINE):
    i = lines.index(headline)            # exact line, no prefix
    return lines[i + 1], lines[i + 2]


def _no_leftovers(d: Path, allowed=()):
    names = sorted(os.listdir(d))
    assert not [n for n in names if n.endswith(".part")], names
    assert sorted(n for n in names if n not in allowed) == []


# ---- free space -------------------------------------------------------------

def test_free_space_refusal_message_and_batch_continues(mixed_dir, tmp_path, monkeypatch,
                                                        capsys):
    out = _out(tmp_path)
    monkeypatch.setattr(kc, "filesystem_type", lambda p: "NTFS")
    monkeypatch.setattr(kc.shutil, "disk_usage",
                        lambda p: Usage(500 * GiB, 500 * GiB - 10 * 2 ** 20, 10 * 2 ** 20))
    rc = kc.main([str(mixed_dir), "-o", str(out), "-q"])
    err = capsys.readouterr().err
    lines = _lines(err)
    assert rc == 1
    assert lines.count(HEADLINE) == 2           # both slides were attempted
    reason, fix = _message_after(lines)
    assert "GB" in reason and "only 0.01 GB is free on" in reason
    assert fix.startswith("  Fix: free up ") and "choose another drive with -o" in fix
    assert "--compression jpeg" not in err     # would not fit in 10 MB either
    assert "Traceback" not in err
    assert "2 failure(s)" in err and "FAILED" in err
    _no_leftovers(out)


def test_traceback_only_with_verbose(kfbf, tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(kc.shutil, "disk_usage", lambda p: Usage(1, 1, 0))
    kc.main([str(kfbf.path), "-o", str(_out(tmp_path)), "-q", "-v"])
    assert "Traceback" in capsys.readouterr().err


# ---- file system ---------------------------------------------------------------

def _scaled_estimate(monkeypatch, kinds, size):
    real = kc.estimate_output

    def fake(slide, kind, writer, copy_tiles=False):
        est = real(slide, kind, writer, copy_tiles)
        if kind in kinds:
            est = dataclasses.replace(est, estimate=size, uncompressed=max(size, est.uncompressed))
        return est
    monkeypatch.setattr(kc, "estimate_output", fake)


def test_fat32_refuses_only_the_big_output(mixed_dir, tmp_path, monkeypatch, capsys):
    out = _out(tmp_path)
    monkeypatch.setattr(kc, "filesystem_type", lambda p: "FAT32")
    _scaled_estimate(monkeypatch, {"ome"}, int(5.2 * GiB))
    rc = kc.main([str(mixed_dir), "-o", str(out), "-q", "--mode", "all"])
    err = capsys.readouterr().err
    lines = _lines(err)
    assert rc == 1 and lines.count(HEADLINE) == 1
    reason, fix = _message_after(lines)
    assert reason == (f"  fluo.ome.tif: would be about 5.20 GB, but the output drive "
                      f"{kc._volume_root(out)} is FAT32, which cannot store files larger "
                      f"than 4.00 GB.")
    assert fix.startswith("  Fix: choose an NTFS or exFAT drive with -o")
    _no_leftovers(out, allowed={"fluo.rgb.tif", "fluo_DAPI.tif", "fluo_FITC.tif",
                                "fluo_TRITC.tif", "bright.tif"})
    assert "5 file(s) written, 1 failure(s)" in err


def test_fat32_warns_when_only_worst_case_exceeds(kfbf, tmp_path, monkeypatch, capsys):
    out = _out(tmp_path)
    monkeypatch.setattr(kc, "filesystem_type", lambda p: "FAT32")
    real = kc.estimate_output
    monkeypatch.setattr(kc, "estimate_output", lambda *a, **k: dataclasses.replace(
        real(*a, **k), uncompressed=6 * GiB))
    assert kc.main([str(kfbf.path), "-o", str(out), "-q"]) == 0
    err = capsys.readouterr().err
    assert "WARNING: fluo.ome.tif: the output drive" in err and "could reach 6.00 GB" in err
    assert HEADLINE not in err
    assert (out / "fluo.ome.tif").exists()


def test_no_preflight_escape_hatch(kfbf, tmp_path, monkeypatch, capsys):
    out = _out(tmp_path)
    monkeypatch.setattr(kc, "filesystem_type", lambda p: "FAT32")
    _scaled_estimate(monkeypatch, {"ome"}, 9 * GiB)
    assert kc.main([str(kfbf.path), "-o", str(out), "-q", "--no-preflight"]) == 0
    assert HEADLINE not in capsys.readouterr().err
    assert (out / "fluo.ome.tif").exists()


def test_filesystem_detection_on_this_machine(tmp_path):
    fs = kc.filesystem_type(tmp_path)
    if os.name == "nt" or sys.platform.startswith("linux"):
        assert isinstance(fs, str) and fs
    assert kc.is_fat("FAT32") and kc.is_fat("vfat") and not kc.is_fat("exFAT")
    assert not kc.is_fat("NTFS") and not kc.is_fat(None)


# ---- format limits (unit) ----------------------------------------------------

def _est(**kw):
    base = dict(estimate=GiB, uncompressed=3 * GiB, exact=False, n_tiles=1000,
                tiles_per_ifd=800, max_dim=50000, tile_size=512, compression="zlib",
                overhead=10 ** 5)
    base.update(kw)
    return kc.SizeEstimate(**base)


@pytest.fixture
def roomy(monkeypatch):
    monkeypatch.setattr(kc, "filesystem_type", lambda p: "NTFS")
    monkeypatch.setattr(kc.shutil, "disk_usage", lambda p: Usage(10 ** 15, 0, 10 ** 15))


def test_classic_tiff_limit(tmp_path, roomy):
    out = tmp_path / "x.tif"
    with pytest.raises(kc.DataTooBigError) as ei:
        kc.preflight_check(_est(estimate=5 * GiB), out, bigtiff=False)
    assert ei.value.lines()[0] == HEADLINE
    assert "classic TIFF files are limited to 4.00 GB" in ei.value.reason
    assert kc.preflight_check(_est(estimate=5 * GiB), out, bigtiff=True) == []
    warn = kc.preflight_check(_est(estimate=GiB, uncompressed=5 * GiB), out, bigtiff=False)
    assert warn and "4 GB classic TIFF limit" in warn[0]


def test_dimension_and_jpeg_tile_limits(tmp_path, roomy):
    with pytest.raises(kc.DataTooBigError, match="32-bit"):
        kc.preflight_check(_est(max_dim=2 ** 32), tmp_path / "x.tif")
    with pytest.raises(kc.DataTooBigError, match="JPEG tiles are limited to 65535"):
        kc.preflight_check(_est(compression="jpeg", tile_size=70000), tmp_path / "x.tif")


def test_jpeg_hint_only_when_it_would_fit(tmp_path, monkeypatch):
    monkeypatch.setattr(kc, "filesystem_type", lambda p: "FAT32")
    monkeypatch.setattr(kc.shutil, "disk_usage", lambda p: Usage(10 ** 15, 0, 10 ** 15))
    with pytest.raises(kc.DataTooBigError) as ei:
        kc.preflight_check(_est(estimate=6 * GiB), tmp_path / "x.tif",
                           jpeg_estimate=int(1.5 * GiB))
    assert ei.value.fix == ("choose an NTFS or exFAT drive with -o, or use --compression "
                            "jpeg (about 1.50 GB).")
    with pytest.raises(kc.DataTooBigError) as ei:
        kc.preflight_check(_est(estimate=12 * GiB), tmp_path / "x.tif",
                           jpeg_estimate=int(4.5 * GiB))
    assert "--compression jpeg" not in ei.value.fix


def test_free_space_margin_numbers(tmp_path, monkeypatch):
    monkeypatch.setattr(kc, "filesystem_type", lambda p: "NTFS")
    monkeypatch.setattr(kc.shutil, "disk_usage", lambda p: Usage(10 * GiB, 8 * GiB, 2 * GiB))
    est = _est(estimate=2 * GiB)       # fits exactly, but not with the margin
    with pytest.raises(kc.DataTooBigError) as ei:
        kc.preflight_check(est, tmp_path / "x.tif")
    margin = kc.free_space_margin(2 * GiB)
    assert f"needs about 2.00 GB (plus a {kc._gb(margin)} safety margin)" in ei.value.reason
    assert f"free up {kc._gb(margin)}" in ei.value.fix
    assert kc.preflight_check(_est(estimate=GiB), tmp_path / "x.tif") == []


def test_estimates_bound_real_outputs(kfbf, svs, tmp_path):
    out = _out(tmp_path)
    for path in (kfbf.path, svs.path):
        with kc.open_slide(path) as s:
            levels = kc.select_export_levels(s.levels)
            w = kc._PyramidWriter(s, 512, "deflate", 2, show_progress=False, levels=levels)
            kind = "brightfield" if s.is_rgb else "ome"
            est = kc.estimate_output(s, kind, w)
            p = (kc.write_brightfield if s.is_rgb else kc.write_ome)(s, out / f"{kind}.tif", w)
            assert p.stat().st_size <= est.estimate <= est.uncompressed
    with kc.open_slide(svs.path) as s:   # tile copy: estimate is (almost) exact
        levels = kc.select_export_levels(s.levels)
        w = kc._PyramidWriter(s, 512, "deflate", 2, show_progress=False, levels=levels)
        est = kc.estimate_output(s, "brightfield", w, copy_tiles=True)
        p = kc.write_brightfield(s, out / "copy.tif", w, copy_tiles=True)
        # exact tile bytes; only the (small) thumbnail page is estimated
        assert est.exact and p.stat().st_size <= est.estimate <= p.stat().st_size + 2 ** 21


# ---- errors while writing ------------------------------------------------------

_REAL_ORDERED_MAP = kc._ordered_map


def _fail_after(monkeypatch, n, exc):
    real = _REAL_ORDERED_MAP

    def failing(fn, items, workers):
        for i, x in enumerate(real(fn, items, workers)):
            if i == n:
                raise exc
            yield x
    monkeypatch.setattr(kc, "_ordered_map", failing)


@pytest.mark.parametrize("which", ["kfbf", "svs"])
def test_enospc_mid_write_leaves_nothing(which, kfbf, svs, tmp_path, monkeypatch, capsys):
    out = _out(tmp_path)
    path = kfbf.path if which == "kfbf" else svs.path
    _fail_after(monkeypatch, 3, OSError(errno.ENOSPC, "No space left on device"))
    errors = []
    assert kc.convert_slide(path, out, show_progress=False, errors=errors) == []
    (_, exc), = errors
    assert isinstance(exc, kc.DataTooBigError) and "ran out of space while writing" in exc.reason
    assert capsys.readouterr().err.splitlines()[0] == HEADLINE
    _no_leftovers(out)
    with pytest.raises(kc.DataTooBigError):             # API without errors= raises
        kc.convert_slide(path, out, show_progress=False)
    _no_leftovers(out)


@pytest.mark.skipif(os.name != "nt", reason="Windows error codes")
@pytest.mark.parametrize("winerror,text", [(112, "ran out of space"), (39, "ran out of space"),
                                           (223, "refused a larger file")])
def test_windows_write_errors(kfbf, tmp_path, monkeypatch, winerror, text):
    out = _out(tmp_path)
    _fail_after(monkeypatch, 2, OSError(0, "disk error", None, winerror))
    errors = []
    kc.convert_slide(kfbf.path, out, show_progress=False, errors=errors)
    assert text in errors[0][1].reason
    _no_leftovers(out)


def test_efbig_and_tifffile_size_errors(kfbf, tmp_path, monkeypatch):
    out = _out(tmp_path)
    _fail_after(monkeypatch, 1, OSError(errno.EFBIG, "File too large"))
    errors = []
    kc.convert_slide(kfbf.path, out, show_progress=False, errors=errors)
    assert "NTFS or exFAT" in errors[0][1].fix
    _fail_after(monkeypatch, 1, ValueError("data too large for standard TIFF file"))
    errors = []
    kc.convert_slide(kfbf.path, out, show_progress=False, errors=errors)
    assert "tifffile refused the data" in errors[0][1].reason
    _no_leftovers(out)


def test_other_errors_are_not_disguised(kfbf, tmp_path, monkeypatch):
    out = _out(tmp_path)
    _fail_after(monkeypatch, 1, RuntimeError("boom"))
    with pytest.raises(RuntimeError, match="boom"):
        kc.convert_slide(kfbf.path, out, show_progress=False, errors=[])
    _no_leftovers(out)


# ---- read-back verification -----------------------------------------------------

def test_truncated_output_fails_verification_and_is_deleted(kfbf, tmp_path, monkeypatch,
                                                            capsys):
    out = _out(tmp_path)
    real = kc.verify_output

    def truncate_then_verify(path, *a, **k):
        with open(path, "r+b") as fh:            # simulate a write cut short
            fh.truncate(path.stat().st_size // 2)
        return real(path, *a, **k)
    monkeypatch.setattr(kc, "verify_output", truncate_then_verify)
    assert kc.main([str(kfbf.path), "-o", str(out), "-q"]) == 1
    err = capsys.readouterr().err
    line = next(l for l in err.splitlines() if l.startswith(kc.VerificationError.PREFIX))
    assert line.startswith("Output failed verification and was deleted: ")
    assert "1 failure(s)" in err and "Traceback" not in err
    _no_leftovers(out)


def test_verification_failure_message(mixed_dir, tmp_path, monkeypatch, capsys):
    out = _out(tmp_path)

    def bad(*a, **k):
        raise kc.VerificationError("level 2 has shape (1, 1), expected (2, 2)")
    monkeypatch.setattr(kc, "verify_output", bad)
    assert kc.main([str(mixed_dir), "-o", str(out), "-q"]) == 1
    err = capsys.readouterr().err
    assert err.splitlines().count("Output failed verification and was deleted: level 2 has "
                                  "shape (1, 1), expected (2, 2)") == 2
    _no_leftovers(out)


def test_verify_output_checks(kfbf, tmp_path):
    out = _out(tmp_path)
    (p,) = kc.convert_slide(kfbf.path, out, mode="ome", show_progress=False)
    shapes = [(520, 700), (260, 350)]
    kc.verify_output(p, shapes, n_pages=3, ome=True)
    with pytest.raises(kc.VerificationError, match="expected 4 image page"):
        kc.verify_output(p, shapes, n_pages=4)
    with pytest.raises(kc.VerificationError, match="pyramid level"):
        kc.verify_output(p, shapes[:1], n_pages=3)
    with pytest.raises(kc.VerificationError, match="level 1 has shape"):
        kc.verify_output(p, [(520, 700), (261, 350)], n_pages=3)
    p.write_bytes(b"II*\0garbage")
    with pytest.raises(kc.VerificationError, match="could not read it back|found 0"):
        kc.verify_output(p, shapes, n_pages=3)


# ---- stale .part ---------------------------------------------------------------

def test_stale_part_files_removed(kfbf, tmp_path, caplog):
    out = _out(tmp_path)
    stale = out / "fluo.ome.tif.part"
    stale.write_bytes(b"left over from a crash")
    other = out / "unrelated.tif.part"
    other.write_bytes(b"not ours")
    with caplog.at_level(logging.INFO, logger="kfbf_converter"):
        kc.convert_slide(kfbf.path, out, show_progress=False)
    assert not stale.exists() and other.exists()
    assert "removed stale partial file" in caplog.text
    assert (out / "fluo.ome.tif").exists()
