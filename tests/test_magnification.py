"""Objective magnification: --info, --magnification, OME-XML Objective."""

import os
import re
from pathlib import Path

import pytest
import tifffile

import kfbf_converter as kc
from conftest import make_kfbf, make_svs

NS = {"o": kc.OME_NS}


def _ome_xml(path: Path) -> str:
    with tifffile.TiffFile(path) as t:
        return t.pages[0].description


def _check_objective(xml: str, expected: str) -> None:
    assert f'<Objective ID="Objective:0:0" NominalMagnification="{expected}"/>' in xml
    # schema order: Instrument before Image; refs inside Image before Pixels
    assert xml.index("<Instrument ") < xml.index("<Image ")
    image = xml[xml.index("<Image "):]
    assert image.index('<InstrumentRef ID="Instrument:0"/>') \
        < image.index('<ObjectiveSettings ID="Objective:0:0"/>') < image.index("<Pixels ")


def _schema():
    """OME 2016-06 XSD, if lxml is installed and OME_XSD points at a local copy
    of ome.xsd (with its xml.xsd import resolvable)."""
    xsd = os.environ.get("OME_XSD")
    if not xsd or not Path(xsd).exists():
        return None
    try:
        from lxml import etree
    except ImportError:
        return None
    return etree.XMLSchema(etree.parse(xsd))


def test_info_shows_magnification(kfbf, svs, capsys):
    assert kc.main([str(kfbf.path), "--info"]) == 0
    assert kc.main([str(svs.path), "--info"]) == 0
    out = capsys.readouterr().out
    assert "Magnification : 40x (header)" in out
    assert "Magnification : 20x (Aperio AppMag)" in out
    assert kc.main([str(svs.path), "--info", "--magnification", "40"]) == 0
    assert "Magnification : 40x (--magnification)" in capsys.readouterr().out


def test_missing_or_implausible_magnification(tmp_path, capsys):
    zero = make_kfbf(tmp_path / "zero.kfbf", magnification=0)
    crazy = make_kfbf(tmp_path / "crazy.kfbf", magnification=100000)
    for info in (zero, crazy):
        with kc.KFBioSlide(info.path) as s:
            assert s.magnification == 0 and s.magnification_source == "missing"
            assert "unknown (set --magnification)" in s.summary()
    with kc.KFBioSlide(crazy.path) as s:
        assert any("implausible magnification" in w for w in s.warnings)
    with kc.KFBioSlide(zero.path, magnification=20) as s:
        assert s.magnification == 20 and s.magnification_source == "--magnification"


def test_cli_rejects_bad_magnification(kfbf, capsys):
    assert kc.main([str(kfbf.path), "--magnification", "0"]) == 2
    assert "--magnification must be between" in capsys.readouterr().err


def test_ome_has_objective(kfbf, tmp_path):
    (p,) = kc.convert_slide(kfbf.path, tmp_path / "out", show_progress=False)
    xml = _ome_xml(p)
    _check_objective(xml, "40.0")
    with tifffile.TiffFile(p) as t:          # still a normal OME-TIFF pyramid
        assert t.is_ome and len(t.series[0].levels) == 2
        assert "DAPI" in t.ome_metadata
    schema = _schema()
    if schema is not None:
        from lxml import etree
        assert schema.validate(etree.fromstring(xml.encode("utf-8"))), schema.error_log


def test_ome_magnification_override_and_zstack(tmp_path):
    info = make_kfbf(tmp_path / "z.kfbf", n_z=2, missing=())
    (p,) = kc.convert_slide(info.path, tmp_path / "out", magnification=20,
                            show_progress=False)
    _check_objective(_ome_xml(p), "20.0")
    with tifffile.TiffFile(p) as t:
        assert len(t.pages) == 6


def test_ome_without_magnification_has_no_objective(tmp_path):
    info = make_kfbf(tmp_path / "nomag.kfbf", magnification=0)
    (p,) = kc.convert_slide(info.path, tmp_path / "out", show_progress=False)
    xml = _ome_xml(p)
    assert "<Instrument" not in xml and "<ObjectiveSettings" not in xml


def test_ome_add_objective_unit():
    xml = ('<?xml version="1.0" encoding="UTF-8"?><OME xmlns="' + kc.OME_NS + '">'
           '<Image ID="Image:0"><Description>x</Description><Pixels ID="Pixels:0"/></Image>'
           '<Image ID="Image:1"><StageLabel Name="s"/><Pixels ID="Pixels:1"/></Image></OME>')
    out = kc.ome_add_objective(xml, 2.5)
    assert 'NominalMagnification="2.5"' in out
    assert out.count("<ObjectiveSettings") == 2
    # after Description, before StageLabel (schema order)
    assert re.search(r"</Description><InstrumentRef[^>]*/><ObjectiveSettings[^>]*/><Pixels", out)
    assert re.search(r"<ObjectiveSettings[^>]*/><StageLabel", out)
    with pytest.raises(ValueError, match="already has"):
        kc.ome_add_objective(out, 40)
    with pytest.raises(ValueError, match="no Image"):
        kc.ome_add_objective('<OME xmlns="x"/>', 40)


def test_non_ascii_ome_survives_rewrite(kfbf, tmp_path):
    """PhysicalSizeXUnit is 'µm' (UTF-8) and must stay intact."""
    (p,) = kc.convert_slide(kfbf.path, tmp_path / "out", show_progress=False)
    assert 'PhysicalSizeXUnit="µm"' in _ome_xml(p)
