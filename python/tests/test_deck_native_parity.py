# SPDX-License-Identifier: Apache-2.0
"""Deck-validation rules that mirror the native reader and the standalone driver.

Each test pins one rule of ``src/CableDyn_DeckDriver.f90`` so that
:class:`cabledyn.DeckFile` neither rejects a deck the solver accepts nor
certifies one the solver refuses.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import subprocess
from pathlib import Path

import numpy as np
import pytest

from cabledyn import DeckFile, DeckFormatError, StudyFormatError, generate_deck_cases, run_study
from cabledyn.builder import DeckModel, OptionSet

ROOT = Path(__file__).resolve().parents[2]
EXAMPLES = ROOT / "examples"

_BAR = "-" * 21
# A held chain on a dynamic clock (the static reference plus dtM/TMax).
_HELD = f"""\
{_BAR} CableDyn Input File {_BAR}
in-memory fixture: one grounded chain line on a dynamic clock
{_BAR} LINE TYPES {_BAR}
Name Diam Mass EA BA EI Cdn Cdt Can Cat
chain 0.252 390.0 1.674e9 -1.0 0.0 1.37 0.64 1.0 0.0
{_BAR} POINTS {_BAR}
ID Type X Y Z Mass Vol CdA Ca
1 Fixed 400.0 0.0 -50.0 0 0 0 0
2 Coupled 0.0 0.0 0.0 0 0 0 0
{_BAR} LINES {_BAR}
ID NodeA NodeB Outputs
1 2 1 -
{_BAR} SECTIONS {_BAR}
LineID LineType Length NumSegs
1 chain 410.0 41
{_BAR} OPTIONS {_BAR}
9.80665 g
1025.0 rhoW
50.0 WtrDpth
0.1 dtM
10.0 TMax
{_BAR} OUTPUTS {_BAR}
FairTen1
{_BAR} need this line {_BAR}
"""
_STATIC = _HELD.replace("0.1 dtM\n10.0 TMax\n", "")
_WAVE_FREQUENCIES = "0 0 0 0\n0.5 0.5 0 0\n0.6 0.3 0.1 0\n"


def _edit(base: str, *edits: tuple[str, str]) -> str:
    text = base
    for old, new in edits:
        assert old in text, old
        text = text.replace(old, new, 1)
    return text


def _options(base: str, *rows: str) -> str:
    return _edit(base, ("50.0 WtrDpth\n", "50.0 WtrDpth\n" + "".join(f"{r}\n" for r in rows)))


def _example(name: str, *edits: tuple[str, str]) -> str:
    return _edit((EXAMPLES / name).read_text(encoding="utf-8"), *edits)


def _section(base: str, heading: str, *rows: str) -> str:
    block = f"{_BAR} {heading} {_BAR}\n" + "".join(f"{row}\n" for row in rows)
    marker = next(line for line in base.splitlines() if "OPTIONS" in line and "---" in line)
    return _edit(base, (marker, block + marker))


def _rejects(text: str, message: str, *, caller_driven: bool = False) -> None:
    with pytest.raises(DeckFormatError, match=message):
        DeckFile.from_text(text, caller_driven=caller_driven)


# --------------------------------------------------------------------------- options


@pytest.mark.parametrize("value", ["0", "-0.1"])
def test_dtwave_must_be_positive(value):
    _rejects(_options(_HELD, f"{value} dtWave"), "OPTION dtWave must be positive")
    assert DeckFile.from_text(_options(_HELD, "0.25 dtWave")).option("dtWave").values == ("0.25",)


@pytest.mark.parametrize(
    ("rows", "message"),
    [
        (("7 WaterKin",), None),
        (("3 WaveKin", "0.1 dtWave"), None),
        (("1 Currents",), None),
        (("airy 2 8 0 waves", "7 WaterKin"), r"both a waves \(or wavetrain\) OPTION and WaveKin 7"),
        (("airy 2 8 0 wavetrain", "3 WaveKin"), "double-counting"),
        (("uniform 0.1 0 0 current", "1 Currents"), "both a current OPTION and Currents 1"),
    ],
)
def test_moordyn_c_kinematics_follow_the_native_rules(rows, message):
    text = _options(_HELD, *rows)
    if message is None:
        DeckFile.from_text(text)
    else:
        _rejects(text, message)


def test_moordyn_c_kinematics_need_the_dynamic_sea_setup():
    _rejects(_options(_STATIC, "7 WaterKin"), "WaveKin 7 waves require dtM, TMax and WtrDpth")
    no_depth = _edit(_HELD, ("50.0 WtrDpth\n", "7 WaterKin\n"))
    _rejects(no_depth.replace("-50.0", "-80.0"), "WaveKin 7 waves require dtM, TMax and WtrDpth")
    _rejects(_options(_STATIC, "1 Currents"), "standalone Currents 1 requires dtM and TMax")


@pytest.mark.parametrize("row", ["7 WaterKin", "1 Currents"])
def test_moordyn_c_kinematics_are_not_certified_for_the_coupled_route(row):
    _rejects(_options(_HELD, row), "not supported on the coupled route", caller_driven=True)


def test_vessel_rao_takes_moordyn_c_component_waves():
    text = _example(
        "lazy_wave_vessel_rao.dat",
        ("jonswap 4.0 10.0 3.3 0.0 waves", "none waves"),
        ("0         WaterKin", "7 WaterKin"),
    )
    assert DeckFile.from_text(text, path=EXAMPLES / "rao.dat").option("WaterKin").values == ("7",)
    without = _example("lazy_wave_vessel_rao.dat", ("jonswap 4.0 10.0 3.3 0.0 waves", "none waves"))
    _rejects(without, "vesselRAO needs deck waves")


@pytest.mark.parametrize(
    "rows",
    [
        ("jonswap 2 8 40 0 waves",),
        ("jonswap 2 8 40 0 waves", "2 WaveSpreading"),
        ("jonswap 2 8 40 0 waves", "300 WaveComponents"),
    ],
)
def test_jonswap_gamma_limit_applies_on_every_path(rows):
    _rejects(_options(_HELD, *rows), r"JONSWAP gamma must be in \[1, 32\.6\)")


def test_spread_jonswap_row_is_validated_as_a_wave_train():
    DeckFile.from_text(_options(_HELD, "jonswap 2 8 3.3 2e6 waves"))
    _rejects(
        _options(_HELD, "jonswap 2 8 3.3 2e6 waves", "2 WaveSpreading"),
        r"\|direction\| <= 1e6",
    )
    DeckFile.from_text(_options(_HELD, "jonswap 2 8 3.3 30 waves", "200 WaveComponents"))


def test_wavetrain_spellings_are_one_option():
    text = _options(_HELD, "airy 2 8 0 wave_train", "airy 1 6 30 wavetrain")
    deck = DeckFile.from_text(text)
    assert deck.option("wave_train").values == ("airy", "1", "6", "30")
    assert deck.option("wavetrain").keyword == "wavetrain"
    options = OptionSet()
    options.add("wave_train", "airy", 2.0, 8.0, 0.0)
    options.add("wavetrain", "airy", 1.0, 6.0, 30.0)
    assert options.remove("wavetrain") == 2


def test_float_edits_use_the_shortest_round_trip_text():
    deck = DeckFile.from_text(_HELD)
    deck.set_option("dtM", 0.1)
    deck.set_option("TMax", 10.0)
    assert "\n0.1 dtM\n10 TMax\n" in deck.text()


# --------------------------------------------------------------------------- seabed, attachments


@pytest.mark.parametrize(
    ("depth", "anchor_z", "accepted"),
    [("200.0", "-200.0015", True), ("200.0", "-200.0021", False), ("1.0", "-1.001", True)],
)
def test_anchor_seabed_tolerance_is_native(depth, anchor_z, accepted):
    text = _edit(
        _HELD, ("50.0 WtrDpth", f"{depth} WtrDpth"), ("400.0 0.0 -50.0", f"400.0 0.0 {anchor_z}")
    )
    if accepted:
        DeckFile.from_text(text)
    else:
        _rejects(text, re.escape(f"z = {float(anchor_z):.10g} m is under"))


@pytest.mark.parametrize(
    ("arc", "accepted"),
    [("5000", False), ("170.215", True), ("70:7:118.12", True), ("100:10:180", False)],
)
def test_attachment_arc_length_lies_on_the_line(arc, accepted):
    text = _example("lazy_wave_buoyancy_modules.dat", ("70.614:5.0:115.614", arc))
    if accepted:
        DeckFile.from_text(text)
    else:
        _rejects(text, "ATTACHMENTS ArcLength lies beyond the end of line 1")


# --------------------------------------------------------------------------- output channels


@pytest.mark.parametrize(
    "channel", ["FairTen01", "AnchIncl001", "TDP01s", "L01N1px", "Ten01N01", "Point01px", "Con01pz"]
)
def test_line_and_point_channel_ids_may_have_leading_zeros(channel):
    DeckFile.from_text(_edit(_HELD, ("FairTen1\n", f"{channel}\n")))


@pytest.mark.parametrize("channel", ["TDP1q", "Point1pq", "Ten9N0"])
def test_malformed_line_and_point_channels_fail_closed(channel):
    _rejects(_edit(_HELD, ("FairTen1\n", f"FairTen1\n{channel}\n")), "unsupported OUTPUT")


@pytest.mark.parametrize("channel", ["Ten1N0", "L1N0px", "BendMom1N0"])
def test_node_zero_channels_name_the_numbering(channel):
    _rejects(_edit(_HELD, ("FairTen1\n", f"FairTen1\n{channel}\n")), "node numbers start at 1")


_R6_ROW = (
    "1   Rigid6  0.0  0.0  -20.0  0.0   0.0    0.0   2.0e4   40.0  0.0  0.0  0.0  8.0   0.5  "
    "3.5e4   3.5e4   3.5e4"
)


def _buoy(*edits: tuple[str, str]) -> str:
    return _example("rigid6_buoy.dat", *edits)


@pytest.mark.parametrize("channel", ["Body01Px", "Body1RVx", "Body1RAz", "Body1Rz", "Body1Mz"])
def test_body_channels_follow_the_native_quantities(channel):
    DeckFile.from_text(_buoy(('"FairTen1"', f'"{channel}"')))


@pytest.mark.parametrize("channel", ["Body1RPx", "Body1TenA", "Body1Sub"])
def test_body_channels_outside_the_native_quantities_fail_closed(channel):
    _rejects(_buoy(('"FairTen1"', f'"{channel}"')), "unsupported OUTPUT channel")


def _point3_buoy(*edits: tuple[str, str]) -> str:
    return _buoy(
        (_R6_ROW, "1 Point3 0.0 0.0 -20.0 0.0 0.0 0.0 2.0e4 40.0 0.0 0.0 0.0 8.0 0.5"),
        ("2     Body1     -0.75      1.299      -2.0      0       0       0      0\n", ""),
        ("3     Body1     -0.75      -1.299     -2.0      0       0       0      0\n", ""),
        ("2     2       5       -", "2     1       5       -"),
        ("3     3       6       -", "3     1       6       -"),
        ('"Point2pz"\n', ""),
        ('"Point3pz"\n', ""),
        *edits,
    )


def test_body_channels_report_rigid6_bodies_only():
    DeckFile.from_text(_point3_buoy())
    _rejects(_point3_buoy(('"FairTen1"', '"Body1Px"')), "Body<N> channels report Rigid6 bodies")


def _spar(*edits: tuple[str, str]) -> str:
    return _example("rod_moored_spar.dat", *edits)


_SPAR_ROD = "1    spar      Free   0.0   0.0   -30.0   0.0   0.0   -20.0   1         p"


@pytest.mark.parametrize(
    ("channel", "accepted"),
    [("Rod1N0px", True), ("Rod01N1pz", True), ("Rod1N2px", False), ("Rod1RVx", True)],
)
def test_rod_node_channels_stay_within_numsegs(channel, accepted):
    text = _spar(('"FairTen1"', f'"{channel}"'))
    if accepted:
        DeckFile.from_text(text)
    else:
        _rejects(text, "rod node beyond NumSegs")


def test_rod_channel_quantities_follow_the_native_vocabulary():
    _rejects(_spar(('"FairTen1"', '"Rod1RPx"')), "unsupported OUTPUT channel")


# --------------------------------------------------------------------------- rods


def test_rod_endpoints_have_an_admissible_magnitude():
    text = _spar((_SPAR_ROD, _SPAR_ROD.replace("0.0   0.0   -30.0", "2.0e6   0.0   -30.0")))
    _rejects(text, "ROD endpoint coordinates must be at most 1e6 m")


def test_rod_numsegs_has_no_upper_bound():
    DeckFile.from_text(_spar((_SPAR_ROD, _SPAR_ROD.replace("   1         p", "   2000000   p"))))


def _zero_length(kind: str, rod_type: str = "spar") -> str:
    return _spar(
        (_SPAR_ROD, f"1    {rod_type}  {kind}   0.0   0.0   -30.0   0.0   0.0   -20.0   0   -"),
        ('"Point2px"\n', ""),
        ('"Point2pz"\n', ""),
    )


def test_zero_length_rod_is_a_point_before_the_rod_checks():
    # native collapse_zero_length_rods runs before the rod-type and motion checks
    DeckFile.from_text(_zero_length("Free", rod_type="nope"))
    DeckFile.from_text(_zero_length("Coupled"))


def test_zero_length_rod_merges_its_end_points():
    # the Rod1B POINT row (id 2) becomes the Rod1A connector: its channels no longer exist
    text = _spar((_SPAR_ROD, "1    spar  Free   0.0   0.0   -30.0   0.0   0.0   -20.0   0   -"))
    _rejects(text, "'Point2px' references an unknown point")


def test_both_rod_end_spellings_name_one_point():
    _rejects(
        _spar(("1     1       3       -", "1     R1A     Rod1A   -")),
        "LINE NodeA and NodeB must differ",
    )
    _rejects(
        _spar(("1     1       3       -", "1     R1A     1       -")),
        "LINE NodeA and NodeB must differ",
    )


def test_body_fixed_rod_ends_resolve_to_body_points():
    text = _buoy(
        (
            f"{_BAR} POINTS",
            f"{_BAR} ROD TYPES {_BAR}\nName Diam Mass Cd Ca CdEnd CaEnd\n"
            "spar 1.0 300.0 0.8 1.0 0.0 0.0\n"
            f"{_BAR} RODS {_BAR}\nID RodType Type XA YA ZA XB YB ZB NumSegs Outputs\n"
            "1 spar Body1 0 0 -2 0 0 -8 2 -\n"
            f"{_BAR} POINTS",
        ),
        ("1     1       4       -", "1     R1B     4       -"),
    )
    deck = DeckFile.from_text(text)
    types: dict[str, str] = {}
    keys, merged = deck._resolve_rod_ends(types, set(), {})
    assert keys["rod1b"] == "rod1b" and types["rod1b"] == "body1" and merged == {}


# --------------------------------------------------------------------------- dispatch gates


_FAILURE = ("FAILURE", "ID Point Lines FailTime FailTen", "1 1 1 10.0 0")


def test_failure_needs_the_dynamic_point_system_route():
    DeckFile.from_text(_section(_HELD, *_FAILURE))
    _rejects(_section(_STATIC, *_FAILURE), r"FAILURE requires a dynamic deck \(dtM and TMax\)")
    _rejects(
        _section(_options(_HELD, "0.5 frictionMu"), *_FAILURE),
        "FAILURE deck does not support seabed friction",
    )
    finite = _edit(_HELD, ("1.674e9 -1.0 0.0", "1.674e9 -1.0 1.0e4"))
    _rejects(_section(finite, *_FAILURE), "FAILURE is not supported on finite-EI decks")
    rod_failure = ("FAILURE", "ID Point Lines FailTime FailTen", "1 3 1 10.0 0")
    _rejects(_section(_spar(), *rod_failure), "ROD decks")


def test_failure_on_rigid6_decks_takes_no_other_dynamic_point_or_motion():
    free = _buoy(
        ("4     Fixed     40.0", "7 Connect 20 0 -60 0 0 0 0\n4     Fixed     40.0"),
        ("3     3       6       -", "3     3       6       -\n4 7 4 -"),
        ("3        poly       86.85    20", "3        poly       86.85    20\n4 poly 70 10"),
    )
    failure = ("FAILURE", "ID Point Lines FailTime FailTen", "1 4 1 10.0 0")
    _rejects(_section(free, *failure), "does not support Connect/Free or Point3 body points")
    moving = _buoy(
        (
            "uniform 0.6 0.0 0.0  current",
            "data/vessel/vessel_surge_heave_pitch_12s_dt005.txt vesselMotion\n"
            "uniform 0.6 0.0 0.0  current",
        )
    )
    _rejects(_section(moving, *failure), "FAILURE deck does not support motionFile coupling")


@pytest.mark.parametrize(
    "body",
    [_R6_ROW, "1 Free 0.0 0.0 -20.0 0.0 0.0 0.0 2.0e4 0|0|0 3.5e4 40.0 8.0 0.5"],
)
def test_multibody_decks_take_no_prescribed_motion_in_either_body_spelling(body):
    text = _buoy(
        (_R6_ROW, body),
        (
            "4     Fixed     40.0",
            "7 Connect 20 0 -60 0 0 0 0\n8 Coupled 0 0 0 0 0 0 0\n4     Fixed     40.0",
        ),
        ("3     3       6       -", "3     3       6       -\n4 8 7 -"),
        ("3        poly       86.85    20", "3        poly       86.85    20\n4 poly 70 10"),
        (
            "uniform 0.6 0.0 0.0  current",
            "data/vessel/vessel_surge_heave_pitch_12s_dt005.txt vesselMotion\n"
            "uniform 0.6 0.0 0.0  current",
        ),
    )
    _rejects(text, "a mixed-topology deck .* does not support motionFile coupling")


def test_connect_point_deck_takes_no_prescribed_motion():
    text = _edit(
        _HELD,
        (
            "1 Fixed 400.0 0.0 -50.0 0 0 0 0",
            "1 Fixed 400.0 0.0 -50.0 0 0 0 0\n3 Connect 200 0 -30 0 0 0 0",
        ),
        ("1 2 1 -", "1 2 3 -\n2 3 1 -"),
        ("1 chain 410.0 41", "1 chain 210.0 21\n2 chain 210.0 21"),
        ("10.0 TMax\n", "10.0 TMax\nmotion.txt motionFile\n"),
    )
    _rejects(text, "a Connect/Free dynamic-point deck does not support motionFile coupling")


def _modes(base: str) -> str:
    marker = next(line for line in base.splitlines() if "OPTIONS" in line and "---" in line)
    return _edit(base, (marker + "\n", marker + "\n5 nModes\n"))


def test_modal_analysis_follows_the_native_deck_rules():
    DeckFile.from_text((EXAMPLES / "chain_modes.dat").read_text(encoding="utf-8"))
    DeckFile.from_text((EXAMPLES / "lazy_wave_modes.dat").read_text(encoding="utf-8"))
    _rejects(
        _example("chain_modes.dat", ("100.0        WtrDpth", "seabed.txt bathymetryFile")),
        "nModes supports a flat WtrDpth seabed only",
    )
    _rejects(_modes(_buoy()), "OPTION nModes needs a deck of lines between Fixed and")
    _rejects(_modes(_section(_HELD, *_FAILURE)), "no BODIES, RODS, Connect/Free points or FAILURE")


def test_modal_analysis_needs_a_fixed_end_b_on_finite_ei_lines():
    text = (
        "title\n"
        f"{_BAR} LINE TYPES {_BAR}\n"
        "chain 0.1 20.0 1.0e9 -1.0 1.0e4 1.2 0.2 1.0 0.0\n"
        f"{_BAR} TURBINES {_BAR}\n"
        "1 0 0 0\n2 600 0 0\n"
        f"{_BAR} POINTS {_BAR}\n"
        "1 Turbine1 0 0 -10 0 0 0 0\n2 Turbine2 0 0 -10 0 0 0 0\n3 Fixed 300 0 -100 0 0 0 0\n"
        f"{_BAR} LINES {_BAR}\n"
        "1 1 3 -\n2 1 2 -\n"
        f"{_BAR} SECTIONS {_BAR}\n"
        "1 chain 380 20\n2 chain 620 30\n"
        f"{_BAR} OPTIONS {_BAR}\n"
        "100 WtrDpth\n0.1 dtM\n1.0 TMax\n5 nModes\n"
        f"{_BAR} need this line {_BAR}\n"
    )
    _rejects(text, "nModes on finite-EI lines needs every End B Fixed")


def test_mixed_aggregate_route_takes_no_modal_analysis():
    # the executable runs a mixed EI=0/finite-EI deck without objects on the aggregate,
    # which writes no modal files
    text = _example("iea15mw_umaine_mixed_cabledyn.dat", ("0.0      TMax", "0.0 TMax\n5 nModes"))
    _rejects(text, "OPTION nModes is not supported for mixed EI=0/finite-EI decks")
    DeckFile.from_text(text, caller_driven=True)
    DeckFile.from_text(text.replace("5 nModes", "0 nModes"))


# --------------------------------------------------------------------------- row syntax


def test_unquoted_ea_text_may_not_mix_quote_characters():
    text = _edit(_HELD, ("1.674e9 -1.0 0.0", "SYROPE:a'b\"c|1|2 1|1 0.0"))
    _rejects(text, "mixes both quote characters")


# --------------------------------------------------------------------------- companion files


def _kinematic_base(folder: Path, *rows: str, files: dict[str, str] | None = None) -> Path:
    folder.mkdir(parents=True, exist_ok=True)
    for name, text in (files or {}).items():
        (folder / name).write_bytes(text.encode())
    path = folder / "base.dat"
    path.write_bytes(_options(_HELD, *rows).encode())
    return path


def test_generated_cases_carry_the_fixed_name_kinematics_files(tmp_path):
    base = _kinematic_base(
        tmp_path / "src",
        "7 WaterKin",
        "1 Currents",
        files={
            "wave_frequencies.txt": _WAVE_FREQUENCIES,
            "current_profile.txt": "-50 0.1 0 0\n0 0.5 0 0\n",
        },
    )
    out = tmp_path / "cases"
    generate_deck_cases(base, out, {"a": {"option.dtM": 0.05}, "b": {"option.WaterKin": 0}})
    assert (out / "wave_frequencies.txt").read_text(encoding="utf-8") == _WAVE_FREQUENCIES
    manifest = json.loads((out / "cases.json").read_text(encoding="utf-8"))
    digest = hashlib.sha256(_WAVE_FREQUENCIES.encode()).hexdigest()
    first, second = manifest["cases"]
    assert first["companion_files"]["wave_frequencies.txt"] == digest
    assert set(first["companion_files"]) == {"wave_frequencies.txt", "current_profile.txt"}
    assert set(second["companion_files"]) == {"current_profile.txt"}
    # run_study refuses a case whose kinematics file changed after generation
    (out / "wave_frequencies.txt").write_text("0 0 0 0\n", encoding="utf-8")
    with pytest.raises(StudyFormatError, match="companion file"):
        run_study(out / "cases.json", executable=tmp_path / "missing.exe")


def test_study_manifest_companion_entries_are_checked(tmp_path):
    base = _kinematic_base(
        tmp_path / "src", "7 WaterKin", files={"wave_frequencies.txt": _WAVE_FREQUENCIES}
    )
    out = tmp_path / "cases"
    generate_deck_cases(base, out, {"a": {}})
    manifest_path = out / "cases.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    for bad, message in (
        ([], "companion_files must be an object"),
        ({"../evil.txt": "0" * 64}, "unknown companion file"),
    ):
        manifest["cases"][0]["companion_files"] = bad
        manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
        with pytest.raises(StudyFormatError, match=message):
            run_study(manifest_path, executable=tmp_path / "missing.exe")


def test_case_generation_fails_closed_without_the_kinematics_file(tmp_path):
    base = _kinematic_base(tmp_path / "src", "7 WaterKin")
    with pytest.raises(DeckFormatError, match=r"reads wave_frequencies\.txt from its folder"):
        generate_deck_cases(base, tmp_path / "cases", {"a": {}})
    assert not (tmp_path / "cases").exists()


def test_case_generation_keeps_a_different_kinematics_file_without_overwrite(tmp_path):
    base = _kinematic_base(
        tmp_path / "src", "7 WaterKin", files={"wave_frequencies.txt": _WAVE_FREQUENCIES}
    )
    out = tmp_path / "cases"
    out.mkdir()
    (out / "wave_frequencies.txt").write_text("0 0 0 0\n", encoding="utf-8")
    with pytest.raises(FileExistsError, match=r"wave_frequencies\.txt"):
        generate_deck_cases(base, out, {"a": {}})
    generate_deck_cases(base, out, {"a": {}}, overwrite=True)
    assert (out / "wave_frequencies.txt").read_text(encoding="utf-8") == _WAVE_FREQUENCIES


def test_case_generation_in_the_source_folder_copies_nothing(tmp_path):
    base = _kinematic_base(
        tmp_path, "7 WaterKin", files={"wave_frequencies.txt": _WAVE_FREQUENCIES}
    )
    (case,) = generate_deck_cases(base, tmp_path, {"a": {}})
    assert case.deck.parent == tmp_path


def test_invalid_case_creates_no_output_folder(tmp_path):
    base = _kinematic_base(tmp_path / "src")
    with pytest.raises(DeckFormatError):
        generate_deck_cases(base, tmp_path / "cases", {"a": {"option.dtM": -1.0}})
    assert not (tmp_path / "cases").exists()


def test_model_save_copies_the_kinematics_files_when_rebasing(tmp_path):
    base = _kinematic_base(
        tmp_path / "src", "3 WaveKin", files={"wave_elevation.txt": "0 0\n1 0.1\n2 0\n3 -0.1\n"}
    )
    model = DeckModel.load(base)
    saved = model.save(tmp_path / "copy" / "deck.dat")
    assert (saved.parent / "wave_elevation.txt").is_file()
    unmoved = model.save(tmp_path / "plain" / "deck.dat", rebase=False)
    assert not (unmoved.parent / "wave_elevation.txt").exists()


def test_rebased_syrope_path_may_contain_a_space(tmp_path):
    source = tmp_path / "src"
    (source / "my data").mkdir(parents=True)
    text = _example(
        "syrope_polyester_mooring.dat",
        ("SYROPE:data/syrope/syrope_settings.dat", "SYROPE:my data/settings.txt"),
    )
    (source / "deck.dat").write_text(text, encoding="utf-8")
    (case,) = generate_deck_cases(source / "deck.dat", tmp_path / "out", {"a": {}})
    rebased = DeckFile.read(case.deck)
    assert rebased.line_types[0].tokens[3].startswith("SYROPE:../src/my data/settings.txt|")


def test_rebased_option_path_with_a_space_fails_closed(tmp_path):
    source = tmp_path / "my src"
    source.mkdir()
    text = _options(_HELD, "motion.txt motionFile")
    (source / "deck.dat").write_text(text, encoding="utf-8")
    with pytest.raises(DeckFormatError, match="which an OPTIONS value cannot hold"):
        generate_deck_cases(source / "deck.dat", tmp_path / "out", {"a": {}})


# --------------------------------------------------------------------------- DeckModel loading


def test_model_loads_a_quoted_line_type_name_with_a_space():
    text = _edit(
        _HELD, ("chain 0.252", '"my chain" 0.252'), ("1 chain 410.0", '1 "my chain" 410.0')
    )
    model = DeckModel.from_text(text)
    assert model.line_types["my chain"].name == "my chain"
    assert '"my chain"' in model.to_text()


def test_external_loads_take_both_column_orders_and_sequential_ids():
    header = "ID Object Force Blin Bquad CSys"
    rows = ("1 Body1 0 1 1 G", "2 Body1 G 1|0|0 0 0", "3 Body1 0|0|5 2 0 L")
    model = DeckModel.from_text(_section(_buoy(), "EXTERNAL LOADS", header, *rows))
    assert [(load.id, load.csys) for load in model.external_loads] == [
        (1, "G"),
        (2, "G"),
        (3, "L"),
    ]
    assert model.external_loads[2].force == (0.0, 0.0, 5.0)
    DeckFile.from_text(model.to_text())
    for bad, message in (
        (("L1 Body1 G 0 1 1",), "ID numbers must be sequential"),
        (("1 Body1 0 1 1 G", "3 Body1 0 1 1 G"), r"expected 2, found '3'"),
        (("1 Body1 G 0 1 L",), "needs the CSys letter"),
        (("1 Body1 0 1 1 -",), "CSys must be G or L for a body"),
    ):
        _rejects(_section(_buoy(), "EXTERNAL LOADS", header, *bad), message)


def test_named_rejections_of_moordyn_only_rows():
    _rejects(_options(_HELD, "0.5 mu_kT"), "CableDyn reads frictionMu")
    _rejects(_options(_HELD, "1.2 StatDynFricScale"), "has no CableDyn equivalent")
    _rejects(_options(_HELD, "airy 2 8 0 waves regular"), "is a positional")
    _rejects(_options(_HELD, "200 WaveComponents", "600 WaveDirections"), "at most 100000")
    DeckFile.from_text(_options(_HELD, "200 WaveComponents", "500 WaveDirections"))
    coupled_pinned = _SPAR_ROD.replace("Free", "CoupledPinned")
    _rejects(_spar((_SPAR_ROD, coupled_pinned)), "CoupledPinned/VesselPinned RODS")
    _rejects(_edit(_HELD, ("0.252 390.0", "0.252 39+1")), "is not a plain number")
    turbine = _edit(_HELD, ("2 Coupled 0.0", "2 Turbine1 0.0"))
    _rejects(turbine, "Turbine<J> points need a TURBINES section")


def test_model_keeps_keyword_end_connection_stiffness():
    source = _example("lozon_gomex80_power_cable.dat")
    marker = next(line for line in source.splitlines() if "OPTIONS" in line and "---" in line)
    text = _edit(
        source,
        (marker, f"{_BAR} END CONNECTIONS {_BAR}\n1 EndA Inf -0.25 0.0 -0.97\n{marker}"),
    )
    model = DeckModel.from_text(text)
    assert model.end_connections[0].stiffness == "Inf"


def test_model_accepts_numpy_numbers_and_reports_argument_errors():
    model = DeckModel.from_text(_HELD)
    model.points[1].x = np.float32(401.5)
    model.add_line_type("wire", diam=np.float64(0.1), mass=np.int64(20), ea=1.0e9)
    assert "401.5" in model.to_text()
    with pytest.raises(TypeError, match="output channel must be a string"):
        model.outputs.remove(1)  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="pass a point id as an integer"):
        model.add_line(9, "2", 1)
    with pytest.raises(ValueError, match="three numbers"):
        model.add_end_connection(1, "A", "Pinned", (0.0, 1.0))


# --------------------------------------------------------------------------- edge rows


@pytest.mark.parametrize(
    ("kind", "resolved"),
    [
        ("free", "free"),
        ("fixed", "fixed"),
        ("pinned", "fixed"),
        ("vessel", "coupled"),
        ("body3pinned", "body3"),
        ("odd", "odd"),
    ],
)
def test_zero_length_rod_connector_takes_the_rod_attachment(monkeypatch, kind, resolved):
    deck = DeckFile.from_text(_spar())
    monkeypatch.setattr(DeckFile, "_rod_rows_by_id", lambda self: {1: (kind, 0)})
    types = {"1": "rod1a", "2": "rod1b"}
    keys, merged = deck._resolve_rod_ends(types, {"1", "2"}, {})
    assert keys == {"rod1a": "1", "rod1b": "1"} and merged == {"2": "1"}
    assert types["1"] == resolved


def test_rod_ids_must_be_integers():
    _rejects(
        _spar((_SPAR_ROD, _SPAR_ROD.replace("1    spar", "x    spar", 1))),
        "rod id must be an integer",
    )


def test_external_load_on_a_point3_body_is_refused_by_name():
    load = ("EXTERNAL LOADS", "ID Object CSys Force Blin Bquad", "1 Body1 G 0 1 1")
    _rejects(_section(_point3_buoy(), *load), "EXTERNAL LOADS apply to Rigid6 bodies only")


@pytest.mark.parametrize(
    ("path", "message"),
    [
        ("", "must not be empty"),
        ("my data/set #1.txt", "comment/header markers"),
        ("a\tb.txt", "quotes or line breaks"),
    ],
)
def test_quoted_syrope_paths_reject_what_a_quote_cannot_hold(path, message):
    from cabledyn.deck_file import _render_quoted_text

    with pytest.raises(ValueError, match=message):
        _render_quoted_text(path)


# --------------------------------------------------------------------------- aggregate routes

_MIXED = _edit(
    _HELD,
    (
        "chain 0.252 390.0 1.674e9 -1.0 0.0 1.37 0.64 1.0 0.0\n",
        "chain 0.252 390.0 1.674e9 -1.0 0.0 1.37 0.64 1.0 0.0\n"
        "cable 0.2 60.0 4.0e8 0.0 1.0e4 1.2 0.1 1.0 0.0\n",
    ),
    (
        "2 Coupled 0.0 0.0 0.0 0 0 0 0",
        "2 Coupled 0.0 0.0 0.0 0 0 0 0\n3 Fixed -400.0 0.0 -50.0 0 0 0 0",
    ),
    ("1 2 1 -", "1 2 1 -\n2 2 3 -"),
    ("1 chain 410.0 41", "1 chain 410.0 41\n2 cable 420.0 42"),
)


def test_standalone_mixed_deck_runs_on_the_still_water_aggregate():
    DeckFile.from_text(_MIXED)
    DeckFile.from_text(_edit(_MIXED, ("2 2 3 -", "2 2 3 r")))
    _rejects(_options(_MIXED, "uniform 0.5 0 0 current"), "runs in still water")
    _rejects(_options(_MIXED, "airy 2 8 0 waves"), "runs in still water")
    _rejects(_edit(_MIXED, ("2 2 3 -", "2 2 3 pt")), "writes no per-line files")
    _rejects(_options(_MIXED, "m.txt motionFile"), "mixed EI=0/finite-EI deck does not support")


def test_coupled_route_takes_no_motion_file_and_needs_a_line():
    _rejects(
        _options(_HELD, "m.txt motionFile"),
        "driven by the host, not a deck motionFile",
        caller_driven=True,
    )
    no_lines = _buoy(
        ("1     1       4       -\n", ""),
        ("2     2       5       -\n", ""),
        ("3     3       6       -\n", ""),
        ("1        poly       86.85    20\n", ""),
        ("2        poly       86.85    20\n", ""),
        ("3        poly       86.85    20\n", ""),
        ('"FairTen1"\n', ""),
        ('"FairTen2"\n', ""),
        ('"FairTen3"\n', ""),
        ('"AnchTen1"\n', ""),
    )
    DeckFile.from_text(no_lines)
    _rejects(no_lines, "needs at least one line", caller_driven=True)


# --------------------------------------------------------------------------- torsion

# A straight finite-EI line clamped at a Coupled and a Fixed point with the optional torsion
# columns of END CONNECTIONS (the pure-torsion deck of test_torsion_deck).
_TORSION = f"""\
{_BAR} CableDyn Input File {_BAR}
in-memory fixture: a straight finite-EI line restrained in torsion at both ends
{_BAR} LINE TYPES {_BAR}
Name Diam Mass EA BA EI GAs GJ Irt Irn Cdn Cdt Can Cat
cab 0.2 32.2 1.0e7 0.0 1.0e6 1.0e8 5.0e4 1.0 1.0 0.0 0.0 0.0 0.0
{_BAR} POINTS {_BAR}
ID Type X Y Z Mass Vol CdA Ca
1 Coupled 0.0 0.0 -50.0 0 0 0 0
2 Fixed 100.0 0.0 -50.0 0 0 0 0
{_BAR} LINES {_BAR}
ID NodeA NodeB Outputs
1 1 2 -
{_BAR} SECTIONS {_BAR}
LineID LineType Length NumSegs
1 cab 100.0 40
{_BAR} END CONNECTIONS {_BAR}
LineID End Stiffness EzX EzY EzZ TorsStiffness NxX NxY NxZ Pretwist
1 A Rigid 1 0 0 Rigid 0 0 1 0
1 B Rigid 1 0 0 Rigid 0 0 1 720
{_BAR} OPTIONS {_BAR}
9.80665 g
1025.0 rhoW
50.0 WtrDpth
0.1 dtM
0.0 TMax
{_BAR} OUTPUTS {_BAR}
Torq1N1
Twist1N41
Twist1
{_BAR} need this line {_BAR}
"""
_TORSION_ROW_A = "1 A Rigid 1 0 0 Rigid 0 0 1 0\n"


def test_torsion_end_connection_columns_are_accepted():
    deck = DeckFile.from_text(_TORSION)
    assert len(deck.end_connections[0].tokens) == 10 + 1
    # ten columns (no pretwist), keywords and a finite stiffness, mixed with six-column rows
    DeckFile.from_text(_edit(_TORSION, (_TORSION_ROW_A, "1 A Rigid 1 0 0 1.0e5 0 0 1\n")))
    DeckFile.from_text(_edit(_TORSION, (_TORSION_ROW_A, "1 A Rigid 1 0 0 Inf 0 1 1 -45\n")))
    six = _edit(
        _TORSION,
        (_TORSION_ROW_A, "1 A Rigid 1 0 0\n"),
        ("Torq1N1\nTwist1N41\nTwist1\n", "FairTen1\n"),
    )
    DeckFile.from_text(six)


@pytest.mark.parametrize(
    ("row", "message"),
    [
        ("1 A Rigid 1 0 0 Rigid 0 0\n", "or 10 or 11 with the torsion columns"),
        ("1 A Rigid 1 0 0 Rigid 0 0 1 0 0\n", "or 10 or 11 with the torsion columns"),
        ("1 A Rigid 1 0 0 Stiff 0 0 1 0\n", "torsional stiffness must be finite and non-negative"),
        ("1 A Rigid 1 0 0 -1.0 0 0 1 0\n", "torsional stiffness must be finite and non-negative"),
        ("1 A Rigid 1 0 0 Rigid 0 0 0 0\n", r"reference normal \(NxX NxY NxZ\) must be non-zero"),
        ("1 A Rigid 1 0 0 Rigid 2 0 0 0\n", "must not be parallel to the direction Ez"),
        ("1 A Rigid 1 0 0 Rigid 0 x 1 0\n", "reference normal column 9 must be numeric"),
        ("1 A Rigid 1 0 0 Rigid 0 0 1 nan\n", "is not a finite number"),
    ],
)
def test_torsion_end_connection_columns_follow_the_native_rules(row, message):
    _rejects(_edit(_TORSION, (_TORSION_ROW_A, row)), message)


@pytest.mark.parametrize(
    ("row", "message"),
    [
        (
            "1 A Rigid 1 0 0 Stiff 0 0 1 0\n",
            r"line 1 End A: END CONNECTIONS torsional stiffness .*, got 'Stiff'",
        ),
        ("1 A Rigid 1 0 0 -1.0 0 0 1 0\n", r"line 1 End A: .*, got '-1\.0'"),
        ("1 A -5 1 0 0 Rigid 0 0 1 0\n", r"line 1 End A: END CONNECTIONS stiffness .*, got '-5'"),
        (
            "1 A Rigid 0 0 0 Rigid 0 0 1 0\n",
            r"line 1 End A: .*direction must be non-zero, got \(0, 0, 0\)",
        ),
        ("1 A Rigid 1 0 0 Rigid 0 0 0 0\n", r"line 1 End A: .*must be non-zero, got \(0, 0, 0\)"),
        (
            "1 A Rigid 1 0 0 Rigid 2 0 0 0\n",
            r"line 1 End A: .*\(2, 0, 0\) must not be parallel to the direction Ez \(1, 0, 0\)",
        ),
        ("1 C Rigid 1 0 0 Rigid 0 0 1 0\n", r"line 1: END CONNECTIONS End must be A or B, got 'C'"),
        ("1 B Rigid 1 0 0 Rigid 0 0 1 0\n", r"line 1 End B: duplicate END CONNECTIONS row"),
        ("7 A Rigid 1 0 0 Rigid 0 0 1 0\n", r"references an undefined line 7"),
        ("0 A Rigid 1 0 0 Rigid 0 0 1 0\n", r"LineID must be positive, got 0"),
        ("x A Rigid 1 0 0 Rigid 0 0 1 0\n", r"LineID must be an integer, got 'x'"),
    ],
)
def test_end_connection_errors_name_the_line_end_and_value(row, message):
    _rejects(_edit(_TORSION, (_TORSION_ROW_A, row)), message)


def test_torsion_line_rules_follow_the_native_checks():
    # torsion has no EI/1.3 default: the 10-column LINE TYPES row has no GJ
    no_gj = _edit(
        _TORSION,
        (
            "Name Diam Mass EA BA EI GAs GJ Irt Irn Cdn Cdt Can Cat\n",
            "Name Diam Mass EA BA EI Cdn Cdt Can Cat\n",
        ),
        (
            "cab 0.2 32.2 1.0e7 0.0 1.0e6 1.0e8 5.0e4 1.0 1.0 0.0 0.0 0.0 0.0",
            "cab 0.2 32.2 1.0e7 0.0 1.0e6 0 0 0 0",
        ),
    )
    _rejects(no_gj, "must give an explicit GJ > 0")
    ei0 = _edit(
        _TORSION,
        ("1.0e7 0.0 1.0e6 1.0e8", "1.0e7 0.0 0.0 1.0e8"),
        ("1 A Rigid 1 0 0 Rigid", "1 A Pinned 1 0 0 Rigid"),
        ("1 B Rigid 1 0 0 Rigid", "1 B Pinned 1 0 0 Rigid"),
        ("Torq1N1\nTwist1N41\nTwist1\n", "FairTen1\n"),
    )
    _rejects(ei0, "require a finite-EI line; line 1 has EI = 0")
    # one restrained end is accepted (it carries no torque), but takes no torsion channel
    one_end = _edit(_TORSION, (_TORSION_ROW_A, "1 A Rigid 1 0 0 Free 0 0 1 0\n"))
    _rejects(one_end, "is not torsionally restrained at both ends")
    DeckFile.from_text(_edit(one_end, ("Torq1N1\nTwist1N41\nTwist1\n", "FairTen1\n")))
    # the dynamic run needs the force blend; modal analysis and the coupled routes refuse torsion
    dynamic = _edit(_TORSION, ("0.0 TMax\n", "10.0 TMax\n"))
    DeckFile.from_text(dynamic)
    _rejects(
        _options(dynamic, "False alpha_force_blend"), "needs the force-blended generalised-alpha"
    )
    _rejects(_options(_TORSION, "3 nModes"), "modal analysis .* with torsion is not yet supported")
    _rejects(_TORSION, "nor in coupled OpenFAST or FAST.Farm runs", caller_driven=True)
    # the coupled routes refuse one restrained end too (any torsion column), after the line rules
    _rejects(one_end, "nor in coupled OpenFAST or FAST.Farm runs", caller_driven=True)
    _rejects(ei0, "require a finite-EI line; line 1 has EI = 0", caller_driven=True)


# The pure-torsion deck dynamic (TMax > 0), for the rules that need a dynamic run.
_TORSION_DYN = _edit(_TORSION, ("0.0 TMax", "1.0 TMax"))
_NO_TORSION_OUTPUTS = ("Torq1N1\nTwist1N41\nTwist1\n", "FairTen1\n")


def _torsion_body(body_row: str) -> str:
    return _edit(
        _TORSION_DYN,
        (
            f"{_BAR} POINTS {_BAR}",
            f"{_BAR} BODIES {_BAR}\n"
            "ID Type X Y Z Roll Pitch Yaw Mass Vol C33 C44 C55 CdA Ca Ixx Iyy Izz\n"
            f"{body_row}\n{_BAR} POINTS {_BAR}",
        ),
        ("1 Coupled 0.0 0.0 -50.0", "1 Body1 0.0 0.0 0.0"),
    )


def _torsion_two_moving(row_a: str) -> str:
    return _edit(
        _TORSION_DYN,
        ("2 Fixed 100.0 0.0 -50.0 0 0 0 0", "2 Free 100.0 0.0 -50.0 100 0 0 0"),
        (_TORSION_ROW_A, row_a),
        ("1 B Rigid 1 0 0 Rigid", "1 B Pinned 1 0 0 Rigid"),
    )


def test_torsion_scope_rules_follow_the_native_checks():
    # End A on a Rigid6 body is supported, on a Point3 body or a rod end it stops by name
    DeckFile.from_text(
        _torsion_body("1 Rigid6 0 0 -50 0 0 0 1.0e5 97.56 0 0 0 0 0 1.0e3 1.0e3 1.0e3")
    )
    _rejects(
        _torsion_body("1 Point3 0 0 -50 0 0 0 1.0e5 97.56 0 0 0 0 0"),
        "a torsional END CONNECTION on a body needs a Rigid6 body",
    )
    rod_end = _edit(
        _TORSION_DYN,
        (
            f"{_BAR} POINTS {_BAR}",
            f"{_BAR} ROD TYPES {_BAR}\nName Diam Mass Cd Ca CdEnd CaEnd\n"
            "spar 1.0 300.0 0.8 1.0 0.0 0.0\n"
            f"{_BAR} RODS {_BAR}\nID RodType Type XA YA ZA XB YB ZB NumSegs Outputs\n"
            f"1 spar Fixed 0 0 -50 0 0 -40 2 -\n{_BAR} POINTS {_BAR}",
        ),
        ("1 1 2 -", "1 R1A 2 -"),
        ("1 Coupled 0.0 0.0 -50.0 0 0 0 0\n", ""),
    )
    _rejects(rod_end, "a torsional END CONNECTION on a rod end is not supported")
    attachments = _edit(
        _TORSION,
        (
            f"{_BAR} OPTIONS {_BAR}",
            f"{_BAR} ATTACHMENTS {_BAR}\nLineID ArcLength Mass Volume CdA Ca\n"
            f"1 50.0 10.0 0.0 0.0 0.0\n{_BAR} OPTIONS {_BAR}",
        ),
    )
    _rejects(attachments, "torsion is not combined with ATTACHMENTS")
    # restrained at both ends, End B must be Fixed
    _rejects(
        _torsion_two_moving("1 A Pinned 1 0 0 Rigid 0 0 1 0\n"),
        "torsional END CONNECTIONS require a finite-EI line with Fixed End B; line 1 has two "
        "moving ends",
    )
    # one restrained end only, or none, is ignored: the same verdict as the six-column rows
    one_end = _edit(_torsion_two_moving("1 A Pinned 1 0 0 Free 0 0 1 0\n"), _NO_TORSION_OUTPUTS)
    six = _edit(
        _torsion_two_moving("1 A Pinned 1 0 0\n"),
        ("1 B Pinned 1 0 0 Rigid 0 0 1 720", "1 B Pinned 1 0 0"),
        _NO_TORSION_OUTPUTS,
    )
    DeckFile.from_text(six)
    DeckFile.from_text(one_end)
    both_free = _edit(
        _TORSION,
        (_TORSION_ROW_A, "1 A Rigid 1 0 0 Free 0 0 1 0\n"),
        ("1 B Rigid 1 0 0 Rigid 0 0 1 720", "1 B Rigid 1 0 0 Free 0 0 1 720"),
    )
    _rejects(both_free, "is not torsionally restrained at both ends")
    DeckFile.from_text(_edit(both_free, _NO_TORSION_OUTPUTS))


def test_torsion_on_a_standalone_mixed_deck_is_refused():
    # a mixed EI = 0 + finite-EI deck without bodies runs on the aggregate, which has no torsion
    gj_cable = _edit(
        _MIXED,
        (
            "Name Diam Mass EA BA EI Cdn Cdt Can Cat\n",
            "Name Diam Mass EA BA EI GAs GJ Irt Irn Cdn Cdt Can Cat\n",
        ),
        (
            "chain 0.252 390.0 1.674e9 -1.0 0.0 1.37 0.64 1.0 0.0\n",
            "chain 0.252 390.0 1.674e9 -1.0 0.0 0 0 0 0 1.37 0.64 1.0 0.0\n",
        ),
        (
            "cable 0.2 60.0 4.0e8 0.0 1.0e4 1.2 0.1 1.0 0.0\n",
            "cable 0.2 60.0 4.0e8 0.0 1.0e4 1.0e8 1.0e4 1.0 1.0 1.2 0.1 1.0 0.0\n",
        ),
    )
    DeckFile.from_text(gj_cable)
    for rows in (
        ("2 A Rigid -1 0 0 Rigid 0 0 1 0", "2 B Rigid -1 0 0 Rigid 0 0 1 90"),
        ("2 A Rigid -1 0 0 Rigid 0 0 1 0", "2 B Rigid -1 0 0 Free 0 0 1 0"),
    ):
        text = _section(
            gj_cable,
            "END CONNECTIONS",
            "LineID End Stiffness EzX EzY EzZ TorsStiffness NxX NxY NxZ Pretwist",
            *rows,
        )
        _rejects(text, "mixes EI = 0 and finite-EI lines without a BODY")
        _rejects(_options(text, "10.0 TMax"), "mixes EI = 0 and finite-EI lines without a BODY")


@pytest.mark.parametrize(
    ("stiffness", "accepted"),
    [
        ("Inf", True),
        ("Infinity", True),
        ("RIGID", True),
        ('"Rigid"', True),
        ("Zero", True),
        ("-0", True),
        ("2.5d4", True),
        ("+Inf", False),
        ("1e-320", False),
        ("1.0q5", False),
        ("Pinned", False),
    ],
)
def test_torsion_stiffness_keywords_and_numbers_follow_the_native_reader(stiffness, accepted):
    text = _edit(_TORSION, (_TORSION_ROW_A, f"1 A Rigid 1 0 0 {stiffness} 0 1 1 -45\n"))
    if stiffness.lower() in {"zero", "-0"}:
        text = _edit(text, _NO_TORSION_OUTPUTS)
    if accepted:
        DeckFile.from_text(text)
    else:
        _rejects(text, "torsional stiffness must be finite and non-negative, Free, or Rigid")


@pytest.mark.parametrize("stiffness", ["1.0q5", "1e-320", "+Inf", "NaN"])
def test_bending_stiffness_numbers_follow_the_native_reader(stiffness):
    _rejects(
        _edit(_TORSION, (_TORSION_ROW_A, f"1 A {stiffness} 1 0 0 Rigid 0 1 1 -45\n")),
        "END CONNECTIONS stiffness must be finite and non-negative, Pinned, or Rigid",
    )


def test_builders_write_or_refuse_the_torsion_columns():
    from cabledyn.deck import DeckWriter

    # DeckModel writes a 10-column row (no pretwist) that reads back the same
    model = DeckModel.from_text(_edit(_TORSION, (_TORSION_ROW_A, ""), _NO_TORSION_OUTPUTS))
    model.add_end_connection(
        1, "A", "Rigid", (1.0, 0.0, 0.0), torsion_stiffness=2.5e4, normal=(0, 1, 1)
    )
    text = model.to_text()
    row = next(r for r in DeckFile.from_text(text).end_connections if r.tokens[1] == "A")
    assert len(row.tokens) == 10 and float(row.tokens[6]) == 2.5e4
    # DeckFile edits the torsion columns of a row in place, keywords unquoted
    deck = DeckFile.from_text(_TORSION)
    deck.set_end_connection(1, "B", torsstiffness="Infinity", pretwist=90.0)
    row_b = next(r for r in deck.end_connections if r.tokens[1] == "B")
    assert row_b.tokens[6] == "Infinity" and float(row_b.tokens[10]) == 90.0
    # the simple writer has no torsion columns and says so
    with pytest.raises(ValueError, match="no torsion columns"):
        DeckWriter().add_end_connection(1, "A", "Rigid", (1.0, 0.0, 0.0), torsion_stiffness="Rigid")


def test_a_free_torsion_end_takes_any_reference_normal():
    for normal in ("0 0 0", "1 0 0"):
        DeckFile.from_text(
            _edit(
                _TORSION,
                (_TORSION_ROW_A, f"1 A Rigid 1 0 0 Free {normal} 0\n"),
                _NO_TORSION_OUTPUTS,
            )
        )


_PARITY_DRIVER = (
    os.environ.get("CABLEDYN_TEST_DRIVER", "").strip()
    or os.environ.get("CABLEDYN_DRIVER", "").strip()
)


@pytest.mark.skipif(
    not _PARITY_DRIVER or not Path(_PARITY_DRIVER).is_file(),
    reason="set CABLEDYN_TEST_DRIVER to a built native driver to run the native torsion parity",
)
@pytest.mark.parametrize(
    "row_a",
    [
        "1 A Rigid 1 0 0 Inf 0 1 1 -45",
        "1 A Rigid 1 0 0 Infinity 0 1 1 -45",
        '1 A Rigid 1 0 0 "Rigid" 0 1 1 -45',
        "1 A Rigid 1 0 0 +Inf 0 1 1 -45",
        "1 A Rigid 1 0 0 1e-320 0 1 1 -45",
        "1 A Rigid 1 0 0 1.0q5 0 1 1 -45",
        "1 A 1.0q5 1 0 0 Rigid 0 1 1 -45",
        "1 A Rigid 1 0 0 Rigid 0 0 1",
        "1 A Rigid 1 0 0 Pinned 0 0 1 0",
        "1 A Rigid 1 0 0 Zero 0 0 1 0",
        "1 A Rigid 1 0 0 2.5d4 +0 .5 1 1d2",
        "1 A Rigid 1 0 0 Rigid 0 0 0 0",
        "1 A Rigid 1 0 0 Rigid 1 0 0.0009 0",
        "1 A Rigid 1 0 0 Rigid 1 0 0.0011 0",
        "1 A Rigid 1 0 0 Rigid 0 0 1 nan",
        "1 A Rigid 1 0 0 Rigid 0 0 1 1e400",
        "1 A Rigid 1 0 0 Free 0 0 0 0",
        "1 A Rigid 1 0 0 Free 1 0 0 -7200",
        "1 A Rigid 1 0 0 Free 0 nan 1 0",
        "1 A Rigid 1 0 0 Free 0 0 1 inf",
    ],
)
def test_torsion_rows_get_the_native_verdict(row_a, tmp_path):
    text = _edit(_TORSION, (_TORSION_ROW_A, f"{row_a}\n"))
    if row_a.split()[6].lower() in {"free", "zero"}:
        # one restrained end carries no torque: judge the row, not the torsion channels
        text = _edit(text, _NO_TORSION_OUTPUTS)
    try:
        DeckFile.from_text(text)
        python_accepts = True
    except DeckFormatError:
        python_accepts = False
    deck = tmp_path / "torsion.dat"
    deck.write_text(text, encoding="utf-8")
    completed = subprocess.run(
        [_PARITY_DRIVER, str(deck), str(tmp_path / "torsion")],
        capture_output=True,
        text=True,
        cwd=tmp_path,
        timeout=600,
        check=False,
    )
    # exit code 1 is an input refusal; 0 a run and 2 a solve failure after validation
    assert python_accepts == (completed.returncode != 1), completed.stdout + completed.stderr


def _native_exit(text: str, tmp_path: Path) -> int:
    deck = tmp_path / "torsion.dat"
    deck.write_text(text, encoding="utf-8")
    completed = subprocess.run(
        [_PARITY_DRIVER, str(deck), str(tmp_path / "torsion")],
        capture_output=True,
        text=True,
        cwd=tmp_path,
        timeout=600,
        check=False,
    )
    return completed.returncode


@pytest.mark.skipif(
    not _PARITY_DRIVER or not Path(_PARITY_DRIVER).is_file(),
    reason="set CABLEDYN_TEST_DRIVER to a built native driver to run the native torsion parity",
)
@pytest.mark.parametrize(
    "channel",
    [
        "TWIST01",
        "torq01n041",
        "Twist1N41",
        "Twist1N42",
        "Torq1N0",
        "Twist1N0",
        "Torq1",
        "Twist",
        "Twist1N",
        "Torq1N1x",
        "Twist1x",
        "Twist2",
    ],
)
def test_torsion_channels_get_the_native_verdict(channel, tmp_path):
    text = _edit(_TORSION, ("Twist1\n", f"{channel}\n"))
    try:
        DeckFile.from_text(text)
        python_accepts = True
    except DeckFormatError:
        python_accepts = False
    # exit code 1 is an input refusal; 0 a run and 2 a solve failure after validation
    assert python_accepts == (_native_exit(text, tmp_path) != 1)


@pytest.mark.parametrize(
    ("channel", "message"),
    [
        ("Torq1N42", "node exceeds the line node count"),
        ("Torq1N0", "node numbers start at 1"),
        ("Torq2N1", "bad torsion channel"),
        ("Twist0", "bad torsion channel"),
        ("TwistN1", "bad torsion channel"),
    ],
)
def test_torsion_channels_follow_the_native_rules(channel, message):
    _rejects(_edit(_TORSION, ("Twist1\n", f"Twist1\n{channel}\n")), message)


def test_torsion_channels_are_one_identity_per_quantity():
    _rejects(_edit(_TORSION, ("Twist1\n", "Twist1\nTORQ01N1\n")), "duplicate")
    DeckFile.from_text(_edit(_TORSION, ("Twist1\n", "Twist1\nTwist1N1\nTorq1N41\n")))


def test_model_round_trips_the_torsion_columns():
    model = DeckModel.from_text(_TORSION)
    row_a, row_b = model.end_connections
    assert (row_a.torsion_stiffness, row_a.normal, row_a.pretwist) == (
        "Rigid",
        (0.0, 0.0, 1.0),
        0.0,
    )
    assert row_b.pretwist == 720.0
    text = model.to_text()
    assert "TorsStiffness" in text and "Pretwist" in text
    again = DeckModel.from_text(text)
    assert [(r.torsion_stiffness, r.normal, r.pretwist) for r in again.end_connections] == [
        (r.torsion_stiffness, r.normal, r.pretwist) for r in model.end_connections
    ]
    DeckFile.from_text(text)
    # a six-column row keeps six columns; the builder checks the torsion arguments
    plain = DeckModel.from_text(
        _edit(
            _TORSION,
            (_TORSION_ROW_A, "1 A Rigid 1 0 0\n"),
            ("1 B Rigid 1 0 0 Rigid 0 0 1 720\n", ""),
            ("Torq1N1\nTwist1N41\nTwist1\n", "FairTen1\n"),
        )
    )
    assert "TorsStiffness" not in plain.to_text()
    with pytest.raises(ValueError, match="must be given together"):
        plain.add_end_connection(1, "B", "Rigid", (1.0, 0.0, 0.0), torsion_stiffness="Rigid")
    with pytest.raises(ValueError, match="pretwist needs"):
        plain.add_end_connection(1, "B", "Rigid", (1.0, 0.0, 0.0), pretwist=10.0)
    with pytest.raises(ValueError, match="normal must be three numbers"):
        plain.add_end_connection(
            1, "B", "Rigid", (1.0, 0.0, 0.0), torsion_stiffness=1.0e4, normal=(0.0, 1.0)
        )
    added = plain.add_end_connection(
        1,
        "B",
        "Rigid",
        (1.0, 0.0, 0.0),
        torsion_stiffness=1.0e4,
        normal=(0.0, 0.0, 1.0),
        pretwist=90.0,
    )
    assert added.torsion_stiffness == 1.0e4
    assert "TorsStiffness" in plain.to_text()
