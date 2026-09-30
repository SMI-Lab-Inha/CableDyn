# SPDX-License-Identifier: Apache-2.0
"""Deck-validation rules that mirror the native reader and the standalone driver.

Each test pins one rule of ``src/CableDyn_DeckDriver.f90`` so that
:class:`cabledyn.DeckFile` neither rejects a deck the solver accepts nor
certifies one the solver refuses.
"""

from __future__ import annotations

import hashlib
import json
import re
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
