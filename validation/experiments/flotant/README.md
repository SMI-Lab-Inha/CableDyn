<!-- SPDX-License-Identifier: Apache-2.0 -->

# FLOTANT OD 48 mm research cable: VIV in tow tests

This directory compares CableDyn with the tests of a power-cable specimen at MARIN for the EU
project FLOTANT: de Wilde, Heijnsdijk & de Ridder (2021), Deliverable D5.5, MARIN report
30995-53-DWB, a public deliverable without a DOI,
<https://flotantproject.eu/wp-content/uploads/2021/11/210430_FLT-D5.5_VIV_test_v0.pdf>. A 6.4 m
cable of 48 mm outer diameter, a simplified 1:3 model of a 66 kV inter-array cable, hung
vertically from the towing carriage with a 200 kg clump weight below it and was towed at 40
speeds from 0.025 to 1.0 m/s. The report gives the cross-flow and in-line A/D_rms at three
stations against the reduced velocity. The same specimen was forced by a hexapod in air and in
water, and those RAOs are used here to set up and check the structural model.

The tow tests are predicted with three published VIV coefficient sets, none of them fitted to
these tests. The structure is identified from the forced-oscillation tests of the same specimen
only, never from the tow data.

In short: on the identified structure the set of Ulveseter et al. (2019), fitted by them to
another riser, gives the height of the first cross-flow peak at mid-length (0.755 against
0.776), but it responds from a lower speed, holds the lock-in 1.2 Ur longer than the test and
stays high over the upper branch. The generic set is 38 % low at the peak and releases at the
measured speed; the wake oscillator is 62 % low. The r.m.s. log errors over the three stations
are 0.73, 0.59 and 1.21, large for all three sets.

## The test (D5.5)

- **Specimen (Table 2-1):** OD 48 mm, 6.4 m, 2.8 kg/m; 9.7 N/m in water (Table 2-2).
- **Set-up in water (Secs. 2.4, 2.5):** from the top, the carriage (or the hexapod), a
  six-component force frame, a Cardan joint, the cable, a Cardan joint, a swivel, a load cell
  and the clump weight. Tension 1970 N at the top and 1436 N at the bottom (Table 6-1).
- **Instrumentation (Sec. 3):** accelerometers at 1/2L, 3/4L and L from the top of the cable
  (L is its lower end). The motions are the accelerations integrated twice, high-passed at
  3.5 rad/s for the tow tests (Table 3-4).
- **A/D (Sec. 3.4.4):** Appendix 3 plots the standard-deviation A/D, A/D_rms = σ/OD (p. 53:
  "The presented results are for standard deviation A/Drms"); the nominal A/D is √2 σ/OD. The
  first peak is A/D_rms 0.78 (nominal 1.10) at Ur 5.7.
- **Ur (Secs. 3.4.3, 6.3):** Ur = V/(OD f1) with f1 = 1.77 Hz.
- **Tow tests (Table 6-2):** one test per speed, 0.025 m/s apart, the towing direction
  alternating from 0.15 m/s on (forwards at odd multiples of 0.025 m/s from there, backwards at
  even ones). Each test is analysed over its steady part, 37.5 s at 1.0 m/s up to 580 s at
  0.025 m/s.

### Caveats of the data

- **The Ur basis.** Ur uses f1 = 1.77 Hz, while Sec. 5.3 calls the 1.3 Hz peak of the in-water
  RAOs the first mode. The origin of 1.77 Hz is not given. On a 1.3 Hz basis every Ur would be
  1.36 times larger. The comparison here is made at the tested speeds, so the basis only
  labels the axis; all Ur values here are on the report's 1.77 Hz basis.
- **The bottom tension.** Table 6-1 gives 1436 N at the bottom, while the in-water FEM table of
  Sec. 5.4 gives 1687–1702 N along the chain above the clump, whose submerged weight it gives as
  1710.5 N. The cable's own wet weight (about 62 N) cannot explain a 534 N difference
  between top and bottom. The model takes the FEM table; its static top tension is in the
  results.
- **Towing direction.** Speed and direction are confounded (one test per speed). Between Ur 5.8
  and 9.2 the backward tests are mostly lower than the forward tests at the neighbouring speeds
  (at 1/2L, 0.30–0.46 against 0.42–0.56). The report attributes small differences to the coil
  curvature the specimen kept. The model is symmetric; the statistics are given over all tests
  and per direction.
- **Displacements from accelerometers.** The accelerometers tilt with the cable, so the
  derived displacement holds a gravity term g θ/ω² (θ the slope of the cable), which the
  report does not remove. The model output is processed the same way (below); without the term
  the cross-flow A/D_rms of the three sets changes by up to 5 % at 1/2L, 8 % at 3/4L and
  25 % at L for Ur ≥ 3, and the r.m.s. log errors by at most 0.07.
- **Not reported:** response frequencies (only the acceleration orbits of one test, Fig. 6-1),
  repeat tests and uncertainties.

## Model

### Structure

The line follows the report's node-by-node FEM input tables (Sec. 4.5 in air, Sec. 5.4 in
water; `fem_air.csv`, `fem_water.csv`): 101 nodes 0.0765 m apart, with the lumped masses, the
weights and the static tensions. Each table element is one CableDyn element (`--refine N`
splits it into N), grouped into runs by the label of its lower node: steel links (EI 1e6 N m²,
rigid next to the cable), Cardan joints and the wire (EI 1 N m², hinges) and the cable.

- **Water:** fresh, 1000 kg/m³. The report does not state it; the 2.8 kg/m and 9.7 N/m of
  Table 2-2 at OD 48 mm imply 1000 kg/m³. The cable's diameter is the OD; the other parts take
  the diameter of the displaced volume of their weight in the table.
- **Axial compliance:** the counterweight (air) or clump (water) keeps the tension nearly
  constant. The model uses an axial spring k_ax = 20 kN/m in series with the table's EA, as
  EA_eff = L/(1/k_ax + L/EA), with the unstretched lengths set so that the static tensions are
  the table's.

What the tests do not give, and how it was identified:

- **Cable bending stiffness and damping, from the in-air RAOs (Appendix 1, 200 kg).** EI is not
  measured: the tables' EI = 3.0 kN m² reproduces the in-air resonance (10.7 rad/s) only with
  an empirical factor λ = 1.10 that MARIN introduced in the tensioned-beam formula to match it
  (Sec. 4.5). The 2.8 % damping ratio of Sec. 4.4 is a single-degree-of-freedom fit, with an
  assumed participation factor of 0.5, to one peak. Here a linear section was identified on a grid
  of EI and the bending damping Bb (the `EI|Bb` form of the LINE TYPES EI column) against the
  surge RAOs at 5 and 10 mm: EI near 1500 N m² puts the resonance at the measured frequency
  (3000 puts it 5–6 % high), and the 10 mm peak needs a larger Bb than the 5 mm one (about 160
  against 110 N m² s): the peak falls with amplitude faster than a linear section can follow,
  as the report notes. The model takes EI 1650 N m² and Bb 130 N m² s, between the two. The ends of the in-air rig are Cardan
  joints, a swivel and a pulley, whose friction is in the measured loss: Bb is an upper bound
  for the cable's own damping.
- **The lower end in water, from the in-water RAOs (Appendix 2, 200 kg, 5 mm, 1.0–3.0 Hz).**
  The measured RAO has two peaks at 1/2L, at 1.3 and 1.8 Hz, with a notch at 1.5–1.6 Hz;
  at 1.8 Hz station L moves little (surge 0.55 against 4.40 at 1/2L). MARIN's FEM table fixes the clump node and makes
  the chain rigid (`--end table`); that model, and the same with the clump as a free point mass
  (`--end free`, a `Free` End B), have one peak near 1.2–1.3 Hz. The model used here
  (`--end identified`) has two elements below the cable that the tests point to, identified on
  grids against the surge and sway RAOs over 1.0–3.0 Hz (run with sea-water density; the check
  below, at 1000 kg/m³, gives the scores of this model):
  - the link chain and its load cell bend with EI 100 N m² (best on surge; 250 N m² on
    sway), neither a hinge (EI 1, which gives a spurious mode with large motion at L at
    2.6–2.7 Hz) nor rigid;
  - the clump (200 kg, its displaced volume, C_d A 0.05 m² assumed, C_a 1) is a point load on
    a short rigid segment held laterally by a spring K = 32.5 kN/m and a damper 20 N s/m (the
    best K on sway, the cross-flow direction of the tow tests; 30 kN/m on surge), modelled as a
    clamped cantilever post, K = 3 EI_p/l³ (l = 0.3 m), from a still anchor.

  This end is phenomenological: an impedance that reproduces the measured split of the first
  mode, not a measured property of the rig. The damper hardly matters (1–80 N s/m changed the
  RAO score by 0.004 in the identification grid). The clump's geometry and the chain are not described in the report.

The RAO check (below) runs the model as identified. The line is axisymmetric, so one
forcing direction stands for surge and sway.

### Hydrodynamics

Uniform current along x over the whole line (the carriage speed), cross-flow (CF) along y.
C_D 1.2 and C_a 1 on the cable and the hardware. Three published VIV coefficient sets on the
cable, `--sets`:

| Set | Source | Values | Fitted to |
| --- | --- | --- | --- |
| `uls19` | Ulveseter et al. (2019), Table 1 and Sec. 4.1 | C_v,y 0.85, f̂0 0.144 ± 0.064, in-line C_v,x 0.75, α 0.15, C_hh 1.5, C_D 1.2, C_M 2 | the NDP 38 m riser in uniform flow: blind here |
| `b1` | The generic set of the Delta Flume riser case | C_l 0.8, f̂0 0.18 ± 0.08; in-line C_x 0.85, α 0.15 (Ulveseter et al. 2018, Table 3); `VIV AMPLITUDE` C_v 0.8 with the drag amplification of Chaplin et al. (2005, Eq. 4); `VIV DAMPING` `0 0 1` | nothing from these tests: blind |
| `wake` | Facchinetti et al. (2004) | `WAKE OSCILLATOR` C_L0 0.3, St 0.2, ε 0.3, A 12; Morison drag | rigid-cylinder lift and lock-in data: blind |

The values were checked against the source articles.

### Run and processing

- dt 2 ms, 60 s, statistics over t > 20 s (frequency resolution 0.025 Hz). The tests' steady
  windows are longer (37.5–580 s, Table 6-2).
- **As the report processes the measurement:** each station's displacement plus the gravity
  term of its tilted accelerometer, g θ/ω² (θ from the two neighbouring nodes), high-passed at
  3.5 rad/s, then A/D_rms = σ/OD. `response_comparison.csv` also holds the plain displacement
  (`..._displacement`).
- Ur = V/(OD f1) with f1 = 1.77 Hz.
- The dominant CF frequency at 1/2L is the spectral peak of the high-passed displacement.
- **Statistics:** the r.m.s. log error, √(mean (ln ratio)²), over the tests from the onset the
  report states, Ur ≥ 3 (30 speeds, 0.275–1.0 m/s); CF and in line, per station and over all
  three. The peak is the largest CF A/D_rms at 1/2L; the release is the first tested Ur above
  the peak where it falls below 70 % of the peak.
- **RAOs:** a harmonic at the forcing frequency fitted over the last 10 s of each forced run
  (5 periods of cosine ramp, then 25 s), of the displacement plus the tilt term, over the
  hexapod amplitude; scored as the r.m.s. log error over the stations, in air over the resonance
  band 1.4–2.0 Hz (1/8L, 1/4L, 1/2L, 3/4L from the bottom), in water over 1.0–3.0 Hz (1/2L,
  3/4L, L from the top).

### Measurements

`extract_d55.py` reads the report's PDF (not distributed; PyMuPDF, NumPy, SciPy and Pillow):

- `fem_air.csv`, `fem_water.csv`: the FEM input tables of Secs. 4.5 and 5.4, from the text.
- `measured_rao.csv`: the 5 and 10 mm RAOs of Appendices 1 and 2 (200 kg, surge and sway).
  The charts are vector graphics: each series is read as plotted coordinates through the tick
  labels (axis residual at most 1.2e-3); the frequencies fall on the 0.1 Hz test grid. The
  axis is labelled "[rad]" and is in rad/s.
- `measured_tow_ad.csv`: Appendix 3, A/D_rms against Ur at the three stations, in line (X) and
  cross-flow (Y), 40 points per chart (39 for Y at L, where an edge marker is clipped). The
  charts are raster images (1386 × 867): the axes are calibrated on the grey grid lines (Ur
  residual 0.004) and the markers found by colour. Every point maps to a tested speed,
  V = Ur OD f1, within 0.0085 m/s of the 0.025 m/s grid; the CF peak at 1/2L reads 0.776 at Ur
  5.58 (the report: 0.78 at 5.7).

The data are digitised from the report; its figures are not reproduced.

## Reproduce

From the repository root, with `path/to/cabledyn.exe` the release executable or a local build:

```powershell
python validation/experiments/flotant/flotant_viv.py `
  path/to/cabledyn.exe validation/experiments/flotant/work --jobs 24 `
  --summary validation/experiments/flotant
```

The default is the committed record: the three sets at the 40 tested speeds (120 runs) and the
90 forced-oscillation runs, about 15 minutes with 24 parallel runs on a 24-core workstation;
each run reproduces its row of the record. `--parts tow` or `--parts rao`, `--sets` and
`--speeds` (comma lists; quote them in PowerShell, `--speeds "0.45,0.475"`) select parts;
`--end` changes the lower end of the tow model; `--dt`, `--refine` and `--tmax` the time step,
the mesh and the run length. A run whose folder holds the same input files and executable is not repeated; `--keep`
keeps the output tables (about 14 MB a tow run).

## Results

`response_comparison.csv` holds every tow run (all stations, both components, the measured
value at the same speed, the response frequency and the mean top tension),
`rao_comparison.csv` every forced run, `summary.json` the statistics and the SHA-256 of the
executable, and `provenance.json` the sources, assumptions and sensitivity runs.

### Structure: the forced-oscillation RAOs

r.m.s. log error of the RAOs over the stations and the band; peaks of the RAO at 1/2L (Hz:
RAO).

| Case | Band | Surge | Sway | CableDyn peaks | Measured, surge | Measured, sway |
| --- | --- | --- | --- | --- | --- | --- |
| Air, 5 mm | 1.4–2.0 Hz | 0.227 | 0.278 | 1.7: 7.39 | 1.7: 8.68 | 1.7: 9.14 |
| Air, 10 mm | 1.4–2.0 Hz | 0.209 | 0.230 | 1.7: 7.23 | 1.7: 5.88 | 1.7: 6.40 |
| Water, 5 mm, MARIN's table end (`table`) | 1.0–3.0 Hz | 1.186 | 1.216 | 1.3: 3.94; 2.5: 1.36 | 1.3: 4.25; 1.8: 4.40 | 1.3: 4.31; 1.9: 2.99 |
| Water, 5 mm, free clump (`free`) | 1.0–3.0 Hz | 1.174 | 1.208 | 1.4: 3.68; 2.6: 1.31 | as above | as above |
| Water, 5 mm, identified end (`identified`) | 1.0–3.0 Hz | **0.459** | **0.743** | 1.2: 3.76; 1.9: 2.45 | as above | as above |

- **In air** the resonance is at the measured 1.7 Hz. With one linear damping the peak is 15 %
  below the 5 mm test and 23 % above the 10 mm test, which brackets the amplitude dependence the
  test shows.
- **In water** only the identified end gives the two peaks and the notch between them. The upper
  peak is low (2.45 against 4.40 in surge; the test frequencies are 0.1 Hz apart about a sharp
  peak), the lower one is 0.1 Hz low, and station L moves too much at the upper peak (1.18
  against 0.55 at 1.8 Hz in surge). The sway fit is worse (0.743), since K and the chain EI are
  the best values of different directions (above). The fixed and the free clump give one peak
  near 1.3–1.4 Hz and a next mode at 2.5–2.6 Hz that the test does not have.

### Tow tests

Cross-flow (CF) and in-line (IL) A/D_rms over the 30 tests with Ur ≥ 3 at the three stations;
the ratio is CableDyn over measured at the same speed.

| Set | Fitted to these tests | CF r.m.s. log error (1/2L, 3/4L, L) | CF ratio, median | IL r.m.s. log error | CF peak at 1/2L (Ur, f) | Release Ur | Top tension (N) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Measured | | 0 | 1 | 0 | 0.776 (5.58) | 5.88 | 1970 |
| `uls19` | no | 0.731 (0.690, 0.605, 0.876) | 1.60 | 0.742 | 0.755 (5.30, 1.25 Hz) | 7.06 | 1962–2076 |
| `b1` | no | **0.594** (0.497, 0.479, 0.769) | 0.95 | 1.017 | 0.483 (3.83, 1.25 Hz) | 5.89 | 1962–2185 |
| `wake` | no | 1.208 (1.190, 1.291, 1.137) | 0.41 | (no in-line VIV) | 0.296 (4.41, 1.30 Hz) | 5.30 | 1962–2031 |

The measured top tension is Table 6-1's; the model's runs from the still-water value to the
mean at 1.0 m/s.

Cross-flow A/D_rms at 1/2L at every test from Ur 2.9, with the dominant cross-flow frequency of
`uls19`:

| V (m/s) | Ur | Direction | Measured | `uls19` | `b1` | `wake` | f (Hz), `uls19` |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 0.250 | 2.94 | backwards | 0.03 | 0.10 | 0.36 | 0.17 | 0.97 |
| 0.275 | 3.24 | forwards | 0.11 | 0.25 | 0.42 | 0.22 | 1.07 |
| 0.300 | 3.53 | backwards | 0.18 | 0.42 | 0.46 | 0.26 | 1.15 |
| 0.325 | 3.83 | forwards | 0.31 | 0.54 | 0.48 | 0.28 | 1.18 |
| 0.350 | 4.12 | backwards | 0.43 | 0.62 | 0.47 | 0.29 | 1.20 |
| 0.375 | 4.41 | forwards | 0.50 | 0.67 | 0.47 | 0.30 | 1.20 |
| 0.400 | 4.71 | backwards | 0.67 | 0.72 | 0.46 | 0.28 | 1.23 |
| 0.425 | 5.00 | forwards | 0.70 | 0.74 | 0.46 | 0.25 | 1.25 |
| 0.450 | 5.30 | backwards | 0.66 | 0.75 | 0.44 | 0.17 | 1.25 |
| 0.475 | 5.59 | forwards | 0.78 | 0.74 | 0.43 | 0.19 | 1.30 |
| 0.500 | 5.88 | backwards | 0.30 | 0.74 | 0.32 | 0.21 | 1.32 |
| 0.525 | 6.18 | forwards | 0.52 | 0.70 | 0.34 | 0.22 | 1.35 |
| 0.550 | 6.47 | backwards | 0.33 | 0.65 | 0.36 | 0.22 | 1.38 |
| 0.575 | 6.77 | forwards | 0.42 | 0.57 | 0.38 | 0.22 | 1.43 |
| 0.600 | 7.06 | backwards | 0.46 | 0.48 | 0.39 | 0.21 | 1.82 |
| 0.625 | 7.36 | forwards | 0.52 | 0.52 | 0.40 | 0.20 | 1.82 |
| 0.650 | 7.65 | backwards | 0.43 | 0.56 | 0.41 | 0.18 | 1.85 |
| 0.675 | 7.95 | forwards | 0.56 | 0.59 | 0.42 | 0.15 | 1.88 |
| 0.700 | 8.24 | backwards | 0.38 | 0.62 | 0.44 | 0.04 | 1.90 |
| 0.725 | 8.53 | forwards | 0.56 | 0.65 | 0.44 | 0.05 | 1.90 |
| 0.750 | 8.83 | backwards | 0.36 | 0.68 | 0.45 | 0.06 | 1.93 |
| 0.775 | 9.12 | forwards | 0.44 | 0.70 | 0.45 | 0.06 | 1.95 |
| 0.800 | 9.42 | backwards | 0.37 | 0.71 | 0.45 | 0.07 | 2.00 |
| 0.825 | 9.71 | forwards | 0.35 | 0.71 | 0.45 | 0.07 | 2.02 |
| 0.850 | 10.01 | backwards | 0.30 | 0.71 | 0.45 | 0.08 | 2.08 |
| 0.875 | 10.30 | forwards | 0.21 | 0.70 | 0.45 | 0.08 | 2.12 |
| 0.900 | 10.59 | backwards | 0.28 | 0.70 | 0.45 | 0.09 | 2.17 |
| 0.925 | 10.89 | forwards | 0.20 | 0.69 | 0.45 | 0.09 | 2.23 |
| 0.950 | 11.18 | backwards | 0.18 | 0.67 | 0.45 | 0.09 | 2.27 |
| 0.975 | 11.48 | forwards | 0.21 | 0.66 | 0.45 | 0.09 | 2.33 |
| 1.000 | 11.77 | backwards | 0.20 | 0.64 | 0.19 | 0.09 | 2.38 |

- **`uls19`** (Ulveseter et al. 2019, fitted to another riser):
  - the first peak is 0.755 at Ur 5.30 against 0.776 at Ur 5.58, and from Ur 4.7 to 5.6 the
    model is within 0.10 of every test (0.72–0.75 against 0.66–0.78);
  - the rise starts at a lower speed: 0.25–0.67 against 0.11–0.50 at Ur 3.24–4.41;
  - the lock-in holds to Ur 6.77 (0.57) and releases at Ur 7.06, 1.2 Ur after the measured
    release at 5.88; the response stays on the lower in-water mode (1.2–1.4 Hz) up to there and
    moves to the upper one (1.8–2.0 Hz) at the release;
  - on the upper branch, Ur 7.0–8.0, it is within 0.03 of the forward tests (0.52–0.59 against
    0.52–0.56); above Ur 9 it holds 0.64–0.71 where the test falls to 0.18–0.44;
  - in line at 1/2L it is about half the measured value during the lock-in (0.08–0.13 against
    0.14–0.23 at Ur 4.1–5.6) and mostly above it after the release;
  - below the compared range it responds at Ur 1.47 (0.11 at 1.3 Hz, f D/U 0.49), where the
    test shows none (0.002).
- **`b1`** (generic): a plateau of 0.42–0.48 from Ur 3.24 to 5.59, 38 % below the measured
  peak, from an early onset (0.36 at Ur 2.94 against 0.03); the release falls at the measured
  Ur (5.89); the upper branch is 0.32–0.45. The in-line response during the lock-in is well
  below the measured one (0.01–0.06 against 0.03–0.23).
- **`wake`** (Facchinetti et al. 2004): 0.22–0.30 up to Ur 5 on the lower mode, 62 % below the
  measured peak; on the upper mode from Ur 5.30 (0.15–0.22) and on a higher one (2.975–3.45 Hz,
  0.04–0.09) from Ur 8.24. The model has no in-line VIV force, so its in-line response is the
  drag's alone.
- **By direction** (1/2L, cross-flow r.m.s. log error, forward / backward tests): `uls19`
  0.647 / 0.729, `b1` 0.564 / 0.419, `wake` 1.214 / 1.164. The forward tests are the larger
  ones on the upper branch, which suits `uls19` and not `b1`.
- **No set predicts these tests.** The height of the first peak is reproduced by `uls19`, its
  speed range by none: `uls19` locks in too early and too long, `b1` releases at the right speed
  at 62 % of the amplitude, and the wake oscillator is too weak.

### Lower end and the tow prediction

`uls19` at all 40 speeds with the other two lower ends (`--end`; `provenance.json`):

| End | CF r.m.s. log error (1/2L) | CF peak at 1/2L (Ur) | Release Ur | RAO r.m.s. log error, surge |
| --- | --- | --- | --- | --- |
| `identified` | 0.731 (0.690) | 0.755 (5.30) | 7.06 | 0.459 |
| `table` (MARIN's fixed clump) | 0.721 (0.433) | 0.805 (5.59) | 8.83 | 1.186 |
| `free` (a `Free` End B) | 0.692 (0.429) | 0.811 (5.59) | 9.12 | 1.174 |

With the fixed or free clump the response at 1/2L drops to 0.19–0.36 from Ur 9.4, which lowers
the error there, while station L moves 0.53–0.67 against the measured 0.14–0.32 and the lock-in
releases 2.9–3.2 Ur late. The tow score over the three stations does not separate the ends;
the identified end is used because it reproduces the RAOs. Every end releases late.

### Convergence

At 0.30, 0.45, 0.50, 0.60 and 0.80 m/s for the three sets (`provenance.json`):

- half the time step (1 ms) changes no A/D_rms by more than 0.006 (cross-flow at most 2.1 %)
  and the dominant frequency by at most one bin;
- two elements per table element change no A/D_rms by more than 0.004 (cross-flow at most
  1.4 %), and the frequency not at all;
- 120 s runs (100 s of statistics) change no A/D_rms by more than 0.005 (cross-flow at most
  1.0 %).

For the RAOs, two elements per table element and half the time step change every r.m.s. log
error by at most 0.001.

## Limitations

- The structure below the cable is an identified impedance, not a measured one, and the
  in-water upper peak and station L are not matched closely; the cable's EI and damping come
  from the in-air tests, whose losses include the end fixtures.
- The tow tests are single runs per speed with the direction alternating, without repeats or
  uncertainties; the forward and backward tests differ on the upper branch.
- The measured response frequencies are not published, so the mode a set locks in to cannot be
  checked against the test except through the amplitudes of the three stations.
- The data are digitised from figures (Appendix 3 from raster images, ±0.004 in Ur).
