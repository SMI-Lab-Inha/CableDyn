<!-- SPDX-License-Identifier: Apache-2.0 -->

# Steel catenary riser: cross-flow VIV in uniform current and under heave

This directory compares CableDyn's cross-flow VIV with the SJTU model tests of a truncated steel
catenary riser: Wang et al. (2014), <https://doi.org/10.1016/j.marstruc.2014.07.002>, as Thorsen
et al. (2017), <https://doi.org/10.1016/j.marstruc.2016.10.007>, use them to verify the
time-domain synchronisation model that CableDyn's cross-flow VIV implements. Both of their cases
are run, each with their published coefficient set for it; nothing is fitted to the measurements.

- **Uniform current, 0.2 m/s (case L3-W):** the dominant strain frequency 8.25 m from the lower
  end is 0.992 Hz against the measured 1.00 Hz (Thorsen et al.'s model 1.02 Hz). The largest r.m.s.
  cross-flow strain is 7.3 % above the measured largest (theirs 0.4 % below), the mean over the 25
  stations 3.6 % above (theirs 3.0 % below), and the r.m.s. log error over the stations 0.338
  (theirs 0.360).
- **Heave-induced VIV in still water (case L3-Wh):** with the phase input of the published model
  (no `VIV FRAME`), the largest r.m.s. strain is 20.8 % above the measured largest (theirs
  5.6 % above) and the station mean 33 % above (theirs 8.7 % above); the dominant frequency is
  1.019 Hz (measured close to 1 Hz). Over the lower half of the span the strain is 1.71 times
  the measured on average, over the upper half 0.93 times; Thorsen et al.'s model has the same
  pattern, 1.44 and 0.71.
- **The heave result depends on the `VIV FRAME` filter time `Tf`.** The heave period, 5.96 s, is
  about 6 VIV periods, below the scale separation the filter needs (a motion period of at least
  about 30 VIV periods; [`VIV FRAME`](../../../doc/driver_format.md)). With `Tf` from 0.35 to 1 s
  the response falls on one of two branches, about 1.02 Hz or 0.88 Hz, and at `Tf` 0.5 s halving
  the time step moves it from the first to the second. The `Tf` runs are a sensitivity, not the
  result.

## Model and data

- **Riser (Thorsen et al. 2017, Table 1; Wang et al. 2014, Table 2):**
  - length 23.71 m, between universal joints (moment free);
  - the top 21.0425 m horizontally and 9 m vertically from the lower end, with a 2.626 m flowline
    on the fake seabed;
  - D 0.024 m, 0.69 kg/m (a solid section of coating, PPR pipe and copper core: no internal water
    is added), EI 10.5 N m², EA 6.66e5 N (Thorsen et al.; Wang et al. give 6.67e5 N).
- **Seabed (Thorsen et al. 2017, Sec. 5.1):** 1 (kN/m)/m in compression (`kBot` = 1000/D),
  friction 0.2.
- **Coefficients (Thorsen et al.'s sets):** C_v 0.7 in the current (their choice for Re ≈ 3000)
  and 0.8 in heave (Re ≈ 4000, from the largest r.m.s. in-plane velocity of their model,
  0.16 m/s); f̂0 0.18, Δf̂ 0.08, cross-flow only; C_D 1.2, C_M 2.
- **Assumed:** the top joint 0.5 m below the surface, in 9.5 m of water (not given).
- **CableDyn:** one finite-EI line, 200 elements, dt 0.005 s (Thorsen et al.'s step; they use 500
  elements).

### Uniform current (Thorsen et al. 2017, Sec. 5.1)

- **Current:** 0.2 m/s along the riser plane, from the lower end towards the top. In the test the
  top and the seabed moved together on tracks through still water.
- **Run:** 150 s, statistics over t > 30 s (frequency resolution 0.0083 Hz).
- **Measured:** the r.m.s. dynamic cross-flow strain at the outer surface at the 25 stations
  (Fig. 17), and the dominant strain frequency 8.25 m from the lower end, 1.00 Hz (Fig. 18).
  Thorsen et al.'s model: 99.6 % of the measured largest r.m.s. strain (their text), 1.02 Hz.

### Heave-induced VIV (Thorsen et al. 2017, Sec. 5.2)

- **Top motion:** the measured x and z motion of the upper end (Fig. 19, top panel), fitted with
  sinusoids of the period the text gives, 5.96 s: x amplitude 0.0485 m, z amplitude 0.1305 m, x
  lagging z by 83°; the residuals of the fits are 3 and 7 mm, about the width of the plotted
  lines. The motion is prescribed about the static position, ramped in over two periods, in still
  water. Thorsen et al. prescribe the measured motion itself; Wang et al. drove the top along
  trajectories they computed for pure heave of the full-length riser.
- **Run:** 200 s, statistics over t > 40 s (frequency resolution 0.0063 Hz).
- **Processing, as Wang et al. (2014, Sec. 3.1) process the measurement:** the cross-flow strain
  band-passed from 4 times the motion frequency (0.67 Hz) to 50 Hz.
- **Measured:** the r.m.s. cross-flow strain at the 25 stations (Fig. 22; the largest, 5.96e-5,
  at s 2.84 m); the dominant frequency 8.25 m from the lower end "close to 1 Hz" (Fig. 23 and
  text). Thorsen et al.'s model: 6 % above the measured largest (their text), also close to 1 Hz.
- **`VIV FRAME`:** without the row the synchronisation phase reads the node's velocity across the
  relative flow, as the published model does. A row with `Tf` > 0 takes the node's velocity less
  its low-pass filtered value, so that the slow motion counts as flow; the runs at `Tf` 0.35,
  0.5, 0.75 and 1 s, and at 0.5 s with half the time step, measure how the result depends on it.

### Analysis (`wang_scr_viv.py`)

- The section axis x of a line not restrained in torsion starts at End A as global +z projected
  normal to the line and is carried along it, so `CurvX` is the in-plane (catenary) curvature
  and `CurvY` the cross-flow one; the script checks this on their static means. Strain = κ D/2,
  the outer-surface strain of Figs. 17 and 22; arc length s from the lower end.
- The r.m.s. strain is compared with the measured values at the stations (linear interpolation
  between nodes); the ratio is CableDyn over measured, and the r.m.s. log error
  √(mean (ln ratio)²). The lower half of the span is stations 1–13 (s < 11.855 m), the upper half
  stations 14–25. Thorsen et al.'s simulated line is scored the same way.
- The dominant frequency is the peak of the Hann-windowed spectrum of the strain at the node
  nearest 8.25 m. In heave, the shares of that spectrum in 0.80–0.95 Hz and 0.95–1.10 Hz tell the
  two response branches apart.

### Data

`extract_thorsen2017.py` reads Thorsen et al.'s PDF (not distributed; PyMuPDF, NumPy and Pillow):

- `thorsen2017_fig17.csv`, `thorsen2017_fig22.csv`: the measured points and Thorsen et al.'s
  simulated line of Figs. 17 and 22, read from the vector graphics. The tick labels are glyph
  outlines, so the axes are calibrated on the tick marks. As a check, the simulated line peaks at
  99.6 % (Fig. 17) and 105.6 % (Fig. 22) of the measured largest value, where the text gives
  99.6 % and 6 % too high.
- `thorsen2017_fig19_top.csv`, `thorsen2017_fig19_fit.json`: the top motion of Fig. 19, a raster
  image, traced by colour (742 and 738 pixel columns over t 30–50 s), with the axis calibrated on
  the dotted grid lines, and the fitted sinusoids.

The data are digitised from the article; its figures are not reproduced.

## Reproduce

From the repository root, with `path/to/cabledyn.exe` the release executable or a local build:

```powershell
python validation/experiments/wang_scr/wang_scr_viv.py `
  path/to/cabledyn.exe validation/experiments/wang_scr/work --jobs 11 `
  --summary validation/experiments/wang_scr
```

The default is the committed record: the 11 runs of the tables below, about 5 minutes with 11
parallel runs on a 24-core workstation. Run into a new work directory with the executable of
`summary.json`, it reproduces `response_comparison.csv` and `summary.json` byte for byte.
`--runs` selects runs by name (a comma list: `current`, `current_dt2`, `current_n400`, `heave`,
`heave_dt2`, `heave_n400`, `heave_tf0.35`, `heave_tf0.5`, `heave_tf0.75`, `heave_tf1`,
`heave_tf0.5_dt2`). A run whose folder holds the same input files and executable is not
repeated; `--keep` keeps the output tables (100–400 MB a run).

## Results

`response_comparison.csv` holds the r.m.s. strain of every run, and of Thorsen et al.'s line, at
the 25 stations; `summary.json` the statistics, the top motion and the SHA-256 of the
executable; `provenance.json` the sources and the build.

### Uniform current (L3-W)

| | Largest r.m.s. strain / measured largest | Station mean ratio | r.m.s. log error | Lower / upper half mean ratio | f at 8.25 m | r.m.s. peaks along the span |
| --- | --- | --- | --- | --- | --- | --- |
| Measured | 1 | 1 | 0 | 1 / 1 | 1.00 Hz | 6 |
| CableDyn | **1.073** (s 14.0 m) | 1.036 | 0.338 | 0.947 / 1.133 | **0.992 Hz** (±0.008) | 7 |
| Thorsen et al. (2017) | 0.996 (s 14.1 m) | 0.970 | 0.360 | 0.879 / 1.068 | 1.02 Hz | 7 |

- **Frequency:** within the 0.008 Hz resolution of the measured value; Thorsen et al.'s model is
  0.02 Hz above it.
- **Strain magnitude:** the largest r.m.s. strain is 7.3 % above the measured largest, and the
  station mean 3.6 % above; Thorsen et al.'s are 0.4 % and 3.0 % below.
- **Spatial pattern:** CableDyn and Thorsen et al.'s model have one r.m.s. peak more along the
  span than the test (7 against 6), so the peaks and troughs fall at different stations and the
  station ratios scatter (0.41–1.71; Thorsen et al.'s 0.34–1.61). Station by station CableDyn's
  r.m.s. strain is 1.07 times Thorsen et al.'s on average (0.87–1.20).

### Heave-induced VIV (L3-Wh)

The primary result is the run without `VIV FRAME`, the published model's phase input:

| | Largest r.m.s. strain / measured largest | Station mean ratio | r.m.s. log error | Lower / upper half mean ratio | f at 8.25 m |
| --- | --- | --- | --- | --- | --- |
| Measured | 1 (s 2.84 m) | 1 | 0 | 1 / 1 | close to 1 Hz |
| CableDyn, no `VIV FRAME` | **1.208** (s 7.2 m) | 1.334 | 0.460 | 1.710 / 0.927 | **1.019 Hz** |
| Thorsen et al. (2017) | 1.056 (s 8.2 m) | 1.087 | 0.406 | 1.438 / 0.708 | close to 1 Hz |

- **Frequency:** 1.019 Hz; 61 % of the band-passed spectrum at 8.25 m lies in 0.95–1.10 Hz and
  1 % in 0.80–0.95 Hz.
- **Strain magnitude:** the largest r.m.s. strain is 20.8 % above the measured largest, at s
  7.2 m instead of 2.84 m; Thorsen et al.'s is 5.6 % above, at s 8.2 m.
- **The lower half reads high:** over stations 1–13 the strain is 1.71 times the measured on
  average, over stations 14–25 0.93 times. Thorsen et al. report the same pattern, "over-predicted
  for the lower part and under-predicted for the upper half"; their line gives 1.44 and 0.71.
  Station 1 (s 1.04 m, next to the lower joint, measured 0.57e-5) contributes most: 5.08 times
  the measured value (Thorsen et al. 3.10); without it the lower-half means are 1.43 and 1.30.
  Over the upper half CableDyn is closer to the measurement than Thorsen et al.'s model (0.93
  against 0.71).

#### Sensitivity to the `VIV FRAME` filter time

| `Tf` (s) | dt (s) | Largest r.m.s. strain / measured largest | Station mean ratio | Lower / upper half mean ratio | f at 8.25 m (Hz) | Spectral share 0.80–0.95 / 0.95–1.10 Hz |
| --- | --- | --- | --- | --- | --- | --- |
| none | 0.005 | 1.208 | 1.334 | 1.710 / 0.927 | 1.019 | 0.01 / 0.61 |
| 0.35 | 0.005 | 1.128 | 1.279 | 1.667 / 0.860 | 1.012 | 0.17 / 0.52 |
| 0.5 | 0.005 | 1.069 | 1.219 | 1.584 / 0.824 | 1.019 | 0.17 / 0.52 |
| 0.5 | 0.0025 | 0.972 | 1.107 | 1.454 / 0.731 | 0.894 | 0.59 / 0.22 |
| 0.75 | 0.005 | 0.967 | 1.029 | 1.400 / 0.628 | 0.881 | 0.69 / 0.14 |
| 1.0 | 0.005 | 0.958 | 1.009 | 1.390 / 0.596 | 0.881 | 0.71 / 0.12 |

- **Two branches:** `Tf` ≤ 0.5 s at dt 0.005 s keeps the response at 1.01–1.02 Hz with the largest
  strain 1.07–1.21 times the measured; `Tf` ≥ 0.75 s moves it to 0.88 Hz, the largest strain
  0.96–0.97 times the measured, at s 4.4 m, and the station mean to 1.01–1.03. At `Tf` 0.5 s,
  halving the time step moves the response from the first branch to the second (largest strain
  1.069 to 0.972, f 1.019 to 0.894 Hz).
- **On the lower branch the strain distribution departs further from the measured one:** the
  r.m.s. log error rises from 0.460 to 0.49–0.52 and the upper-half mean falls to 0.60–0.63, so
  the station mean near 1 comes from a high lower half and a low upper half.
- **Why no `Tf` separates the scales:** the first-order filter keeps the fraction
  ωTf/√(1 + (ωTf)²) of the VIV velocity in the phase input and leaves the same fraction of the
  motion's velocity, at the motion frequency, in it. At the VIV frequency (1.02 Hz) and the heave
  frequency (0.168 Hz):

  | `Tf` (s) | VIV velocity kept | Heave velocity left in the phase input |
  | --- | --- | --- |
  | 0.35 | 0.913 | 0.346 |
  | 0.5 | 0.955 | 0.466 |
  | 0.75 | 0.979 | 0.620 |
  | 1.0 | 0.988 | 0.726 |

  Keeping 95 % of the VIV velocity while leaving at most 10 % of the motion's needs a motion
  period of at least about 30 VIV periods; here it is about 6. Every `Tf` that keeps the VIV
  velocity leaves 35–73 % of the heave velocity in the phase input, against 100 % without the row.

### Convergence

Re-run on this build (`summary.json`, runs `*_dt2` and `*_n400`):

- **Uniform current:** half the time step and 400 elements change the largest r.m.s. strain by
  +0.13 % and +0.24 %, the station mean by at most 0.05 %, and leave the frequency and the number
  of peaks unchanged.
- **Heave, no `VIV FRAME`:** half the time step and 400 elements change the largest r.m.s.
  strain by −1.0 % and −0.7 %, the station mean by −2.9 % and −2.3 %, the upper-half mean by
  −5.8 % and −4.2 %, and leave the frequency (1.019 Hz) and the branch unchanged.
- **Heave with `VIV FRAME` `Tf` 0.5 s:** not converged in the time step; it changes branch (above).

## Limitations

- The heave case tests the model in a regime where its phase input is not defined uniquely: the
  published model reads the node's velocity in the earth frame, and a filtered frame cannot
  separate a motion only about 6 VIV periods long. Over the runs here the largest strain moves
  by up to 21 % and the frequency by up to 14 % with that choice.
- The cause of the lower-half excess in heave is not identified. It is larger than in Thorsen et
  al.'s model; the depth of the top joint, the seabed contact and the static procedure (Thorsen
  et al. reach the static shape by a prescribed top displacement) differ from theirs or are
  assumed.
- The top motion is a two-sinusoid fit to a raster figure; the measured motion's higher
  harmonics (residual r.m.s. 3–7 mm, about the plotted line width) are not reproduced.
- Thorsen et al. do not state which case of the Wang et al. (2014) test matrix their heave case
  is; the comparison uses the top motion and strains they plot.
- The measured strains were taken at a diameter of 19.5 mm and scaled to the outer surface by
  Thorsen et al.; the uncertainties of the measurements are not published.
