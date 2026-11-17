<!-- SPDX-License-Identifier: Apache-2.0 -->

# Southampton untensioned catenary riser: VIV with high curvature and touchdown

This directory compares CableDyn with the tests of Chaplin & King (2018),
<https://doi.org/10.1016/j.jfluidstructs.2018.01.008>: a 5.36 m, 56 mm model riser of very low
bending stiffness, pinned on the floor of a 2 m deep flume and 0.2 m above the surface 4.30 m
away, tensioned only by its weight and the drag, at two mass ratios in currents of 0.2-1.4 m/s in
both directions. Much of it lies on the floor, and the touchdown moves with the flow: a highly
curved, untensioned catenary with touchdown in VIV, the configuration of a free-hanging dynamic
cable near the seabed. Its frequencies test what a local time-domain VIV model does where the
normal velocity varies threefold along the span: the article measures one cross-flow frequency
along the whole riser, close to f = 0.11 U/d with U the flow speed, well below the Strouhal number
of about 0.2 of a stationary cylinder. No model parameter is fitted to these tests.

In short:

- **All three VIV models reproduce the uniform frequency along the span** in positive current
  (within 10 % over the suspended stations in 26 of 26 runs with the generic set, 22 with the
  wake oscillator, 20 with Ulveseter et al.'s set), **but at their own Strouhal number, not at
  the measured 0.11:** f d/U 0.121 (Ulveseter et al. 2019, 9 % high), 0.138 (generic set, 24 %
  high) and 0.17 (wake oscillator, St 0.2, 53 % high). The frequency follows the flow at the
  near-vertical top, as measured, and the lower part follows it.
- **In negative current** the heavy riser's measured frequency stays at 1.05-1.16 Hz from
  1.1 m/s, where the straight part's normal velocity is constant; Ulveseter et al.'s set also
  levels off, about 35 % higher at 1.45-1.5 Hz (f d/(U sin ψ) 0.130 against 0.119 measured), the
  generic set gives 1.45-1.85 Hz and the wake oscillator 2.1-2.45 Hz.
- **The mean top tension needs the drag amplification:** with it (the generic set) the light
  riser's tension in positive current is within 9.8 % r.m.s. and rises by 88.8 N to 1.4 m/s
  against 88.6 N measured; with C_d 1.2 it rises by 49.8 N without VIV and by 60.1 and 51.4 N
  with the two sets that keep C_d. The same set gives the heavy riser's straight-line inclination
  in negative current to 5.0° r.m.s. (the drag coefficient it implies by the article's Eq. 9,
  2.55, against 2.47 measured), and overshoots its tension rise above 1.2 m/s.
- **The still-water tension of the heavy riser is 37 % above the measured one** (81.2 against
  59.1 N; light riser 15.1 against 13.4 N) with the article's length, span and weights, a
  geometry the article's numbers do not close; the tension is therefore also compared as the
  increment from still water.

## The tests (Chaplin & King 2018)

- **Riser (Table 1, Sec. 2):** a 4 mm wire rope through nylon bobbins (30 mm long, 10 mm apart) in
  a silicone sheath, d 56 mm, 5.36 m between pinned ends; EI 1.3 N m² and EA 790 kN estimated from
  the parts; structural damping about 3 % of critical in air (free decays, hanging vertically);
  submerged weight 2.9 N/m (m* 1.12) and 20.9 N/m (m* 1.84, lead plugs; the twelve instrumented
  bobbins have none), m* the mass with the internal water over the displaced mass.
- **Set-up:** a recirculating flume, 2.0 m deep; the lower pin on the floor, the upper pin 0.2 m
  above the surface and 4.30 m away horizontally. Positive current: the upstream pin on the floor
  (the profile concave to the flow); negative: the upstream pin at the surface, obtained in the
  laboratory by swapping the pin elevations (Sec. 3). U from 0.2 to 1.4 m/s in 0.05 m/s steps,
  Re up to 70 000.
- **Measured:** tension load cells at both ends; twelve 3-axis accelerometers 0.4 m apart; video
  of up to three bobbins. The article reports the mean top tension with ±1 standard deviation
  (Figs. 10 and 13), spectral peak frequencies of the cross-flow accelerations (Figs. 7-9, 15, 16),
  the mean profile and the inclination of the straight part in negative current (Figs. 5, 11, 12),
  and no amplitude statistics.
- **Findings the case tests:** in positive current the cross-flow frequency is uniform along the
  riser and close to f = 0.11 U/d, U within 3 % of the normal velocity at the near-vertical top
  (Fig. 9), so that the local Strouhal number on the normal velocity runs from 0.11 at the top
  to above 0.3 lower down (Fig. 8); in negative current the riser takes a straight inclined
  profile on which the normal drag balances the normal weight (their Eq. 9; C_d 2.0-2.8 bracket
  the measured inclinations), U sin ψ tends to a constant, and the frequency follows
  f = 0.13 U sin ψ/d (Fig. 15(b)) except where the response locks in.

### Caveats of the data

- **Still-water tension.** With the stated length, span, depth and weights the riser hangs as a
  catenary whose top tension is 81.2 N (m* 1.84) and 15.1 N (m* 1.12), with the 0.2 m above the
  surface in air. The article's still-water points are 59.1 and 61.8 N (Figs. 10 and 13, the
  two rigs) and 13.4 and 11.5 N. A catenary of 5.36 m between the stated pins cannot hang that
  slack; the touchdown in Fig. 5 is also further from the floor pin than the model's. A laid part
  that is not straight (slack held by the floor friction) or a geometry other than Table 1's would
  explain both; the article does not say. The comparison is made on the absolute tension and on
  its increment from still water, each figure against its own still-water point.
- **The reference lines W h of Fig. 13** (41.8 and 5.8 N, Table 1's weight times the water
  depth) leave out the in-air weight of the part above the surface, m g × 0.2 m, 9.0 and 5.4 N
  in the model; the measured tensions in negative current sit on these lines.
- **Frequencies** are spectral peaks of accelerometer signals, which gravity contaminates at the
  vibration frequency (Sec. 3); the article shows that the peaks agree with those of video
  displacements (Fig. 4). Several peaks are given at some speeds (Fig. 15).
- **Not given:** the accelerometer positions (only the 0.4 m spacing; the profiles of Fig. 5 put
  the first 0.5 m from the floor pin), the floor friction, the flume's water density, the spectral
  method, and the scatter of the frequencies.

## Model

- **Line:** one finite-EI line of 54 elements (0.099 m), End A the upper pin at (4.30, 0, 0.2) m,
  End B the lower pin on the floor at (0, 0, -2.0) m, both pinned. d 56 mm, EI 1.3 N m², EA
  790 kN, mass ρA + W/g (2.759 and 4.594 kg/m; Table 1's m* 1.84 gives 4.532 kg/m, 1.4 % less),
  the same mass above the surface. Fresh water, 1000 kg/m³ (not stated; m* 1.12 and 2.9 N/m agree
  at it). C_d 1.2, C_a 1.
- **Structural damping:** Rayleigh damping (`EI|-zeta|f1|f2`), 3 % at 0.3 and 3.0 Hz, the range of
  the measured cross-flow frequencies; run without it as a sensitivity.
- **Floor:** the default penalty contact and a friction coefficient 0.5 (a silicone sheath on the
  flume floor; not given). Without friction the laid part slides sideways under the suspended
  part's motion; frictionless runs are a sensitivity.
- **Current:** uniform over the depth along the riser's plane, +x positive, -x negative, ramped
  linearly from still water over 10 s (a static solve in the full current finds no equilibrium
  for some light-riser cases in negative current, where the riser lifts off the floor; every run
  starts from the same still-water state). The part above the surface is dry.
- **Runs:** dt 2.5 ms, 110 s, statistics over 30-110 s (80 s); without VIV 40 s, statistics over
  30-40 s. U = 0.2-1.4 m/s in 0.1 m/s steps for both mass ratios and directions (52 runs per set).

### Coefficient sets

The sets of the FLOTANT and NDP cases, none fitted to these tests, and a run
without VIV:

| Set | Source | Values | Fitted to |
| --- | --- | --- | --- |
| `uls19` | Ulveseter et al. (2019), Table 1 | C_l 0.85, f̂0 0.144 ± 0.064; in-line C_x 0.75, α 0.15, C_hh 1.5 | the NDP riser in uniform flow: blind here |
| `b1` | The generic set of the Delta Flume case | C_l 0.8, f̂0 0.18 ± 0.08; in-line C_x 0.85, α 0.15 (Ulveseter et al. 2018, Table 3); `VIV AMPLITUDE` C_v 0.8 with the drag amplification of Chaplin et al. (2005, Eq. 4); `VIV DAMPING` `0 0 1` | nothing from these tests: blind |
| `wake` | Facchinetti et al. (2004) | `WAKE OSCILLATOR` C_L0 0.3, St 0.2, ε 0.3, A 12; cross-flow only | rigid-cylinder data: blind |
| `none` | | no VIV force, C_d 1.2 | |

Only `b1` amplifies the drag with the response; `uls19` and `wake` keep C_d 1.2, so their mean
tension differs from `none` only by the dynamics.

### Processing

- **Tension:** mean and standard deviation of the tension at the upper pin (`FairTen1`); the
  measured value of each figure at the same speed; the increment from still water, each against
  its own still-water value.
- **Frequency:** the peak of the averaged spectrum (20 s Hann segments, half overlapping,
  0.05 Hz resolution) of the cross-flow (y) displacement within 0.3-6 Hz (the lowest measured
  frequency is 0.34 Hz; below it some runs wander slowly sideways), the dominant frequency the
  article shows to agree in its acceleration and video-displacement spectra (Fig. 4); the peak
  of the same spectrum weighted by ω⁴, the acceleration spectrum, within 0.15-6 Hz, is in
  `response_comparison.csv` (`model_f_acc_…`). Compared points, by their mean elevation: the two
  dots of the inset of Fig. 9(a) (0.65 and 0.37 of the height between the pins; the mean of the
  two) and the dot of Fig. 9(b) (0.19) in positive current, mid-depth (Fig. 15, "around
  mid-depth") in negative current. The model's frequency is set against the nearest measured
  peak at the same speed and mass ratio (`nearest_measured_…`) and against the median measured
  f d/U. The spread along the span is (max - min)/median over the twelve accelerometer stations,
  taken 0.48 m from each pin and 0.4 m apart, that are off the floor and under water.
- **Inclination:** of the mean profile to the horizontal at mid-depth; the drag coefficient it
  implies by the article's Eq. 9, C_d = 2 W cos ψ / (ρ d U² sin² ψ), for the model and the
  measurement, as a comparison only: nothing is fitted to it.

### Measurements

`extract_chaplin2018.py` reads the article's PDF (not distributed; PyMuPDF, NumPy, SciPy, Pillow)
into `measured_chaplin2018.csv`. The figures are greyscale raster images:

- each panel is calibrated on its frame and the tick marks on both opposite sides of it (largest
  residual 0.0010 m/s, 0.0036 Hz, 0.07°, 0.14 N);
- the markers are found by a circle Hough transform, classed solid or open by their fill and
  snapped to the 0.05 m/s test grid (largest snap 0.012 m/s); peaks off the grid or between the
  classes are dropped and listed by the script; the error bars are read from their caps;
- **checks:** the line f d/U = 0.11 drawn in Fig. 9 reads 0.105 and 0.108; the lines W h in
  Fig. 13 read 41.9 and 5.65 N against Table 1's 41.8 and 5.8 N;
- **points read:** Fig. 10, 26 open and 24 solid (every marker; 46 standard deviations); Fig. 13,
  25 and 26 (one open marker lost where the series cross; 44 standard deviations); Fig. 12(a),
  20 and 25; Fig. 15(a), 26 and 47; Fig. 9(a), 12 and 40, Fig. 9(b), 22 and 34. Markers hidden
  under others are lost, which in Fig. 9(a), where both series and both accelerometers fall on
  one line, is most of the open ones; `--check DIR` draws what was read over each figure.

The data are digitised from the article; its figures are not reproduced.

## Reproduce

From the repository root, with `path/to/cabledyn.exe` the release executable or a local build:

```powershell
python validation/experiments/chaplin_catenary/extract_chaplin2018.py CHAPLIN2018.pdf
python validation/experiments/chaplin_catenary/chaplin_catenary_viv.py `
  path/to/cabledyn.exe validation/experiments/chaplin_catenary/work --jobs 4 `
  --summary validation/experiments/chaplin_catenary
```

The default is the committed record: 2 still-water runs, 208 primary runs, 12 convergence runs and
72 sensitivity runs, about an hour with 4 parallel runs. `--runs` selects runs by name; a run
whose folder holds the same input files and executable is not repeated; `--keep` keeps the output
tables (about 110 MB a run). `response_comparison.csv` holds every run, `summary.json` the scores
and the SHA-256 of the executable, `provenance.json` the sources and assumptions.

## Results

`response_comparison.csv` holds every run against the measured values at its speed;
`summary.json` the scores below.

### Still water

| m* | CableDyn (N) | Measured, Fig. 10 (N) | Measured, Fig. 13 (N) |
| --- | --- | --- | --- |
| 1.12 | 15.1 | 13.4 | 11.5 |
| 1.84 | 81.2 | 59.1 | 61.8 |

The model's touchdown is 2.18 m (m* 1.84) and 1.99 m (m* 1.12) from the floor pin; Fig. 5 puts
the heavy riser's touchdown at the lowest speeds near 2.6 m (by eye). See the caveats above.

### Mean top tension

r.m.s. relative error of the mean tension over the 13 speeds (12 where Fig. 10 has no marker),
r.m.s. error of its increment from still water (N), and the median ratio of the standard
deviations (CableDyn over measured, where the error bar was read):

| Set | Positive, m* 1.12 | Positive, m* 1.84 | Negative, m* 1.12 | Negative, m* 1.84 |
| --- | --- | --- | --- | --- |
| `none` (C_d 1.2) | 0.245; 19.6 N | 0.135; 18.4 N | 0.365; 11.5 N | 0.385; 22.8 N |
| `uls19` | 0.178; 14.5 N; sd 0.81 | 0.150; 12.2 N; sd 1.14 | 0.322; 6.9 N; sd 1.43 | 0.296; 14.3 N; sd 1.61 |
| `b1` | **0.098; 1.8 N**; sd 0.59 | 0.238; **4.9 N**; sd 0.69 | 0.301; **3.9 N**; sd 1.22 | 0.248; 13.7 N; sd 1.21 |
| `wake` | 0.234; 18.8 N; sd 0.19 | 0.139; 17.7 N; sd 0.42 | 0.351; 10.5 N; sd 0.42 | 0.316; 16.4 N; sd 0.76 |

Tension at 1.4 m/s (N; the increment from still water in brackets):

| | Positive, m* 1.12 | Positive, m* 1.84 | Negative, m* 1.12 | Negative, m* 1.84 |
| --- | --- | --- | --- | --- |
| Measured | 102.0 (+88.6) | 154.9 (+95.8) | 58.3 (+46.8) | 49.5 (−12.3) |
| `none` | 64.9 (+49.8) | 145.6 (+64.4) | 39.7 (+24.6) | 41.9 (−39.3) |
| `uls19` | 75.2 (+60.1) | 158.5 (+77.3) | 48.2 (+33.1) | 50.8 (−30.4) |
| `b1` | 103.8 (+88.8) | 188.8 (+107.6) | 70.7 (+55.7) | 68.4 (−12.8) |
| `wake` | 66.5 (+51.4) | 146.4 (+65.2) | 41.5 (+26.4) | 49.7 (−31.5) |

- **Positive current:** the tension rises with the drag on the suspended part. Only the generic
  set, whose drag grows with the response (Chaplin et al. 2005, Eq. 4), follows the light
  riser's rise: within 5.6 N at every speed and within 1.8 N from 1.0 m/s (80.2 against
  80.2 N at 1.2 m/s). The sets that keep C_d 1.2 fall 13-36 N short at 1.0-1.4 m/s, close to
  the run without VIV: Ulveseter et al.'s response adds 5-10 N to it, the wake oscillator's
  1-2 N. For the heavy riser the still-water offset dominates the absolute error; the generic
  set's increment is 12 % above the measured at 1.4 m/s, the others 19-33 % below.
- **Negative current:** the heavy riser's tension falls from the still-water value as the drag
  lifts it into a straight line whose top tension the article puts at W h (41.8 N): measured
  59 → 40 N at 0.2-0.65 m/s, then 41-49 N to 1.4 m/s. Every set reproduces the fall, at higher
  speeds: the generic set reaches 40 N at 0.7 m/s, Ulveseter et al.'s and the wake oscillator
  near 0.9 m/s. Above 1.2 m/s the generic set rises to 68 N where the measured stays below
  50 N. The light riser lifts off the floor entirely from 0.5-0.6 m/s in the model, as measured
  (lifted at 0.6 m/s, Sec. 4.2); its tension rise to 1.4 m/s is 55.7 N with the generic set
  against 46.8 N measured, 33.1 N with Ulveseter et al.'s.
- **Fluctuation:** the generic set's standard deviation is 0.59-0.69 of the measured in positive
  current; Ulveseter et al.'s is 0.81-1.14; the wake oscillator's 0.19-0.42.

### Cross-flow frequency

Median f d/U of CableDyn at the points compared and of the measured points read (all speeds),
and the median ratio of CableDyn's frequency to the nearest measured point at the same speed:

| Set | Positive, upper (Fig. 9(a)) | Positive, lower (Fig. 9(b)) | Negative, mid-depth (Fig. 15(a)) | Uniform along the span (positive / negative runs) |
| --- | --- | --- | --- | --- |
| Measured, m* 1.12 / 1.84 | 0.111 / 0.113 | 0.090 / 0.102 | 0.064 / 0.061 | uniform (Fig. 7) |
| `uls19` | 0.121; ratio 1.11 / 1.07 | 0.121; 1.37 / 1.08 | 0.090 / 0.103; 1.47 / 1.32 | 20 of 26 / 12 of 26 |
| `b1` | 0.138; 1.24 / 1.25 | 0.138; 1.52 / 1.28 | 0.128 / 0.098; 1.84 / 1.45 | 26 of 26 / 15 of 26 |
| `wake` | 0.173 / 0.168; 1.31 / 1.40 | 0.168; 1.83 / 1.47 | 0.074 / 0.128; 0.98 / 2.03 | 22 of 26 / 22 of 26 |

Uniform: the spread of the frequency over the suspended stations under water within 10 % of its
median.

- **Positive current: one frequency along the span.** In every set the cross-flow response of
  the whole suspended span, the curved and nearly horizontal lower part included, is at one
  frequency set by the flow at the near-vertical top; the local normal velocity there is within
  3 % of U, and further down the same frequency means a local f d/(U sin ψ) rising towards the
  touchdown, as Fig. 8 measures. The lower part does not respond at the frequency of its own
  normal velocity: its synchronisation or wake follows the motion the span imposes on it. The
  level is each model's own: Ulveseter et al.'s set, whose band is centred at 0.144, settles at
  0.121, 9 % above the measured 0.11 (7 of 12 heavy-riser speeds within 10 % of a measured
  peak); the generic set (centre 0.18) at 0.138; the wake oscillator at its Strouhal number,
  0.17.
- **Not reproduced in positive current:** the heavy riser's second branch at the lower
  accelerometer, 0.67-0.83 Hz at 0.6-1.0 m/s (about half the dominant frequency, Fig. 9(b)),
  and the light riser's frequency at the same point levelling at 1.57-1.81 Hz from 0.9 m/s;
  CableDyn keeps the dominant frequency there.
- **Negative current, heavy riser:** the measured frequency at mid-depth rises to about 1.05 Hz
  at 0.5-0.6 m/s, dips (the article attributes it to lock-in while the touchdown moves) and
  stays at 1.05-1.16 Hz from 1.1 m/s, where U sin ψ is constant. Ulveseter et al.'s set has the
  same plateau at 1.45-1.5 Hz from 1.0 m/s, f d/(U sin ψ) 0.130 on the model's own inclination
  (measured 0.119 on the measured one, the article's line 0.13); the generic set 1.45-1.85 Hz;
  the wake oscillator 2.1-2.45 Hz.
- **Negative current, light riser:** the measured spectra have a low branch at 0.34-0.45 Hz to
  0.85 m/s and a branch rising from 0.5 Hz at 0.6 m/s to 1.75 Hz at 1.3 m/s, which the article
  attributes to lock-in under the rising tension; Ulveseter et al.'s and the generic set have
  0.5-0.85 Hz to 0.8 m/s and then mostly 2.2-3.6 Hz at mid-depth, above the measured branch; the
  wake oscillator stays at 0.9-1.15 Hz.
- **Slow sideways wandering:** with Ulveseter et al.'s set the light riser at 0.4-0.5 m/s, where
  it lifts off the floor, and the heavy riser at 1.0-1.3 m/s in negative current drift sideways
  at under 0.3 Hz (A/D_rms 2-7.3 with the drift, 0.5-0.7 without); the article does not report
  it.
- **The acceleration spectrum** (`model_f_acc_…`), the article's measured quantity, peaks at a
  higher harmonic of the displacement's peak in 2-5 of the 13 runs of a point with Ulveseter et
  al.'s set (its higher-harmonic force, C_hh 1.5) and in 8 of 13 with the wake oscillator in
  light negative current, at most 2 with the generic set; the measured acceleration spectra
  peak where the displacement's do (Fig. 4).

### Inclination of the straight part in negative current

r.m.s. difference of the inclination at mid-depth from Fig. 12(a), and the median drag
coefficient the inclination implies by Eq. 9 (measured: 1.84 for m* 1.12 and 2.47 for m* 1.84):

| Set | m* 1.12 | m* 1.84 |
| --- | --- | --- |
| `none` (C_d 1.2) | 18.6°; C_d 0.80 | 10.1°; C_d 1.33 |
| `uls19` | 15.0°; 0.90 | 10.1°; 1.36 |
| `b1` | 14.7°; 0.94 | **5.0°; 2.55** |
| `wake` | 15.7°; 0.92 | 10.3°; 1.36 |

- For the heavy riser the generic set's drag amplification gives the measured inclination and
  drag: the article's C_d 2.0-2.8 follows from the response, not from an input. The sets with
  C_d 1.2 differ by 10° r.m.s., steeper at the higher speeds (34° against 25° at 1.2 m/s).
- The light riser's model inclination at mid-depth is not a tension-free straight line: the run
  without VIV, C_d 1.2 by input, reads 0.80 through Eq. 9. Its comparison is a profile check,
  not a drag check; every set differs by 15-19° r.m.s., steeper on average.

### Sensitivity and convergence

- **Convergence** (`uls19`, `b1`, `wake` at m* 1.84, +0.8 m/s and m* 1.12, −1.2 m/s): half the
  time step changes the mean tension by at most 0.6 % except one run (+5.3 %, `b1`, light,
  −1.2 m/s, whose frequency moves by −11 %); the standard deviation by at most 5.5 % and the
  largest A/D_rms by at most 7.9 %. Twice the elements change the mean
  tension by at most 2.9 %, its standard deviation by at most 17 % (`uls19`, light, −1.2 m/s)
  and the frequency by at most one bin.
- **Without structural damping** (36 runs at 0.4, 0.8 and 1.2 m/s): the mean tension changes by
  −1.3 to +15.5 % (above 5 % in five runs, all in negative current); the standard deviation of
  the tension by up to +296 % (the wake oscillator, light, +0.8 m/s), so the 3 % damping
  matters for the fluctuation and not for the mean; 12 of the 54 compared frequencies change.
- **Frictionless floor** (36 runs): the mean tension changes by −3 to +40 % (above 5 % in five
  runs, the largest the wake oscillator, light, −0.4 m/s, where the riser lifts off); the
  standard deviation by up to +412 %, the laid part sliding sideways; 34 of the 54 frequencies
  change, by −15 % to +150 %. The friction coefficient, not given, is the largest modelling
  uncertainty of the case after the still-water geometry.
- The default command, run into a new work directory, reproduces `response_comparison.csv` and
  `summary.json` byte for byte.

## Limitations

- The still-water top tension of the heavy riser is not reproduced (+37 %), so the absolute
  tension carries that offset; the article's numbers do not close the geometry.
- The frequencies are read from raster figures in which many markers overlap: in Fig. 9(a) most
  of the light riser's points are hidden (12 read), and the accelerometer positions are assumed.
- The floor friction and the representation of the 3 % structural damping are assumptions, and
  the tension fluctuation depends on both.
- The current is uniform to the floor; the flume's boundary layer is not modelled, and the laid
  part sees the full current.
- One 80 s record per speed; the article does not give the scatter of its frequencies.
