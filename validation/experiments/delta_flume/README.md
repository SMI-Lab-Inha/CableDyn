<!-- SPDX-License-Identifier: Apache-2.0 -->

# Delta Flume vertical riser in a stepped current: cross-flow and in-line VIV

This directory compares CableDyn's VIV with the nine blind-trial cases of the Delta Flume riser
tests. Chaplin et al. (2005a), <https://doi.org/10.1016/j.jfluidstructs.2005.04.010>, describe the
riser and the tests, and Chaplin et al. (2005b),
<https://doi.org/10.1016/j.jfluidstructs.2005.05.016>, the blind predictions made for them. The
lower 5.94 m of the 13.12 m riser is in a uniform current of 0.16–0.95 m/s, the rest in still
water. The VIV coefficients are a generic set; nothing is fitted to these measurements. The
still-water natural frequencies of the same riser are case L3-C
([`../chaplin`](../chaplin/README.md)).

## Model and data

- **Riser:** as L3-C: 13.12 m between pinned joints, D 28 mm, EI 29.9 N m², mass ratio 3.0 with
  the internal water (1.847 kg/m), fresh water, C_a 1. From the project's mechanical-properties
  page (University of Southampton, archived; `provenance.json`): EA 5.88 MN, structural damping
  0.33 % of critical, top spring 11.4 kN/m.
- **Top spring:** in series with the riser, one effective axial stiffness
  EA_eff = L/(1/k + L/EA) = 1.459e5 N. The unstretched length of each case is iterated until the
  static top tension in the current is the case's mean top tension. A rigid top would turn the
  current's bow into tension (1799 N instead of 1002 N in case 9) and stiffen the riser.
- **Current:** the towing speed U over the lower 5.94 m (the bottom of the carriage's vacuum
  tank), still water above, with a 2 cm linear step between them.
- **Cases (slide 3 of the presentation below):**

  | Case | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
  | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
  | U (m/s) | 0.16 | 0.21 | 0.31 | 0.40 | 0.54 | 0.60 | 0.70 | 0.85 | 0.95 |
  | Mean top tension (N) | 405 | 407 | 457 | 506 | 598 | 670 | 743 | 923 | 1002 |
  | Std of top tension (N) | 2 | 3 | 3 | 8 | 15 | 14 | 8 | 14 | 65 |

- **VIV (all cases):**
  - synchronisation lift C_l 0.8 in the band f̂0 0.18 ± 0.08;
  - in-line force C_x 0.85, α 0.15 (Ulveseter et al. 2018);
  - C_D 1.2;
  - a `VIV AMPLITUDE` table with C_v 0.8 constant and the drag amplification of Chaplin et al.
    (2005a, Eq. 4), CdAmp = 1 + 2.3 s + 1.7 s², s = (A/D)/√2, at A/D 0–1.5;
  - `VIV DAMPING` `0 0 1`: the amplification acts on the in-line drag only, and the cross-flow drag
    stays Morison's.
- **Measurements:** the experimenters' processed results, from their presentation of the blind
  trial (Chaplin et al., Trondheim, October 2004). The slides are vector graphics, so the values are
  the plotted coordinates mapped through the tick labels, not a pixel digitisation:
  - `measured_summary.csv`, from slides 17–19: maximum cross-flow displacement, maximum mean and
    maximum dynamic in-line displacement (in D), and the overall standard deviation of the in-line
    and cross-flow curvature (κD). Its `U` column is the speed read from each point; it matches
    the tabulated speed within 0.006 m/s, which checks the axis mapping.
  - `envelopes/`: the measured envelopes along the riser of each case (slides 7–15), with their
    columns, checks and the assumed curvature scale in [`envelopes/README.md`](envelopes/README.md).
  - The dominant cross-flow mode numbers 2, 3, 4, 5, 6, 6, 7, 7, 7 (slide 23).
- **Assumed:** the 11.4 kN/m top spring (the page also lists 38.1 kN/m, without saying which tests
  used which); the 0.33 % structural damping as bending damping of a pure bending mode at 2 Hz (the
  page gives no frequency).
- **CableDyn:** 64 elements; dt the largest of 5, 2.5, 2, 1 and 0.5 ms at or below 1/200 of the
  shedding period D/(0.18 U) (2.5 ms in cases 1–3, 1 ms in 4–7, 0.5 ms in 8–9); `VIVrmsWindow`
  0.25 s; 40 s, statistics over t > 15 s.
- **Metrics:** x is in-line and y cross-flow, both from the line between the joints.
  - dominant cross-flow mode: the pinned-pinned sine with the largest r.m.s. projection of y;
  - maximum cross-flow amplitude: the largest half peak-to-peak of y along the riser, against the
    slide-18 maximum;
  - maximum mean in-line displacement: the largest time mean of x;
  - dominant cross-flow frequency (`response_comparison.csv` only): the spectral peak of y a
    quarter of the length above the lower joint, at a resolution of 0.04 Hz;
  - overall standard deviation of curvature: the r.m.s. over the nodes of the standard deviation
    of the curvature (second differences of the node positions), times D.

## Reproduce

From the repository root, with `path/to/cabledyn.exe` the release executable or a local build:

```powershell
python validation/experiments/delta_flume/delta_flume_viv.py `
  path/to/cabledyn.exe validation/experiments/delta_flume/work `
  --summary validation/experiments/delta_flume
```

The nine cases run in parallel (`--jobs`), in about two minutes on a 24-core workstation. `--cases`
(for example `5,9`) selects cases; `--segments`, `--dt-scale` and `--spring` change the mesh, the
time step and the top spring.

## Results

CableDyn / measured; the measured cross-flow amplitude is the slide-18 maximum, with the half
peak-to-peak of the measured envelope in brackets.

| Case | U (m/s) | Mode | A_CF/D | Mean IL/D | Std κD, CF | Std κD, IL | Mean top tension (N) |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 0.16 | 2 / 2 | 0.653 / 0.757 (0.763) | 1.19 / 1.06 | 5.7e-5 / 9.3e-5 | 3.5e-5 / 5.0e-5 | 409 / 405 |
| 2 | 0.21 | 3 / 3 | 0.602 / 0.579 (0.561) | 1.60 / 1.19 | 1.11e-4 / 0.89e-4 | 7.2e-5 / 3.4e-5 | 413 / 407 |
| 3 | 0.31 | 4 / 4 | 0.628 / 0.779 (0.809) | 3.27 / 3.45 | 2.05e-4 / 2.70e-4 | 0.92e-4 / 1.95e-4 | 476 / 457 |
| 4 | 0.40 | **4** / 5 | 0.700 / 0.674 (0.595) | 4.78 / 4.97 | 2.30e-4 / 3.18e-4 | 1.36e-4 / 2.69e-4 | 541 / 506 |
| 5 | 0.54 | 6 / 6 | 0.658 / 0.843 (0.889) | 6.53 / 6.69 | 4.45e-4 / 4.81e-4 | 2.00e-4 / 2.31e-4 | 660 / 598 |
| 6 | 0.60 | 6 / 6 | 0.679 / 0.709 (0.779) | 6.94 / 8.13 | 4.55e-4 / 4.77e-4 | 2.32e-4 / 4.57e-4 | 738 / 670 |
| 7 | 0.70 | **6** / 7 | 0.738 / 0.777 (0.717) | 8.28 / 8.90 | 4.63e-4 / 6.47e-4 | 2.69e-4 / 3.03e-4 | 835 / 743 |
| 8 | 0.85 | **6** / 7 | 0.769 / 0.907 (0.890) | 9.66 / 11.18 | 4.94e-4 / 6.89e-4 | 2.68e-4 / 4.48e-4 | 1046 / 923 |
| 9 | 0.95 | **8** / 7 | 0.595 / 1.074 (0.996) | 10.15 / 13.51 | 6.87e-4 / 7.36e-4 | 3.67e-4 / 6.24e-4 | 1128 / 1002 |

- **Dominant cross-flow mode:** right in 5 of 9 cases and within one in the other four (one low in
  cases 4, 7 and 8, one high in case 9). Four of the eleven blind predictions found 7 or 8 of the 9
  (slide 24).
- **Mean in-line displacement:** within ±17 % in 7 of 9 cases; +34 % in case 2, where the
  measured maximum is 1.2 D, and −25 % in case 9. The blind predictions under-predicted the drag
  amplification, and with it the mean in-line displacement, by 20–40 % (slide 26).
- **Cross-flow amplitude:** −45 % to +4 % against the slide-18 maxima (−40 % to +18 % against
  the envelopes), within ±12 %, the scatter between the two measured values, in cases 2, 4, 6
  and 7. The −45 % of case 9 goes with the response in mode 8.
- **Curvature:** the cross-flow standard deviation is −39 % to +26 %, low in every case but 2
  (−5 % to −28 % in cases 3–9); the in-line one −53 % to +111 %. The blind predictions' r.m.s.
  curvatures were 20 % to 180 % of the measured (slide 26).
- **Top tension:** the static top tension in the current is the measured mean, but the drag
  amplification bows the riser further: the mean top tension of the run is 1–13 % above the
  measured mean (7–13 % in cases 4–9), so the run is stiffer than the test at the higher
  speeds.
- **Convergence (`provenance.json`):**
  - cases 1–4, 7 and 8: half the time step leaves the mode and changes the mean in-line
    displacement and the cross-flow curvature standard deviation by at most 0.6 %, the peak
    amplitude by at most 2.7 % and the in-line curvature by at most 1.2 %; 128 elements leave the
    mode and change them by at most 2.7 %, 5.9 % and 3.1 %, except the in-line curvature of case
    8 (+19 %);
  - case 5: half the time step and 128 elements change its in-line curvature by +22 % and +25 %
    and the other metrics by at most 6.4 %; an 80 s run raises its peak amplitude by 2.6 %;
  - case 6: its peak amplitude depends on the step and the mesh, −11 % at half the time step
    (kept over an 80 s run: 0.604 against 0.679) and −13 % with 128 elements, while the mean
    in-line displacement and the cross-flow curvature move by at most 2 % and the in-line
    curvature by at most 12 %;
  - case 9 has two responses: at half the time step it settles in mode 7 (5.80 Hz, A/D 0.751,
    mean in-line displacement 10.68 D), with the base step and with 128 elements in mode 8
    (6.60 and 6.64 Hz), which an 80 s run keeps.
- **Top spring 38.1 kN/m (`--spring 38.1e3`):** the mean top tension of the run is up to 14 %
  higher and the mean in-line displacement 2–17 % lower, within ±17 % in 4 of 9 cases; the mode
  is right in 5 of 9 (case 9 in mode 7, case 5 in mode 5, case 8 in mode 5).

`response_comparison.csv` holds the per-case values with their errors, and `summary.json` the
statistics and the SHA-256 of the executable.
