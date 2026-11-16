<!-- SPDX-License-Identifier: Apache-2.0 -->

# NDP 38 m riser in uniform flow: high-mode cross-flow VIV

This directory compares CableDyn's VIV with the bare-riser, uniform-flow tests of the Norwegian
Deepwater Programme's high-mode riser: Trim et al. (2005),
<https://doi.org/10.1016/j.jfluidstructs.2005.07.014>. A 38 m riser of 27 mm diameter was towed
horizontally at 22 speeds from 0.3 to 2.4 m/s and responded at up to 12.5 Hz cross-flow. The
same riser model runs with four coefficient sets. One is generic and not fitted to these tests.
The other three were fitted by their authors to NDP data, two of them to these very tests, so
for those the comparison checks CableDyn's implementation of the model, not its predictive
power.

## Model and data

- **Riser (Fontaine et al. 2011, Table 1; Thorsen et al. 2016, Table 1):** 38.00 m between
  universal joints (pinned), D 27 mm, EI 598.8 N m², EA 8.19e6 N, 0.933 kg/m flooded, as tested.
- **Water:** fresh, 1000 kg/m³. The mass ratio 1.62 of Fontaine et al. is 0.933 kg/m over the
  displaced water at 1006 kg/m³; 1000 kg/m³ gives 1.63 and sea water 1.59.
- **Tension (Trim et al. 2005, Sec. 2.1):** clump weights on inclined pendulum arms held the
  riser at about 4–6 kN over the test speeds. The tension of each test is not published
  (Ulveseter et al. 2019 give the range 4067–6285 N). The rig is modelled as an axial spring
  k_ax in series with the riser, EA_eff = L/(1/k_ax + L/EA), with T0 = 4 kN in still water.
  - k_ax = 3.5 kN/m comes from the reported range, not from the response: at 2.4 m/s the Morison
    drag (C_D 1.2) of about 93 N/m sags a 6 kN string by wL²/(8T) ≈ 2.8 m, which takes
    8s²/(3L) ≈ 0.55 m more length, and 2 kN over 0.55 m is about 3.6 kN/m.
  - The resulting mean tension of each set is in the table below.
- **Assumed:** the riser 5 m below the surface (not given; the current is uniform), in 10 m of
  water.
- **Current:** uniform, normal to the riser and horizontal, so cross-flow (CF) is vertical and
  in-line (IL) horizontal.
- **CableDyn:** one finite-EI line, 200 elements (0.19 m), dt 1 ms, 25 s, statistics over
  t > 10 s (frequency resolution 0.067 Hz); `VIVrmsWindow` 0.5 s.
- **Analysis (`ndp_riser_viv.py`):**
  - `CurvX` at every node is the CF curvature and `CurvY` the IL one: the section axis x starts
    at End A as global +z projected normal to the line and is carried along it, which keeps it
    within 2° of vertical on the bowed riser. Strain = κ D/2.
  - The dominant CF frequency is the peak of the CF curvature spectra summed over the nodes,
    searched for in f D/U = 0.05–0.30. This covers the synchronisation band of every set
    (0.065–0.26) and excludes a slow swing of the riser below it and the third harmonic of the
    Ulveseter et al. set above it. The IL frequency is searched for in 0.10–0.60.
  - The fundamental CF strain is band-passed to 0.5–1.5 times the dominant CF frequency, as
    Thorsen et al. (2016) filter the measured strain. The compared value is its largest r.m.s.
    along the riser.
  - σ/D is the spatial mean, over every fourth node, of the standard deviation of the
    displacement band-passed in the same way, as Trim et al. reconstruct the vibration from
    modes.

### Coefficient sets

`--sets` selects them. C_a = 1 (C_M = 2) in all.

| Set | Source | Values | Fitted to |
| --- | --- | --- | --- |
| `b1` | The generic set of the [Delta Flume case](../delta_flume/README.md) | C_l 0.8, f̂0 0.18 ± 0.08, C_D 1.2; in-line C_x 0.85, α 0.15 (Ulveseter et al. 2018, Table 3); `VIV AMPLITUDE` C_v 0.8 with the drag amplification of Chaplin et al. (2005, Eq. 4); `VIV DAMPING` `0 0 1` (the amplification on the in-line drag only) | nothing from the NDP tests: blind |
| `and24b` | Andersen et al. (2024b), Table 3, cross-flow average | C_v,y 1.28, f̂0,y 0.131, Δf̂y 0.066, C_D 1.16; Morison drag | the 22 sheared-flow tests of this riser: semi-blind here. Their in-line set (C_v,x 1.09, f̂0,x 0.321, Δf̂x 0.161) belongs to an independent in-line synchronisation, not CableDyn's `VIV` row, and is left out |
| `uls19` | Ulveseter et al. (2019), Table 1 and Sec. 4.1 | C_D 1.2, C_v,y 0.85, f̂0 0.144 ± 0.064, in-line C_v,x 0.75, α 0.15, C_hh 1.5; Morison drag | these 22 uniform-flow tests ("a visual best fit with cross-flow strain measurements"): an implementation check |
| `th16` | Thorsen et al. (2016), Fig. 7 and Eq. 1 | `VIV AMPLITUDE` with their fitted C_v(y0/D) (`thorsen2016_fig7_cv.csv`); `VIV DAMPING` `0.31 0.89 0`, their cross-flow damping without Morison's linear term; band 0.18 ± 0.08 for the 0.10–0.26 range of their Fig. 4; C_D 1.2 in line | these 22 uniform and the 22 sheared-flow tests: an implementation check. Their model is cross-flow only |
| `th16c` | `th16` with Morison's linear cross-flow damping kept (`VIV DAMPING` `0.31 0.89 1`) | as `th16` | a variant, not a published set |

The published values were checked against the source articles. CableDyn's synchronisation
frequency is f̂0 + Δf̂ sin(...), so `th16` matches the range of Thorsen et al.'s synchronisation
function (their Fig. 4), not its shape.

### Measurements

- `thorsen2016_fig8_strain.csv`: the maximum r.m.s. CF strain of each test, at the fundamental
  frequency and in total, with Thorsen et al.'s own simulation using their fitted curve (their
  Fig. 8; the uniform panel is compared, the sheared one is kept). The figure is vector
  graphics: `extract_thorsen2016.py` maps the plotted coordinates through the plot frames, whose
  ends are the tick labels 0 and 2.5 m/s and 0 and 5e-4. Each series has 22 points, one per tested
  speed. `thorsen2016_fig7_cv.csv` is their fitted C_v curve from the same page (0.35 at 0, a
  maximum of 1.00 at y0/D 0.46, zero from 0.67).
- `measured_trim2005_fig16.csv`: the spatial mean of σ/D, CF and IL, at all 22 speeds (Trim et
  al. 2005, Fig. 16). `measured_trim2005_fig18.csv`: the peak frequency of the dominant
  displacement mode, CF and IL, at five speeds (Fig. 18). The figures are raster images and their
  markers overlap, so the values were read by eye on four-times enlarged crops with a guide line
  at every speed and every 0.01 or 1 Hz (`digitise_trim2005.py` calibrates the axes on the
  panel frames and grid lines and draws the guides). Uncertainty: ±0.01 in σ/D, about ±0.02 where
  markers of neighbouring speeds overlap; ±0.3 Hz, ±0.5 Hz for CF at 2.0 m/s. The measured CF
  frequency follows f D/U ≈ 0.14 and the IL one twice that.

The two extraction scripts need the articles (not distributed), PyMuPDF and Pillow; the
comparison needs only the CSV files.

## Reproduce

From the repository root, with `path/to/cabledyn.exe` the release executable or a local build:

```powershell
python validation/experiments/ndp_riser/ndp_riser_viv.py `
  path/to/cabledyn.exe validation/experiments/ndp_riser/work
```

The default runs the five sets at 0.5, 1.0, 1.5, 2.0 and 2.4 m/s, 25 runs in parallel, in
about two minutes on a 24-core workstation, and prints the statistics over those speeds; each
run reproduces its row of the record. The committed record is the full sweep, 110 runs in about
twelve minutes with 24 parallel runs:

```powershell
python validation/experiments/ndp_riser/ndp_riser_viv.py `
  path/to/cabledyn.exe validation/experiments/ndp_riser/work `
  --speeds all --jobs 24 --summary validation/experiments/ndp_riser
```

`--sets` and `--speeds` take comma lists (quote them in PowerShell, `--speeds "1.0,2.4"`);
`--segments`, `--dt`, `--tension` and `--kax` change the mesh, the time step and the rig. A run
whose folder holds the same deck, executable and statistics window is not repeated; `--keep`
keeps the output tables (about 150 MB a run) so that the analysis can be redone without them.

## Results

All 22 speeds; the strain ratio is CableDyn's fundamental over the measured fundamental, and
the r.m.s. log error is √(mean (ln ratio)²).

| Set | Fitted to these tests | Strain ratio, median (range) | r.m.s. log error | U ≥ 0.5 m/s | f D/U, CF | σ/D, CF | σ/D, IL | Mean tension (kN) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Measured | | 1 | 0 | 0 | 0.143 | 0.436 | 0.149 | 4–6 |
| `b1` | no | 1.26 (0.99–1.57) | 0.277 | 0.277 | 0.192 | 0.429 | 0.041 | 4.0–10.0 |
| `and24b` | sheared tests | 1.00 (0.67–1.27) | 0.186 | 0.181 | 0.134 | 0.687 | 0.088 | 4.0–7.5 |
| `uls19` | yes | 0.89 (0.64–1.33) | **0.178** | **0.156** | **0.145** | 0.461 | **0.146** | 4.0–6.9 |
| `th16` | yes | 1.40 (0.91–1.86) | 0.363 | 0.377 | 0.193 | 0.390 | 0.067 | 4.0–6.5 |
| `th16c` | variant | 0.73 (0.51–1.06) | 0.365 | 0.376 | 0.153 | 0.307 | 0.025 | 4.0–6.4 |
| Thorsen et al.'s own simulation (their Fig. 8a) | yes | 1.12 (0.75–1.51) | 0.187 | 0.182 | | | | |

f D/U is the mean at the five speeds of Fig. 18; σ/D the mean over the 22 speeds.
`response_comparison.csv` holds every run, `summary.json` the statistics and the SHA-256 of the
executable.

- **Implementation check:** with the set Ulveseter et al. fitted to these tests, CableDyn's
  three-dimensional model matches every metric: the strain with an r.m.s. log error of 0.178
  (0.184 at half the time step), comparable to Thorsen et al.'s own cross-flow-only fit (0.187), f
  D/U 0.145 against the measured 0.143, and σ/D within 6 % cross-flow and 3 % in line, with the
  higher-harmonic force (C_hh 1.5) of their set.
- **Semi-blind:** the cross-flow set Andersen et al. fitted to the sheared-flow tests predicts the
  uniform-flow strain as well as Thorsen et al.'s fit (0.186 against 0.187) at f D/U 0.134, but
  with σ/D 57 % high: a lower mode at a larger amplitude gives the same curvature.
- **Blind:** the generic set gives σ/D within 2 % and the strain 26 % high in the median (r.m.s.
  log error 0.277), and responds at f D/U 0.19, a third above the measured frequency. Its
  amplified drag raises the mean tension to 10.0 kN at 2.4 m/s, above the reported 6 kN.
- **The damping form sets the frequency:** with Thorsen et al.'s curve and band, dropping
  Morison's linear cross-flow damping (`th16`, as they do) gives f D/U 0.193; keeping it
  (`th16c`) 0.153, consistent with the higher modes being less damped without it. `th16` also
  develops a slow cross-flow swing, at 0.6–0.8 Hz up to 1.9 m/s (raw σ/D up to 2.0, excluded from
  the band-passed values), the known consequence of dropping the linear term while an in-line
  force acts (theory guide, Cross-flow damping). Neither reproduces the strain of their own fit
  (0.363 and 0.365 against 0.187): their curve was fitted in a cross-flow-only model run at the
  tension of each test, and both differ here.
- **One metric does not discriminate:** `and24b` matches the strain with the wrong amplitude,
  `b1` the amplitude with the wrong frequency, `th16c` the frequency with the wrong strain. Only
  `uls19` matches strain, amplitude and frequency together.
- **In-line:** `b1` and `uls19` also swing slowly in line at 0.3–0.9 Hz, more at the higher
  speeds (raw IL σ/D up to 3.1 and 2.6). The band-passed values exclude the swing; the measured
  σ/D, a modal reconstruction of the vibration, could not show it, and whether the tests had
  one is not reported.
- **Sensitivity to the water density:** with sea water (1025 kg/m³) the r.m.s. log errors from
  0.5 m/s are 0.299 (`b1`), 0.181 (`and24b`), 0.154 (`uls19`), 0.434 (`th16`) and 0.353
  (`th16c`). The change of `th16` (0.377 to 0.434) shows how much the mode chosen at single
  speeds moves the statistic.

### Convergence

- **Single speeds** (1.0 and 2.4 m/s, every set; `provenance.json`):
  - at 1.0 m/s, half the time step and 400 elements change the fundamental strain by at most
    0.8 % and σ/D by at most 3 % for `b1`, `and24b`, `th16` and `th16c`, at the same frequency
    within one bin; `uls19` moves by +1.9 % and +18.5 % in strain and responds at 5.60 Hz instead
    of 5.80 Hz with both (σ/D −2 %);
  - at 2.4 m/s the response is not converged speed by speed: the strain moves by −14 % to +12 %
    with half the time step and by +0.6 % to +20 % with 400 elements, and the dominant frequency
    moves by more than one bin in 7 of the 10 pairs;
  - 55 s runs (45 s of statistics) of `b1`, `and24b` and `th16c` at 2.4 m/s change the strain of
    the 25 s runs by −0.3 %, −1.2 % and +12.5 %, and at 55 s half the time step and 400 elements
    still move it by −0.7 % to +12 % and −7.8 % to +21 %: a longer window does not remove the
    spread. It comes from the mix of neighbouring high modes the response settles in, which a
    small change of the model moves; it is of the size of the speed-to-speed scatter of the
    ratios.
- **The statistics over the 22 speeds hold:** the full sweep at half the time step gives
  r.m.s. log errors of 0.280 (`b1`), 0.187 (`and24b`), 0.184 (`uls19`), 0.375 (`th16`) and
  0.376 (`th16c`), f D/U within 0.001 and σ/D within 1.4 % of the record. Differences of a few
  hundredths in the r.m.s. log error between sets are within this spread.

## Limitations

- The tension of each test is not published; the rig spring gives the reported range, but each
  run's tension differs from the test's by an unknown amount, up to the 10 kN of `b1`.
- The σ/D and frequency data are read by eye from raster figures.
- The 3D model carries in-line drag and slow swings that the 2D cross-flow model of Thorsen et
  al. does not, so their fitted curve is not expected to transfer exactly.
