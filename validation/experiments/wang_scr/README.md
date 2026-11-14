<!-- SPDX-License-Identifier: Apache-2.0 -->

# Steel catenary riser in uniform current: cross-flow VIV

This directory compares CableDyn's cross-flow VIV with the SJTU model test of a truncated steel
catenary riser in a 0.2 m/s uniform current: Wang et al. (2014),
<https://doi.org/10.1016/j.marstruc.2014.07.002>, as Thorsen et al. (2017),
<https://doi.org/10.1016/j.marstruc.2016.10.007>, use it to verify the time-domain synchronisation
model that CableDyn's cross-flow VIV implements. The model runs with their coefficient set for this
case; nothing is fitted to the measurements.

## Model and data

- **Riser (Thorsen et al. 2017, Table 1):**
  - length 23.71 m, between universal joints (pinned);
  - the top 21.0425 m horizontally and 9 m vertically from the lower end, with a 2.626 m flowline on
    the seabed;
  - D 0.024 m, 0.69 kg/m (a solid section: no internal water is added), EI 10.5 N m², EA 6.66e5 N.
- **Seabed:** 1 (kN/m)/m in compression (`kBot` = 1000/D), friction 0.2.
- **Current:** 0.2 m/s along the riser plane, from the lower end towards the top.
- **VIV (their set for this case):** C_v 0.7 (their choice for Re ≈ 3000), f̂0 0.18, Δf̂ 0.08,
  cross-flow only; C_D 1.2, C_M 2.
- **Measurements:** `measured_rms_strain.csv` holds the 25 measured r.m.s. dynamic cross-flow strains
  at the outer surface (Thorsen et al. 2017, Fig. 17), read from the vector figure. The axes are
  calibrated on its tick marks; as a check, the simulated maximum of the same figure reads 99.6 % of
  the measured one, as their text states. The measured dominant strain frequency 8.25 m from the
  lower end is 1.00 Hz (their Fig. 18).
- **Assumed:** the top joint 0.5 m below the surface (not given; the current is uniform).
- **CableDyn:** 200 elements, dt 0.005 s, 150 s, statistics over t > 30 s. The cross-flow curvature
  is the section component with no static mean (the in-plane one carries the catenary curvature);
  strain = κ D/2.

## Reproduce

From the repository root, with `path/to/cabledyn.exe` the release executable or a local build:

```powershell
python validation/experiments/wang_scr/wang_scr_uniform_current.py `
  path/to/cabledyn.exe validation/experiments/wang_scr/work `
  --summary validation/experiments/wang_scr
```

`--segments` and `--dt` change the mesh and the time step.

## Results

| | CableDyn | Measured | Thorsen et al. (2017) | Davies et al. (2026) |
| --- | --- | --- | --- | --- |
| Dominant cross-flow strain frequency, 8.25 m from the lower end | **0.992 Hz** (±0.008) | 1.00 Hz | 1.02 Hz | 1.04 Hz |
| Maximum r.m.s. strain / measured maximum | **1.073** | 1 | 0.996 | - |
| Mean ratio at the 25 stations | 1.036 | 1 | - | - |
| r.m.s. peaks along the span | 7 | 6 | 7 | - |

- **Frequency:** within the 0.008 Hz resolution of the measured value.
- **Strain magnitude:** the largest r.m.s. strain is 7 % above the measured largest, and the mean
  over the stations 4 % above.
- **Spatial pattern:** the predicted response has one r.m.s. peak more along the span than the test
  (7 against 6), as in Thorsen et al. (2017). So the peaks and troughs fall at different stations,
  and the station-by-station ratio scatters (`rms_strain_comparison.csv`, r.m.s. log error 0.34).
- **Convergence:** 400 elements and half the time step change the maximum ratio by 0.2 % and 0.1 %
  and leave the frequency unchanged (`provenance.json`).
