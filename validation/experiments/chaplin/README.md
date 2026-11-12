<!-- SPDX-License-Identifier: Apache-2.0 -->

# Delta Flume vertical riser: still-water natural frequencies

This directory compares CableDyn's modal analysis with the still-water natural frequencies of the
13.12 m vertical tension riser that Chaplin et al. (2005) tested in the Delta Flume,
<https://doi.org/10.1016/j.jfluidstructs.2005.04.010>. The riser is the subject of the paper's
stepped-current VIV measurements. Its frequencies at five top tensions (Table 2) check the riser
model, the mass, added mass and tension along the line, which a VIV comparison then builds on. No
model parameter is fitted to the measurements.

## Model and data

- From Section 2 of the paper: length 13.12 m, diameter 28 mm, bending stiffness 29.9 N m², mass
  ratio 3.0 with the internal water (so 1.847 kg/m in fresh water, whose submerged weight
  12.08 N/m matches the paper's 12.1 N/m), universal joints at both ends (pinned in bending). The
  riser was submerged over its whole length, the upper part in the water-filled tank.
- `natural_frequencies.csv` holds the top tensions and the first seven measured frequencies of
  Table 2, typed from the article, next to CableDyn's.
- Assumed: fresh water (1000 kg/m³), the normal added-mass coefficient 1, and an axial stiffness
  of 5.5e6 N, which the paper does not give. The unstretched length of each series gives its
  Table 2 top tension, with the tension falling by the submerged weight down the riser.
- CableDyn: a finite-EI line between two stationary points, 64 elements, `nModes` about the
  static state (TMax 0). The two transverse planes of the straight riser give each frequency
  twice; the script checks that they agree and keeps one.

## Reproduce

From the repository root, with `path/to/cabledyn.exe` the release executable or a local build:

```powershell
python validation/experiments/chaplin/chaplin_natural_frequencies.py `
  path/to/cabledyn.exe validation/experiments/chaplin/work `
  --summary validation/experiments/chaplin
```

`--segments`, `--ca` and `--ea` change the mesh, the added-mass coefficient and the axial
stiffness.

## Results

| Statistic (34 frequencies, series C mode 2 left out) | Value |
| --- | ---: |
| Mean error | +1.87 % |
| Mean absolute error | 2.28 % |
| Largest absolute error | 4.88 % (series A, mode 1) |

- CableDyn is 2–3 % high in most modes. The difference corresponds to a slightly larger added
  mass: with `--ca 1.1` the mean error is +0.62 % and with `--ca 1.2` it is −0.59 %. The default
  keeps the potential-flow value 1.
- Series C mode 2 is measured at 1.423 Hz, where CableDyn gives 1.619 Hz (+13.8 %). Its ratio
  to the first mode, 1.82, departs from the near-harmonic ratios of the other four series
  (1.97–2.16), while the neighbouring C modes agree with CableDyn within 3.3 %. It is reported
  in the CSV and counted in the all-35 statistic (mean absolute error 2.61 %).
- 128 elements change no frequency by more than 3e-6, and ten times the axial stiffness by more
  than 1.2e-4 (relative).

`natural_frequencies_summary.json` records the statistics and the SHA-256 of the executable.
