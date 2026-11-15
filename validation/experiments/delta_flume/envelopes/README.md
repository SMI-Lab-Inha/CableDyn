<!-- SPDX-License-Identifier: Apache-2.0 -->

# Measured response envelopes, blind-trial cases 1–9

`case<N>_envelopes.csv` holds the measured envelopes along the riser of case N, read from the
measurement curves (the top row of slides 7–15) of the experimenters' presentation of the blind
trial of Chaplin et al. (2005b), <https://doi.org/10.1016/j.jfluidstructs.2005.05.016>. The slides
are vector graphics: each measured curve is a polyline of 102 nearly equally spaced vertices from
one end of the riser to the other, so the values are the plotted vertices, linearly resampled to
101 equally spaced stations, not a pixel digitisation.

| Column | Meaning |
| --- | --- |
| `z_over_L` | height above the lower joint over the 13.12 m length: 0 at the lower joint (in the current), 1 at the top joint |
| `x_min_D`, `x_mean_D`, `x_max_D` | in-line displacement (along the current) over the diameter: minimum, mean and maximum; zero is the panel edge that the end stations touch |
| `y_min_D`, `y_max_D` | cross-flow displacement over the diameter: minimum and maximum, about the panel centre |
| `kx_rms_D`, `ky_rms_D` | in-line and cross-flow r.m.s. curvature times the diameter (see below) |

- **Displacement scales:** the in-line panel's `10` tick and the cross-flow panel's `-1` and `1`
  ticks. As a check, the largest `x_mean_D` of each case is within 0.05 D (4 %) of the measured
  maximum mean in-line displacement of slide 17 (`../measured_summary.csv`).
- **Cross-flow scatter:** half the largest peak-to-peak of the cross-flow envelope differs from
  the measured maximum cross-flow displacement of slide 18 by −12 % to +10 %. The experimenters
  processed quasi-steady windows of each record, and slide 18 summarises them; this difference is
  the scale of the measurement scatter in a peak amplitude.
- **Curvature scale (assumed):** the curvature panels carry no scale. Each shape is scaled so that
  its r.m.s. over the 101 stations equals the overall standard deviation of curvature of slide 19
  for the case, assuming that number is the r.m.s. over the span of the r.m.s. curvature. The
  shapes and the positions of their peaks do not depend on the assumption; the magnitudes do.
