<!-- SPDX-License-Identifier: Apache-2.0 -->

# Same-host computational performance

This page lists the wall times and timing definitions behind the lazy-wave performance
comparison. The protocol is in [README.md](README.md).

All runs used one thread of an Intel Core i9-13900K. Each value is the median of three measured
runs after one warm-up; the range of the three runs is in brackets. The values come from
`performance_summary.csv`.

| Water depth | CableDyn | OrcaFlex | Refined MoorDyn-F |
| ---: | ---: | ---: | ---: |
| 80 m | 1.55 s (1.55–1.84), 240 elements, 0.05 s | 2.29 s (2.26–2.31), 179 segments, 0.05 s | 249.07 s (249.04–363.90), 208 segments, 0.125 ms |
| 200 m | 4.64 s (4.62–4.69), 480 elements, 0.05 s | 4.56 s (4.51–4.74), 350 segments, 0.05 s | 144.61 s (142.70–214.11), 224 segments, 0.125 ms |
| 800 m | 7.37 s (7.37–7.54), 960 elements, 0.05 s | 23.17 s (23.00–23.42), 1227 segments, 0.05 s | 243.08 s (240.15–249.30), 232 segments, 0.5 ms |

The MoorDyn-F models are refined fourfold, and their time step is the largest at which every
reported response changes by less than 0.01 % when the step is halved. That accuracy
requirement is much stricter than the one normally used in design work, and it sets the
explicit time step. The figures are specific to this host, these models and these settings.

## Timing definitions

- CableDyn and OrcaFlex totals comprise static initialisation and the 60 s dynamic march.
- The MoorDyn-F total is the complete dynamic-driver call, including its dynamic-relaxation
  initialisation. The separately timed MoorDyn-F static-only call is a diagnostic and is not
  added a second time.

CableDyn used 34, 16 and 8 line-search backtracking trials over the 1,200 timed intervals at
80, 200 and 800 m, respectively. No timed interval required recovery subdivision.
