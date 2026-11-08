<!-- SPDX-License-Identifier: Apache-2.0 -->

# Same-host wall times

This page lists the wall times behind the L5 rows of [`VALIDATION.md`](../../VALIDATION.md) with
their ranges and the timing definitions. The cases, protocol and builds are in
[README.md](README.md); the values come from `performance_summary.csv` and
`wallclock_summary.csv`.

All runs used one thread on an AMD Ryzen Threadripper PRO 7965WX under Windows 11 Pro on
2026-11-08. Each value is the median of three timed runs after one warm-up; the range of the
three runs is in brackets. `IFX` is the v0.4.1 release driver, `GNU` a GNU Release build of `dev`.

## L5-DPC: pinned-touchdown lazy-wave spans, 60 s

| Water depth | CableDyn (IFX) | CableDyn (GNU) | OrcaFlex 11.6d | Refined MoorDyn-F |
| ---: | ---: | ---: | ---: | ---: |
| 80 m | 1.44 s (1.43–1.45) | 2.06 s (2.04–2.06) | 2.69 s (2.68–2.71) | 326.6 s (323.8–343.3) |
| 200 m | 2.56 s (2.54–2.56) | 3.40 s (3.38–3.43) | 5.31 s (5.18–5.64) | 191.9 s (191.7–193.5) |
| 800 m | 5.29 s (5.28–5.48) | 7.59 s (7.47–7.59) | 27.81 s (27.44–28.15) | 320.7 s (319.3–326.9) |
| Mesh, step | 240 / 480 / 960 elements, 0.05 s | as IFX | 179 / 350 / 1227 segments, 0.05 s | 208 / 224 / 232 segments; 0.125 / 0.125 / 0.5 ms |

CableDyn's static initialisation takes 0.11 / 0.20 / 0.34 s (IFX) of these totals; it takes
3.04 / 2.75 / 3.27 linear solves per step, and no interval needs recovery subdivision. OrcaFlex's
statics take 0.32 / 0.30 / 2.65 s. The MoorDyn-F total includes its dynamic relaxation to
`threshIC` 1e-5. Halving the MoorDyn-F step moves the scored curvature and tension by at most
0.0005 %, inside the 0.01 % that selects the step; that accuracy requirement is much stricter
than the one normally used in design work, and it sets the explicit step.

## L5-V: VolturnUS-S chain, one hour

| | CableDyn (IFX) | CableDyn (GNU) | OrcaFlex 11.6d |
| --- | ---: | ---: | ---: |
| Wall time | 4.44 s (4.44–4.45) | 6.18 s (6.16–6.19) | 19.72 s (19.33–19.74) |
| Model | 50 elements, dt 0.1 s, 3600 s | as IFX | 50 segments, dt 0.1 s, 200 s build-up and 3600 s |

CableDyn takes 2.03 Newton iterations per step over the 36,000 steps.

## L5-wc: WD0050 chain, 80 s

| | CableDyn (IFX) | CableDyn (GNU) | MoorDyn-C v2.7.2 |
| --- | ---: | ---: | ---: |
| Wall time | 0.200 s (0.198–0.204) | 0.277 s (0.275–0.406) | 0.482 s (0.476–0.483) |
| Model | 41 elements, dt 0.05 s | as IFX | 41 segments, RK2, dtM 1 ms, coupling step 0.1 s |

MoorDyn-C's stepping loop takes 0.452 s of its total. Over the scored 20 s the fairlead tensions
of the two codes differ by 1.35 % in mean and 1.71 % in maximum.

## Timing definitions

- CableDyn, MoorDyn-F and MoorDyn-C: the wall time of the whole process, from start to exit,
  including input, initial condition, dynamics and output.
- OrcaFlex: `CalculateStatics` plus `RunSimulation`, timed inside a fresh Python process; model
  building, licence checkout and result extraction are outside the timed interval.
- The static and dynamic split of CableDyn comes from a separate run with `CABLEDYN_PROFILE`.
- Peak working set is that of the whole process; for OrcaFlex it includes the Python
  interpreter.

The figures are specific to this host, these models and these settings.
