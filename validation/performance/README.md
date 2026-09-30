<!-- SPDX-License-Identifier: Apache-2.0 -->

# Lazy-wave cable performance comparison

This directory records the same-host wall-time comparison of CableDyn, OrcaFlex and MoorDyn-F
on the three spatially resolved Lozon lazy-wave cables. Use it to check the timing figures in
[`VALIDATION.md`](../../VALIDATION.md).

| File | Contents |
| --- | --- |
| `performance_summary.csv` | Per-site, per-solver mesh size, time step, wall and CPU times, memory and Newton statistics |
| `performance_metrics.json` | The same rows in JSON |
| [`PERFORMANCE_AUDIT.md`](PERFORMANCE_AUDIT.md) | Wall times with their ranges, and timing definitions |
| [`ORCAFLEX_MESH_CONVERGENCE.md`](ORCAFLEX_MESH_CONVERGENCE.md) | 800 m OrcaFlex segment-length study; data in `orcaflex_mesh_convergence.csv` |

## Protocol

- Models: matched 80, 200 and 800 m suspended spans with a pinned touchdown point, no grounded
  tail and no seabed contact.
- Timed work: static initialisation and a 60 s dynamic march.
- Host: one Intel Core i9-13900K processor thread; median of three runs after one warm-up, with
  the range retained.
- CableDyn and OrcaFlex use implicit 0.05 s steps.
- MoorDyn-F finite-bending models are refined fourfold and their explicit step is halved until
  every reported response changes by less than 0.01%. The retained steps are 0.125 ms at 80 and
  200 m and 0.5 ms at 800 m. These are below the shortest-segment axial transit time because
  the resolved bending modes of the contact-free cable, not seabed penalty stiffness, set the
  explicit stability limit.
- The 800 m OrcaFlex segment length, 0.75 m, is the converged value from
  [`ORCAFLEX_MESH_CONVERGENCE.md`](ORCAFLEX_MESH_CONVERGENCE.md).

## Provenance

The rows record the host (one thread of an Intel Core i9-13900K), the mesh, the time step and the
timing statistics of each run. The OrcaFlex segment-length study used OrcaFlex 11.6d
([`ORCAFLEX_MESH_CONVERGENCE.md`](ORCAFLEX_MESH_CONVERGENCE.md)). The executables of the three
codes, their build identifiers and the measurement dates were not recorded with the rows, and
the timing models and driver scripts are not included in the repository. Treat the figures as
indicative of this host and these settings; the release timings, with the code versions and
builds used, are in [`PERFORMANCE_0_1_0.md`](../PERFORMANCE_0_1_0.md).

## Result

Total times at 80, 200 and 800 m are 1.55, 4.64 and 7.37 s for CableDyn, 2.29, 4.56 and
23.17 s for OrcaFlex, and 249.07, 144.61 and 243.08 s for the refined MoorDyn-F models. The
MoorDyn-F times reflect the explicit step that the fourfold-refined bending mesh requires under
this protocol. The figures are specific to this host, these models and these settings, and do
not apply to the separate tension-only chain benchmark in `VALIDATION.md`.
