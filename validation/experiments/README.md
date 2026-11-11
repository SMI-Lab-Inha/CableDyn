<!-- SPDX-License-Identifier: Apache-2.0 -->

# Experimental comparisons

This directory lets reviewers reproduce CableDyn's comparisons with three published physical
experiments. Each subdirectory gives the data source, the translation of the published
parameters into a CableDyn model, the assessed quantity, the reproduction commands and the
committed results.

| Directory | Experiment | Model |
| --- | --- | --- |
| [`holcombe/`](holcombe/README.md) | Holcombe et al. (2025), 1:70 static lazy-wave cable | Finite-EI power cable |
| [`bergdahl/`](bergdahl/README.md) | Bergdahl et al. (2016), 33 m dynamically scaled chain | Tension-only mooring chain |
| [`chaplin/`](chaplin/README.md) | Chaplin et al. (2005), 13.12 m vertical tension riser, still-water natural frequencies | Finite-EI riser, modal analysis |

The results are summarised as cases L3-H, L3-B and L3-C in [`VALIDATION.md`](../../VALIDATION.md).
