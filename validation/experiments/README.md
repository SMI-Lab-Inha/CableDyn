<!-- SPDX-License-Identifier: Apache-2.0 -->

# Experimental comparisons

This directory lets reviewers reproduce CableDyn's comparisons with eight published physical
experiments. Each subdirectory gives the data source, the translation of the published
parameters into a CableDyn model, the assessed quantity, the reproduction commands and the
committed results.

| Directory | Experiment | Model |
| --- | --- | --- |
| [`holcombe/`](holcombe/README.md) | Holcombe et al. (2025), 1:70 static lazy-wave cable | Finite-EI power cable |
| [`bergdahl/`](bergdahl/README.md) | Bergdahl et al. (2016), 33 m dynamically scaled chain | Tension-only mooring chain |
| [`chaplin/`](chaplin/README.md) | Chaplin et al. (2005), 13.12 m vertical tension riser, still-water natural frequencies | Finite-EI riser, modal analysis |
| [`delta_flume/`](delta_flume/README.md) | Chaplin et al. (2005), 13.12 m vertical tension riser, stepped current, blind-trial cases 1–9 | Finite-EI riser, synchronisation VIV, cross-flow and in-line |
| [`wang_scr/`](wang_scr/README.md) | Wang et al. (2014), truncated steel catenary riser in uniform current and under heave motion of its upper end, cross-flow VIV | Finite-EI riser, synchronisation VIV |
| [`ndp_riser/`](ndp_riser/README.md) | Trim et al. (2005), NDP 38 m high-mode riser in uniform flow, cross-flow VIV with a generic and three published coefficient sets | Finite-EI riser, synchronisation VIV |
| [`flotant/`](flotant/README.md) | de Wilde et al. (2021), FLOTANT OD 48 mm power-cable specimen with a clump weight, tow tests and forced-oscillation RAOs, cross-flow and in-line VIV with three published coefficient sets | Finite-EI cable, synchronisation and wake-oscillator VIV |
| [`fu_oscillatory/`](fu_oscillatory/README.md) | Fu et al. (2014), SJTU 4 m flexible cylinder moved harmonically through still water, cross-flow VIV on the acceleration and deceleration paths with three published coefficient sets | Finite-EI cylinder, synchronisation and wake-oscillator VIV, `VIV FRAME` |
| [`chaplin_catenary/`](chaplin_catenary/README.md) | Chaplin & King (2018), 5.36 m untensioned catenary riser with high curvature and floor touchdown in current of both directions, two mass ratios: mean top tension, cross-flow frequency, inclination, with three published coefficient sets | Finite-EI riser with floor contact, synchronisation and wake-oscillator VIV |

The results are summarised as cases L3-H, L3-B, L3-C, L3-D, L3-W, L3-Wh, L3-N, L3-F and L3-K in [`VALIDATION.md`](../../VALIDATION.md).
