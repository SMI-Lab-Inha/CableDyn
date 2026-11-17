<!-- SPDX-License-Identifier: Apache-2.0 -->

# SJTU flexible cylinder in oscillatory flow: cross-flow VIV hysteresis

This directory compares CableDyn's cross-flow VIV of a flexible cylinder moved harmonically through
still water with the tests of Fu et al. (2014), <https://doi.org/10.1115/1.4025759>. The relative
flow reverses twice per motion cycle and the shedding frequency sweeps from zero to its largest value
and back, so the case tests what a steady current cannot: the VIV model through flow reversal, and
the response on the accelerating and decelerating halves of each half-cycle. Fu et al. measure
different amplitudes on the two halves (hysteresis) and give the paths for three cases. No model
parameter is fitted to these tests.

## Model and data

- **Cylinder** (Table 1): length 4 m, diameter 24 mm, 0.69 kg/m, EI 10.5 N m², EA 6.67e5 N,
  pretension 500 N; their calculated first natural frequency in water f1 2.68 Hz (with C_a 1)
  defines U_R. The ends are universal joints (pinned here), 40 elements, dt 1 ms (the axial period in
  8 steps), C_d 1.2 and C_a 1.
- **Motion:** the whole cylinder moved across its axis, y = A_m sin(2πt/T), A_m 0.68 m (KC 178), with
  a one-period cosine ramp, 8 cycles after it. The three KC 178 cases are the periods of the captions
  of Figs. 15-17, T 16.5, 10.2 and 8.45 s, U_R,max = 2πA_m/(T f1 D) 4.03, 6.51 and 7.86 (the
  article's 4, 6.5 and 7.9).
- **Measured** (Figs. 15-17, station CF4): A/D against the instantaneous U_R = |V|/(f1 D) on the
  acceleration and deceleration paths, averaged by Fu et al. over 10 middle cycles of 15 (their
  Eq. 20). `extract_fu2014.py` reads the paths from the figures, embedded raster images, calibrated
  on the plot frame and traced by colour, at U_R steps of 0.25 (`fu2014_hysteresis.csv`). The traced
  gaps between the paths, 0.164 at U_R 3.5 (U_R,max 6.5) and 0.222 at U_R 4 (7.9), match the
  figures' own annotations, 0.15D and 0.2D.
- **Model output:** the cross-flow (vertical) displacement at mid-span, band-passed 0.5-25 Hz, its
  Hilbert envelope over D, averaged in U_R bins of ±0.125 over the half-cycles from two periods on,
  on the two paths. CF4 is the fourth of seven cross-flow stations at 0.5 m spacing (Table 2), taken
  as mid-span.
- **Lock-in share:** the share of each half-cycle the envelope spends above √2/2 of its largest value,
  the article's definition (Sec. 4.6), against the shares it gives for KC > 100: about 17, 34 and 50 %.

## Coefficient sets

The sets of the FLOTANT and NDP cases, none fitted to these tests:

- **uls19:** Ulveseter et al. (2019), Table 1 (Cl 0.85, cF 0.144, dF 0.064; VIV ClX 0.75, alpha 0.15,
  Chh 1.5), fitted by them to the NDP riser;
- **b1:** the generic set of the Delta Flume case (Cl 0.8, cF 0.18, dF 0.08; ClX 0.85, alpha 0.15,
  Ulveseter et al. 2018, Table 3; VIV AMPLITUDE with the drag amplification of Chaplin et al. 2005,
  Eq. 4; VIV DAMPING 0 0 1);
- **wake:** the wake oscillator of Facchinetti et al. (2004) (CL0 0.3, St 0.2, eps 0.3, A 12),
  cross-flow only.

The cylinder moves through still water, so every set runs with a `VIV FRAME` row: without it the
wake oscillator never sheds (the ambient flow is zero) and the amplitude estimate and damping split
of b1 act across a flow that is not there. Its filter time cannot meet both of the frame's
conditions here: the motion period is 23-44 periods of f1, around the 30 the conditions need.
The primary runs take `Tf` 0.2 s (ωTf 3.4 at f1; ω_m Tf 0.08, 0.12 and 0.15 at the three periods), and
the record reports the runs without the row and with `Tf` 0.1 and 0.3 s.

## Results (primary runs, `Tf` 0.2 s)

| U_R,max | Set | A/D r.m.s. error | Mean gap dec. − acc. (measured) | Largest A/D (measured) | Lock-in share (measured) |
|---|---|---|---|---|---|
| 4 | uls19 | 0.079 | −0.005 (+0.009) | 0.124 (0.286) | 0.08 (0.17) |
| 4 | b1 | 0.138 | +0.123 (+0.009) | 0.472 (0.286) | 0.31 (0.17) |
| 4 | wake | 0.088 | +0.106 (+0.009) | 0.280 (0.286) | 0.27 (0.17) |
| 6.5 | uls19 | 0.340 | +0.426 (+0.089) | 1.007 (0.475) | 0.34 (0.34) |
| 6.5 | b1 | 0.339 | +0.364 (+0.089) | 0.874 (0.475) | 0.53 (0.34) |
| 6.5 | wake | 0.177 | +0.271 (+0.089) | 0.609 (0.475) | 0.46 (0.34) |
| 7.9 | uls19 | 0.338 | +0.498 (+0.123) | 1.043 (0.450) | 0.34 (0.50) |
| 7.9 | b1 | 0.345 | +0.380 (+0.123) | 0.937 (0.450) | 0.60 (0.50) |
| 7.9 | wake | 0.132 | +0.255 (+0.123) | 0.570 (0.450) | 0.50 (0.50) |

The largest measured A/D is the largest on the U_R grid; the measured deceleration line ends at 0.32 at
U_R 4.12 in Fig. 15.

- **Through the flow reversals** every run converges, the synchronisation phases and the wake
  variable running continuously through U → 0.
- **U_R,max 4:** the wake oscillator gives the measured largest amplitude (0.280 against 0.286); uls19
  barely responds (its band centre, 0.144, is reached only near the velocity peak); b1 is 65 % high.
- **U_R,max 6.5 and 7.9: every set holds the lock-in on the deceleration path.** The deceleration path
  stays near its peak while the flow slows, so the mean gap between the paths is 2-5 times the
  measured one and the largest amplitude 1.3-2.3 times the measured. The measured hysteresis is the
  same in kind (the deceleration path above the acceleration path), smaller in size. The wake
  oscillator comes closest.
- The same persistence shows in the FLOTANT case, whose lock-in releases late for the set that
  reaches its peak.

## Sensitivity and convergence

- **Without `VIV FRAME`:** the wake oscillator does not shed (largest A/D 0.000-0.021). The
  synchronisation sets respond less: b1 0.297, 0.769 and 0.760 (A/D r.m.s. error 0.034 at U_R,max 4,
  the best of any run there; 0.229 and 0.212 at 6.5 and 7.9), uls19 0.080, 0.738 and 0.824. The
  deceleration hold remains (mean gap +0.25 to +0.35 at 6.5 and 7.9).
- **`Tf` 0.1 and 0.3 s:** b1's largest A/D 0.404-0.511, 0.710-0.922 and 0.738-0.978; uls19
  0.110-0.139, 0.945-0.969 and 1.096-1.141; the wake oscillator 0.273-0.284, 0.607-0.621 and
  0.579-0.633. The synchronisation sets depend on the filter time because their amplitude estimate and
  damping split take the filtered frame; the wake oscillator's wake frequency barely does.
- **Convergence (U_R,max 6.5):** half the time step changes the largest A/D by +0.03 % (b1), −0.1 %
  (wake) and −4.1 % (uls19); 80 elements by +0.1 %, +0.07 % and −4.1 %.
- The default command, run into a new work directory, reproduces `response_comparison.csv` and
  `summary.json` byte for byte.

## Limitations

- The measured A/D is reconstructed from strain by a sinusoidal modal fit with constant tension (Fu
  et al., Sec. 3.2), while the tension varies with the motion; the model's A/D is the displacement.
- The submergence of the cylinder, the band-pass cut-offs and the structural damping are not given;
  the model takes 2 m, 0.5-25 Hz and no structural damping.
- The paths are one station (CF4) and three cases; Fu et al. do not show the scatter between cycles.
- The lock-in shares are "about" values from the text, for KC > 100.

## Reproduction

```text
python extract_fu2014.py FU2014.pdf
python fu_oscillatory_viv.py path/to/cabledyn.exe WORKDIR --jobs 8 --summary .
```

The executable of the record is named by its SHA-256 in `summary.json`; `provenance.json` gives the
sources, the data used and the build.
