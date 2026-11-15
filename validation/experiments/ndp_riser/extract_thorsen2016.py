#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
"""Read Figs. 7 and 8 of Thorsen et al. (2016) from the vector graphics of the article's page 8.

Fig. 7 is their excitation curve C_v(y0/D), fitted to the NDP riser tests; Fig. 8 the maximum
r.m.s. cross-flow strain of the 22 uniform-flow (a) and 22 sheared-flow (b) tests, measured at
the fundamental frequency (blue squares) and in total (red circles), with their simulation using
the fitted curve (blue line). The figures are vector paths, so the values are the plotted
coordinates mapped through the axes, not a pixel digitisation. The tick labels are drawn glyphs,
not text, so the axes come from the plot frames:

- Fig. 7: y0/D = 0 and 1 at x = 191.7 and 381.1 pt (the frame, where the curve starts and
  ends), C_v = 0 and 1.5 at y = 161.25 and 59.7 pt (the bottom and top ticks);
- Fig. 8: x from 0 to 2.5 m/s across the frame, y from 0 to 5e-4 (a) and 0 to 3.5e-4 (b) up it,
  the end tick labels of each axis on the rendered page.

The legend symbols, in the upper left of each panel, are dropped.

    python extract_thorsen2016.py path/to/article.pdf [--out DIR]

writes DIR/thorsen2016_fig7_cv.csv (the curve sampled every 0.05 in y0/D) and
DIR/thorsen2016_fig8_strain.csv (panel, series, U, strain in units of 1e-4). Needs PyMuPDF.

Source: M.J. Thorsen, S. Saevik, C.M. Larsen (2016), Time domain simulation of vortex-induced
vibrations in stationary and oscillating flows, Journal of Fluids and Structures 61, 1-19,
https://doi.org/10.1016/j.jfluidstructs.2015.11.006
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PAGE = 7  # zero-based: the article's page 8 holds Figs. 7 and 8


def fig7(draws) -> np.ndarray:
    """The C_v curve: the polyline with the most vertices above y = 175 pt."""

    def fx(px: float) -> float:
        return (px - 191.7) / (381.1 - 191.7)

    def fy(py: float) -> float:
        return 1.5 * (161.25 - py) / (161.25 - 59.7)

    best: list[tuple[float, float]] = []
    for d in draws:
        if d["rect"].y1 > 175 or d.get("color") is None or d.get("width") is None:
            continue
        pts: list[tuple[float, float]] = []
        for it in d["items"]:
            if it[0] == "l":
                if not pts:
                    pts.append((it[1].x, it[1].y))
                pts.append((it[2].x, it[2].y))
        if len(pts) > 20 and len(pts) > len(best):
            best = pts
    arr = np.array([(fx(x), fy(y)) for x, y in best])
    return arr[np.argsort(arr[:, 0])]


def fig8(draws) -> list[tuple[str, str, float, float]]:
    # the white plot frames from top to bottom: Fig. 7, Fig. 8(a), Fig. 8(b)
    boxes = sorted(
        {
            (
                round(d["rect"].x0, 1),
                round(d["rect"].y0, 1),
                round(d["rect"].x1, 1),
                round(d["rect"].y1, 1),
            )
            for d in draws
            if d.get("fill") == (1.0, 1.0, 1.0) and d["rect"].width > 100
        },
        key=lambda b: b[1],
    )
    panels = {"uniform": (boxes[1], 5.0), "sheared": (boxes[2], 3.5)}
    rows = []
    for name, (box, ymax) in panels.items():
        x0, y0, x1, y1 = box
        ax, bx = 2.5 / (x1 - x0), -2.5 * x0 / (x1 - x0)
        ay, by = -ymax / (y1 - y0), ymax * y1 / (y1 - y0)
        print(f"Fig. 8 {name}: frame {box}, y 0..{ymax}e-4")
        series = defaultdict(list)
        for d in draws:
            r = d["rect"]
            col = d.get("color")
            if col is None or not (y0 - 5 < r.y0 < y1 + 5 and x0 - 5 < r.x0 < x1 + 5):
                continue
            cx, cy = 0.5 * (r.x0 + r.x1), 0.5 * (r.y0 + r.y1)
            kinds = {it[0] for it in d["items"]}
            if 2 < r.width < 8 and 2 < r.height < 8:
                if col[2] > 0.8 and col[0] < 0.3 and "re" in kinds:
                    series["exp_fundamental"].append((ax * cx + bx, ay * cy + by))
                elif col[0] > 0.8 and col[2] < 0.3 and "c" in kinds:
                    series["exp_total"].append((ax * cx + bx, ay * cy + by))
            elif col[2] > 0.8 and col[0] < 0.3 and "l" in kinds and r.width > 20:
                for it in d["items"]:
                    if it[0] == "l":
                        for p in (it[1], it[2]):
                            series["sim"].append((ax * p.x + bx, ay * p.y + by))
        for s, pts in series.items():
            # drop the legend symbols (upper left, above every data point there)
            pts = [(u, v) for u, v in pts if not (u < 1.35 and v > 0.5 * ymax)]
            pts = sorted({(round(u, 3), round(v, 4)) for u, v in pts})
            if s == "sim":  # polyline vertices, one per test speed
                pts = sorted({round(u, 1): v for u, v in pts}.items())
            rows += [(name, s, u, v) for u, v in pts]
            print(f"  {s}: {len(pts)} points")
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("pdf")
    ap.add_argument("--out", default=str(HERE))
    a = ap.parse_args()
    import pymupdf  # only this script needs it

    draws = pymupdf.open(a.pdf)[PAGE].get_drawings()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    cv = fig7(draws)
    i = int(np.argmax(cv[:, 1]))
    zero = cv[cv[:, 1] <= 1e-3][:, 0]
    print(
        f"Fig. 7: {len(cv)} vertices; Cv(0) {cv[0, 1]:.4f}, maximum {cv[i, 1]:.4f} at y0/D "
        f"{cv[i, 0]:.4f}, zero from y0/D {zero.min():.4f}"
    )
    grid = np.arange(0.0, 1.0001, 0.05)
    table = np.clip(np.interp(grid, cv[:, 0], cv[:, 1], right=0.0), 0.0, None)
    # the curve runs on along the axis, a line width above it, from its zero crossing: the
    # fitted curve is zero there (their Sec. 3: C_v is positive or zero)
    table[grid >= zero.min()] = 0.0
    with open(out / "thorsen2016_fig7_cv.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["y0_over_D", "Cv"])
        w.writerows((f"{g:.2f}", f"{c:.4f}") for g, c in zip(grid, table, strict=True))

    rows = fig8(draws)
    with open(
        out / "thorsen2016_fig8_strain.csv", "w", newline="", encoding="utf-8"
    ) as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["panel", "series", "U", "strain_1e-4"])
        w.writerows(rows)


if __name__ == "__main__":
    main()
