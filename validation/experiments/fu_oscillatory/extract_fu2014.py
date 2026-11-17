#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
"""Read the measured hysteresis paths of Fu et al. (2014), Figs. 15-17, from the article's PDF.

The three figures are embedded raster images (1503 px wide): the cross-flow A/D at CF4 against the
instantaneous reduced velocity U_R = |V|/(f1 D) on the acceleration path (black solid line) and the
deceleration path (red dashed line), for KC 178 and U_R,max 4, 6.5 and 7.9. Each image is calibrated
on its plot frame, the outermost dark lines that span the plot, with the axis limits printed at the
frame's corners (U_R 0 to 5, 7 and 8; A/D 0 to 0.4, 0 to 0.6 and 0.1 to 0.6). The inset, the legend
and the title box are masked, and so are the tick marks inside the frame; the line is followed column by column from its
left end (the run of its colour nearest the previous column's centre), its centre the run's median row,
and columns where an arrowhead widens the line are bridged linearly. The
paths are written at U_R steps of 0.25 between the ends of each line.

    python extract_fu2014.py FU2014.pdf [--out fu2014_hysteresis.csv] [--check DIR]

Source:
  S. Fu, J. Wang, R. Baarholm, J. Wu, C.M. Larsen (2014), Features of vortex-induced vibration in
  oscillatory flow, Journal of Offshore Mechanics and Arctic Engineering 136(1), 011801,
  https://doi.org/10.1115/1.4025759 (Figs. 15-17, p. 011801-7 and 011801-8)
"""

from __future__ import annotations

import argparse
import csv
import io
import json
from pathlib import Path

import numpy as np
import pymupdf
from PIL import Image

HERE = Path(__file__).resolve().parent

# (page index, image index on the page, U_R,max, T [s], x range, y range, masked boxes in pixels
# (x0, y0, x1, y1): the inset, the legend and the title box)
FIGURES = (
    (
        "15",
        6,
        1,
        4.0,
        16.5,
        (0.0, 5.0),
        (0.0, 0.4),
        ((160, 115, 780, 340), (1160, 10, 1500, 110), (650, 15, 955, 110)),
    ),
    (
        "16",
        6,
        2,
        6.5,
        10.2,
        (0.0, 7.0),
        (0.0, 0.6),
        ((180, 80, 750, 310), (1160, 10, 1500, 110), (630, 15, 960, 105)),
    ),
    (
        "17",
        7,
        0,
        7.9,
        8.45,
        (0.0, 8.0),
        (0.1, 0.6),
        ((165, 80, 760, 310), (1155, 10, 1500, 110), (635, 15, 960, 105)),
    ),
)


def frame(dark: np.ndarray) -> tuple[float, float, float, float]:
    """The plot frame: centre rows of the top and bottom lines, centre columns of the side lines."""
    h, w = dark.shape
    rows = np.where(dark.sum(axis=1) > 0.8 * w)[0]
    cols = np.where(dark.sum(axis=0) > 0.8 * h)[0]

    def centre(idx: np.ndarray, first: bool) -> float:
        groups = np.split(idx, np.where(np.diff(idx) > 1)[0] + 1)
        g = groups[0] if first else groups[-1]
        return float(g.mean())

    return (
        centre(cols, True),
        centre(cols, False),
        centre(rows, True),
        centre(rows, False),
    )


def trace(mask: np.ndarray, x0: float, x1: float) -> tuple[np.ndarray, np.ndarray]:
    """The line's centre row in each column inside the frame, following the line from its left end.

    In each column the coloured pixels split into contiguous runs; the run nearest the previous
    column's centre is the line (a stray mark elsewhere in the column is ignored), and its median row
    is the centre. Columns where an arrowhead widens the run are dropped and bridged linearly.
    """
    cols, rows, widths = [], [], []
    prev = None
    for c in range(int(np.ceil(x0)) + 3, int(np.floor(x1)) - 2):
        r = np.where(mask[:, c])[0]
        if not len(r):
            continue
        runs = np.split(r, np.where(np.diff(r) > 2)[0] + 1)
        if prev is None:
            run = max(runs, key=len)
        else:
            run = min(runs, key=lambda q: abs(float(np.median(q)) - prev))
            if (
                abs(float(np.median(run)) - prev) > 40
            ):  # nothing of this line in the column
                continue
        prev = float(np.median(run))
        cols.append(c)
        rows.append(prev)
        widths.append(int(run.max() - run.min() + 1))
    cols_a, rows_a, w = np.array(cols), np.array(rows), np.array(widths)
    keep = w <= 2.5 * np.median(w) + 2
    return cols_a[keep], rows_a[keep]


def extract(pdf: Path, check: Path | None) -> list[dict]:
    doc = pymupdf.open(pdf)
    out = []
    for fig, page, index, urmax, period, (ux0, ux1), (ay0, ay1), boxes in FIGURES:
        xref = doc[page].get_images(full=True)[index][0]
        img = np.asarray(
            Image.open(io.BytesIO(doc.extract_image(xref)["image"])).convert("RGB")
        ).astype(int)
        r, g, b = img[..., 0], img[..., 1], img[..., 2]
        dark = img.max(axis=2) < 100
        left, right, top, bottom = frame(dark)
        # the solid line is anti-aliased: any neutral grey darker than the blue-grey grid
        black = (img.max(axis=2) < 150) & (np.abs(r - b) < 30) & (np.abs(r - g) < 30)
        red = (r > 170) & (g < 110) & (b < 110)
        for x0, y0, x1, y1 in boxes:
            black[y0:y1, x0:x1] = False
            red[y0:y1, x0:x1] = False
        # exclude the frame lines and the tick marks drawn inward from them (at most 11 px on the
        # horizontal axis, 20 px on the vertical axes)
        for m in (black, red):
            m[:, : int(left) + 21] = False
            m[:, int(right) - 20 :] = False
            m[: int(top) + 12, :] = False
            m[int(bottom) - 11 :, :] = False

        def ur(c):
            return ux0 + (np.asarray(c, float) - left) / (right - left) * (ux1 - ux0)

        def ad(rw):
            return ay1 - (np.asarray(rw, float) - top) / (bottom - top) * (ay1 - ay0)

        for path, m in (("acceleration", black), ("deceleration", red)):
            c, rw = trace(m, left, right)
            u, a = ur(c), ad(rw)
            grid = np.arange(np.ceil(u.min() / 0.25) * 0.25, u.max() + 1e-9, 0.25)
            vals = np.interp(grid, u, a)
            for gu, gv in zip(grid, vals):
                out.append(
                    {
                        "figure": fig,
                        "U_R_max": urmax,
                        "T_s": period,
                        "KC": 178,
                        "path": path,
                        "U_R": round(float(gu), 4),
                        "A_D": round(float(gv), 4),
                    }
                )
            out.append(
                {
                    "figure": fig,
                    "U_R_max": urmax,
                    "T_s": period,
                    "KC": 178,
                    "path": path,
                    "U_R": round(float(u.max()), 4),
                    "A_D": round(float(a[np.argmax(u)]), 4),
                    "end": 1,
                }
            )
        if check is not None:
            check.mkdir(parents=True, exist_ok=True)
            vis = img.astype(np.uint8).copy()
            for path, m, colour in (
                ("acceleration", black, (0, 200, 0)),
                ("deceleration", red, (0, 120, 255)),
            ):
                c, rw = trace(m, left, right)
                for cc, rr in zip(c, rw):
                    vis[int(round(rr)) - 1 : int(round(rr)) + 2, cc] = colour
            Image.fromarray(vis).save(check / f"fig{fig}_check.png")
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("pdf", type=Path)
    ap.add_argument("--out", type=Path, default=HERE / "fu2014_hysteresis.csv")
    ap.add_argument(
        "--check", type=Path, help="write the traced lines over each figure into DIR"
    )
    a = ap.parse_args()
    rows = extract(a.pdf, a.check)
    with open(a.out, "w", encoding="utf-8", newline="\n") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=["figure", "U_R_max", "T_s", "KC", "path", "U_R", "A_D"],
            lineterminator="\n",
        )
        w.writeheader()
        for row in rows:
            if not row.get("end"):
                w.writerow({k: row[k] for k in w.fieldnames})
    ends = [
        {k: row[k] for k in ("figure", "path", "U_R", "A_D")}
        for row in rows
        if row.get("end")
    ]
    print(json.dumps(ends, indent=1))


if __name__ == "__main__":
    main()
