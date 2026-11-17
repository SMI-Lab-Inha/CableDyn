#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
"""Read the Wang et al. (2014) riser data of Thorsen, Saevik & Larsen (2017) from the article's PDF.

The PDF is not distributed; this script documents how the committed CSV and JSON files were read
from it (PyMuPDF, NumPy and Pillow):

- Fig. 17 (uniform current, PDF page 14) and Fig. 22 (heave-induced VIV, PDF page 16), vector
  graphics: the measured r.m.s. dynamic cross-flow strains (4 x 4 pt stroked squares, the legend
  square removed) and the simulated line of Thorsen et al. (two polylines of 200-300 straight
  segments each, joined at s 9.5 m; their vertices).
  The tick labels are glyph outlines, not text, so the axes are calibrated on the tick marks:
  x 0, 5, ..., 25 m; y 0, 1, ..., 6 x 1e-5 (Fig. 17) and 0, 2, ..., 8 x 1e-5 (Fig. 22).
  -> thorsen2017_fig17.csv, thorsen2017_fig22.csv (series measured or thorsen2017)
- Fig. 19 top panel (PDF page 15), a raster image: the measured x (red) and z (blue)
  displacement of the riser's upper end, t 30-50 s. The panel's frame is found from its dark
  pixels, the time axis spans the frame, the displacement axis is calibrated on the dotted grid
  lines at -0.1, 0 and 0.1 m, and each colour gives one median row per pixel column. A sinusoid
  of the period the text gives (5.96 s) is fitted to each component by least squares.
  -> thorsen2017_fig19_top.csv, thorsen2017_fig19_fit.json

    python extract_thorsen2017.py THORSEN2017.pdf [--check fig19_check.png]

Source: M.J. Thorsen, S. Saevik, C.M. Larsen (2017), Non-linear time domain analysis of
cross-flow vortex-induced vibrations, Marine Structures 51, 134-151,
https://doi.org/10.1016/j.marstruc.2016.10.007
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
PERIOD = 5.96  # s, the top motion period of Sec. 5.2


def strain_figure(page, y_top: float, name: str) -> None:
    """Measured squares and the simulated line of an r.m.s. strain figure in the lower part of
    a page; y_top is the strain at the highest tick mark."""
    h = page.rect.height
    dr = [g for g in page.get_drawings() if g["rect"].y0 > 0.6 * h]
    strokes = [g for g in dr if g["type"] == "s"]
    frame = max(
        (g for g in strokes if len(g["items"]) == 2), key=lambda g: g["rect"].width
    )["rect"]
    xt = sorted(
        {
            round(g["rect"].x0, 2)
            for g in strokes
            if g["rect"].width < 0.1
            and 1.5 < g["rect"].height < 2.5
            and abs(g["rect"].y1 - frame.y1) < 3
        }
    )
    yt = sorted(
        {
            round(g["rect"].y0, 2)
            for g in strokes
            if g["rect"].height < 0.1
            and 1.5 < g["rect"].width < 2.5
            and abs(g["rect"].x0 - frame.x0) < 3
        },
        reverse=True,
    )
    assert len(xt) == 6, xt
    assert len(yt) >= 5, yt

    def fx(x: float) -> float:
        return (x - xt[0]) / (xt[-1] - xt[0]) * 25.0

    def fy(y: float) -> float:
        return (yt[0] - y) / (yt[0] - yt[-1]) * y_top

    squares = [
        g["rect"]
        for g in strokes
        if abs(g["rect"].width - 4.0) < 0.4 and abs(g["rect"].height - 4.0) < 0.4
    ]
    # the legend square, in the upper right corner of the frame
    squares = [
        r
        for r in squares
        if not (
            r.x0 > frame.x0 + 0.7 * frame.width and r.y1 < frame.y0 + 0.3 * frame.height
        )
    ]
    measured = sorted((fx((r.x0 + r.x1) / 2), fy((r.y0 + r.y1) / 2)) for r in squares)
    # the simulated line is drawn as two polylines of straight segments, joined at x 254.7 pt
    pts = set()
    for line in strokes:
        if len(line["items"]) >= 20 and all(it[0] == "l" for it in line["items"]):
            for it in line["items"]:
                pts.update(((it[1].x, it[1].y), (it[2].x, it[2].y)))
    sim = sorted((fx(x), fy(y)) for x, y in pts)
    with open(HERE / name, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(["series", "s_m", "rms_strain"])
        for series, data in (("measured", measured), ("thorsen2017", sim)):
            w.writerows((series, round(a, 4), f"{b:.4e}") for a, b in data)
    peak = max(measured, key=lambda p: p[1])
    print(
        f"{name}: {len(measured)} measured points, max {peak[1]:.3e} at s {peak[0]:.2f} m; "
        f"simulated line {len(sim)} points, max {max(v for _, v in sim):.3e} "
        f"({100 * max(v for _, v in sim) / peak[1]:.1f} % of the measured max)"
    )


def top_motion(doc, check: Path | None) -> None:
    page = doc[14]
    big = max(page.get_image_info(xrefs=True), key=lambda i: i["width"] * i["height"])
    img = Image.open(io.BytesIO(doc.extract_image(big["xref"])["image"])).convert("RGB")
    a = np.asarray(img).astype(int)
    r, g, b = a[..., 0], a[..., 1], a[..., 2]
    red = (r > 180) & (g < 110) & (b < 110)
    blue = (b > 180) & (r < 110) & (g < 140)
    dark = (r < 90) & (g < 90) & (b < 90)
    # the top panel lies above the strain contour plot: the first row where coloured pixels
    # fill more than 30 % of the width
    colour = (np.abs(r - g) > 60) | (np.abs(g - b) > 60) | (np.abs(r - b) > 60)
    top_end = int(np.flatnonzero(colour.mean(axis=1) > 0.3)[0]) - 8
    red[top_end:] = False
    blue[top_end:] = False
    rows = np.flatnonzero(red.any(axis=1))
    r0, r1 = rows.min(), rows.max()
    # the panel frame: dark columns with long vertical runs over the red rows
    band = dark[max(0, r0 - 40) : r1 + 40]
    cand = np.flatnonzero(band.sum(axis=0) > 0.6 * band.shape[0])
    x0, x1 = cand.min(), cand.max()
    hor = np.flatnonzero(dark[:, x0:x1].sum(axis=1) > 0.6 * (x1 - x0))
    hor = hor[(hor > r0 - 60) & (hor < r1 + 60)]
    y_top, y_bot = hor.min(), hor.max()
    # dotted grid lines: rows inside the frame with many dark pixels
    darkish = (r < 140) & (g < 140) & (b < 140) & ~red & ~blue
    inner = darkish[y_top + 4 : y_bot - 3, x0 + 4 : x1 - 3].sum(axis=1)
    grid = np.flatnonzero(inner > 0.04 * (x1 - x0)) + y_top + 4
    groups = np.split(grid, np.flatnonzero(np.diff(grid) > 2) + 1)
    gy = sorted(float(gr.mean()) for gr in groups if gr.size)
    assert len(gy) == 3, gy  # 0.1, 0 and -0.1 m from the top

    def disp(row: float) -> float:
        return 0.1 - (row - gy[0]) / (gy[2] - gy[0]) * 0.2

    def time(col: int) -> float:
        return 30.0 + (col - x0) / (x1 - x0) * 20.0

    w = 2 * np.pi / PERIOD
    fit, traces = {}, {}
    for name, mask in (("x", red), ("z", blue)):
        cols, vals = [], []
        for c in range(x0 + 2, x1 - 1):
            rr = np.flatnonzero(mask[y_top:y_bot, c])
            if rr.size:
                cols.append(time(c))
                vals.append(disp(y_top + float(np.median(rr))))
        t, y = np.array(cols), np.array(vals)
        basis = np.column_stack([np.sin(w * t), np.cos(w * t), np.ones_like(t)])
        coef, *_ = np.linalg.lstsq(basis, y, rcond=None)
        resid = y - basis @ coef
        fit[name] = {
            "amp_m": float(np.hypot(coef[0], coef[1])),
            "phase_rad": float(np.arctan2(coef[1], coef[0])),
            "mean_m": float(coef[2]),
            "rms_resid_m": float(resid.std()),
            "n": int(t.size),
        }
        traces[name] = (t, y)
        print(
            f"{name}: amplitude {fit[name]['amp_m']:.4f} m, phase "
            f"{np.degrees(fit[name]['phase_rad']):.2f} deg (A sin(w t + phase) + mean), "
            f"residual r.m.s. {resid.std():.4f} m over {t.size} columns"
        )
    (HERE / "thorsen2017_fig19_fit.json").write_text(
        json.dumps({**fit, "period_s": PERIOD}, indent=1) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    with open(
        HERE / "thorsen2017_fig19_top.csv", "w", newline="", encoding="utf-8"
    ) as fh:
        fh.write("component,t_s,disp_m\n")
        for name, (t, y) in traces.items():
            fh.writelines(
                f"{name},{tt:.4f},{yy:.5f}\n" for tt, yy in zip(t, y, strict=True)
            )
    if check is not None:  # the fitted sinusoids drawn over the image
        chk = img.copy()
        px = chk.load()
        for v in fit.values():
            for c in range(x0, x1):
                yy = v["amp_m"] * np.sin(w * time(c) + v["phase_rad"]) + v["mean_m"]
                row = round(gy[0] + (0.1 - yy) / 0.2 * (gy[2] - gy[0]))
                if 0 <= row < chk.height:
                    px[c, row] = (0, 0, 0)
        box = (max(0, x0 - 60), max(0, y_top - 20), min(chk.width, x1 + 20), y_bot + 40)
        chk.crop(box).save(check)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("pdf")
    ap.add_argument("--check", type=Path, default=None, help="PNG of the Fig. 19 fit")
    a = ap.parse_args()
    doc = pymupdf.open(a.pdf)
    strain_figure(doc[13], 6.0e-5, "thorsen2017_fig17.csv")
    strain_figure(doc[15], 8.0e-5, "thorsen2017_fig22.csv")
    top_motion(doc, a.check)


if __name__ == "__main__":
    main()
