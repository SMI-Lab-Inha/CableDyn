#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
"""Reading aids for the bare-riser, uniform-flow panels of Trim et al. (2005), Figs. 16 and 18.

Both figures are raster images embedded in the article (962 px wide). The panel 'bare riser,
uniform flow' is the top-left one of each:

- Fig. 16 (article page 21): the spatial mean of the temporal standard deviation of the
  displacement over the diameter, sigma/D, cross-flow (blue diamonds) and in-line (red stars);
- Fig. 18 (page 23): the peak frequency of the dominant displacement mode [Hz].

The markers (about 10 px) are wider than the speed step (0.1 m/s, about 5.7 px) and overlap, so
an automatic marker search does not separate them. This script calibrates the axes on the
panel frame and its dotted grid lines (x: the 0.5 m/s grid lines and the frame's right edge at
2.5 m/s; y: the frame bottom at zero and the horizontal grid lines) and writes a four-times
enlarged crop of each panel with a guide line at every tested speed and every 0.01 (Fig. 16) or
1 Hz (Fig. 18). measured_trim2005_fig16.csv and measured_trim2005_fig18.csv were read by eye on
these crops: the centre of the marker on each speed's guide.

    python digitise_trim2005.py path/to/article.pdf [--out DIR]

writes DIR/trim2005_fig16_guides.png and DIR/trim2005_fig18_guides.png and prints the axis
calibration. Needs PyMuPDF and Pillow.

Source: A.D. Trim, H. Braaten, H. Lie, M.A. Tognarelli (2005), Experimental investigation of
vortex-induced vibration of long marine risers, Journal of Fluids and Structures 21(3), 335-361,
https://doi.org/10.1016/j.jfluidstructs.2005.07.014
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
# figure: (zero-based page, search region x0 y0 x1 y1 [px] around the panel frame, y tick values
# from the frame bottom up, guide step in y)
FIGS = {
    "fig16": (20, (100, 60, 262, 345), [0.0, 0.1, 0.2, 0.3, 0.4, 0.5], 0.01),
    "fig18": (22, (60, 60, 250, 345), [0.0, 10.0, 20.0, 30.0], 1.0),
}
SPEEDS = np.round(np.arange(0.3, 2.45, 0.1), 2)  # the 22 tested speeds [m/s]


def frame(rgb: np.ndarray, reg) -> tuple[int, int, int, int]:
    """The panel frame: the dark columns and rows that span most of the region."""
    x0, y0, x1, y1 = reg
    dark = rgb[y0:y1, x0:x1].sum(axis=2) < 250
    cols = np.where(dark.sum(axis=0) > 0.45 * dark.shape[0])[0]
    rows = np.where(dark.sum(axis=1) > 0.3 * dark.shape[1])[0]
    return (
        int(x0 + cols.min()),
        int(y0 + rows.min()),
        int(x0 + cols.max()),
        int(y0 + rows.max()),
    )


def grid(rgb: np.ndarray, box, axis: int) -> list[float]:
    """The dotted grey grid lines inside the frame; axis 0: vertical lines, 1: horizontal."""
    x0, y0, x1, y1 = box
    sub = rgb[y0 + 3 : y1 - 2, x0 + 3 : x1 - 2].astype(int)
    grey = (sub.sum(axis=2) < 450) & (np.ptp(sub, axis=2) < 60)
    prof = grey.sum(axis=axis)
    groups: list[list[int]] = []
    for i in np.where(prof > 0.12 * grey.shape[axis])[0]:
        if groups and i - groups[-1][-1] <= 2:
            groups[-1].append(int(i))
        else:
            groups.append([int(i)])
    off = (x0 if axis == 0 else y0) + 3
    return [float(np.mean(g)) + off for g in groups]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("pdf")
    ap.add_argument("--out", default=str(HERE))
    a = ap.parse_args()
    import pymupdf
    from PIL import Image, ImageDraw

    doc = pymupdf.open(a.pdf)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    for fig, (page, reg, yt, step) in FIGS.items():
        pix = pymupdf.Pixmap(doc, doc[page].get_images(full=True)[0][0])
        if pix.n - pix.alpha > 3:
            pix = pymupdf.Pixmap(pymupdf.csRGB, pix)
        img = Image.frombytes("RGB", (pix.width, pix.height), pix.samples)
        rgb = np.asarray(img).astype(int)
        box = frame(rgb, reg)
        gx = grid(rgb, box, 0)
        gy = grid(rgb, box, 1)
        # the vertical grid lines are at multiples of 0.5 m/s; the right frame edge is 2.5 m/s
        sp = np.median(np.diff(gx))
        gx = [g for g in gx if g < box[2] - 3] + [float(box[2])]
        ticks = [2.5 - 0.5 * round((box[2] - g) / sp) for g in gx]
        ax, bx = np.polyfit(ticks, gx, 1)  # pixel x = ax U + bx
        ypix = [box[3]] + sorted((y for y in gy if y < box[3] - 5), reverse=True)
        n = min(len(ypix), len(yt))
        ay, by = np.polyfit(yt[:n], ypix[:n], 1)  # pixel y = ay value + by
        print(
            f"{fig}: frame {box}; x ticks {ticks} at {[round(g, 1) for g in gx]} px "
            f"({ax:.2f} px per m/s); y ticks {yt[:n]} at {[round(y, 1) for y in ypix[:n]]} px "
            f"({-ay:.2f} px per unit)"
        )
        x0, y0, x1, y1 = box[0] - 7, box[1] - 10, box[2] + 7, box[3] + 8
        k = 4
        crop = img.crop((x0, y0, x1, y1)).resize(
            ((x1 - x0) * k, (y1 - y0) * k), Image.NEAREST
        )
        dr = ImageDraw.Draw(crop)
        for u in SPEEDS:
            px = (ax * u + bx - x0) * k
            dr.line([(px, 0), (px, crop.height)], fill=(0, 190, 0), width=1)
            dr.text((px + 1, 2), f"{u:.1f}", fill=(0, 0, 0))
        for j in range(int(round(yt[-1] / step)) + 1):
            py = (ay * j * step + by - y0) * k
            major = j % 10 == 0 if step < 0.1 else j % 5 == 0
            dr.line(
                [(0, py), (crop.width, py)],
                fill=(200, 100, 0) if major else (255, 190, 0),
                width=1,
            )
        crop.save(out / f"trim2005_{fig}_guides.png")


if __name__ == "__main__":
    main()
