#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
"""Read the measured points of Chaplin & King (2018), Figs. 9, 10, 12(a), 13 and 15(a), from the PDF.

The figures are embedded greyscale raster images. Each panel is calibrated on its plot frame and
its tick marks: the frame's left and bottom lines are the zero of each axis (printed there), the
marks found on both opposite sides of the frame give the tick spacing (a mark on one side only,
or a reference line drawn across the frame, is not a tick), and the axis is the least-squares
line through the origin and the ticks; its largest residual is reported.

The markers are circles of one radius per kind (open, m* 1.12; solid, m* 1.84). Each is a peak of
the share of dark pixels among the dark and white ones on a circle about the point (a circle
Hough transform: grey pixels, an error bar drawn over the outline, count neither way), classed
by the share of dark pixels inside it (less its centre column, where a frame line can pass) as
solid or open; one between the classes is reported and dropped, and so is a peak off the
0.05 m/s grid of the tested speeds by more than a quarter step (two overlapping markers make
such a peak). A marker hidden under another is lost. The speeds are snapped to the grid.

The error bars of Figs. 10 and 13 (one standard deviation of the top tension each way) end in
caps, grey rows at least a marker radius wide through the marker's column: the bar is the
nearest pair of caps about equally far above and below the marker, its half-length the
standard deviation and its centre the mean (the bar is drawn about it). A bar hidden by its
marker (shorter than the radius), or among the overlapping bars of the two series, is not read
(`sd` blank). Two lines drawn in the figures check the reading: f d/U = 0.11 in Fig. 9 (read
where no marker or tick is near it) and the total submerged weight times the water depth, W h,
in Fig. 13 (Table 1: 41.8 and 5.8 N).

    python extract_chaplin2018.py CHAPLIN2018.pdf [--out measured_chaplin2018.csv] [--check DIR]

Source:
  J.R. Chaplin, R. King (2018), Laboratory measurements of the vortex-induced vibrations of an
  untensioned catenary riser with high curvature, Journal of Fluids and Structures 79, 26-38,
  https://doi.org/10.1016/j.jfluidstructs.2018.01.008 (Figs. 9, 10, 12, 13 and 15)
"""

from __future__ import annotations

import argparse
import csv
import io
import json
from pathlib import Path

import numpy as np
import pymupdf
from PIL import Image, ImageDraw
from scipy import ndimage

HERE = Path(__file__).resolve().parent
U_GRID = 0.05  # m/s, the speed step of the tests
DIAMETER = 0.056  # m (Table 1), for the slope of the line f d/U = 0.11 in Fig. 9
DEPTH = 2.0  # m (Table 1), for the lines W h in Fig. 13

# Panel: figure label, page index (from 0), image index on the page, panel (0 the left frame, 1
# the right), quantity, current direction, the point on the riser, x and y tick steps, the grey
# level below which a pixel is dark, marker radii (open, solid) in pixels, masked boxes
# (x0, y0, x1, y1) in pixels (insets, legends, panel labels, notes), and the grey levels of the
# error bars or None.
PANELS = (
    {
        "figure": "9a",
        "page": 8,
        "image": 1,
        "panel": 0,
        "quantity": "f_cf_hz",
        "current": "positive",
        "where": "upper (two points of the inset)",
        "x_step": 0.25,
        "y_step": 0.5,
        "black": 128,
        "radius": (7.5, 7.5),
        "masks": ((150, 18, 340, 125), (265, 312, 470, 378), (455, 325, 516, 380)),
        "bars": None,
    },
    {
        "figure": "9b",
        "page": 8,
        "image": 1,
        "panel": 1,
        "quantity": "f_cf_hz",
        "current": "positive",
        "where": "lower (the point of the inset)",
        "x_step": 0.25,
        "y_step": 0.5,
        "black": 128,
        "radius": (7.5, 7.5),
        "masks": ((585, 18, 790, 125), (655, 318, 880, 378), (890, 325, 951, 380)),
        "bars": None,
    },
    {
        "figure": "10",
        "page": 8,
        "image": 2,
        "panel": 0,
        "quantity": "mean_top_tension_N",
        "current": "positive",
        "where": "top",
        "x_step": 0.25,
        "y_step": 50.0,
        "black": 70,
        "radius": (9.5, 9.5),
        "masks": (),
        "bars": (75, 140),
    },
    {
        "figure": "12a",
        "page": 9,
        "image": 3,
        "panel": 0,
        "quantity": "psi_deg",
        "current": "negative",
        "where": "straight part",
        "x_step": 0.25,
        "y_step": 15.0,
        "black": 128,
        "radius": (7.0, 7.0),
        "masks": ((330, 8, 556, 120), (142, 270, 200, 328)),
        "bars": None,
    },
    {
        "figure": "13",
        "page": 10,
        "image": 0,
        "panel": 0,
        "quantity": "mean_top_tension_N",
        "current": "negative",
        "where": "top",
        "x_step": 0.25,
        "y_step": 25.0,
        "black": 110,
        "radius": (10.0, 10.0),
        "masks": ((390, 40, 720, 160),),
        "bars": (140, 225),
    },
    {
        "figure": "15a",
        "page": 11,
        "image": 0,
        "panel": 0,
        "quantity": "f_cf_hz",
        "current": "negative",
        "where": "around mid-depth",
        "x_step": 0.4,
        "y_step": 0.5,
        "black": 128,
        "radius": (7.5, 4.5),
        "masks": ((165, 17, 225, 70),),
        "bars": None,
    },
)


def image(doc, page: int, index: int) -> np.ndarray:
    xref = doc[page].get_images(full=True)[index][0]
    return np.asarray(
        Image.open(io.BytesIO(doc.extract_image(xref)["image"])).convert("L")
    ).astype(int)


def groups(idx: np.ndarray) -> list[float]:
    if not len(idx):
        return []
    return [float(g.mean()) for g in np.split(idx, np.where(np.diff(idx) > 1)[0] + 1)]


def frames(dark: np.ndarray) -> list[tuple[float, float, float, float]]:
    """The plot frames of the image, left to right: (left, right, top, bottom) centre lines."""
    h, w = dark.shape
    rows = groups(np.where(dark.sum(axis=1) > 0.35 * w)[0])
    top, bottom = rows[0], rows[-1]
    t, b = int(top) + 3, int(bottom) - 2
    cols = groups(np.where(dark[t:b].sum(axis=0) > 0.9 * (b - t))[0])
    return [(cols[i], cols[i + 1], top, bottom) for i in range(0, len(cols) - 1, 2)]


def axis(
    dark: np.ndarray, frame: tuple, horizontal: bool, step: float
) -> tuple[float, float, float, int]:
    """Origin pixel and pixels per unit of one axis from the frame and its tick marks, the
    largest residual of the ticks (in axis units) and their number."""
    left, right, top, bottom = (int(round(v)) for v in frame)
    if horizontal:
        near = groups(
            np.where(dark[bottom - 13 : bottom - 3, left + 4 : right - 3].sum(0) >= 7)[
                0
            ]
            + left
            + 4
        )
        far = groups(
            np.where(dark[top + 3 : top + 13, left + 4 : right - 3].sum(0) >= 7)[0]
            + left
            + 4
        )
        origin = frame[0]
    else:
        near = groups(
            np.where(dark[top + 4 : bottom - 3, left + 3 : left + 13].sum(1) >= 7)[0]
            + top
            + 4
        )
        far = groups(
            np.where(dark[top + 4 : bottom - 3, right - 13 : right - 3].sum(1) >= 7)[0]
            + top
            + 4
        )
        origin = frame[3]
    ticks = [p for p in near if any(abs(p - q) <= 1.5 for q in far)]
    # the spacing that puts marks on the most whole steps from the origin, less the steps
    # without a mark (a reference line drawn across the frame is not on a step)
    cands = [abs(p - origin) / k for p in ticks for k in range(1, 11)]

    def score(s: float) -> int:
        steps = [abs(p - origin) / s for p in ticks]
        on = {round(k) for k in steps if abs(k - round(k)) < 0.1 and round(k) > 0}
        return 2 * len(on) - max(on, default=0)

    spacing = max(cands, key=lambda s: (score(s), s))
    pix, val = [origin], [0.0]
    for p in ticks:
        k = (p - origin) / spacing * (1 if horizontal else -1)
        if abs(k - round(k)) < 0.15 and round(k) > 0:
            pix.append(p)
            val.append(round(k) * step)
    slope, icept = np.polyfit(val, pix, 1)
    resid = max(abs((p - icept) / slope - v) for p, v in zip(pix, val, strict=True))
    return icept, slope, float(resid), len(pix) - 1


def ring_votes(black: np.ndarray, white: np.ndarray, radius: float) -> np.ndarray:
    """The share of dark pixels among the dark and white ones on the circle of the radius about
    each pixel (a grey pixel, an error bar drawn over the outline, is no evidence either way);
    zero where under half the circle is dark or white."""
    h, w = black.shape
    pad = int(radius) + 3
    pb, pw = np.pad(black.astype(float), pad), np.pad(white.astype(float), pad)
    ang = np.linspace(0.0, 2.0 * np.pi, 96, endpoint=False)
    offs = {(round(radius * np.sin(a)), round(radius * np.cos(a))) for a in ang}
    nb, nw = np.zeros((h, w)), np.zeros((h, w))
    for dy, dx in offs:
        nb += pb[pad + dy : pad + dy + h, pad + dx : pad + dx + w]
        nw += pw[pad + dy : pad + dy + h, pad + dx : pad + dx + w]
    seen = nb + nw
    return np.where(seen >= 0.5 * len(offs), nb / np.maximum(seen, 1.0), 0.0)


def disc_fill(black: np.ndarray, white: np.ndarray, radius: float) -> np.ndarray:
    """The share of dark pixels among the dark and white ones in the disc of the radius, less
    its centre column on a large disc (where a frame line passes through a marker)."""
    r = int(np.ceil(radius))
    yy, xx = np.mgrid[-r : r + 1, -r : r + 1]
    k = (
        (xx * xx + yy * yy <= radius * radius)
        & (np.abs(xx) > (2 if radius > 4.0 else -1))
    ).astype(float)
    nb = ndimage.convolve(black.astype(float), k, mode="constant")
    nw = ndimage.convolve(white.astype(float), k, mode="constant")
    return nb / np.maximum(nb + nw, 1.0)


def markers(
    black: np.ndarray, white: np.ndarray, r_open: float, r_solid: float, vote: float
) -> list[tuple]:
    """(x, y, kind, ring share, fill share) of the circles found: a peak of the share of dark
    pixels on the circle (the largest of four radii, the outline being anti-aliased), classed
    by the share of dark pixels inside it."""
    found = []
    for kind, r in (("open", r_open), ("solid", r_solid)):
        v = np.max(
            [ring_votes(black, white, r - d) for d in (0.0, 0.5, 1.0, 1.5)], axis=0
        )
        fill = disc_fill(black, white, max(r - 4.0, 1.5) if kind == "open" else r - 2.0)
        size = int(2 * round(0.7 * r) + 1)
        peak = (v == ndimage.maximum_filter(v, size=size)) & (v >= vote)
        for y, x in zip(*np.where(peak), strict=True):
            f = float(fill[y, x])
            if kind == "open" and r_open == r_solid:
                k = "solid" if f > 0.7 else ("open" if f < 0.3 else "unclassified")
            elif kind == "open":
                k = "open" if f < 0.3 else None
            else:
                k = "solid" if f > 0.7 else None
            if k:
                found.append((float(x), float(y), k, float(v[y, x]), f))
        if r_open == r_solid:
            break
    # one marker per cluster of peaks (a plateau or an anti-aliased edge gives several)
    found.sort(key=lambda m: -m[3])
    out: list[tuple] = []
    for m in found:
        r = r_open if m[2] == "open" else r_solid
        if all((m[0] - o[0]) ** 2 + (m[1] - o[1]) ** 2 > (0.9 * r) ** 2 for o in out):
            out.append(m)
    return out


def error_bar(
    grey: np.ndarray, x: float, y: float, r: float
) -> tuple[float, float] | None:
    """The centre row and half-length in pixels of the error bar of the marker at (x, y), or
    None when the bar is hidden by the marker.

    The bar ends in caps, grey rows at least a marker radius wide through the marker's column;
    the bar is the nearest pair of caps about equally far (within half a radius) above and below
    the marker's centre. The bar is drawn about the mean, so its centre refines the marker's."""
    c = int(round(x))
    caps = []
    for row in range(grey.shape[0]):
        if not grey[row, c - 1 : c + 2].any():
            continue
        a = b = c
        while a > 0 and grey[row, a - 1]:
            a -= 1
        while b < grey.shape[1] - 1 and grey[row, b + 1]:
            b += 1
        if b - a + 1 >= r:
            caps.append(row)
    centres = groups(np.array(caps, dtype=int))
    pairs = [
        (0.5 * (d - u), 0.5 * (u + d))
        for u in centres
        if u < y
        for d in centres
        if d > y and abs((y - u) - (d - y)) <= 0.5 * r
    ]
    return min(pairs, key=lambda q: q[0])[::-1] if pairs else None


def extract(pdf: Path, check: Path | None) -> tuple[list[dict], dict]:
    doc = pymupdf.open(pdf)
    rows, report = [], {}
    for p in PANELS:
        g = image(doc, p["page"], p["image"])
        dark = g < 110
        fr = frames(dark)[p["panel"]]
        x0, sx, rx, nx = axis(dark, fr, True, p["x_step"])
        y0, sy, ry, ny = axis(dark, fr, False, p["y_step"])
        black = g < p["black"]
        left, right, top, bottom = (int(round(v)) for v in fr)
        inner = np.zeros_like(black)
        inner[top + 3 : bottom - 2, left - 12 : right + 12] = True
        if p["figure"] not in ("10", "13"):  # no marker on the side lines
            inner[:, : left + 3] = False
            inner[:, right - 2 :] = False
        for a, b, c, d in p["masks"]:
            inner[b:d, a:c] = False
        black &= inner
        white = (g > 215) & inner
        # reference lines drawn across the frame are no evidence for or against a marker
        h, w = g.shape
        for q in groups(np.where(dark.sum(axis=1) > 0.35 * w)[0])[1:-1]:
            black[int(q) - 2 : int(q) + 3] = False
            white[int(q) - 2 : int(q) + 3] = False
        found = markers(black, white, *p["radius"], p.get("vote", 0.8))
        grey = None
        if p["bars"]:
            lo, hi = p["bars"]
            grey = (g >= lo) & (g <= hi) & inner
        snaps, dropped, off_grid, shifts, panel_rows = [], [], [], [], []
        for x, y, kind, vote, fill in sorted(found):
            u = (x - x0) / sx
            val = (y - y0) / sy
            ug = round(u / U_GRID) * U_GRID
            if kind == "unclassified":
                dropped.append(
                    {"U": round(u, 3), "value": round(val, 3), "fill": round(fill, 2)}
                )
                continue
            if abs(u - ug) > 0.25 * U_GRID:  # between two tested speeds: not a marker
                off_grid.append(
                    {"U": round(u, 3), "value": round(val, 3), "kind": kind}
                )
                continue
            snaps.append(abs(u - ug))
            row = {
                "figure": p["figure"],
                "quantity": p["quantity"],
                "current": p["current"],
                "where": p["where"],
                "m_star": "1.84" if kind == "solid" else "1.12",
                "U_mps": f"{ug:.2f}",
                "value": round(float(val), 4),
                "sd": "",
            }
            bar = None if grey is None else error_bar(grey, x, y, p["radius"][0])
            if bar is not None:
                shifts.append(abs(bar[0] - y))
                row["value"] = round(float((bar[0] - y0) / sy), 4)
                row["sd"] = round(abs(bar[1] / sy), 3)
            panel_rows.append(row)
        panel_rows.sort(key=lambda q: (q["m_star"], float(q["U_mps"]), q["value"]))
        rows += panel_rows
        anchor = {}
        if p["figure"] in ("9a", "9b"):
            # the line f = 0.11 U/d drawn in the panel, read where no marker or tick mark is
            # near it (below the markers, and in (b) above them at the higher speeds), as the
            # slope through the origin
            us, fs = [], []
            for c in range(left + 4, right - 3):
                u = (c - x0) / sx
                if not (
                    0.07 <= u <= 0.15 or (p["figure"] == "9b" and 1.1 <= u <= 1.45)
                ):
                    continue
                rr = np.where(dark[top + 14 : bottom - 14, c])[0] + top + 14
                if len(rr) and np.ptp(rr) < 6:
                    us.append(u)
                    fs.append((float(np.mean(rr)) - y0) / sy)
            us_a, fs_a = np.array(us), np.array(fs)
            anchor = {
                "line_fd_over_U": round(
                    float(us_a @ fs_a / (us_a @ us_a)) * DIAMETER, 4
                ),
                "line_columns": len(us),
                "stated": 0.11,
            }
        if p["figure"] == "13":
            # the lines of the total submerged weight times the water depth, W h
            h, w = g.shape
            lines = groups(np.where(dark.sum(axis=1) > 0.35 * w)[0])[1:-1]
            anchor = {
                "W_h_lines_N": [round(float((q - y0) / sy), 3) for q in lines],
                "table_1_W_h_N": [round(wt * DEPTH, 3) for wt in (20.9, 2.9)],
            }
        report[p["figure"]] = {
            **({"anchor": anchor} if anchor else {}),
            "x_ticks": nx,
            "x_residual": round(rx, 4),
            "y_ticks": ny,
            "y_residual": round(ry, 4),
            "points": {
                m: int(sum(q["m_star"] == m for q in panel_rows))
                for m in ("1.12", "1.84")
            },
            "sd_read": {
                m: int(sum(q["m_star"] == m and q["sd"] != "" for q in panel_rows))
                for m in ("1.12", "1.84")
            }
            if grey is not None
            else None,
            "unclassified": dropped,
            "off_grid": off_grid,
            "largest_U_snap": round(max(snaps), 4) if snaps else None,
            "largest_bar_centre_shift": round(max(shifts) / abs(sy), 4)
            if shifts
            else None,
        }
        if check is not None:
            check.mkdir(parents=True, exist_ok=True)
            vis = Image.fromarray(g.astype(np.uint8)).convert("RGB")
            dr = ImageDraw.Draw(vis)
            for x, y, kind, _v, _f in found:
                r = p["radius"][0] if kind == "open" else p["radius"][1]
                col = {"open": (0, 90, 255), "solid": (230, 0, 0)}.get(
                    kind, (255, 160, 0)
                )
                dr.ellipse((x - r, y - r, x + r, y + r), outline=col, width=2)
                bar = None if grey is None else error_bar(grey, x, y, p["radius"][0])
                if bar is not None:
                    yc, half = bar
                    dr.line(
                        (x + 3, yc - half, x + 3, yc + half), fill=(0, 180, 0), width=1
                    )
            for k in range(0, 40):
                xx = x0 + k * p["x_step"] * sx
                if left - 1 <= xx <= right + 1:
                    dr.line((xx, bottom - 20, xx, bottom), fill=(255, 0, 255))
                yy = y0 + k * p["y_step"] * sy
                if top - 1 <= yy <= bottom + 1:
                    dr.line((left, yy, left + 20, yy), fill=(255, 0, 255))
            vis.save(check / f"fig{p['figure']}_check.png")
    return rows, report


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("pdf", type=Path)
    ap.add_argument("--out", type=Path, default=HERE / "measured_chaplin2018.csv")
    ap.add_argument(
        "--check", type=Path, help="write the markers found over each figure into DIR"
    )
    a = ap.parse_args()
    rows, report = extract(a.pdf, a.check)
    with open(a.out, "w", encoding="utf-8", newline="\n") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(json.dumps(report, indent=1))


if __name__ == "__main__":
    main()
