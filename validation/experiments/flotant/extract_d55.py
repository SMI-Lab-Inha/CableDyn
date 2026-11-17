#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 Jae Hoon Seo, SMI Lab, Inha University
"""Extract the data of the FLOTANT D5.5 report used by flotant_viv.py from the report's PDF.

    python extract_d55.py path/to/210430_FLT-D5.5_VIV_test_v0.pdf [OUTDIR]

Writes, by default next to this script:
  fem_air.csv, fem_water.csv  the node-by-node FEM input tables (Sec. 4.5, in air, and Sec. 5.4,
                              in water): element label, position, accelerometer station, X Y Z,
                              constraint flags, lumped mass [kg], weight [N] (submerged in water),
                              static tension [N], EA [N], EI [N m^2]
  measured_rao.csv            the RAOs of Appendix 1 (in air, 200 kg, surge and sway) and Appendix
                              2 (in water, 200 kg, surge and sway) at 5 and 10 mm
  measured_tow_ad.csv         A/D_rms against Ur of the 200 kg tow tests, Appendix 3

The RAO charts are vector graphics: a series is a polyline with filled markers of one colour,
named by the legend entry beside the same marker; data coordinates follow from the tick-label
words by a linear map per chart (the frequency axis is labelled "[rad]" and is in rad/s). The
Appendix 3 charts are embedded raster images (1386 x 867 px) with the tick labels inside the
raster: their axes are calibrated on the light-grey gridlines (Ur -10, -5, 0, 5, 10; A/D from 0
in the steps of PLOTS) and the markers are the tab:blue discs. Each digitised Ur is checked
against the test speeds of Table 6-2 (multiples of 0.025 m/s, V = |Ur| OD f1, OD 0.048 m,
f1 1.77 Hz).

Needs PyMuPDF, NumPy, SciPy and Pillow; the report is not distributed with CableDyn.
"""

from __future__ import annotations

import csv
import io
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pymupdf
from PIL import Image
from scipy import ndimage
from scipy.cluster.vq import kmeans2

HERE = Path(__file__).resolve().parent
OD, F1 = 0.048, 1.77

# ---------------------------------------------------------------- FEM input tables
FEM_PAGES = {"air": (58, 59), "water": (49, 50)}
ROW = re.compile(
    r"(?P<node>\d+)\s+(?P<el>[A-Za-z]+)\s+(?P<pos>-?\d+\.\d+)\s+(?:(?P<acc>\d\.\d+L?|L|Hex|Top)\s+)?"
    r"(?P<x>-?\d+\.\d+)\s+(?P<y>-?\d+\.\d+)\s+(?P<z>-?\d+\.\d+)\s+(?P<c1>[01])(?=\s)\s+(?P<c2>[01])(?=\s)\s+"
    r"(?P<c3>[01])(?=\s)\s+(?P<c4>[01])(?=\s)"
    r"(?:\s+(?P<m>-?[\d.]+)\s+(?P<w>-?[\d.]+)\s+(?P<t>-?[\d.]+)\s+(?P<ea>[\d.E+-]+)\s+(?P<ei>[\d.E+-]+))?"
)


def fem_tables(doc, out: Path) -> None:
    for medium, pages in FEM_PAGES.items():
        text = re.sub(r"\s+", " ", " ".join(doc[p - 1].get_text() for p in pages))
        rows = []
        for m in ROW.finditer(text):
            d = m.groupdict()
            if rows and int(d["node"]) != rows[-1]["node"] + 1:
                continue
            if not rows and int(d["node"]) != 0:
                continue
            rows.append(
                {
                    "node": int(d["node"]),
                    "element": d["el"],
                    "pos": float(d["pos"]),
                    "acc": d["acc"] or "",
                    "x": float(d["x"]),
                    "y": float(d["y"]),
                    "z": float(d["z"]),
                    "fixed": "".join(d[k] for k in ("c1", "c2", "c3", "c4")),
                    "mass": float(d["m"]) if d["m"] else None,
                    "weight": float(d["w"]) if d["w"] else None,
                    "tension": float(d["t"]) if d["t"] else None,
                    "EA": float(d["ea"]) if d["ea"] and d["ea"] not in "-" else None,
                    "EI": float(d["ei"]) if d["ei"] and d["ei"] not in "-" else None,
                }
            )
        with open(out / f"fem_{medium}.csv", "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
        cab = [r for r in rows if r["element"] == "Cable"]
        print(
            f"fem_{medium}.csv: {len(rows)} nodes, {len(cab)} cable nodes, total lumped mass "
            f"{sum(r['mass'] or 0 for r in rows):.2f} kg"
        )


# ---------------------------------------------------------------- RAOs (Appendices 1 and 2)
RAO_PAGES = {63: ("air", "X"), 64: ("air", "Y"), 67: ("water", "X"), 68: ("water", "Y")}
NUM = re.compile(r"^-?\d+,\d+$")


def num(s: str) -> float:
    return float(s.replace(",", "."))


def linear_fit(pairs):
    n = len(pairs)
    mx = sum(p for p, _ in pairs) / n
    mv = sum(v for _, v in pairs) / n
    a = sum((p - mx) * (v - mv) for p, v in pairs) / sum(
        (p - mx) ** 2 for p, _ in pairs
    )
    res = max(abs(mv + a * (p - mx) - v) for p, v in pairs)
    return (lambda q: mv + a * (q - mx)), res


def axis_map(words, box):
    """Linear maps (x -> data, y -> data) from the tick labels around the chart box."""
    x0, y0, _, y1 = box
    yt = [
        ((w[1] + w[3]) / 2, num(w[4]))
        for w in words  # y labels carry two decimals, x labels one
        if re.match(r"^-?\d+,\d\d$", w[4])
        and w[2] < x0 + 30
        and y0 - 8 < (w[1] + w[3]) / 2 < y1 + 8
    ]
    # x ticks: the row of numbers on one baseline with the most entries in the chart's band
    rows = defaultdict(list)
    for w in words:
        if NUM.match(w[4]) and w[0] > x0 - 10 and y0 - 8 < (w[1] + w[3]) / 2 < y1 + 20:
            rows[round(w[1])].append(w)
    xt = [((w[0] + w[2]) / 2, num(w[4])) for w in max(rows.values(), key=len)]
    fx, rx = linear_fit(xt)
    fy, ry = linear_fit(yt)
    return fx, fy, rx, ry


def raos(doc, out: Path) -> None:
    rows_all = []
    for pg, (medium, motion) in RAO_PAGES.items():
        page = doc[pg - 1]
        words = page.get_text("words")
        ticks = sorted(
            (w for w in words if NUM.match(w[4]) and w[0] < 110), key=lambda w: w[1]
        )
        mid = page.rect.height / 2  # the RAO chart is the upper one of the page
        upper = [w for w in ticks if w[1] < mid]
        ys = [(w[1] + w[3]) / 2 for w in upper]
        xs = [w for w in words if NUM.match(w[4]) and w[1] < mid and w[0] > 110]
        box = (min(w[2] for w in upper) + 2, min(ys), max(w[2] for w in xs), max(ys))
        fx, fy, rx, ry = axis_map(words, box)
        markers: dict[tuple, list] = defaultdict(list)
        for g in page.get_drawings():
            if g.get("type") != "fs" or g.get("fill") is None:
                continue
            r = g["rect"]
            if r.width > 8 or r.height > 8:
                continue
            markers[tuple(round(c, 3) for c in g["fill"])].append(
                ((r.x0 + r.x1) / 2, (r.y0 + r.y1) / 2)
            )
        legend: dict[tuple, str] = {}
        for col, pts in markers.items():
            for cx, cy in pts:
                if box[3] + 25 < cy < box[3] + 140 and col not in legend:
                    row = sorted(
                        (
                            w
                            for w in words
                            if abs((w[1] + w[3]) / 2 - cy) < 3 and cx < w[0] < cx + 120
                        ),
                        key=lambda w: w[0],
                    )
                    legend[col] = " ".join(w[4] for w in row[:5])
        n = 0
        for col, pts in markers.items():
            m = re.match(
                r"([XY])_(\w+?)\s*\[m\]\s*-\s*amp\s*(\d+)mm", legend.get(col, "")
            )
            if not m or int(m.group(3)) not in (5, 10):
                continue
            for cx, cy in pts:
                if box[0] - 3 <= cx <= box[2] + 3 and box[1] - 3 <= cy <= box[3] + 3:
                    omega = fx(cx)
                    rows_all.append(
                        (
                            medium,
                            motion,
                            m.group(2),
                            int(m.group(3)),
                            round(omega / (2 * np.pi), 1),
                            round(omega, 4),
                            round(fy(cy), 4),
                            pg,
                        )
                    )
                    n += 1
        print(
            f"p{pg} ({medium}, {motion}): {n} points, axis residual {rx:.2g} rad/s, {ry:.2g}"
        )
    rows_all.sort()
    with open(out / "measured_rao.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(
            [
                "medium",
                "motion",
                "station",
                "amp_mm",
                "f_hz",
                "omega_rad_s",
                "rao",
                "page",
            ]
        )
        w.writerows(rows_all)


# ---------------------------------------------------------------- tow tests (Appendix 3)
PLOTS = {  # (page, image index): component, station, A/D gridline step
    (70, 2): ("X", "1_2L", 0.05),
    (70, 3): ("Y", "1_2L", 0.1),
    (71, 2): ("X", "3_4L", 0.05),
    (71, 3): ("Y", "3_4L", 0.1),
    (72, 2): ("X", "L", 0.02),
    (72, 3): ("Y", "L", 0.1),
}
UR_GRID = [-10.0, -5.0, 0.0, 5.0, 10.0]


def gridlines(mask: np.ndarray, axis: int, min_frac: float) -> list[float]:
    """Centres of the runs of rows (axis=1) or columns (axis=0) that are mostly gridline grey."""
    frac = mask.mean(axis=axis)
    idx = np.flatnonzero(frac > min_frac * frac.max())
    runs, cur = [], [idx[0]]
    for i in idx[1:]:
        if i == cur[-1] + 1:
            cur.append(i)
        else:
            runs.append(cur)
            cur = [i]
    runs.append(cur)
    return [float(np.mean(r)) for r in runs]


def digitise(image: Image.Image, step: float):
    a = np.asarray(image.convert("RGB")).astype(int)
    grey = np.all(np.abs(a - 230) <= 3, axis=2)
    blue = (
        (np.abs(a[..., 0] - 31) < 25)
        & (np.abs(a[..., 1] - 119) < 25)
        & (np.abs(a[..., 2] - 180) < 25)
    )
    cols = gridlines(grey, axis=0, min_frac=0.6)
    rows = sorted(
        gridlines(grey, axis=1, min_frac=0.6), reverse=True
    )  # bottom (A/D 0) first
    if len(cols) != len(UR_GRID):
        raise SystemExit(f"found {len(cols)} vertical gridlines {cols}")
    yvals = [step * k for k in range(len(rows))]
    ax, bx = np.polyfit(cols, UR_GRID, 1)
    ay, by = np.polyfit(rows, yvals, 1)
    resx = max(abs(ax * c + bx - u) for c, u in zip(cols, UR_GRID))
    lab, n = ndimage.label(blue)
    sizes = ndimage.sum(blue, lab, range(1, n + 1))
    med = float(np.median([s for s in sizes if s > 60]))
    pts = []
    for i in [i + 1 for i, s in enumerate(sizes) if s > 0.4 * med]:
        yy, xx = np.nonzero(lab == i)
        cx, cy = xx.mean(), yy.mean()
        if cy < 0.08 * a.shape[0] and cx > 0.85 * a.shape[1]:
            continue  # the legend marker
        k = round(len(xx) / med)
        if (
            k >= 2
        ):  # overlapping markers (the slowest speeds): split the blob into k discs
            xy = np.column_stack([xx, yy]).astype(float)
            init = np.column_stack(
                [np.linspace(xx.min() + 7, xx.max() - 7, k), np.full(k, cy)]
            )
            cent, _ = kmeans2(xy, init, minit="matrix", seed=1)
            pts += [(ax * cent[j, 0] + bx, ay * cent[j, 1] + by) for j in range(k)]
            continue
        pts.append((ax * cx + bx, ay * cy + by))
    return sorted(pts), resx


def tow(doc, out: Path) -> None:
    rows = []
    for (pg, idx), (comp, st, step) in PLOTS.items():
        xref = doc[pg - 1].get_images(full=True)[idx][0]
        image = Image.open(io.BytesIO(pymupdf.Pixmap(doc, xref).tobytes("png")))
        pts, resx = digitise(image, step)
        print(
            f"p{pg} image {idx} ({comp}, {st}): {len(pts)} points, Ur axis residual {resx:.3f}"
        )
        for ur, ad in pts:
            v = abs(ur) * OD * F1
            rows.append(
                (
                    st,
                    comp,
                    "forwards" if ur > 0 else "backwards",
                    round(abs(ur), 3),
                    round(v, 4),
                    f"{round(v / 0.025) * 0.025:.3f}",
                    round(ad, 4),
                )
            )
    dv = max(abs(r[4] - float(r[5])) for r in rows)
    print(f"speed check: |V - nearest multiple of 0.025 m/s| at most {dv:.4f} m/s")
    with open(out / "measured_tow_ad.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(
            [
                "station",
                "component",
                "direction",
                "Ur",
                "V_digitised",
                "V_test",
                "AD_rms",
            ]
        )
        w.writerows(sorted(rows, key=lambda r: (r[0], r[1], float(r[5]))))


def main() -> None:
    doc = pymupdf.open(sys.argv[1])
    out = Path(sys.argv[2]) if len(sys.argv) > 2 else HERE
    out.mkdir(parents=True, exist_ok=True)
    fem_tables(doc, out)
    raos(doc, out)
    tow(doc, out)


if __name__ == "__main__":
    main()
