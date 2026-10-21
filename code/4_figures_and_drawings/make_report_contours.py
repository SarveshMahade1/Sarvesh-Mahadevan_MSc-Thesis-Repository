# -*- coding: utf-8 -*-
"""
CONTOUR FIGURES FOR THE REVISED REPORT (October 2026)

Rebuilds the contour figures of the report from the single Abaqus captures
(ABAQUS_capture_all.py output, <AC>_<case>_f<fill>_<field>_front.png) so that
they meet the review comments:

  - no overhead titles: each frame carries only a panel label (a), (b), ...,
    and everything else is explained in the figure caption of the report;
  - legible legends: the Abaqus legend, a small raster, is replaced by a
    vector colour bar with the SAME discrete levels and colours. The levels
    are read from the legend by OCR and the colours are sampled from its
    swatches, so the bar is a faithful copy, not a re-scaling;
  - values in report units: stress in MPa, temperature in K, flux in W/m2.

Inputs : a folder of captures, CAPS/<AC>/<AC>_<case>_f<ff>_<field>_front.png
Outputs: <out>.pdf and <out>.png in the figure folder.

    python make_report_contours.py <captures folder> <figure folder>
"""
import os
import re
import subprocess
import sys

import numpy as np
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

CAPS = sys.argv[1] if len(sys.argv) > 1 else 'caps'
OUT = sys.argv[2] if len(sys.argv) > 2 else '.'
TEXT_W = 6.30          # report text width in inches (A4, 2.5 cm margins)
NAME = {'E190': 'E190', 'A320': 'A320', 'A350': 'A350-900'}
UNIT = {'mises': ('MPa', 1e-6), 'nt11': ('K', 1.0), 'hfl': (r'W/m$^2$', 1.0)}
SHORT = {'mises': r'$\sigma_{\mathrm{vM}}$ [MPa]', 'nt11': r'$T$ [K]', 'hfl': r'$q$ [W/m$^2$]'}
LABEL = {'mises': 'von Mises stress', 'nt11': 'Wall temperature', 'hfl': 'In-plane heat flux'}

plt.rcParams.update({'font.family': 'serif', 'font.serif': ['Liberation Serif', 'TeX Gyre Termes', 'Times'], 'pdf.fonttype': 42,
                     'mathtext.fontset': 'stix', 'font.size': 8, 'axes.linewidth': 0.5})


def runs(mask_cols):
    out, start = [], None
    for x, v in enumerate(mask_cols):
        if v and start is None:
            start = x
        if not v and start is not None:
            out.append((start, x))
            start = None
    if start is not None:
        out.append((start, len(mask_cols)))
    return [r for r in out if r[1] - r[0] > 3]


def split(path):
    """Return (legend array, tank array) of one capture."""
    im = np.asarray(Image.open(path).convert('RGB')).astype(int)
    nw = im.min(2) < 235
    rr = runs(nw.sum(0) > 0)
    leg = im[:, rr[0][0]:rr[0][1]]
    x0, x1 = rr[1][0], rr[-1][1]
    tank = im[:, x0:x1]
    rows = np.where((tank.min(2) < 235).sum(1) > 0)[0]
    tank = tank[max(rows[0] - 4, 0):rows[-1] + 5]
    return leg, tank


def legend_levels(leg, tmp='/tmp/_leg.png'):
    """Levels (OCR) and band colours (sampled) of an Abaqus legend."""
    sat = (leg.max(2) - leg.min(2)) > 60
    cols = np.where(sat.sum(0) > 20)[0]
    sx0, sx1 = cols[0], cols[-1]
    xs = (sx0 + sx1) // 2
    col = leg[:, xs]
    satc = sat[:, xs]
    ys = np.where(satc)[0]
    y0, y1 = ys[0], ys[-1]
    # bands: consecutive rows of (nearly) the same colour
    bands, cur, start = [], None, y0
    for y in range(y0, y1 + 2):
        c = tuple(col[y]) if y <= y1 else None
        if cur is None:
            cur, start = c, y
            continue
        if c is None or max(abs(a - b) for a, b in zip(c, cur)) > 25:
            if y - start >= 3:
                bands.append((start, y, cur))
            cur, start = c, y
    colours = [tuple(v / 255.0 for v in b[2]) for b in bands]       # top to bottom
    # OCR of the labels to the right of the swatches
    txt = leg[:, sx1 + 2:]
    Image.fromarray(txt.astype(np.uint8)).resize((txt.shape[1] * 3, txt.shape[0] * 3),
                                                 Image.LANCZOS).save(tmp)
    s = subprocess.run(['tesseract', tmp, '-', '--psm', '6'], capture_output=True, text=True).stdout
    s = s.replace('/7', '7').replace('/', '7').replace('O', '0').replace('o', '0')
    raw = sorted({float(v.replace(' ', '')) for v in
                  re.findall(r'[+-]\s?\d\.\d{3}e[+-]\d{2}', s)}, reverse=True)
    # Abaqus spaces the in-range labels uniformly; a label above the uniform
    # spacing is the maximum of an "above range" band, drawn in grey.
    n = len(colours)
    over = None
    if len(raw) >= 3:
        d1, d2 = raw[0] - raw[1], raw[1] - raw[2]
        if abs(d1 - d2) > 0.05 * abs(d2):
            g = leg[max(y0 - 6, 0), xs]
            over = tuple(v / 255.0 for v in g) if (g.max() - g.min()) < 25 and g.max() < 235 else (0.75, 0.75, 0.75)
            raw = raw[1:]
    top, bot = raw[0], raw[-1]
    vals = list(np.linspace(bot, top, n + 1))
    # every OCR label must sit on the uniform grid, otherwise stop
    stepv = (top - bot) / n
    for v in raw:
        k = round((v - bot) / stepv)
        if abs(v - (bot + k * stepv)) > 0.02 * stepv + 1e-9 * abs(top):
            raise RuntimeError('legend label %g off the uniform grid' % v)
    return vals, colours, over


def colourbar(ax, vals, colours, over, field):
    unit, sc = UNIT[field]
    v = np.array(sorted(vals)) * sc
    cols = list(reversed(colours))                                   # bottom to top
    n = len(v) - 1
    if len(cols) != n:
        raise RuntimeError('legend parse: %d levels, %d colours' % (len(v), len(cols)))
    cmap = ListedColormap(cols)
    if over is not None:
        cmap.set_over(over)
    norm = BoundaryNorm(v, n)
    cb = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=ax,
                      extend='max' if over is not None else 'neither', extendfrac=0.08)
    step = 2 if n > 8 else 1
    ticks = v[::step]
    cb.set_ticks(ticks)
    fmt = '%.0f' if (field != 'nt11' and v[-1] - v[0] > 20) else '%.1f' if field == 'nt11' else '%.0f'
    cb.set_ticklabels([fmt % t for t in ticks])
    cb.ax.tick_params(labelsize=7, length=2, width=0.4)
    cb.outline.set_linewidth(0.4)
    cb.set_label(SHORT[field], fontsize=7, labelpad=2)


def figure(rows, out, width=0.80):
    """rows = [(label, path, legend_path, field), ...]; consecutive rows sharing a
    legend_path are drawn against one colour bar spanning them."""
    tanks = [(lab, split(p)[1], lp, fld) for lab, p, lp, fld in rows]
    tank_w = TEXT_W * width
    hs = [tank_w * t.shape[0] / t.shape[1] for _, t, _, _ in tanks]
    gap = 0.06
    H = sum(hs) + gap * (len(hs) - 1) + 0.05
    FW = tank_w + 1.0
    fig = plt.figure(figsize=(FW, H))
    y = H
    groups = {}
    for (lab, t, lp, fld), h in zip(tanks, hs):
        ax = fig.add_axes([0.0, (y - h) / H, tank_w / FW, h / H])
        ax.imshow(t.astype(np.uint8), interpolation='lanczos')
        ax.axis('off')
        ax.text(0.012, 0.97, lab, transform=ax.transAxes, ha='left', va='top', fontsize=8,
                bbox=dict(boxstyle='square,pad=0.15', fc='white', ec='none'))
        groups.setdefault(lp, [fld, []])[1].append((y - h, y))
        y -= h + gap
    x_cb = (tank_w + 0.12) / FW
    for lp, (fld, spans) in groups.items():
        vals, cols, over = legend_levels(split(lp)[0])
        lo, hi = min(s[0] for s in spans), max(s[1] for s in spans)
        cax = fig.add_axes([x_cb, (lo + 0.06 * (hi - lo)) / H, 0.13 / FW, 0.88 * (hi - lo) / H])
        colourbar(cax, vals, cols, over, fld)
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUT, out + '.' + ext), dpi=300)
    plt.close(fig)
    print('wrote', out)


def cap(ac, case, f, field):
    return os.path.join(CAPS, ac, '%s_%s_f%02d_%s_front.png' % (ac, case, f, field))


# (R_mid, L_barrel) in m, as AIRCRAFT_DATA in the master script; head height R/1.6
GEO = {'E190': (1.3015, 13.9675), 'A320': (1.7080, 10.0120), 'A350': (2.5771, 33.8240)}


def models_figure(out):
    """(a)-(c) the three inner vessels to one length scale, as silhouettes of
    the MESH captures, with a tick at the mid ring; (d)-(f) the central part
    of the ZOOM capture of each tank at the mid ring (LC7, fill 0.90), where
    the individual elements can be seen. No overhead titles: the caption of
    the report carries the dimensions, seeds and element counts."""
    Lt = {a: L + 2 * R / 1.6 for a, (R, L) in GEO.items()}
    Lmax = max(Lt.values())
    W = TEXT_W
    scale = (W - 0.3) / Lmax                         # inch per metre
    hs = [2 * GEO[a][0] * scale for a in ('E190', 'A320', 'A350')]
    gap, sb, zw = 0.10, 0.30, (W - 0.2) / 3.0
    crops = []
    for a in ('E190', 'A320', 'A350'):
        t = split(os.path.join(CAPS, a, '%s_ZOOM_LC7_Emergency_f90_mises_front_mesh.png' % a))[1]
        w = t.shape[1]
        crops.append(t[:, int(0.25 * w):int(0.75 * w)])
    zh = zw * crops[0].shape[0] / crops[0].shape[1]
    H = sum(hs) + gap * 3 + sb + 0.12 + zh + 0.05
    fig = plt.figure(figsize=(W, H))
    y = H - 0.02
    for k, (a, h) in enumerate(zip(('E190', 'A320', 'A350'), hs)):
        R, L = GEO[a]
        im = np.asarray(Image.open(os.path.join(CAPS, a, '%s_MESH_front.png' % a)).convert('RGB')).astype(int)
        m = im.min(2) < 235
        rr, cc = np.where(m.any(1))[0], np.where(m.any(0))[0]
        m = m[rr[0]:rr[-1] + 1, cc[0]:cc[-1] + 1]
        rgb = np.full(m.shape + (3,), 255, np.uint8)
        rgb[m] = (175, 175, 175)
        ax = fig.add_axes([0.25 / W, (y - h) / H, Lt[a] * scale / W, h / H])
        ax.imshow(rgb, interpolation='lanczos', aspect='auto', extent=(0, Lt[a], -R, R))
        xc = R / 1.6 + L / 2.0                          # mid ring at mid barrel
        ax.plot([xc, xc], [-R, R], color='k', lw=0.8, ls='--')
        ax.set_xlim(0, Lt[a]); ax.set_ylim(-R, R); ax.axis('off')
        fig.text(0.0, (y - 0.02) / H, '(%s)' % 'abc'[k], ha='left', va='top', fontsize=8)
        y -= h + gap
    ax = fig.add_axes([0.25 / W, (y - sb) / H, Lmax * scale / W, sb / H])
    ax.set_xlim(0, Lmax); ax.set_ylim(0, 1); ax.axis('off')
    ax.plot([0, 5], [0.8, 0.8], color='k', lw=1.0)
    for xx in (0, 5):
        ax.plot([xx, xx], [0.65, 0.95], color='k', lw=0.8)
    ax.text(2.5, 0.15, '5 m', ha='center', va='bottom', fontsize=7.5)
    y -= sb + 0.12
    for k, c in enumerate(crops):
        ax = fig.add_axes([(0.0 + k * (zw + 0.1)) / W, (y - zh) / H, zw / W, zh / H])
        ax.imshow(c.astype(np.uint8), interpolation='lanczos', aspect='auto')
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_linewidth(0.4)
        ax.text(0.02, 0.96, '(%s)' % 'def'[k], transform=ax.transAxes, ha='left', va='top', fontsize=8,
                bbox=dict(boxstyle='square,pad=0.15', fc='white', ec='none'))
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUT, out + '.' + ext), dpi=300)
    plt.close(fig)
    print('wrote', out)


def main():
    os.makedirs(OUT, exist_ok=True)
    models_figure('fig4_models')
    # E190 reference case: three fields, three panels, one bar each
    rows = []
    for k, (field, case) in enumerate((('nt11', 'Baseline'), ('hfl', 'Baseline'),
                                        ('mises', 'LC7_Emergency'))):
        p = cap('E190', case, 90, field)
        rows.append(('(%s)' % 'abc'[k], p, p, field))
    figure(rows, 'fig5_e190_fields', width=0.74)
    # three aircraft, two fills: one bar per aircraft
    for case, field, out in (('Baseline', 'nt11', 'fig5_thermal_three'),
                             ('LC7_Emergency', 'mises', 'fig5_emergency_three'),
                             ('Baseline', 'hfl', 'fig5_heatflux_three')):
        rows, k = [], 0
        for ac in ('E190', 'A320', 'A350'):
            lp = cap(ac, case, 90, field)
            for f in (10, 90):
                rows.append(('(%s)' % 'abcdef'[k], cap(ac, case, f, field), lp, field))
                k += 1
        figure(rows, out, width=0.62)
    # appendix: fill sweep 0.10, 0.50, 0.90 per aircraft and field
    for ac in ('E190', 'A320', 'A350'):
        for case, field, tag in (('Baseline', 'nt11', 'temperature'), ('Baseline', 'hfl', 'heatflux'),
                                 ('LC7_Emergency', 'mises', 'emergency')):
            lp = cap(ac, case, 90, field)
            rows = [('(%s)' % 'abc'[i], cap(ac, case, f, field), lp, field)
                    for i, f in enumerate((10, 50, 90))]
            figure(rows, 'appx_contours_%s_%s' % (ac, tag), width=0.55)


if __name__ == '__main__':
    main()
