# -*- coding: utf-8 -*-
"""
================================================================================
MULTI-AIRCRAFT CONTOUR PANELS FOR CHAPTERS 4, 5 AND THE APPENDIX
TU Delft MSc thesis - Sarvesh Mahadevan
================================================================================

Assembles the Abaqus renders in abaqus_captures_v2 into report figures that
show the SAME case for all three aircraft, so the reader sees that one
procedure produced three comparable results.

  - one legend per aircraft, shared by every frame of that aircraft in the
    panel (the contour window is fixed per aircraft and load case, so the
    frames genuinely share it);
  - legend and tank are cut from separate renders (see compose_side_legend.py),
    so a legend can never cover a tank;
  - the mesh figure is drawn to ONE length scale, so the three tanks can be
    compared by eye.

Inputs : abaqus_captures_v2\\_raw\\<AC>\\<name>__L.png / __N.png  (contours)
         abaqus_captures_v2\\<AC>\\<AC>_MESH_front.png            (meshes)
         abaqus_captures_v2\\capture_log.txt                      (element counts)
Outputs: latex\\figures\\37_... to 40_... and appx_contours_<AC>_*.{pdf,png}

RUN (from the thesis root):   python python_scripts\\make_contour_panels.py
================================================================================
"""
import os
import re
import sys

import numpy as np
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import os as _os
import sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import repo_paths as RP  # noqa: E402  (folder layout of the repository)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)
import thesis_figstyle as tfs                     # noqa: E402
import compose_side_legend as csl                 # noqa: E402

CAP = RP.CAPTURES   # 28 Sep: sized walls
FIGDIR = RP.FIGURES
NAVY, GREY = '#000000', '#000000'    # 28 Sep: all figure text black, as in the report

ORDER = ['E190', 'A320', 'A350']                  # table order of the report
NAME = {'A320': 'A320', 'E190': 'E190', 'A350': 'A350-900'}
# (R_mid, L_barrel, cap height) in m, as AIRCRAFT_DATA in the master script
GEO = {'A320': (1.7080, 10.0120), 'E190': (1.3015, 13.9675),
       'A350': (2.5771, 33.8240)}
MAX_PX = 1800                                     # downsample for file size


def _small(img):
    if img.width > MAX_PX:
        img = img.resize((MAX_PX, int(round(img.height * MAX_PX / img.width))),
                         Image.LANCZOS)
    return np.asarray(img.convert('RGB'))


def frame(ac, name):
    """(legend, tank) for one render, e.g. frame('A320', 'LC7_Emergency_f90_mises')."""
    base = os.path.join(CAP, '_raw', ac, '%s_%s_front' % (ac, name))
    legend, tank, clipped = csl.pieces(base + '__L.png', base + '__N.png')
    if clipped:
        print('  WARNING: %s %s clipped by the camera' % (ac, name))
    return legend, tank


def panel(blocks, out, width_in=tfs.FIG_W, legend_frac=0.19, gap_in=0.08,
          label_h=0.20, block_gap=0.16):
    """
    blocks = [(title, legend_img, [(sublabel, tank_img), ...]), ...]
    Each block: the legend on the left, the tanks stacked on the right.
    """
    lw = width_in * legend_frac
    tw = width_in - lw - gap_in
    # block heights in inches
    geo = []
    for title, leg, tanks in blocks:
        th = [tw * t.height / t.width for _, t in tanks]
        stack_h = sum(th) + label_h * len(tanks)
        # One legend size for every block, so the legends read alike; a
        # block shorter than its legend is made as tall as the legend.
        lwe = lw
        leg_h = lw * leg.height / leg.width
        stack_h = max(stack_h, leg_h + label_h)
        geo.append((th, stack_h, lwe, leg_h))
    total_h = sum(g[1] for g in geo) + block_gap * (len(blocks) - 1) + 0.05

    fig = plt.figure(figsize=(width_in, total_h))
    y = total_h
    for (title, leg, tanks), (th, stack_h, lwe, leg_h) in zip(blocks, geo):
        # legend, vertically centred on the block, below the block title
        ax = fig.add_axes([0.0, (y - label_h - leg_h * 0.97) / total_h,
                           lwe / width_in, leg_h * 0.97 / total_h])
        ax.imshow(_small(leg), interpolation='lanczos')
        ax.axis('off')
        fig.text(0.0, (y - 0.02) / total_h, title, ha='left', va='top',
                 fontsize=9.5, fontweight='bold', color=NAVY)
        yy = y
        for (sub, t), h in zip(tanks, th):
            yy -= label_h
            fig.text((lw + gap_in) / width_in, (yy + 0.03) / total_h, sub,
                     ha='left', va='bottom', fontsize=8.5, color=GREY)
            ax = fig.add_axes([(lw + gap_in) / width_in, (yy - h) / total_h,
                               tw / width_in, h / total_h])
            ax.imshow(_small(t), interpolation='lanczos')
            ax.axis('off')
            yy -= h
        y -= stack_h + block_gap
    tfs.save(fig, os.path.join(FIGDIR, out))


def models_figure(out, width_in=455.24 / 72.27):
    """
    (a) the three inner vessels drawn to one length scale, as silhouettes;
    (b) the production mesh at the mid ring of each, zoomed so the element
        size can be seen, with the LC7 von Mises field it carries.
    Element counts and seed sizes are read from final_numbers.json.
    """
    import json
    D = json.load(open(os.path.join(RP.PROCESSED, 'final_numbers.json')))
    Lmax = max(L + 2 * R / 1.6 for R, L in GEO.values())
    sil = []
    for ac in ORDER:
        img = Image.open(os.path.join(CAP, ac, '%s_MESH_front.png' % ac)).convert('RGB')
        a = np.asarray(img).astype(np.int16)
        x0, y0, x1, y1 = csl.bbox(a.min(axis=2) < csl.WHITE_LEVEL)
        a = a[y0:y1, x0:x1]
        m = a.min(axis=2) < csl.WHITE_LEVEL
        rgb = np.full(a.shape, 255, np.uint8)
        rgb[m] = (170, 170, 170)
        sil.append((ac, rgb))
    scale = (width_in - 0.1) / Lmax
    heights = [2 * GEO[ac][0] * scale for ac, _ in sil]
    lab = 0.24
    zoom_h = (width_in - 0.2) / 3.0
    WIN = 0.074          # fraction of the tank length visible in the zoomed frame (square crop, zoom 5)
    total_h = 0.2 + sum(heights) + lab * len(sil) + 0.3 + 0.2 + lab + zoom_h + 0.05
    fig = plt.figure(figsize=(width_in, total_h))
    y = total_h
    fig.text(0.0, (y - 0.02) / total_h, '(a) The three vessels to one scale; the box marks the window enlarged in (b)', ha='left', va='top',
             fontsize=10, fontweight='bold')
    y -= 0.2
    for (ac, im), h in zip(sil, heights):
        R, L = GEO[ac]
        Lt = L + 2 * R / 1.6
        y -= lab
        d = D[ac]
        txt = '%s:  $L$ = %.2f m,  $2R$ = %.2f m,  $t$ = %.3f mm,  %s elements' % (
            NAME[ac], Lt, 2 * R, d['t_wall_mm'], format(int(d['production']['Baseline|0.10']['n_el']), ','))
        fig.text(0.05 / width_in, (y + 0.03) / total_h, txt, ha='left', va='bottom', fontsize=9.5)
        ax = fig.add_axes([0.05 / width_in, (y - h) / total_h, Lt * scale / width_in, h / total_h])
        ax.imshow(im, interpolation='lanczos', aspect='auto', extent=(0, Lt, -R, R))
        # window shown enlarged in (b): centred on the mid ring, width WIN of the tank length
        xc = R / 1.6 + L / 2.0
        wz_m = WIN * Lt
        ax.add_patch(plt.Rectangle((xc - wz_m / 2, -min(R, wz_m / 2)), wz_m, 2 * min(R, wz_m / 2),
                                   fill=False, ec='k', lw=1.1))
        ax.set_xlim(0, Lt); ax.set_ylim(-R, R)
        ax.axis('off')
        y -= h
    ax = fig.add_axes([0.05 / width_in, (y - 0.3) / total_h, (width_in - 0.1) / width_in, 0.25 / total_h])
    ax.set_xlim(0, Lmax); ax.set_ylim(0, 1); ax.axis('off')
    ax.plot([0, 5], [0.75, 0.75], color='k', lw=1.5)
    for xx in (0, 5):
        ax.plot([xx, xx], [0.55, 0.95], color='k', lw=1.0)
    ax.text(2.5, 0.05, '5 m', ha='center', va='bottom', fontsize=8.5)
    y -= 0.3 + 0.2
    fig.text(0.0, (y + 0.15) / total_h, '(b) Production mesh at the mid ring, LC7 Emergency, $\\phi_f$ = 0.90',
             ha='left', va='top', fontsize=10, fontweight='bold')
    y -= lab
    wz = (width_in - 0.2) / 3.0
    zoom_h = wz
    for k, ac in enumerate(ORDER):
        img = Image.open(os.path.join(CAP, '_raw', ac, '%s_ZOOM_LC7_Emergency_f90_mises_front_mesh__N.png' % ac)).convert('RGB')
        a = np.asarray(img)
        H_, W_ = a.shape[:2]
        # central square crop: the mid ring sits at the centre of the zoomed view
        side = min(H_, W_)
        c0 = (W_ - side) // 2
        crop = a[(H_ - side) // 2:(H_ + side) // 2, c0:c0 + side]
        ax = fig.add_axes([(k * (wz + 0.1)) / width_in, (y - zoom_h) / total_h, wz / width_in, zoom_h / total_h])
        ax.imshow(crop, interpolation='lanczos', aspect='auto')
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(True); sp.set_linewidth(0.5)
        fig.text((k * (wz + 0.1)) / width_in, (y + 0.02) / total_h,
                 '%s, seed %.1f mm' % (NAME[ac], D[ac]['seed_mm']), ha='left', va='bottom', fontsize=9.5)
    tfs.save(fig, os.path.join(FIGDIR, out))


def two_fill_panel(case, field, fills, out, label_fmt):
    blocks = []
    for ac in ORDER:
        # legend from the highest fill: it carries the out-of-range band, if any
        leg, _ = frame(ac, '%s_f%02d_%s' % (case, fills[-1], field))
        tanks = [(label_fmt % (f / 100.0), frame(ac, '%s_f%02d_%s' % (case, f, field))[1])
                 for f in fills]
        blocks.append((NAME[ac], leg, tanks))
    panel(blocks, out)


def main():
    tfs.apply_style()
    lab = r'$\phi_f$ = %.2f'
    # Chapter 4: the three models and their meshes
    models_figure('fig4_models')
    # Chapter 5: the three fields of the reference case, each with its legend
    blocks = []
    for name, title in (('Baseline_f90_nt11', 'Wall temperature NT11 [K], Baseline, $\\phi_f$ = 0.90'),
                        ('Baseline_f90_hfl', 'Heat flux HFL [W/m$^2$], Baseline, $\\phi_f$ = 0.90'),
                        ('LC7_Emergency_f90_mises', 'Von Mises stress [Pa], LC7 Emergency, $\\phi_f$ = 0.90')):
        leg, tank = frame('E190', name)
        blocks.append((title, leg, [('', tank)]))
    panel(blocks, 'fig5_e190_fields', label_h=0.24)
    # one observation, three aircraft
    two_fill_panel('Baseline', 'nt11', (10, 90), 'fig5_thermal_three', lab)
    two_fill_panel('LC7_Emergency', 'mises', (10, 90), 'fig5_emergency_three', lab)
    two_fill_panel('Baseline', 'hfl', (10, 90), 'fig5_heatflux_three', lab)
    # Appendix: fill sweep per aircraft
    for ac in ORDER:
        for case, field, tag in (('Baseline', 'nt11', 'temperature'),
                                 ('Baseline', 'hfl', 'heatflux'),
                                 ('LC7_Emergency', 'mises', 'emergency')):
            leg, _ = frame(ac, '%s_f90_%s' % (case, field))
            tanks = [(lab % (f / 100.0), frame(ac, '%s_f%02d_%s' % (case, f, field))[1])
                     for f in (10, 50, 90)]
            panel([(NAME[ac], leg, tanks)], 'appx_contours_%s_%s' % (ac, tag))
    return 0


if __name__ == '__main__':
    sys.exit(main())
