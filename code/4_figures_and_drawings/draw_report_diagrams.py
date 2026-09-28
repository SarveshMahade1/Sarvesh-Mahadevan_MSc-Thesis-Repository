# -*- coding: utf-8 -*-
"""
THE THREE DRAWN DIAGRAMS KEPT IN THE REPORT (28 Sep 2026)
=========================================================

Drawn at their printed width (0.8 of the text width) in the report typeface,
with every number taken from the results files:

  fig3_resistance_network  the wall stack as thermal resistances in series,
                           E190, with the share of each layer computed here
  fig3_wall_stack          the E190 wall stack at true relative thickness
  fig4_free_surface        the E190 free surface at phi_f = 0.90, vertical
                           and LC7 resultant, from the same volume-sampling
                           rule as the model

    python python_scripts/draw_report_diagrams.py
"""
import json
import math
import os
import sys

import os as _os
import sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import repo_paths as RP  # noqa: E402  (folder layout of the repository)
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np                                      # noqa: E402
import thesis_figstyle as tfs                           # noqa: E402
import matplotlib.pyplot as plt                         # noqa: E402
from matplotlib.patches import Rectangle, Circle, Wedge, PathPatch, FancyArrowPatch  # noqa: E402
from matplotlib.path import Path                        # noqa: E402
import boiloff_holdtime_sensitivity as bh               # noqa: E402

tfs.apply_style()
FIG = RP.FIGURES
D = json.load(open(os.path.join(RP.PROCESSED, 'final_numbers.json')))
I = json.load(open(os.path.join(RP.PROCESSED, 'insulation.json')))
W = tfs.FIG_W
# black and white (28 Sep 2026): grey levels only
C_WALL, C_MLI, C_JACK, C_FILM = '#404040', '#f2f2f2', '#a6a6a6', '#d9d9d9'
C_LIQ = '#e6e6e6'


def bare(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


# ---------------------------------------------------------------------------
def resistance_network():
    e = D['E190']
    R, Lt = e['R_m'], e['L_total_m']
    t_w = e['t_wall_mm'] / 1e3
    t_mli = I['E190']['mli_mm'] / 1e3
    t_j = 3.45e-3
    r0, r1 = R - t_w / 2, R + t_w / 2
    r2 = r1 + t_mli
    r3 = r2 + t_j
    Rw = math.log(r1 / r0) / (2 * math.pi * bh.K_INNER * Lt)
    Rm = math.log(r2 / r1) / (2 * math.pi * bh.K_EFF_MLI * Lt)
    Rj = math.log(r3 / r2) / (2 * math.pi * bh.K_OUTER * Lt)
    Rf = 1.0 / (bh.H_CONV * 2 * math.pi * r3 * Lt)
    tot = Rw + Rm + Rj + Rf
    share = [100 * x / tot for x in (Rw, Rm, Rj, Rf)]

    fig, ax = plt.subplots(figsize=(W, 1.8))
    bare(ax)
    names = ['inner\nwall', 'vacuum\nand MLI', 'outer\njacket', 'outer\nfilm']
    cols = [C_WALL, C_MLI, C_JACK, C_FILM]
    tcol = ['white', 'black', 'black', 'black']
    x0, bw, y0, bh_ = 0.9, 1.0, 1.35, 0.45
    for k in range(4):
        xx = x0 + k * bw
        ax.add_patch(Rectangle((xx, y0), bw, bh_, fc=cols[k], ec='black', lw=0.8))
        ax.text(xx + bw / 2, y0 + bh_ / 2, names[k], ha='center', va='center', color=tcol[k], fontsize=8, linespacing=1.0)
        # zigzag resistor
        zx = np.linspace(xx + 0.2, xx + bw - 0.2, 13)
        zy = 0.55 + 0.09 * np.array([0 if i in (0, 12) else (1 if i % 2 else -1) for i in range(13)])
        ax.plot(zx, zy, color='black', lw=0.9)
        ax.plot([xx, xx + 0.2], [0.55, 0.55], color='black', lw=0.9)
        ax.plot([xx + bw - 0.2, xx + bw], [0.55, 0.55], color='black', lw=0.9)
        ax.plot([xx + bw / 2] * 2, [y0, 0.72], color='0.5', lw=0.5, ls=':')
        lab = ['$R_{\\mathrm{wall}}$', '$R_{\\mathrm{MLI}}$', '$R_{\\mathrm{jacket}}$', '$R_{\\mathrm{film}}$'][k]
        ax.text(xx + bw / 2, 0.80, lab, ha='center', va='bottom', fontsize=8.5)
        s = share[k]
        txt = ('%.1f %%' % s) if s > 1 else ('%.1e %%' % s).replace('e-0', r'$\times10^{-') + '}$' if s < 0.01 else ('%.2f %%' % s)
        ax.text(xx + bw / 2, 0.36, ('%.1f %%' % s) if s >= 0.1 else '< 0.01 %', ha='center', va='top', fontsize=8.5)
    ax.plot([0.35, x0], [0.55, 0.55], color='black', lw=0.9)
    ax.plot([x0 + 4 * bw, x0 + 4 * bw + 0.55], [0.55, 0.55], color='black', lw=0.9)
    ax.text(0.30, 0.55, '$T_{\\mathrm{cold}}$ = 20 K', ha='right', va='center', fontsize=8.5)
    ax.text(x0 + 4 * bw + 0.6, 0.55, '$T_{\\mathrm{warm}}$ = 293 K', ha='left', va='center', fontsize=8.5)
    ax.text(x0 + 2 * bw, -0.12, r'$Q = (T_{\mathrm{warm}} - T_{\mathrm{cold}})\,/\,\sum R_i$ = %.1f W; share of $\sum R_i$ below each resistor' % I['E190']['Q_W'],
            ha='center', va='bottom', fontsize=8.5)
    ax.set_xlim(-0.6, x0 + 4 * bw + 1.7)
    ax.set_ylim(-0.15, y0 + bh_ + 0.05)
    fig.tight_layout(pad=0.1)
    tfs.save(fig, os.path.join(FIG, 'fig3_resistance_network'))
    return share


# ---------------------------------------------------------------------------
def wall_stack():
    e = D['E190']
    tw, tm, tj = e['t_wall_mm'], I['E190']['mli_mm'], 3.45
    fig = plt.figure(figsize=(W, 2.0))
    # (a) section, gap exaggerated
    ax = fig.add_axes([0.0, 0.02, 0.30, 0.84])
    bare(ax)
    ax.set_aspect('equal')
    ax.add_patch(Circle((0, 0), 1.0, fc=C_JACK, ec='black', lw=0.8))
    ax.add_patch(Circle((0, 0), 0.93, fc=C_MLI, ec='black', lw=0.6))
    ax.add_patch(Circle((0, 0), 0.80, fc=C_WALL, ec='black', lw=0.6))
    ax.add_patch(Circle((0, 0), 0.74, fc=C_LIQ, ec='black', lw=0.6))
    ax.text(0, 0.08, 'LH$_2$', ha='center', va='center', fontsize=8.5)
    ax.text(0, -0.18, '20 K', ha='center', va='center', fontsize=8.5)
    ax.add_patch(Wedge((0, 0), 1.04, -12, 12, width=0.34, fc='none', ec='black', lw=1.0))
    ax.set_xlim(-1.1, 1.25)
    ax.set_ylim(-1.1, 1.1)
    ax.set_title('(a) Section, gap exaggerated', loc='left', fontsize=9)
    # (b) true relative thickness
    ax = fig.add_axes([0.36, 0.02, 0.64, 0.84])
    bare(ax)
    tot = tw + tm + tj
    x = 0.0
    for thk, col, nm in ((tw, C_WALL, 'inner wall, Al 2219-T87'), (tm, C_MLI, 'vacuum and MLI'), (tj, C_JACK, 'outer jacket')):
        ax.add_patch(Rectangle((x, 0.3), thk, 0.5, fc=col, ec='black', lw=0.7))
        ax.annotate('', xy=(x, 0.18), xytext=(x + thk, 0.18),
                    arrowprops=dict(arrowstyle='<->', lw=0.6, shrinkA=0, shrinkB=0))
        ax.text(x + thk / 2, 0.12, '%.2f mm' % thk if thk < 2 else '%.2f mm' % thk, ha='center', va='top', fontsize=8.5)
        if nm.startswith('inner'):
            ax.annotate(nm, xy=(x + thk / 2, 0.8), xytext=(x - 1.0, 1.02), ha='left', va='bottom', fontsize=8.5,
                        arrowprops=dict(arrowstyle='-', lw=0.5))
        else:
            ax.text(x + thk / 2, 0.86, nm, ha='center', va='bottom', fontsize=8.5)
        x += thk
    ax.text(-0.4, 0.55, '20 K', ha='right', va='center', fontsize=8.5)
    ax.text(tot + 0.4, 0.55, '293 K', ha='left', va='center', fontsize=8.5)
    ax.set_xlim(-2.8, tot + 3.0)
    ax.set_ylim(-0.1, 1.28)
    ax.set_title('(b) E190 wall stack at true relative thickness', loc='left', fontsize=9)
    tfs.save(fig, os.path.join(FIG, 'fig3_wall_stack'))


# ---------------------------------------------------------------------------
def free_surface():
    e = D['E190']
    R, L, h = e['R_m'], e['L_barrel_m'], e['h_cap_m']

    def outline():
        th = np.linspace(-np.pi / 2, np.pi / 2, 90)
        aft = [(-h * np.cos(a), R * np.sin(a)) for a in th[::-1]]
        fwd = [(L + h * np.cos(a), R * np.sin(a)) for a in th]
        pts = [(-h * 0, -R)] + [(L, -R)] + fwd[1:] + [(0, R)] + aft[1:]
        # build closed contour: bottom line, fwd cap, top line, aft cap
        bottom = [(0, -R), (L, -R)]
        fwdc = [(L + h * np.cos(a), R * np.sin(a)) for a in np.linspace(-np.pi / 2, np.pi / 2, 90)]
        top = [(L, R), (0, R)]
        aftc = [(-h * np.cos(a), R * np.sin(a)) for a in np.linspace(np.pi / 2, 3 * np.pi / 2, 90)]
        aftc = [(-h * abs(np.cos(a)), R * np.sin(a)) for a in np.linspace(np.pi / 2, -np.pi / 2, 90)]
        return bottom + fwdc + top + aftc

    def s_surf(gx, gy, gz, fill, n=60):
        pts = []
        xs = np.linspace(-h, L + h, 3 * n)
        ys = np.linspace(-R, R, n)
        for x in xs:
            rl = R * math.sqrt(max(0.0, 1 - (x / h) ** 2)) if x < 0 else (
                R * math.sqrt(max(0.0, 1 - ((x - L) / h) ** 2)) if x > L else R)
            for y in ys:
                for z in ys:
                    if y * y + z * z <= rl * rl:
                        pts.append(x * gx + y * gy + z * gz)
        pts.sort()
        return pts[int(round((1 - fill) * len(pts)))]

    fig, axs = plt.subplots(2, 1, figsize=(W, 2.9))
    poly = outline()
    for ax, lc, key in ((axs[0], (0.0, 0.0, 1.0), 'Baseline|0.90'), (axs[1], (9.0, 3.0, 6.0), 'LC7_Emergency|0.90')):
        nx, ny, nz = lc
        m = math.sqrt(nx * nx + ny * ny + nz * nz)
        g = (nx / m, -nz / m, -ny / m)
        ss = s_surf(g[0], g[1], g[2], 0.9)
        path = Path(poly + [poly[0]], closed=True)
        ax.add_patch(PathPatch(path, fc='white', ec='black', lw=1.0, zorder=3, fill=False))
        # liquid in the vertical mid-plane z = 0: x gx + y gy >= s_surf
        X, Y = np.meshgrid(np.linspace(-h, L + h, 800), np.linspace(-R, R, 200))
        inside = path.contains_points(np.c_[X.ravel(), Y.ravel()]).reshape(X.shape)
        wet = (X * g[0] + Y * g[1] >= ss) & inside
        ax.contourf(X, Y, wet.astype(float), levels=[0.5, 1.5], colors=[C_LIQ], zorder=1)
        ax.contour(X, Y, (X * g[0] + Y * g[1] - ss) * inside, levels=[0.0], colors=['black'], linewidths=1.0, zorder=2)
        hd = D['E190']['heads'][key]
        cx, cy = 0.62 * L, 0.55
        ax.add_patch(FancyArrowPatch((cx, cy), (cx + 1.2 * g[0], cy + 1.2 * g[1]), arrowstyle='-|>',
                                     mutation_scale=9, color='black', lw=1.3, zorder=6))
        ax.text(cx + 1.2 * g[0] + 0.3, cy + 1.2 * g[1] + 0.1, r'$\hat{\mathbf{g}}$', color='black', fontsize=9,
                va='center', zorder=6)
        ax.set_aspect('equal')
        ax.set_xlim(-h - 0.3, L + h + 0.3)
        ax.set_ylim(-R - 0.25, R + 0.25)
        bare(ax)
        if key.startswith('Baseline'):
            ax.set_title('(a) Vertical resultant, 1.0 g: depth %.2f m, head %.1f %% of $p_{\\mathrm{gauge}}$'
                         % (hd['depth_m'], 100 * hd['head_over_p']), loc='left', fontsize=9)
        else:
            ax.set_title('(b) LC7 resultant, 11.22 g tilted %.1f$^\\circ$: depth %.2f m, head %.1f %% of $p_{\\mathrm{gauge}}$'
                         % (hd['tilt_deg'], hd['depth_m'], 100 * hd['head_over_p']), loc='left', fontsize=9)
    fig.tight_layout(pad=0.2, h_pad=0.6)
    tfs.save(fig, os.path.join(FIG, 'fig4_free_surface'))


if __name__ == '__main__':
    sh = resistance_network()
    print('resistance shares wall, MLI, jacket, film [%]:', ['%.3g' % s for s in sh])
    wall_stack()
    free_surface()
