# -*- coding: utf-8 -*-
"""
TANK DESIGN DETAIL DRAWINGS, BLACK AND WHITE (28 Sep 2026)
Drafts for Sarvesh to redraw. E190 values, from results/processed/final_numbers.json,
results/processed/insulation.json and the jacket trade (Chapter 3).

    12_tank_general_arrangement   side elevation of the E190 tank: inner vessel,
                                  jacket, the three support rings and main dimensions
    13_tank_details               (a) wall build-up at a support ring (detail A)
                                  (b) semi-ellipsoidal end cap geometry

Output: design_drawings/12_*.{pdf,png}, 13_*.{pdf,png}
"""
import json
import os
import sys
import numpy as np

import os as _os
import sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import repo_paths as RP  # noqa: E402  (folder layout of the repository)
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import thesis_figstyle as tfs                              # noqa: E402
import matplotlib.pyplot as plt                            # noqa: E402
from matplotlib.patches import Polygon, Rectangle, Circle, FancyArrowPatch   # noqa: E402

tfs.apply_style()
OUT = RP.DRAWINGS
os.makedirs(OUT, exist_ok=True)
FW = 455.24 / 72.27                                        # full text width, inches
FN = json.load(open(os.path.join(RP.PROCESSED, 'final_numbers.json')))
INS = json.load(open(os.path.join(RP.PROCESSED, 'insulation.json')))
E = FN['E190']
R, L, H = E['R_m'], E['L_barrel_m'], E['h_cap_m']
T_MM, MLI_MM, JKT_MM = E['t_wall_mm'], INS['E190']['mli_mm'], 3.45
K = 'black'


def dim(ax, p0, p1, text, off=(0, 0), fs=8, rot=0):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle='<|-|>', mutation_scale=7, color=K, lw=0.7))
    ax.text((p0[0] + p1[0]) / 2 + off[0], (p0[1] + p1[1]) / 2 + off[1], text,
            ha='center', va='center', fontsize=fs, rotation=rot,
            bbox=dict(fc='white', ec='none', pad=0.6))


def outline(x0, r, h, n=80):
    """Meridian of barrel + two semi-ellipsoidal caps, returned as a closed polygon."""
    th = np.linspace(-np.pi / 2, np.pi / 2, n)
    fwd = np.c_[x0[1] + h * np.cos(th), r * np.sin(th)]
    aft = np.c_[x0[0] - h * np.cos(th[::-1]), r * np.sin(th[::-1])]
    return np.r_[fwd, aft]


def save(fig, name):
    tfs.save(fig, os.path.join(OUT, name))


# =========================================================================
# 12  GENERAL ARRANGEMENT
# =========================================================================
fig, ax = plt.subplots(figsize=(FW, 3.6))
ax.set_aspect('equal')
ax.axis('off')
G = 0.16                                   # drawn gap, exaggerated
Rj = R + G
# jacket
ax.add_patch(Polygon(outline((0, L), Rj, H + G), closed=True, fc='white', ec=K, lw=1.3))
# jacket stiffening rings (every 0.25 m on the barrel)
for x in np.arange(0.25, L, 0.25):
    ax.plot([x, x], [Rj, Rj + 0.07], color=K, lw=0.45)
    ax.plot([x, x], [-Rj, -Rj - 0.07], color=K, lw=0.45)
# gap with MLI (light hatch)
ax.add_patch(Polygon(outline((0, L), Rj - 0.01, H + G - 0.01), closed=True, fc='none', ec='none', hatch='.....', lw=0))
# inner vessel
ax.add_patch(Polygon(outline((0, L), R, H), closed=True, fc='white', ec=K, lw=1.6))
# liquid at 90 % fill, level surface (Baseline): depth 2.19 m from the bottom
depth = FN['E190']['heads']['Baseline|0.90']['depth_m']
ys = -R + depth
ax.plot([-H * np.sqrt(max(0, 1 - (ys / R) ** 2)), L + H * np.sqrt(max(0, 1 - (ys / R) ** 2))], [ys, ys],
        color=K, lw=0.8, ls='--')
ax.text(L * 0.30, -0.35, r'LH$_2$, 20.3 K', ha='center', fontsize=8.5)
ax.text(L * 0.30, (ys + R) / 2, 'ullage', ha='center', va='center', fontsize=7.5)
# support rings: section lines through the vessel and symbols below
for x, name, anchor in ((0.0, 'aft ring\nsliding', False), (L / 2, 'mid ring\naxial anchor', True),
                        (L, 'forward ring\nsliding', False)):
    ax.plot([x, x], [-Rj - 0.12, Rj + 0.12], color=K, lw=1.1, ls='-.')
    yb = -Rj - 0.12
    ax.add_patch(Polygon([(x, yb), (x - 0.28, yb - 0.38), (x + 0.28, yb - 0.38)], closed=True, fc='white', ec=K, lw=0.9))
    if anchor:
        ax.add_patch(Rectangle((x - 0.42, yb - 0.50), 0.84, 0.12, fc='none', ec='none', hatch='////'))
        ax.plot([x - 0.42, x + 0.42], [yb - 0.38, yb - 0.38], color=K, lw=0.9)
    else:
        for dx in (-0.14, 0.14):
            ax.add_patch(Circle((x + dx, yb - 0.45), 0.07, fc='white', ec=K, lw=0.8))
        ax.plot([x - 0.42, x + 0.42], [yb - 0.52, yb - 0.52], color=K, lw=0.9)
    ax.text(x, yb - 0.68, name, ha='center', va='top', fontsize=7.5)
# axis
ax.plot([-H - G - 0.6, L + H + G + 0.9], [0, 0], color=K, lw=0.6, ls=(0, (8, 3, 1, 3)))
ax.add_patch(FancyArrowPatch((L + H + G + 0.3, 0), (L + H + G + 1.1, 0), arrowstyle='-|>', mutation_scale=8, color=K, lw=0.8))
ax.text(L + H + G + 1.15, 0.08, 'x, forward', fontsize=7.5, va='bottom')
# dimensions
yt = Rj + 0.45
dim(ax, (0, yt), (L, yt), 'barrel  L = %.2f m' % L, fs=8)
dim(ax, (-H, yt), (0, yt), '', fs=7.5)
dim(ax, (L, yt), (L + H, yt), '', fs=7.5)
ax.text(-H - 0.1, yt, 'h = %.2f m' % H, ha='right', va='center', fontsize=7.5)
ax.text(L + H + 0.1, yt, 'h = %.2f m' % H, ha='left', va='center', fontsize=7.5)
yt2 = Rj + 1.05
dim(ax, (-H, yt2), (L + H, yt2), 'inner vessel  %.2f m' % (L + 2 * H), fs=8)
for x in (-H, L + H):
    ax.plot([x, x], [0.2, yt2 + 0.1], color=K, lw=0.4)
for x in (0, L):
    ax.plot([x, x], [Rj + 0.2, yt + 0.1], color=K, lw=0.4)
xd = L * 0.82
dim(ax, (xd, -R), (xd, R), '2R = %.2f m' % (2 * R), off=(0.62, 0), fs=7.5, rot=90)
# labels with leaders
def lead(xy, xyt, text, ha='left'):
    ax.annotate(text, xy=xy, xytext=xyt, fontsize=7.5, ha=ha, va='center',
                arrowprops=dict(arrowstyle='-', color=K, lw=0.5))
lead((L * 0.25, -R), (L * 0.25, -Rj - 1.05), 'inner vessel, Al 2219-T87\nt = %.3f mm' % T_MM, ha='center')
lead((L * 0.75, -Rj - 0.07), (L * 0.75, -Rj - 1.05), 'outer jacket %.2f mm\nstiffening rings at 0.25 m' % JKT_MM, ha='center')
lead((-H - 0.07, 0.55), (-H - 1.2, 1.1), 'vacuum gap with\nMLI %.1f mm' % MLI_MM, ha='right')
lead((-H, -0.4), (-H - 1.2, -1.0), 'semi-ellipsoidal\ncap, a/b = 1.6', ha='right')
# detail A marker at the mid ring, top
ax.add_patch(Circle((L / 2, R + G / 2), 0.42, fc='none', ec=K, lw=0.9))
ax.text(L / 2 + 0.5, R - 0.25, 'A', fontsize=9, fontweight='bold', va='center')
ax.set_xlim(-H - 2.6, L + H + 2.3)
ax.set_ylim(-Rj - 1.75, yt2 + 0.35)
fig.tight_layout(pad=0.2)
save(fig, '12_tank_general_arrangement')

# =========================================================================
# 13  DETAILS: (a) wall at a support ring, (b) end cap geometry
# =========================================================================
fig, (a1, a2) = plt.subplots(1, 2, figsize=(FW, 3.0), gridspec_kw=dict(width_ratios=[1.15, 1]))
for a in (a1, a2):
    a.set_aspect('equal')
    a.axis('off')

# (a) radial build-up, not to scale; x is axial, y is radial outward
yw0, yw1 = 0.0, 0.25          # inner wall
yg1 = 1.55                    # gap top
yj1 = 1.85                    # jacket top
X0, X1 = -2.2, 2.2
a1.add_patch(Rectangle((X0, -0.9), X1 - X0, 0.9, fc='none', ec='none', hatch='---'))
a1.text(X0 + 0.15, -0.45, r'LH$_2$ 20.3 K', fontsize=8, va='center', bbox=dict(fc='white', ec='none', pad=1))
a1.add_patch(Rectangle((X0, yw0), X1 - X0, yw1 - yw0, fc='0.35', ec=K, lw=0.9))
a1.add_patch(Rectangle((X0, yw1), X1 - X0, yg1 - yw1, fc='white', ec='none', hatch='.....'))
a1.add_patch(Rectangle((X0, yg1), X1 - X0, yj1 - yg1, fc='0.75', ec=K, lw=0.9))
for xr in (-1.6, 1.6):        # jacket stiffening rings
    a1.add_patch(Rectangle((xr - 0.07, yj1), 0.14, 0.35, fc='0.75', ec=K, lw=0.8))
# support ring across the gap, idealised
a1.add_patch(Rectangle((-0.45, yw1), 0.9, yg1 - yw1, fc='white', ec=K, lw=0.9, ls='--'))
a1.plot([-0.45, 0.45], [yw1, yw1], color=K, lw=2.0)
a1.text(0.0, (yw1 + yg1) / 2 + 0.2, 'support\nring', ha='center', va='center', fontsize=7,
        bbox=dict(fc='white', ec='none', pad=0.5))
a1.plot([0, 0], [-0.9, yj1 + 0.55], color=K, lw=0.7, ls='-.')
a1.text(0, yj1 + 0.6, 'ring station', ha='center', va='bottom', fontsize=7.5)
# labels right
def rl(y, text):
    a1.annotate(text, xy=(X1 - 0.05, y), xytext=(X1 + 0.35, y), fontsize=7.5, va='center',
                arrowprops=dict(arrowstyle='-', color=K, lw=0.5))
rl((yw0 + yw1) / 2, 'inner wall, %.3f mm' % T_MM)
rl((yw1 + yg1) / 2, 'vacuum gap, MLI %.1f mm' % MLI_MM)
rl((yg1 + yj1) / 2, 'outer jacket, %.2f mm' % JKT_MM)
rl(yj1 + 0.2, 'jacket stiffening ring')
a1.text(X0 + 0.15, yj1 + 0.55, 'bay, 293 K', fontsize=8)
a1.add_patch(FancyArrowPatch((-1.0, yj1 + 0.45), (-1.0, yw1 + 0.05), arrowstyle='-|>', mutation_scale=8, color=K, lw=0.8))
a1.text(-1.1, (yj1 + yw1) / 2 - 0.1, r'$q_{\mathrm{MLI}}$ =' + '\n' + r'%.2f W/m$^2$' % E['q_Wm2'], fontsize=7, va='center', ha='right',
        bbox=dict(fc='white', ec='none', pad=0.5))
a1.text(X0, -1.15, '(a) Detail A: wall at a support ring, not to scale', fontsize=8.5, va='top')
a1.set_xlim(X0 - 0.1, X1 + 3.4)
a1.set_ylim(-1.5, yj1 + 1.0)

# (b) end cap geometry, to scale (E190)
th = np.linspace(0, np.pi / 2, 100)
a2.plot(H * np.cos(th), R * np.sin(th), color=K, lw=1.6)
a2.plot(H * np.cos(th), -R * np.sin(th), color=K, lw=1.6)
a2.plot([0, -1.2], [R, R], color=K, lw=1.6)
a2.plot([0, -1.2], [-R, -R], color=K, lw=1.6)
a2.plot([0, 0], [-R - 0.25, R + 0.25], color=K, lw=0.9, ls='-.')
a2.text(0.0, R + 0.3, 'cap equator = cap-to-barrel\njunction = forward ring station', ha='center', va='bottom', fontsize=7)
a2.plot([-1.3, H + 0.35], [0, 0], color=K, lw=0.6, ls=(0, (8, 3, 1, 3)))
dim(a2, (0.0, -R - 0.3), (H, -R - 0.3), '', fs=7)
a2.plot([H, H], [-0.1, -R - 0.4], color=K, lw=0.4)
a2.text(H / 2, -R - 0.42, 'b = h = R/1.6 = %.3f m' % H, ha='center', va='top', fontsize=7)
dim(a2, (-0.55, 0), (-0.55, R), 'a = R\n= %.3f m' % R, off=(-0.33, 0), fs=7)
a2.plot(H, 0, 'o', color=K, ms=3)
a2.text(H + 0.06, 0.12, 'pole', fontsize=7.5)
a2.text(-0.6, -R + 0.25, 'barrel', fontsize=7.5)
a2.text(-1.3, -R - 0.8, '(b) Semi-ellipsoidal end cap, a/b = 1.6, to scale', fontsize=8.5, va='top')
a2.set_xlim(-1.5, H + 0.7)
a2.set_ylim(-R - 1.15, R + 0.8)
fig.tight_layout(pad=0.3, w_pad=0.8)
save(fig, '13_tank_details')
print('written 12 and 13 to', OUT)
