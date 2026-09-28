# -*- coding: utf-8 -*-
"""
Where the two cuts are made in a pressurised cylinder, and what each exposes.
Draft for Sarvesh to redraw. Output: design_drawings/01_pressure_vessel_cuts.{pdf,png}
"""
import os
import sys
import numpy as np

import os as _os
import sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
import repo_paths as RP  # noqa: E402  (folder layout of the repository)
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import thesis_figstyle as tfs                          # noqa: E402
import matplotlib.pyplot as plt                        # noqa: E402
from matplotlib.patches import Ellipse, Polygon, FancyArrowPatch   # noqa: E402

tfs.apply_style()
OUT = RP.FIGURES
WALL, CUT, P = 'black', 'black', 'black'
FILL = 'white'


def arr(ax, p0, p1, col, lw=1.0, ms=11):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle='-|>', mutation_scale=ms*0.8, color=col, lw=lw))


def cylinder(ax, x0, x1, r, e):
    """Side view of a cylinder along x with elliptical ends (e = end-ellipse half-width)."""
    ax.add_patch(Polygon([(x0, -r), (x1, -r), (x1, r), (x0, r)], closed=True, fc=FILL, ec='none'))
    ax.plot([x0, x1], [r, r], color=WALL, lw=1.0)
    ax.plot([x0, x1], [-r, -r], color=WALL, lw=1.0)
    ax.add_patch(Ellipse((x1, 0), 2 * e, 2 * r, fc=FILL, ec=WALL, lw=1.0))
    t = np.linspace(np.pi / 2, 3 * np.pi / 2, 60)
    ax.plot(x0 + e * np.cos(t), r * np.sin(t), color=WALL, lw=1.0)


fig, axs = plt.subplots(2, 2, figsize=(tfs.FIG_W, 3.9), gridspec_kw=dict(height_ratios=[1, 1.05]))
for a in axs.ravel():
    a.set_aspect('equal')
    a.axis('off')

# ---------------- (a) longitudinal cut: plane CONTAINS the axis ----------------
ax = axs[0, 0]
cylinder(ax, 0, 6, 1.4, 0.45)
ax.plot([-0.9, 7.3], [0, 0], color='k', lw=0.8, ls='-.')
arr(ax, (6.6, 0), (7.4, 0), 'k', lw=1.0, ms=9)
ax.text(7.45, 0.12, 'axis $x$', fontsize=7.5, va='bottom')
# cutting plane seen edge-on: horizontal plane through the axis
ax.add_patch(Polygon([(-0.6, -0.12), (6.6, -0.12), (6.6, 0.12), (-0.6, 0.12)], fc='none', hatch='////', ec=CUT, lw=0.8))
ax.text(3.0, -0.35, 'plane contains the axis', color=CUT, fontsize=7.5, ha='center', va='top')
ax.set_title('(a) Longitudinal cut', loc='left', fontsize=8.5)
ax.set_xlim(-1.2, 9.4)
ax.set_ylim(-2.0, 2.0)

# ---------------- (b) circumferential cut: plane PERPENDICULAR to the axis ----------------
ax = axs[0, 1]
cylinder(ax, 0, 6, 1.4, 0.45)
ax.plot([-0.9, 7.3], [0, 0], color='k', lw=0.8, ls='-.')
arr(ax, (6.6, 0), (7.4, 0), 'k', lw=1.0, ms=9)
ax.text(7.45, 0.12, 'axis $x$', fontsize=7.5, va='bottom')
ax.add_patch(Ellipse((3.0, 0), 0.9, 3.1, fc='white', hatch='////', ec=CUT, lw=0.8))
ax.text(3.0, -1.65, 'plane normal to the axis', color=CUT, fontsize=7.5, ha='center', va='top')
ax.set_title('(b) Circumferential cut', loc='left', fontsize=8.5)
ax.set_xlim(-1.2, 9.4)
ax.set_ylim(-2.0, 2.0)

# ---------------- (c) free body of the upper half: end view ----------------
ax = axs[1, 0]
R, t = 1.5, 0.16
th = np.linspace(0, np.pi, 120)
ax.fill_between(R * np.cos(th), R * np.sin(th), (R - t) * np.sin(th) * 0 + (R - t) * np.sin(th), color='none')
xo, yo = R * np.cos(th), R * np.sin(th)
xi, yi = (R - t) * np.cos(th), (R - t) * np.sin(th)
ax.add_patch(Polygon(np.r_[np.c_[xo, yo], np.c_[xi[::-1], yi[::-1]]], fc='white', ec=WALL, lw=0.9))
for a_ in np.linspace(0.35, np.pi - 0.35, 6):
    arr(ax, (0.55 * np.cos(a_), 0.55 * np.sin(a_)), (1.22 * np.cos(a_), 1.22 * np.sin(a_)), P, lw=0.8, ms=7)
ax.text(0, 0.28, '$p$', color=P, fontsize=9.5, ha='center')
for s in (-1, 1):
    arr(ax, (s * (R - t / 2), -0.05), (s * (R - t / 2), -0.85), WALL, lw=1.0)
    ax.text(s * (R - t / 2) + s * 0.25, -0.65, r'$\sigma_\theta$', fontsize=9.5, ha='center', color=WALL)
ax.text(0, -1.25, r'top half, end view: $p\,(2R)L = 2\,\sigma_\theta\,tL$', ha='center', fontsize=7.5)
ax.text(0, -1.62, r'$\sigma_\theta = pR/t$ (hoop)', ha='center', fontsize=8, fontweight='bold')
ax.set_xlim(-2.4, 2.4)
ax.set_ylim(-1.9, 1.8)

# ---------------- (d) free body of the forward part: side view ----------------
ax = axs[1, 1]
ax.add_patch(Polygon([(0, -1.0), (2.6, -1.0), (2.6, 1.0), (0, 1.0)], fc=FILL, ec='none'))
ax.plot([0, 2.6], [1.0, 1.0], color=WALL, lw=1.0)
ax.plot([0, 2.6], [-1.0, -1.0], color=WALL, lw=1.0)
t_ = np.linspace(-np.pi / 2, np.pi / 2, 60)
ax.plot(2.6 + 0.6 * np.cos(t_), np.sin(t_), color=WALL, lw=1.0)
ax.plot([0, 0], [-1.0, 1.0], color=CUT, lw=0.9, ls='--')
for y in (-1.0, 1.0):
    arr(ax, (0, y), (-0.85, y), WALL, lw=1.0)
ax.text(-0.95, 0.0, r'$\sigma_x$', fontsize=9.5, ha='right', va='center', color=WALL)
for y in (-0.5, 0.0, 0.5):
    arr(ax, (1.7, y), (2.9, y), P, lw=0.8, ms=7)
ax.text(1.4, 0.0, '$p$', color=P, fontsize=9.5, ha='right', va='center')
ax.text(1.0, -1.45, r'forward part, side view: $p\,\pi R^2 = \sigma_x\,2\pi R t$', ha='center', fontsize=7.5)
ax.text(1.0, -1.82, r'$\sigma_x = pR/2t$ (axial)', ha='center', fontsize=8, fontweight='bold')
ax.set_xlim(-2.0, 4.2)
ax.set_ylim(-2.1, 1.7)

fig.tight_layout(h_pad=0.4, w_pad=0.6)
tfs.save(fig, os.path.join(OUT, '01_pressure_vessel_cuts'))
