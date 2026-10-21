# -*- coding: utf-8 -*-
"""
DIGITALIZED WORKFLOW DIAGRAM (October 2026)

Draws the flow of the workflow as one figure for Chapter 4 of the report:
the centre column is the sequence of stages, the boxes on the left are the
inputs of each stage, and the shaded boxes on the right are its outputs and
checks. The loop on the left is the finite element sizing of the wall.

    python draw_workflow.py [output folder]
"""
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Polygon

OUT = sys.argv[1] if len(sys.argv) > 1 else '.'
plt.rcParams.update({'font.family': 'serif', 'font.serif': ['Liberation Serif', 'TeX Gyre Termes', 'Times'], 'pdf.fonttype': 42,
                     'mathtext.fontset': 'stix', 'font.size': 7.5})

W, H = 6.3, 6.9
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W)
ax.set_ylim(0, H)
ax.axis('off')

XC, XL, XR = 3.15, 0.95, 5.35          # centre, left and right column centres
WC, WS, BH = 2.1, 1.6, 0.36            # box widths and height
GREY = '#e6e6e6'


def box(x, y, w, text, shade=False, bold=False):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - BH / 2), w, BH, boxstyle='square,pad=0',
                                fc=GREY if shade else 'white', ec='k', lw=0.6))
    ax.text(x, y, text, ha='center', va='center', fontsize=7.2,
            fontweight='bold' if bold else 'normal', linespacing=1.15)


def arrow(x0, y0, x1, y1):
    ax.annotate('', xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle='-|>', lw=0.6, color='k', shrinkA=0, shrinkB=0,
                                mutation_scale=7))


stages = [
    ('Requirements and design basis\nmission, CS-25 load factors, pressures, hold time', 'Reference aircraft\nkerosene mass, fuselage diameter', None),
    ('Tank geometry and insulation\n$R$, $L$, heat leak $Q$, MLI, jacket', 'Energy equivalence (LHV)\nmission fill level 0.90', 'Analytical sizing\nof insulation and jacket'),
    ('Parametric shell model\nelement size $0.5\\sqrt{Rt}$', 'Wall thickness $t$\nfrom the sizing loop', None),
    ('Free surface\nplane normal to $\\hat{\\mathbf{g}}$, shared', 'Load cases and fill levels\n4 load cases, 5 fill levels', None),
    ('Steady-state heat transfer\n(Abaqus)', 'Coefficients $h_{\\mathrm{L}}$, $h_{\\mathrm{V}}$\nflux $q_{\\mathrm{MLI}} = Q/A$', 'Temperature, heat flux\ncheck $\\Delta T_{\\mathrm{w}} = q_{\\mathrm{MLI}}/h_{\\mathrm{V}}$'),
    ('Linear static analysis\n(Abaqus)', 'Internal and hydrostatic pressure\ninertia, ring constraints, $T$', 'Far-field stress, reactions\ndisplacement, equilibrium'),
    ('Reserve factors\n$\\mathrm{RF}_{\\mathrm{y}}$ at 1.5 bar, $\\mathrm{RF}_{\\mathrm{b}}$ at 0.187 bar', None, None),
]
ys = [H - 0.35 - i * 0.72 for i in range(len(stages))]
for (c, l, r), y in zip(stages, ys):
    box(XC, y, WC, c, bold=False)
    if l:
        box(XL, y, WS, l)
        arrow(XL + WS / 2, y, XC - WC / 2, y)
    if r:
        box(XR, y, WS, r, shade=True)
        arrow(XC + WC / 2, y, XR - WS / 2, y)
for y0, y1 in zip(ys[:-1], ys[1:]):
    arrow(XC, y0 - BH / 2, XC, y1 + BH / 2)

# decision diamond
yd = ys[-1] - 0.78
dw, dh = 1.5, 0.5
ax.add_patch(Polygon([(XC, yd + dh / 2), (XC + dw / 2, yd), (XC, yd - dh / 2), (XC - dw / 2, yd)],
                     closed=True, fc='white', ec='k', lw=0.6))
ax.text(XC, yd, 'Lower RF\nbetween 1.00 and 1.06?', ha='center', va='center', fontsize=7.0)
arrow(XC, ys[-1] - BH / 2, XC, yd + dh / 2)
# loop back to the wall-thickness box
xl = 0.08
ax.plot([XC - dw / 2, xl, xl, XL - WS / 2], [yd, yd, ys[2], ys[2]], color='k', lw=0.6)
arrow(xl + 0.01, ys[2], XL - WS / 2, ys[2])
ax.text(XC - dw / 2 - 0.08, yd + 0.06, 'no: updated $t$', ha='right', va='bottom', fontsize=7.0, style='italic')

# production
yp = yd - 0.78
box(XC, yp, WC, 'Production and stability runs\n120 cases at the sized walls')
arrow(XC, yd - dh / 2, XC, yp + BH / 2)
ax.text(XC + 0.06, yd - dh / 2 - 0.05, 'yes', ha='left', va='top', fontsize=7.0, style='italic')
box(XL, yp, WS, 'Mesh convergence\nthree element sizes')
arrow(XL + WS / 2, yp, XC - WC / 2, yp)
box(XR, yp, WS, 'Verification\nanalytical solutions, benchmark', shade=True)
arrow(XC + WC / 2, yp, XR - WS / 2, yp)
yf = yp - 0.72
box(XC, yf, WC, 'Comparison of the three tanks\nidentical rules and extraction', shade=True)
arrow(XC, yp - BH / 2, XC, yf + BH / 2)

ax.set_ylim(yf - 0.3, H)
fig.set_size_inches(W, H - (yf - 0.3))
ax.set_position([0, 0, 1, 1])
for ext in ('pdf', 'png'):
    fig.savefig(os.path.join(OUT, 'fig4_workflow.' + ext), dpi=300)
print('written fig4_workflow')
