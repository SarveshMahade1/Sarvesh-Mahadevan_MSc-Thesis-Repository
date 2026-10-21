"""
================================================================================
SHARED FIGURE STYLE AND OUTPUT — TU Delft MSc thesis
Author: Sarvesh Mahadevan

Import this from every plotting script so all figures share one look and one
output path. Two things it fixes.

--------------------------------------------------------------------------------
1. VECTOR OUTPUT
--------------------------------------------------------------------------------
Every figure was previously written as a 200 dpi PNG. A raster figure is fixed
at the resolution it was written: text edges are soft on screen, lines alias
when the figure is scaled, and zooming in on a printed PDF shows pixels.

These figures are line art — axes, curves, markers, text. Line art belongs in a
vector format, where there is no resolution at all: the curve is stored as a
curve and the renderer draws it at whatever resolution the output device has.
pdflatex embeds PDF natively, so this is also the cheaper path.

save() therefore writes a PDF as the primary artefact and a 400 dpi PNG
alongside it for previewing. In LaTeX, reference the figure WITHOUT an
extension:

    \\includegraphics[width=\\textwidth]{36_scale_effects}

graphicx then prefers the PDF automatically.

--------------------------------------------------------------------------------
2. TEXT BELONGS IN THE CAPTION, NOT IN THE IMAGE
--------------------------------------------------------------------------------
Each script previously drew a bold suptitle and a grey explanatory paragraph
inside the figure. That text is already in the LaTeX caption, so it appeared
twice on the page, and it crowded the panels that carry the actual result.

A figure should carry labelling only: axis labels, units, tick values, legend,
and any annotation that points at a specific feature. Everything else is the
caption's job. apply_style() sets type sizes on the assumption that the figure
is scaled to roughly text width, so the smallest text on the page is still
legible at print size.
================================================================================
"""

import os

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

NAVY, RED, GREEN, GREY, BLUE = '#1F3864', '#C00000', '#2E7D32', '#595959', '#2E75B6'

# Printed width of a figure included at 0.8\\textwidth (textwidth 455.24 pt).
FIG_W = 0.8 * 455.24 / 72.27

# One colour per aircraft, used in every figure of the report.
AC_COLOR = {'E190': '#1F3864', 'A320': '#2E7D32', 'A350': '#C00000'}
AC_LABEL = {'E190': 'E190', 'A320': 'A320', 'A350': 'A350-900'}
# Load cases are named by their physical condition (October 2026); the
# internal keys keep the original identifiers of the result files.
LC_LABEL = {'Baseline': 'Baseline condition', 'LC1_Maneuver': 'Symmetric manoeuvre',
            'LC3_Combined': 'Combined manoeuvre', 'LC7_Emergency': 'Emergency landing'}


def apply_style():
    """Type sizes and line weights for a figure reproduced at text width."""
    plt.rcParams.update({
        # --- type ---------------------------------------------------------
        # 28 Sep 2026: the report is set in Times (newtx), so the figures use
        # TeX Gyre Termes (a Times design) with STIX mathematics, and sizes
        # chosen for a figure drawn at its printed width (FIG_W below), so
        # the type in a figure prints at the same size as the caption text.
        'font.family': 'serif',
        'font.serif': ['Liberation Serif', 'TeX Gyre Termes', 'Times'],   # Times-metric TrueType: embeds cleanly in PDF
        'mathtext.fontset': 'stix',
        'font.size': 9.5,
        'axes.labelsize': 9.5,
        'axes.titlesize': 9.5,
        'xtick.labelsize': 8.5,
        'ytick.labelsize': 8.5,
        'legend.fontsize': 8.5,

        # --- line art -----------------------------------------------------
        'axes.linewidth': 0.9,
        'lines.linewidth': 1.4,
        'lines.markersize': 4.5,
        'grid.linewidth': 0.6,
        'grid.alpha': 0.30,
        'xtick.major.width': 0.9,
        'ytick.major.width': 0.9,

        # --- frame: drop the top and right spines -------------------------
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.edgecolor': '#333333',

        # --- vector output hygiene ----------------------------------------
        # Keep text as text in the PDF rather than outlines, so it stays
        # selectable, searchable and crisp at any zoom.
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
        'savefig.transparent': False,
        'figure.facecolor': 'white',
        'savefig.facecolor': 'white',
    })


def save(fig, path, png_dpi=400):
    """
    Write `fig` as PDF (primary) and PNG (preview) next to each other.

    `path` may carry any extension or none; it is stripped and both are
    written. bbox_inches='tight' crops the whitespace left behind once the
    suptitle and the explanatory paragraph are removed, so no per-figure
    retuning of subplots_adjust(top=...) is needed.
    """
    stem = os.path.splitext(path)[0]
    for ext, kw in (('.pdf', {}), ('.png', {'dpi': png_dpi})):
        fig.savefig(stem + ext, bbox_inches='tight', pad_inches=0.04, **kw)
    print('written %s.pdf + %s.png' % (stem, os.path.basename(stem)))
    plt.close(fig)


# ---------------------------------------------------------------------------
# Panel conventions (October 2026). Every result graph of the report follows
# the same rules: one quantity and one aircraft (or load case) per panel, an
# individual y-scale per panel so that small changes remain visible, the value
# written next to every plotted point, and a note in the upper left corner
# giving the change over the plotted range and, where one exists, the value
# predicted by the governing analytical relation. Panels with different scales
# are therefore compared by their numbers, not by their slopes.
# ---------------------------------------------------------------------------
def label_points(ax, x, y, col, fmt, side='auto', fs=7, skip=()):
    """Write the value of every point next to it. side='auto' places the label
    below and to the right of a point on a rising segment, above a point on a
    falling segment, and above and to the left of the last point; 'above' and
    'below' force the position."""
    n = len(x)
    sides = side if isinstance(side, (list, tuple)) else [side] * n
    for k, (xi, yi) in enumerate(zip(x, y)):
        if k in skip:
            continue
        side = sides[k]
        if side == 'peak':      # above a local maximum, below a local minimum
            nb = [y[j] for j in (k - 1, k + 1) if 0 <= j < n]
            side = 'above' if yi >= sum(nb) / len(nb) else 'below'
        if side == 'above':
            off, ha, va = (0, 5), 'center', 'bottom'
        elif side == 'below':
            off, ha, va = (0, -5), 'center', 'top'
        elif side == 'ur':
            off, ha, va = (4, 3), 'left', 'bottom'
        elif side == 'lr':
            off, ha, va = (4, -4), 'left', 'top'
        elif side == 'ul':
            off, ha, va = (-4, 3), 'right', 'bottom'
        elif side == 'll':
            off, ha, va = (-4, -3), 'right', 'top'
        elif side == 'auto':
            # keep the label on the side of the point away from both
            # adjoining segments of the curve
            tol = 1e-9 + 0.02 * (max(y) - min(y))
            if k < n - 1:
                d = y[k + 1] - yi
                if abs(d) <= tol:
                    off, ha, va = (0, -5), 'center', 'top'          # flat: below
                elif d < 0:
                    off, ha, va = (4, 3), 'left', 'bottom'          # falling: above right
                else:
                    off, ha, va = (4, -4), 'left', 'top'            # rising: below right
            else:
                if y[k - 1] > yi + tol:
                    off, ha, va = (-4, -3), 'right', 'top'          # came down: below left
                else:
                    off, ha, va = (-5, 4), 'right', 'bottom'        # came up: above left
        ax.annotate(fmt % yi, (xi, yi), textcoords='offset points', xytext=off,
                    ha=ha, va=va, fontsize=fs, color=col)


def note(ax, text, loc='ul', fs=7):
    """Boxed note inside the panel: 'ul' upper left, 'ur' upper right,
    'll' lower left, 'lr' lower right."""
    x = 0.03 if loc[1] == 'l' else 0.97
    y = 0.96 if loc[0] == 'u' else 0.04
    ax.text(x, y, text, transform=ax.transAxes, ha='left' if loc[1] == 'l' else 'right',
            va='top' if loc[0] == 'u' else 'bottom', fontsize=fs, linespacing=1.3,
            bbox=dict(fc='white', ec='0.7', lw=0.5, pad=2.0), zorder=10)


def scale(ax, lo, hi, minspan, top=0.8, bottom=0.3):
    """Individual y-scale: the data occupy the middle of the panel, with free
    space above for the note and the labels and below for the labels."""
    span = max(hi - lo, minspan)
    mid = 0.5 * (hi + lo)
    lo, hi = mid - span / 2, mid + span / 2
    ax.set_ylim(lo - bottom * span, hi + top * span)
