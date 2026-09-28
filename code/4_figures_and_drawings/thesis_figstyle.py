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
LC_LABEL = {'Baseline': 'Baseline', 'LC1_Maneuver': 'LC1 Manoeuvre',
            'LC3_Combined': 'LC3 Combined', 'LC7_Emergency': 'LC7 Emergency'}


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
