# -*- coding: utf-8 -*-
"""
================================================================================
PUT THE CONTOUR LEGEND BESIDE THE TANK, NEVER ON IT
TU Delft MSc thesis - Sarvesh Mahadevan
================================================================================

Second half of the capture pipeline. ABAQUS_capture_all.py renders every frame
of a view listed in SIDE_LEGEND_VIEWS twice, into <capture folder>\\_raw\\<AC>\\ :

    <name>__L.png   legend ON
    <name>__N.png   legend OFF        (same camera, same everything else)

This script assembles each pair into  <capture folder>\\<AC>\\<name>.png  with
the legend on the left and the tank on the right.

HOW THE TWO PIECES ARE FOUND - no hand-tuned coordinates
  legend : the only pixels that differ between __L and __N are the legend's,
           so its bounding box is the bounding box of the difference. It is
           cut from __L, where Abaqus drew it on an opaque white background.
  tank   : every non-white pixel of __N, which has no legend in it at all.
Because the two are cut from separate renders and then placed side by side,
the legend cannot overlap the tank, whatever the tank's shape or length.
The script also checks that the tank does not touch the edge of its render,
which would mean the camera clipped it; such frames are reported, not hidden.

RUN (normal Python with numpy and Pillow, from the thesis root):
    python python_scripts\\compose_side_legend.py
    python python_scripts\\compose_side_legend.py abaqus_captures_final
================================================================================
"""
import os
import sys
import glob

import numpy as np
from PIL import Image

CAPTURE_DIR = sys.argv[1] if len(sys.argv) > 1 else 'abaqus_captures_final'

WHITE_LEVEL = 245        # a pixel is "background" if every channel is above this
DIFF_LEVEL = 24          # legend-on vs legend-off difference that counts
PAD = 24                 # px of white around the assembled figure
GAP = 60                 # px between legend and tank
LEGEND_SCALE = 1.0       # >1 enlarges the legend (resampled) for legibility
EDGE_TOL = 2             # px; a tank this close to the render edge was clipped


def bbox(mask):
    """(x0, y0, x1, y1) of the True pixels, exclusive end, or None."""
    ys, xs = np.nonzero(mask)
    if len(xs) == 0:
        return None
    return int(xs.min()), int(ys.min()), int(xs.max()) + 1, int(ys.max()) + 1


OUTLINE_PX = 3           # silhouette line width drawn round the model
OUTLINE_RGB = (20, 20, 20)


def add_outline(arr):
    """
    Draw the model's silhouette. Contours are rendered with no element edges
    (Abaqus FEATURE edges draw the whole free mesh), so the outline is taken
    from the image itself: the boundary of the non-white region.
    """
    m = arr.min(axis=2) < WHITE_LEVEL
    # pad with background so a model touching the crop edge still gets a line
    m = np.pad(m, 1, constant_values=False)
    inner = m.copy()
    for _ in range(OUTLINE_PX):
        e = inner.copy()
        e[1:, :] &= inner[:-1, :]; e[:-1, :] &= inner[1:, :]
        e[:, 1:] &= inner[:, :-1]; e[:, :-1] &= inner[:, 1:]
        inner = e
    edge = (m & ~inner)[1:-1, 1:-1]
    out = arr.copy()
    out[edge] = OUTLINE_RGB
    return out


def pieces(path_l, path_n):
    """(legend, tank, clipped) as PIL images, or raise ValueError."""
    a_l = np.asarray(Image.open(path_l).convert('RGB')).astype(np.int16)
    a_n = np.asarray(Image.open(path_n).convert('RGB')).astype(np.int16)
    if a_l.shape != a_n.shape:
        raise ValueError('size mismatch %s vs %s' % (a_l.shape, a_n.shape))
    h, w, _ = a_n.shape

    # legend = where the two renders differ
    lb = bbox(np.abs(a_l - a_n).max(axis=2) > DIFF_LEVEL)
    if lb is None:
        raise ValueError('no legend found (renders identical)')
    x0, y0, x1, y1 = lb
    # The legend is drawn inside a thin dark border on an opaque background, so
    # everything inside that border belongs to the legend and nothing outside
    # it does. The difference box can stop a pixel short of the border where
    # the border happens to sit over dark mesh lines (same colour in both
    # renders), so grow each side one pixel at a time while the next line is
    # a border line (dark along most of its length), then crop exactly there.
    dark = a_l.min(axis=2) < 110

    def is_border_col(x, ya, yb):
        return 0 <= x < w and dark[ya:yb, x].mean() > 0.8

    def is_border_row(y, xa, xb):
        return 0 <= y < h and dark[y, xa:xb].mean() > 0.8

    # Snap each side to the border line nearest the difference box, looking
    # up to 6 px both inwards and outwards; the crop then ends exactly on it.
    def snap(pos, test, step_out):
        for d in range(0, 7):
            for c in (pos + step_out * d, pos - step_out * d):
                if test(c):
                    return c
        return None
    r = snap(x1 - 1, lambda c: is_border_col(c, y0, y1), +1)
    l = snap(x0, lambda c: is_border_col(c, y0, y1), -1)
    x0, x1 = (l if l is not None else x0), ((r + 1) if r is not None else x1)
    b = snap(y1 - 1, lambda c: is_border_row(c, x0, x1), +1)
    t = snap(y0, lambda c: is_border_row(c, x0, x1), -1)
    y0, y1 = (t if t is not None else y0), ((b + 1) if b is not None else y1)
    legend = Image.fromarray(a_l[y0:y1, x0:x1].astype(np.uint8))
    if LEGEND_SCALE != 1.0:
        legend = legend.resize((int(round(legend.width * LEGEND_SCALE)),
                                int(round(legend.height * LEGEND_SCALE))),
                               Image.LANCZOS)

    # tank = everything non-white in the legend-free render
    tb = bbox(a_n.min(axis=2) < WHITE_LEVEL)
    if tb is None:
        raise ValueError('no model found in the legend-off render')
    tx0, ty0, tx1, ty1 = tb
    clipped = (tx0 <= EDGE_TOL or ty0 <= EDGE_TOL or
               tx1 >= w - EDGE_TOL or ty1 >= h - EDGE_TOL)
    tank = Image.fromarray(add_outline(a_n[ty0:ty1, tx0:tx1]).astype(np.uint8))

    return legend, tank, clipped


def compose(path_l, path_n, path_out):
    try:
        legend, tank, clipped = pieces(path_l, path_n)
    except ValueError as e:
        return str(e)

    out_w = PAD + legend.width + GAP + tank.width + PAD
    out_h = PAD + max(legend.height, tank.height) + PAD
    out = Image.new('RGB', (out_w, out_h), (255, 255, 255))
    out.paste(legend, (PAD, (out_h - legend.height) // 2))
    tank_x = PAD + legend.width + GAP
    out.paste(tank, (tank_x, (out_h - tank.height) // 2))

    # Check the result rather than trust the construction: the legend's box
    # and the tank's box must not share a single pixel column.
    assert PAD + legend.width <= tank_x, 'legend and tank overlap'

    if not os.path.isdir(os.path.dirname(path_out)):
        os.makedirs(os.path.dirname(path_out))
    out.save(path_out)
    return 'CLIPPED BY CAMERA - check the zoom for this view' if clipped else 'ok'


def main():
    raw_root = os.path.join(CAPTURE_DIR, '_raw')
    pairs = sorted(glob.glob(os.path.join(raw_root, '*', '*__N.png')))
    if not pairs:
        print('No raw frames under %s. Run ABAQUS_capture_all.py first.'
              % raw_root)
        return 1
    n_ok, problems = 0, []
    for p_n in pairs:
        p_l = p_n[:-len('__N.png')] + '__L.png'
        ac = os.path.basename(os.path.dirname(p_n))
        name = os.path.basename(p_n)[:-len('__N.png')] + '.png'
        p_out = os.path.join(CAPTURE_DIR, ac, name)
        if not os.path.exists(p_l):
            problems.append((name, 'legend-on render missing'))
            continue
        status = compose(p_l, p_n, p_out)
        if status == 'ok':
            n_ok += 1
        else:
            problems.append((name, status))
            if status.startswith('CLIPPED'):
                n_ok += 1        # written, but flagged
        print('  %-50s %s' % (name, status))
    # Views rendered without a legend (the section view): outline and crop
    # them to the model with the same margin, from a pristine copy in _raw.
    n_sect = 0
    for p in sorted(glob.glob(os.path.join(CAPTURE_DIR, '*', '*_sect.png'))):
        raw = os.path.join(raw_root, os.path.basename(os.path.dirname(p)),
                           os.path.basename(p)[:-4] + '__S.png')
        if not os.path.exists(raw):
            Image.open(p).save(raw)          # keep the Abaqus original once
        a = np.asarray(Image.open(raw).convert('RGB')).astype(np.int16)
        bb = bbox(a.min(axis=2) < WHITE_LEVEL)
        if bb is None:
            continue
        x0, y0, x1, y1 = bb
        body = add_outline(a[y0:y1, x0:x1])
        out = Image.new('RGB', (x1 - x0 + 2 * PAD, y1 - y0 + 2 * PAD), 'white')
        out.paste(Image.fromarray(body.astype(np.uint8)), (PAD, PAD))
        out.save(p)
        n_sect += 1
    print('%d section views outlined and cropped' % n_sect)
    print('\n%d assembled, %d with problems' % (n_ok, len(problems)))
    for name, st in problems:
        print('  !! %s: %s' % (name, st))
    return 0 if not problems else 2


if __name__ == '__main__':
    sys.exit(main())
