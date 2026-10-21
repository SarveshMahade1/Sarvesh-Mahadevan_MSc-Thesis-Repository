# -*- coding: utf-8 -*-
"""
ALL DATA GRAPHS OF THE REPORT, FROM THE FINAL RESULTS (28 Sep 2026, v2)
======================================================================

Reads results/processed/final_numbers.json and results/processed/hand_calcs.json
(written by final_numbers.py and hand_calcs.py) and draws every graph used in
Chapters 3 to 5. Nothing is typed in by hand.

Version 2 (readability):
  * every graph is drawn at the full text width and included at \\textwidth,
    so the type prints at the size set here (10 pt labels, 9 pt ticks);
  * every line, marker and reference line in a figure is named in a legend
    placed below the panels, never on top of the data;
  * no free text is written inside the axes; what a text label used to say
    is now either a legend entry or in the caption.

    python code/3_results_processing/final_numbers.py
    python code/3_results_processing/hand_calcs.py
    python code/4_figures_and_drawings/plot_final_figures.py
"""
import csv
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
import numpy as np                                             # noqa: E402
import thesis_figstyle as tfs                                  # noqa: E402
import matplotlib.pyplot as plt                                # noqa: E402
from matplotlib.lines import Line2D                            # noqa: E402

tfs.apply_style()
plt.rcParams.update({'font.size': 10, 'axes.labelsize': 10, 'axes.titlesize': 10,
                     'xtick.labelsize': 9, 'ytick.labelsize': 9, 'legend.fontsize': 9,
                     'lines.linewidth': 1.6, 'lines.markersize': 5.5})
FIG = RP.FIGURES
D = json.load(open(os.path.join(RP.PROCESSED, 'final_numbers.json')))
H = json.load(open(os.path.join(RP.PROCESSED, 'hand_calcs.json')))
AC = ['E190', 'A320', 'A350']
LCS = ['Baseline', 'LC1_Maneuver', 'LC3_Combined', 'LC7_Emergency']
FILLS = [0.1, 0.3, 0.5, 0.75, 0.9]
W = 455.24 / 72.27            # full text width, inches
C, L = tfs.AC_COLOR, tfs.AC_LABEL
LCL = tfs.LC_LABEL
LC_STYLE = {'Baseline': dict(ls='-', marker='o', ms=7, mfc='none', mew=1.2, zorder=4),
            'LC1_Maneuver': dict(ls='--', marker='s', ms=3.5),
            'LC3_Combined': dict(ls='-.', marker='^'), 'LC7_Emergency': dict(ls='-', marker='D')}
LC_COL = {'Baseline': '#7f7f7f', 'LC1_Maneuver': '#000000', 'LC3_Combined': '#2E75B6',
          'LC7_Emergency': '#C00000'}
ALLOW = 262.0
REF = dict(color='k', lw=0.9, ls=':')


def pr(ac, lc, f, k, stage='production'):
    return D[ac][stage]['%s|%.2f' % (lc, f)][k]


def series(ac, lc, k, stage='production'):
    return [pr(ac, lc, f, k, stage) for f in FILLS]


def title(ax, s):
    # Review, Oct 2026: no overhead titles; the panel letter only, the rest
    # is said in the caption of the report.
    ax.set_title(s.split(')')[0] + ')', loc='left', fontsize=10)


def ac_handles():
    return [Line2D([], [], color=C[a], marker='o', lw=1.6) for a in AC], [L[a] for a in AC]


def legend_below(fig, handles, labels, ncol, y=0.0):
    fig.legend(handles, labels, loc='lower center', ncol=ncol, frameon=False,
               bbox_to_anchor=(0.5, y), handlelength=2.6, columnspacing=1.6)


def save(fig, name, bottom):
    fig.tight_layout(rect=(0, bottom, 1, 1), w_pad=1.4)
    tfs.save(fig, os.path.join(FIG, name))


# ---------------------------------------------------------------------------
# Chapter 3: the finite element sizing of the wall
# ---------------------------------------------------------------------------
def fig_sizing_history():
    rows = {'E190': D['E190']['sizing_iterations'], 'A320': D['A320']['sizing_iterations']}
    a350 = []
    rd = csv.reader(open(os.path.join(RP.SIZING, 'A350', 'A350_sizing_results_run1_oscillation.csv')))
    hd = next(rd)
    for r in rd:
        if len(r) == len(hd) + 1:
            r = [r[0] + ';' + r[1]] + r[2:]
        d = dict(zip(hd, r))
        a350.append(dict(stage=d['stage'], t_wall_mm=float(d['t_wall_mm']), p_gauge_bar=float(d['p_gauge_bar']),
                         fos=float(d['fos']), buckling_fos=float(d['buckling_fos'])))
    rows['A350'] = a350
    fig, axs = plt.subplots(2, 3, figsize=(W, 4.7), gridspec_kw=dict(height_ratios=[1.3, 1]))
    for k, ac in enumerate(AC):
        rr = [r for r in rows[ac] if r['stage'].startswith('iteration')]
        it = lambda r: int(r['stage'].split(';')[0].split()[1])
        y = [(it(r), r['fos']) for r in rr if r['p_gauge_bar'] > 1.0]
        # reserve factor RF = allowable / applied, 1.5 inside both allowables
        b = [(it(r), r['buckling_fos'] / 1.5) for r in rr if r['p_gauge_bar'] < 1.0 and math.isfinite(r['buckling_fos'])]
        t = [(it(r), r['t_wall_mm']) for r in rr if r['p_gauge_bar'] > 1.0]
        ax = axs[0, k]
        xs, ys = [i for i, _ in y], [v for _, v in y]
        ax.plot(xs, ys, '-o', color='#1F3864')
        y0b = b[0][1] if b else -1.0
        tfs.label_points(ax, xs, ys, '#1F3864', '%.2f',
                         side=([('ll' if ys[0] < y0b else 'ul') if b else 'ur']
                               + ['below' if b else 'peak'] * (len(xs) - 1)))
        allv = list(ys)
        if b:
            xb, yb = [i for i, _ in b], [v for _, v in b]
            ax.plot(xb, yb, '--s', color='#C00000')
            tfs.label_points(ax, xb, yb, '#C00000', '%.2f', side=[('ul' if yb[0] > ys[0] else 'll'), 'ur'] + ['above'] * (len(xb) - 2))
            allv += yb
        else:
            tfs.note(ax, 'Barrel in axial tension at 0.187 bar:\nno buckling reserve factor', loc='ur')
        ax.axhline(1.0, color='k', lw=0.9, ls=':')
        tfs.scale(ax, min(allv + [1.0]), max(allv), 0.5, top=0.25 if b else 0.75, bottom=0.2)
        ax.grid(True)
        ax.set_title('(%s) %s' % ('abc'[k], L[ac]), loc='left', fontsize=10)
        ax2 = axs[1, k]
        xt, yt = [i for i, _ in t], [v for _, v in t]
        ax2.plot(xt, yt, '-^', color='#555555')
        tfs.label_points(ax2, xt, yt, '#555555', '%.2f', side=['ur'] + ['peak'] * (len(xt) - 1))
        tfs.scale(ax2, min(yt), max(yt), 0.2, top=0.3, bottom=0.15)
        ax2.set_title('(%s) %s' % ('def'[k], L[ac]), loc='left', fontsize=10)
        for a_ in (ax, ax2):        # every panel carries both axis labels
            a_.set_xlabel('Sizing iteration [-]')
            a_.set_xticks(sorted(set(i for i, _ in t)))
            a_.set_xlim(min(xt) - 0.9, max(xt) + 0.5)
            a_.tick_params(labelbottom=True)
        ax.set_ylabel('Reserve factor RF [-]')
        ax2.set_ylabel('Trial wall thickness $t$ [mm]')
        ax2.grid(True)
    hs =[Line2D([], [], color='#1F3864', marker='o'), Line2D([], [], color='#C00000', ls='--', marker='s'),
          Line2D([], [], color='#555555', marker='^'), Line2D([], [], color='k', ls=':')]
    legend_below(fig, hs, ['Yield reserve factor at 1.5 bar', 'Buckling reserve factor at 0.187 bar', 'Trial wall thickness',
                           'Requirement RF = 1'], 2)
    save(fig, 'fig3_sizing_history', 0.12)


# ---------------------------------------------------------------------------
# Chapter 5
# ---------------------------------------------------------------------------
def fig_e190_matrix():
    fig, axs = plt.subplots(1, 2, figsize=(W, 2.9))
    for lc in LCS:
        axs[0].plot(FILLS, series('E190', lc, 'mises_MPa'), color=LC_COL[lc], **LC_STYLE[lc])
        axs[1].plot(FILLS, series('E190', lc, 'dT'), color=LC_COL[lc], **LC_STYLE[lc])
    axs[0].axhline(ALLOW, **REF)
    axs[0].set_ylim(150, 275)
    axs[0].set_ylabel(r'Far-field $\sigma_{\mathrm{vM}}$ [MPa]')
    axs[1].axhline(H['E190']['thermal']['dT_asymptote_K'], color='#555555', lw=1.0, ls=(0, (5, 2)))
    axs[1].set_ylabel(r'$\Delta T_{\mathrm{w}}$ [K]')
    axs[1].set_ylim(3.35, 3.78)
    for ax, s in zip(axs, ['(a)', '(b)']):
        ax.set_xlabel(r'Fill level $\phi_{\mathrm{f}}$ [-]')
        ax.set_xticks(FILLS)
        ax.grid(True)
        title(ax, s)
    hs = [Line2D([], [], color=LC_COL[lc], **LC_STYLE[lc]) for lc in LCS] + \
         [Line2D([], [], **REF), Line2D([], [], color='#555555', lw=1.0, ls=(0, (5, 2)))]
    legend_below(fig, hs, [LCL[lc] for lc in LCS] + [r'Allowable 262 MPa (a)', r'Asymptote $q_{\mathrm{MLI}}/h_{\mathrm{V}}$ (b)'], 3)
    save(fig, 'fig5_e190_matrix', 0.2)


def fig_head_mechanism():
    fig, axs = plt.subplots(1, 2, figsize=(W, 2.9))
    for ac in AC:
        for lc, ls in (('Baseline', ':'), ('LC1_Maneuver', '--'), ('LC3_Combined', '-.'), ('LC7_Emergency', '-')):
            dep = [D[ac]['heads']['%s|%.2f' % (lc, f)]['depth_m'] for f in FILLS]
            hp = [100 * D[ac]['heads']['%s|%.2f' % (lc, f)]['head_over_p'] for f in FILLS]
            kw = dict(color=C[ac], ls=ls, marker='o', ms=3.5)
            axs[0].plot(FILLS, dep, **kw)
            axs[1].plot(FILLS, hp, **kw)
    axs[0].set_ylabel(r'Liquid depth along $\hat{\mathbf{g}}$ [m]')
    axs[1].set_ylabel(r'Hydrostatic pressure ratio $p_{\mathrm{h}}/p$ [%]')
    axs[1].set_yscale('log')
    for ax, s in zip(axs, ['(a) Depth of liquid', '(b) Hydrostatic head']):
        ax.set_xlabel(r'Fill level $\phi_{\mathrm{f}}$ [-]')
        ax.set_xticks(FILLS)
        ax.grid(True, which='both')
        title(ax, s)
    hs, ls_ = ac_handles()
    hs += [Line2D([], [], color='k', ls=l) for l in (':', '--', '-.', '-')]
    ls_ += [r'Baseline, 1.0 g', r'LC1 Manoeuvre, 2.5 g', r'LC3 Combined, 2.92 g at 31.0$^\circ$',
            r'LC7 Emergency, 11.22 g at 57.7$^\circ$']
    legend_below(fig, hs, ls_, 4)
    save(fig, 'fig5_head_mechanism', 0.2)


def fig_stress_vs_head():
    fig, ax = plt.subplots(figsize=(W, 3.6))
    mk = {'Baseline': 'o', 'LC1_Maneuver': 's', 'LC3_Combined': '^', 'LC7_Emergency': 'D'}
    for ac in AC:
        pRt = D[ac]['hoop_pR_t_design_MPa']
        for lc in LCS:
            x = [1 + D[ac]['heads']['%s|%.2f' % (lc, f)]['head_over_p'] for f in FILLS]
            y = [pr(ac, lc, f, 'mises_MPa') / pRt for f in FILLS]
            ax.plot(x, y, ls='none', marker=mk[lc], color=C[ac],
                    mfc='none' if lc in ('Baseline', 'LC1_Maneuver') else C[ac])
    xx = np.linspace(1.0, 2.5, 50)
    ax.plot(xx, 0.866 * xx, **REF)
    ax.set_xlabel(r'$1 + p_{\mathrm{h}}/p_{\mathrm{gauge}}$ at the deepest point of the liquid [-]')
    ax.set_ylabel(r'Far-field $\sigma_{\mathrm{vM}}\,/\,(pR/t)$ [-]')
    ax.grid(True)
    hs, ls_ = ac_handles()
    hs = [Line2D([], [], color=C[a], marker='o', ls='none') for a in AC]
    hs += [Line2D([], [], color='k', marker=mk[l], ls='none', mfc='none' if l in ('Baseline', 'LC1_Maneuver') else 'k')
           for l in LCS] + [Line2D([], [], **REF)]
    ls_ += [LCL[l] for l in LCS] + [r'Barrel membrane, $0.866\,(1 + p_{\mathrm{h}}/p)$']
    legend_below(fig, hs, ls_, 4)
    save(fig, 'fig5_stress_vs_head', 0.16)


def fig_sensitivity_vs_tilt():
    # Review, Oct 2026: the change with fill is plotted against the quantity the
    # mechanism says controls it, the hydrostatic head relative to the gauge
    # pressure at full fill, p_{\mathrm{h}}/p = rho |n| g d / p (magnitude AND depth),
    # and, for LC7, against the liquid depth along g_hat, which is set by the
    # absolute length of the tank. Slenderness ranks the three tanks the same
    # way, so the data cannot separate length from slenderness.
    mk = {'Baseline': 'o', 'LC1_Maneuver': 's', 'LC3_Combined': '^', 'LC7_Emergency': 'D'}
    fig, axs = plt.subplots(1, 2, figsize=(W, 2.9))
    for ac in AC:
        s_ = D[ac]['fill_sensitivity_pct']
        for lc in LCS:
            h = 100 * D[ac]['heads']['%s|0.90' % lc]['head_over_p']
            axs[0].plot([h], [s_[lc]], ls='none', marker=mk[lc], color=C[ac], ms=6,
                        mfc='none' if lc in ('Baseline', 'LC1_Maneuver') else C[ac])
        axs[1].plot([D[ac]['heads']['LC7_Emergency|0.90']['depth_m']], [s_['LC7_Emergency']], 'D',
                    color=C[ac], ms=7)
    axs[0].set_xscale('log')
    axs[0].set_xlabel(r'Hydrostatic pressure ratio $p_{\mathrm{h}}/p$ at $\phi_{\mathrm{f}} = 0.90$ [%]')
    axs[0].set_ylabel(r'Fill-level sensitivity $S_{\phi}$ [%]')
    axs[1].set_xlabel(r'Liquid depth along $\hat{\mathbf{g}}$, LC7, $\phi_{\mathrm{f}} = 0.90$ [m]')
    axs[1].set_ylabel(r'Fill-level sensitivity $S_{\phi}$, LC7 [%]')
    axs[1].set_xlim(0, 30)
    axs[1].set_ylim(0, 115)
    for ax, t in zip(axs, ['(a)', '(b)']):
        ax.grid(True, which='both')
        title(ax, t)
    hs, ls_ = ac_handles()
    hs = [Line2D([], [], color=C[a], marker='o', ls='none') for a in AC]
    hs += [Line2D([], [], color='k', marker=mk[l], ls='none', mfc='none' if l in ('Baseline', 'LC1_Maneuver') else 'k')
           for l in LCS]
    ls_ += [LCL[l] for l in LCS]
    legend_below(fig, hs, ls_, 4)
    save(fig, 'fig5_sensitivity', 0.2)


def fig_three_aircraft():
    # Review, Oct 2026: absolute far-field stress in MPa, with the allowable drawn.
    fig, axs = plt.subplots(1, 2, figsize=(W, 2.9), sharey=True)
    for ax, lc, t in zip(axs, ('Baseline', 'LC7_Emergency'), ('(a)', '(b)')):
        for ac in AC:
            ax.plot(FILLS, [pr(ac, lc, f, 'mises_MPa') for f in FILLS], '-o', color=C[ac])
        ax.axhline(ALLOW, **REF)
        title(ax, t)
        ax.set_xlabel(r'Fill level $\phi_{\mathrm{f}}$ [-]')
        ax.set_xticks(FILLS)
        ax.grid(True)
    axs[0].set_ylabel(r'Far-field $\sigma_{\mathrm{vM}}$ [MPa]')
    axs[0].set_ylim(0.0, 290)
    hs, ls_ = ac_handles()
    legend_below(fig, hs + [Line2D([], [], **REF)], ls_ + [r'Allowable $F_{\mathrm{ty}}/1.5$ = 262 MPa'], 4)
    save(fig, 'fig5_three_aircraft', 0.13)


def fig_thermal():
    fig, axs = plt.subplots(1, 2, figsize=(W, 2.9))
    for ac in AC:
        for lc, ls, mk in (('Baseline', '-', 'o'), ('LC7_Emergency', '--', 'D')):
            axs[0].plot(FILLS, series(ac, lc, 'dT'), ls=ls, marker=mk, color=C[ac], ms=4)
        axs[0].axhline(H[ac]['thermal']['dT_asymptote_K'], color=C[ac], lw=0.9, ls=':')
        axs[1].plot([H[ac]['thermal']['q_Wm2']], [pr(ac, 'Baseline', 0.1, 'dT')], 'o', color=C[ac], ms=7)
    qq = np.linspace(4, 12, 20)
    axs[1].plot(qq, qq / 1.46, **REF)
    axs[0].set_xlabel(r'Fill level $\phi_{\mathrm{f}}$ [-]')
    axs[0].set_ylabel(r'$\Delta T_{\mathrm{w}}$ [K]')
    axs[0].set_xticks(FILLS)
    axs[1].set_xlabel(r'Insulation flux $q_{\mathrm{MLI}}$ [W/m$^2$]')
    axs[1].set_ylabel(r'$\Delta T_{\mathrm{w}}$ at $\phi_{\mathrm{f}} = 0.10$ [K]')
    for ax, s in zip(axs, ['(a) Against fill', '(b) Against insulation flux']):
        ax.grid(True)
        title(ax, s)
    hs, ls_ = ac_handles()
    hs += [Line2D([], [], color='k', ls='-', marker='o', ms=4), Line2D([], [], color='k', ls='--', marker='D', ms=4),
           Line2D([], [], color='k', lw=0.9, ls=':')]
    ls_ += ['Baseline', 'LC7 Emergency', r'$\Delta T_{\mathrm{w}} = q_{\mathrm{MLI}}/h_{\mathrm{V}}$']
    legend_below(fig, hs, ls_, 3)
    save(fig, 'fig5_thermal', 0.2)


def fig_reactions():
    fig, axs = plt.subplots(1, 3, figsize=(W, 2.9), sharey=True)
    for ax, ac in zip(axs, AC):
        for lc, ls in (('Baseline', '--'), ('LC1_Maneuver', '-')):
            tot = series(ac, lc, 'sumRF_N')
            aft = series(ac, lc, 'RF_aft_kN')
            fwd = series(ac, lc, 'RF_fwd_kN')
            mid = [t / 1e3 - a - b for t, a, b in zip(tot, aft, fwd)]
            ax.plot(FILLS, [100 * m / (t / 1e3) for m, t in zip(mid, tot)], ls=ls, marker='o', color=C[ac], ms=4)
        ax.axhline(62.5, **REF)
        ax.axhline(0.0, color='k', lw=0.6)
        title(ax, '(%s) %s' % ('abc'[AC.index(ac)], L[ac]))
        ax.set_xlabel(r'Fill level $\phi_{\mathrm{f}}$ [-]')
        ax.set_xticks([0.1, 0.5, 0.9])
        ax.grid(True)
        ax.set_ylim(-100, 90)
    axs[0].set_ylabel('Mid-ring share of reaction [%]')
    hs = [Line2D([], [], color='k', ls='--', marker='o', ms=4), Line2D([], [], color='k', ls='-', marker='o', ms=4),
          Line2D([], [], **REF)]
    legend_below(fig, hs, ['Baseline (1.0 g)', 'LC1 Manoeuvre (2.5 g)', 'Uniform beam on three supports, 62.5 %'], 3)
    save(fig, 'fig5_ring_share', 0.13)


def fig_margins():
    # Review, Oct 2026: one convention, reserve factor RF = allowable / applied,
    # with the project margin of 1.5 inside both allowables; requirement RF >= 1.
    fig, axs = plt.subplots(1, 2, figsize=(W, 3.0))
    for ac in AC:
        axs[0].plot(FILLS, series(ac, 'LC7_Emergency', 'fos'), '-o', color=C[ac])
        ax_ = series(ac, 'LC7_Emergency', 'axial_mem_min_MPa', 'stability')
        rf = [(f, D[ac]['sigma_cr_MPa'] / 1.5 / -a) for f, a in zip(FILLS, ax_) if a < 0]
        if rf:
            axs[1].plot([f for f, _ in rf], [r for _, r in rf], '-D', color=C[ac])
    for ax in axs:
        ax.axhline(1.0, **REF)
    axs[0].set_ylabel(r'Yield RF $=\sigma_{\mathrm{allow}}/\sigma_{\mathrm{vM}}$ [-]')
    axs[0].set_ylim(0.8, 4.0)
    axs[1].set_ylabel(r'Buckling RF $=\sigma_{\mathrm{cr}}/(1.5\,|\sigma_{x,\mathrm{min}}|)$ [-]')
    axs[1].set_yscale('log')
    axs[1].set_ylim(0.8, 12)
    for ax, t in zip(axs, ['(a)', '(b)']):
        ax.set_xlabel(r'Fill level $\phi_{\mathrm{f}}$ [-]')
        ax.set_xticks(FILLS)
        ax.grid(True, which='both')
        title(ax, t)
    hs, ls_ = ac_handles()
    hs += [Line2D([], [], **REF)]
    ls_ += ['Requirement RF = 1']
    legend_below(fig, hs, ls_, 4)
    save(fig, 'fig5_margins', 0.2)


def fig_convergence():
    fig, axs = plt.subplots(1, 2, figsize=(W, 2.9))
    for ac in AC:
        cd, cm = D[ac]['convergence']['pdesign'], D[ac]['convergence']['pmin']
        n = [c['n_el'] for c in cd]
        axs[0].plot(n, [100 * (c['mises_MPa'] / cd[-1]['mises_MPa'] - 1) for c in cd], '-o', color=C[ac])
        axs[0].plot([n[1]], [100 * (cd[1]['mises_MPa'] / cd[-1]['mises_MPa'] - 1)], 'o', ms=10, mfc='none', color=C[ac])
        if cm[-1]['axial_mem_min_MPa'] < 0:
            axs[1].plot(n, [c['buckling_fos'] / 1.5 for c in cm], '-D', color=C[ac])
            axs[1].plot([n[1]], [cm[1]['buckling_fos'] / 1.5], 'D', ms=10, mfc='none', color=C[ac])
    axs[0].axhline(0, color='k', lw=0.6)
    axs[0].set_ylabel(r'$\sigma_{\mathrm{vM}}$ relative to finest mesh [%]')
    axs[1].axhline(1.0, **REF)
    axs[1].set_ylabel('Buckling RF [-]')
    axs[1].set_ylim(0.95, 1.17)
    for ax, t in zip(axs, ['(a)', '(b)']):
        ax.set_xscale('log')
        ax.set_xlabel('Elements [-]')
        ax.grid(True, which='both')
        title(ax, t)
    hs, ls_ = ac_handles()
    hs += [Line2D([], [], color='k', marker='o', ms=10, mfc='none', ls='none'), Line2D([], [], **REF)]
    ls_ += ['Production mesh', 'Requirement RF = 1']
    legend_below(fig, hs, ls_, 5)
    save(fig, 'fig5_convergence', 0.13)


def fig_ovality():
    fig, ax = plt.subplots(figsize=(W, 2.9))
    for ac in AC:
        ax.plot(FILLS, series(ac, 'Baseline', 'ovality_mm'), '--o', color=C[ac], ms=4)
        ax.plot(FILLS, series(ac, 'LC7_Emergency', 'ovality_mm'), '-D', color=C[ac], ms=4)
    ax.set_xlabel(r'Fill level $\phi_{\mathrm{f}}$ [-]')
    ax.set_ylabel('Section distortion [mm]')
    ax.set_xticks(FILLS)
    ax.grid(True)
    hs, ls_ = ac_handles()
    hs += [Line2D([], [], color='k', ls='--', marker='o', ms=4), Line2D([], [], color='k', ls='-', marker='D', ms=4)]
    ls_ += ['Baseline', 'LC7 Emergency']
    legend_below(fig, hs, ls_, 5)
    save(fig, 'fig5_ovality', 0.13)


# ---------------------------------------------------------------------------
# Absolute-axis versions (October 2026): the same data as fig5_stress_vs_head
# and fig5_sensitivity, but in MPa and kPa, one line per load case connecting
# the five fill levels in order, so that each point can be read physically.
# ---------------------------------------------------------------------------
MK = {'Baseline': 'o', 'LC1_Maneuver': 's', 'LC3_Combined': '^', 'LC7_Emergency': 'D'}
RC = json.load(open(os.path.join(RP.PROCESSED, 'review_checks.json')))


def local_pressure_kPa(ac, lc, f):
    return 150.0 + 100.0 * D[ac]['heads']['%s|%.2f' % (lc, f)]['head_bar']


def fig_stress_vs_pressure_abs():
    fig, axs = plt.subplots(1, 3, figsize=(W, 3.0), sharey=True)
    for ax, ac, t in zip(axs, AC, ['(a)', '(b)', '(c)']):
        Rt = D[ac]['R_over_t']
        xx = np.linspace(140, 380, 50)
        ax.plot(xx, 0.866 * Rt * xx * 1e-3, color='k', lw=0.9, ls=':')
        a = RC[ac]['peak_fit']['slope_a']; b = RC[ac]['peak_fit']['intercept_b_MPa']
        ax.plot(xx, a * Rt * xx * 1e-3 + b, color=C[ac], lw=0.9, ls='--')
        ax.axhline(ALLOW, color='k', lw=0.9, ls='-.')
        for lc in LCS:
            x = [local_pressure_kPa(ac, lc, f) for f in FILLS]
            y = [pr(ac, lc, f, 'mises_MPa') for f in FILLS]
            ax.plot(x, y, ls='-', lw=0.8, marker=MK[lc], color=C[ac], ms=5,
                    mfc='white' if lc in ('Baseline', 'LC1_Maneuver') else C[ac])
        ax.set_xlim(140, 380)
        ax.set_ylim(0, 300)
        ax.set_xlabel('Pressure at the deepest\npoint of the liquid [kPa]')
        ax.grid(True)
        title(ax, t)
    axs[0].set_ylabel(r'Far-field $\sigma_{\mathrm{vM}}$ [MPa]')
    hs = [Line2D([], [], color=C[a], marker='o', ls='-', lw=0.8) for a in AC]
    ls_ = [L[a] for a in AC]
    hs += [Line2D([], [], color='k', marker=MK[l], ls='none', mfc='white' if l in ('Baseline', 'LC1_Maneuver') else 'k')
           for l in LCS]
    ls_ += [LCL[l] for l in LCS]
    hs += [Line2D([], [], color='k', lw=0.9, ls=':'), Line2D([], [], color='grey', lw=0.9, ls='--'),
           Line2D([], [], color='k', lw=0.9, ls='-.')]
    ls_ += [r'Membrane theory $0.866\,(p+p_{\mathrm{h}})R/t$', 'Linear fit of the twenty cases',
            r'Allowable $F_{\mathrm{ty}}/1.5$']
    legend_below(fig, hs, ls_, 4)
    save(fig, 'fig5_stress_vs_pressure_abs', 0.24)


def fig_stress_increase_abs():
    fig, ax = plt.subplots(figsize=(W, 3.3))
    xmax = 0
    for ac in AC:
        Rt = D[ac]['R_over_t']
        a = RC[ac]['peak_fit']['slope_a']
        for lc in LCS:
            dph = 100.0 * (D[ac]['heads']['%s|0.90' % lc]['head_bar'] - D[ac]['heads']['%s|0.10' % lc]['head_bar'])
            ds = pr(ac, lc, 0.9, 'mises_MPa') - pr(ac, lc, 0.1, 'mises_MPa')
            xmax = max(xmax, dph)
            ax.plot([dph], [ds], ls='none', marker=MK[lc], color=C[ac], ms=6,
                    mfc='white' if lc in ('Baseline', 'LC1_Maneuver') else C[ac])
        xx = np.linspace(0, 220, 50)
        ax.plot(xx, a * Rt * xx * 1e-3, color=C[ac], lw=0.9, ls='--')
    ax.set_xlim(0, 220)
    ax.set_ylim(-5, 90)
    ax.set_xlabel(r'Increase of the hydrostatic pressure at the deepest point, $\phi_{\mathrm{f}}$ 0.10 to 0.90 [kPa]')
    ax.set_ylabel(r'Increase of far-field $\sigma_{\mathrm{vM}}$ [MPa]')
    ax.grid(True)
    hs = [Line2D([], [], color=C[a], marker='o', ls='--', lw=0.9) for a in AC]
    ls_ = [L[a] + r', slope $a\,R/t$' for a in AC]
    hs += [Line2D([], [], color='k', marker=MK[l], ls='none', mfc='white' if l in ('Baseline', 'LC1_Maneuver') else 'k')
           for l in LCS]
    ls_ += [LCL[l] for l in LCS]
    legend_below(fig, hs, ls_, 4)
    save(fig, 'fig5_stress_increase_abs', 0.2)



# ---------------------------------------------------------------------------
# One column per load case (October 2026): Baseline, LC3, LC7 shown apart so
# that each load case can be read on its own.
# ---------------------------------------------------------------------------
SPLIT = [('Baseline', '(a)', '(d)'), ('LC3_Combined', '(b)', '(e)'), ('LC7_Emergency', '(c)', '(f)')]


def fig_cause_effect_by_case():
    """Row 1: hydrostatic pressure at the deepest point of the liquid against
    the fill level (the cause). Row 2: far-field von Mises stress against the
    fill level (the effect). Columns: Baseline, LC3, LC7. Shared y-axes per
    row, so the three load cases are compared on the same scale."""
    fig, axs = plt.subplots(2, 3, figsize=(W, 5.0), sharex=True, sharey='row')
    for k, (lc, t1, t2) in enumerate(SPLIT):
        for ac in AC:
            ph = [100.0 * D[ac]['heads']['%s|%.2f' % (lc, f)]['head_bar'] for f in FILLS]
            sg = [pr(ac, lc, f, 'mises_MPa') for f in FILLS]
            axs[0, k].plot(FILLS, ph, '-o', color=C[ac], ms=4.5)
            axs[1, k].plot(FILLS, sg, '-o', color=C[ac], ms=4.5)
        axs[1, k].axhline(ALLOW, **REF)
        axs[1, k].set_xlabel(r'Fill level $\phi_{\mathrm{f}}$ [-]')
        axs[1, k].set_xticks(FILLS)
        axs[1, k].set_xticklabels(['0.10', '0.30', '0.50', '0.75', '0.90'], rotation=0, fontsize=8)
        for ax, t in ((axs[0, k], t1), (axs[1, k], t2)):
            ax.grid(True)
            title(ax, t)
    axs[0, 0].set_ylim(0, 230)
    axs[1, 0].set_ylim(0, 290)
    axs[0, 0].set_ylabel('Hydrostatic pressure at the\ndeepest point, $p_{\\mathrm{h}}$ [kPa]')
    axs[1, 0].set_ylabel(r'Far-field $\sigma_{\mathrm{vM}}$ [MPa]')
    hs, ls_ = ac_handles()
    hs += [Line2D([], [], **REF)]
    ls_ += [r'Allowable $F_{\mathrm{ty}}/1.5$ = 262 MPa']
    legend_below(fig, hs, ls_, 4)
    save(fig, 'fig5_cause_effect_by_case', 0.08)


def fig_stress_vs_pressure_by_case():
    """Far-field von Mises stress against the total pressure at the deepest
    point of the liquid, p + p_h, one panel per load case, with the straight
    line fitted to the twenty cases of each tank (slope a R/t)."""
    fig, axs = plt.subplots(1, 3, figsize=(W, 3.0), sharex=True, sharey=True)
    xx = np.linspace(140, 380, 50)
    for ax, (lc, t, _) in zip(axs, SPLIT):
        for ac in AC:
            Rt = D[ac]['R_over_t']
            a = RC[ac]['peak_fit']['slope_a']; b = RC[ac]['peak_fit']['intercept_b_MPa']
            ax.plot(xx, a * Rt * xx * 1e-3 + b, color=C[ac], lw=0.8, ls='--')
            x = [local_pressure_kPa(ac, lc, f) for f in FILLS]
            y = [pr(ac, lc, f, 'mises_MPa') for f in FILLS]
            ax.plot(x, y, '-o', color=C[ac], ms=4.5, lw=1.2)
        ax.axhline(ALLOW, **REF)
        ax.set_xlim(140, 380)
        ax.set_ylim(0, 290)
        ax.set_xlabel('$p + p_{\\mathrm{h}}$ at the deepest\npoint of the liquid [kPa]')
        ax.grid(True)
        title(ax, t)
    axs[0].set_ylabel(r'Far-field $\sigma_{\mathrm{vM}}$ [MPa]')
    hs = [Line2D([], [], color=C[a], marker='o', lw=1.2) for a in AC]
    ls_ = [L[a] + ', five fill levels' for a in AC]
    hs += [Line2D([], [], color='grey', lw=0.8, ls='--'), Line2D([], [], **REF)]
    ls_ += [r'Fitted line, slope $a\,R/t$', r'Allowable $F_{\mathrm{ty}}/1.5$']
    legend_below(fig, hs, ls_, 3)
    save(fig, 'fig5_stress_vs_pressure_by_case', 0.24)



# ---------------------------------------------------------------------------
# One figure per load case (October 2026): each figure has its own scales so
# that small values remain visible, and the values at the lowest and highest
# fill level are written next to the points so that load cases can still be
# compared across figures. (a) cause, (b) effect, (c) the relation between them.
# ---------------------------------------------------------------------------
CASE_NAME = {'Baseline': 'baseline', 'LC3_Combined': 'lc3', 'LC7_Emergency': 'lc7'}
LETTERS = 'abcdef'


def _label_points(ax, x, y, col, fmt):
    """Write the value of every point next to it: below and to the right of
    the point (under a rising line), the last point above and to the left."""
    n = len(x)
    for k, (xi, yi) in enumerate(zip(x, y)):
        if k < n - 1 and y[k + 1] < yi:          # falling segment: label above
            off, ha, va = (0, 5), 'center', 'bottom'
        elif k < n - 1:
            off, ha, va = (4, -4), 'left', 'top'
        else:
            off, ha, va = (-5, 4), 'right', 'bottom'
        ax.annotate(fmt % yi, (xi, yi), textcoords='offset points', xytext=off,
                    ha=ha, va=va, fontsize=7, color=col)


def _note(ax, text):
    ax.text(0.03, 0.96, text, transform=ax.transAxes, ha='left', va='top', fontsize=7,
            linespacing=1.3, bbox=dict(fc='white', ec='0.7', lw=0.5, pad=2.0))


def _scale(ax, lo, hi, minspan):
    """Own y-scale for each panel: the data fill the lower part, the upper
    part is kept free for the value labels and the note."""
    span = max(hi - lo, minspan)
    mid = 0.5 * (hi + lo)
    lo, hi = mid - span / 2, mid + span / 2
    ax.set_ylim(lo - 0.3 * span, hi + 0.8 * span)


def fill_axis(ax, fs=None):
    """Fill level on the abscissa of every panel, with its label."""
    ax.set_xlabel(r'Fill level $\phi_{\mathrm{f}}$ [-]')
    ax.set_xticks(FILLS)
    ax.set_xticklabels(['0.10', '0.30', '0.50', '0.75', '0.90'], fontsize=fs or 8, rotation=90)
    ax.set_xlim(0.0, 1.0)


def fig_case(lc):
    """One figure per load case. Columns: one aircraft each. Top row: the
    hydrostatic pressure at the deepest point of the liquid, p_h, against the
    fill level (the cause). Bottom row: the far-field von Mises stress against
    the fill level (the effect). Every panel has its own y-scale, so that small
    changes are visible; every point carries its value, and the note gives the
    increase from fill level 0.10 to 0.90 together with the increase predicted
    by the pressure-vessel relation a (R/t) dp_h, so that panels with different
    scales can be compared by their numbers."""
    fig, axs = plt.subplots(2, 3, figsize=(W, 4.6))
    for k, ac in enumerate(AC):
        ph = [100.0 * D[ac]['heads']['%s|%.2f' % (lc, f)]['head_bar'] for f in FILLS]
        sg = [pr(ac, lc, f, 'mises_MPa') for f in FILLS]
        Rt = D[ac]['R_over_t']
        a = RC[ac]['peak_fit']['slope_a']
        dph, dsg = ph[-1] - ph[0], sg[-1] - sg[0]
        top, bot = axs[0, k], axs[1, k]
        top.plot(FILLS, ph, '-o', color=C[ac], ms=4)
        bot.plot(FILLS, sg, '-o', color=C[ac], ms=4)
        _label_points(top, FILLS, ph, C[ac], '%.1f')
        _label_points(bot, FILLS, sg, C[ac], '%.1f')
        _scale(top, min(ph), max(ph), 0.5)
        _scale(bot, min(sg), max(sg), 1.0)
        hd = D[ac]['heads']['%s|0.90' % lc]
        _note(top, 'Rise 0.10 to 0.90: %+.1f kPa' % dph + '\n'
              + 'Depth %.2f m, load factor %.2f' % (hd['depth_m'], hd['n_mag']) + '\n'
              + '%.1f %% of internal pressure' % (100.0 * ph[-1] / 150.0))
        _note(bot, 'Rise 0.10 to 0.90:\n%+.1f MPa (%+.1f %%)' % (dsg, 100.0 * dsg / sg[0]))
        top.set_title('(%s) %s' % (LETTERS[k], L[ac]), loc='left', fontsize=10)
        bot.set_title('(%s) %s' % (LETTERS[k + 3], L[ac]), loc='left', fontsize=10)
        for ax in (top, bot):
            fill_axis(ax)
            ax.grid(True)
        top.set_ylabel('Hydrostatic pressure at the\ndeepest point, $p_{\\mathrm{h}}$ [kPa]')
        bot.set_ylabel('Far-field von Mises\nstress, $\\sigma_{\\mathrm{vM}}$ [MPa]')
    fig.tight_layout(h_pad=1.2, w_pad=1.0)
    tfs.save(fig, os.path.join(FIG, 'fig5_case_%s' % CASE_NAME[lc]))


# ---------------------------------------------------------------------------
# The remaining result graphs in the same panel convention (October 2026)
# ---------------------------------------------------------------------------
LC_SHORT = {'Baseline': 'Baseline', 'LC1_Maneuver': 'Symmetric\nmanoeuvre', 'LC3_Combined': 'Combined\nmanoeuvre',
            'LC7_Emergency': 'Emergency\nlanding'}


def ph_kPa(ac, lc, f):
    return 100.0 * D[ac]['heads']['%s|%.2f' % (lc, f)]['head_bar']


def fig_e190_matrix_panels():
    """E190, one column per load case. Top: far-field von Mises stress against
    the fill level. Bottom: wall temperature difference against the fill level
    with the asymptote q_MLI/h_V."""
    ac = 'E190'
    asym = D[ac]['dT_q_over_hdry_K']
    fig, axs = plt.subplots(4, 2, figsize=(W, 6.9))
    for k, lc in enumerate(LCS):
        sg = series(ac, lc, 'mises_MPa')
        dt = series(ac, lc, 'dT')
        top, bot = axs[k, 0], axs[k, 1]
        top.plot(FILLS, sg, '-o', color=C[ac], ms=4)
        tfs.label_points(top, FILLS, sg, C[ac], '%.1f', fs=7)
        tfs.scale(top, min(sg), max(sg), 1.0, top=1.0)
        tfs.note(top, 'Rise 0.10 to 0.90: %+.1f %%' % (100.0 * (sg[-1] / sg[0] - 1.0)) + '\n'
                 + 'Liquid adds %.1f %% to the pressure at 0.90' % (ph_kPa(ac, lc, 0.9) / 1.5), fs=7)
        bot.axhline(asym, color='k', lw=0.8, ls=':')
        bot.plot(FILLS, dt, '-o', color=C[ac], ms=4)
        tfs.label_points(bot, FILLS, dt, C[ac], '%.3f', fs=7,
                         side=['below', 'below', 'below', 'll', 'below'])
        tfs.scale(bot, min(dt + [asym]), max(dt + [asym]), 0.06, top=1.0, bottom=0.4)
        tfs.note(bot, 'Change 0.10 to 0.90: %+.1f %%' % (100.0 * (dt[-1] / dt[0] - 1.0)) + '\n'
                 + 'Limit $q_{\\mathrm{MLI}}/h_{\\mathrm{V}}$ = %.2f K (dotted)' % asym, fs=7)
        nm = LC_SHORT[lc].replace('\n', ' ')
        top.set_title('(%s) %s: stress' % ('aceg'[k], nm), loc='left', fontsize=10)
        bot.set_title('(%s) %s: wall temperature' % ('bdfh'[k], nm), loc='left', fontsize=10)
        for ax in (top, bot):
            fill_axis(ax)
            ax.grid(True)
            ax.tick_params(labelsize=8)
        top.set_ylabel('Far-field von Mises\nstress, $\\sigma_{\\mathrm{vM}}$ [MPa]')
        bot.set_ylabel('Wall temperature\ndifference, $\\Delta T_{\\mathrm{w}}$ [K]')
    fig.tight_layout(h_pad=0.8, w_pad=1.5)
    tfs.save(fig, os.path.join(FIG, 'fig5_e190_matrix'))


SVP_NAME = {'Baseline': 'baseline', 'LC1_Maneuver': 'symmetric', 'LC3_Combined': 'combined',
            'LC7_Emergency': 'emergency'}


def fig_stress_vs_pressure_case(lc):
    """One figure per load case, one panel per aircraft, each with its own
    scales: far-field von Mises stress of the five fill levels against the
    total pressure difference at the deepest point of the liquid, p + p_h.
    The dashed line of an aircraft is the same in all four figures: the
    straight line a (p + p_h) R/t + b fitted to all twenty cases of that
    aircraft. No notes inside the panels; a and b are given in the legend."""
    fig, axs = plt.subplots(1, 3, figsize=(W, 2.75))
    for c, ac in enumerate(AC):
        Rt = D[ac]['R_over_t']
        a = RC[ac]['peak_fit']['slope_a']; b = RC[ac]['peak_fit']['intercept_b_MPa']
        ax = axs[c]
        x = [local_pressure_kPa(ac, lc, f) for f in FILLS]
        y = [pr(ac, lc, f, 'mises_MPa') for f in FILLS]
        xr = max(x) - min(x)
        span_x = max(xr, 2.0)
        x0, x1 = min(x) - 0.25 * span_x, max(x) + 0.35 * span_x
        xx = np.linspace(x0, x1, 20)
        yl = a * Rt * xx * 1e-3 + b
        ax.plot(xx, yl, color=C[ac], lw=0.9, ls='--', zorder=1)
        ax.plot(x, y, 'o', color=C[ac], ms=4.5, zorder=3)
        # label above a point that lies above the fitted line, below otherwise;
        # closely spaced points alternate above and below
        sides = ['above' if yi >= a * Rt * xi * 1e-3 + b else 'below' for xi, yi in zip(x, y)]
        if xr < 6.0:
            sides = ['above' if k % 2 == 0 else 'below' for k in range(len(x))]
        tfs.label_points(ax, x, y, C[ac], '%.1f', side=sides, fs=7)
        lo, hi = min(y + [yl.min()]), max(y + [yl.max()])
        sp = max(hi - lo, 1.0)
        ax.set_ylim(lo - 0.2 * sp, hi + 0.2 * sp)
        ax.set_xlim(x0, x1)
        ax.set_title('(%s) %s' % ('abc'[c], L[ac]), loc='left', fontsize=10)
        ax.set_xlabel('Total pressure at the deepest\npoint, $p + p_{\\mathrm{h}}$ [kPa]', fontsize=9)
        ax.set_ylabel('Far-field von Mises\nstress, $\\sigma_{\\mathrm{vM}}$ [MPa]', fontsize=9)
        ax.tick_params(labelsize=8)
        ax.grid(True)
    hs = [Line2D([], [], color='0.3', marker='o', ls='none'), Line2D([], [], color='0.3', lw=0.9, ls='--')]
    fit = ', '.join('%s $a$ = %.3f, $b$ = %.1f MPa' % (L[ac], RC[ac]['peak_fit']['slope_a'],
                                                      RC[ac]['peak_fit']['intercept_b_MPa']) for ac in AC)
    legend_below(fig, hs, ['Finite element result at fill levels 0.10, 0.30, 0.50, 0.75 and 0.90',
                           'Line $a\\,(p + p_{\\mathrm{h}})\\,R/t + b$ fitted to all twenty cases of each aircraft:\n'
                           + fit], 1)
    save(fig, 'fig5_svp_%s' % SVP_NAME[lc], 0.24)


def fig_sensitivity_bars():
    """One panel per load case: fill-level sensitivity S_phi of the three tanks
    from the finite element results (solid bars) and the relative increase of
    the local pressure, (p_h,0.90 - p_h,0.10)/(p + p_h,0.10) (hatched bars),
    which the pressure-vessel relation predicts it to equal."""
    fig, axs = plt.subplots(2, 2, figsize=(W, 4.6))
    xk = np.arange(3)
    for k, lc in enumerate(LCS):
        ax = axs.flat[k]
        fe = [100.0 * (pr(ac, lc, 0.9, 'mises_MPa') / pr(ac, lc, 0.1, 'mises_MPa') - 1.0) for ac in AC]
        pd = [100.0 * (ph_kPa(ac, lc, 0.9) - ph_kPa(ac, lc, 0.1)) / (150.0 + ph_kPa(ac, lc, 0.1)) for ac in AC]
        cols = [C[ac] for ac in AC]
        b1 = ax.bar(xk - 0.19, fe, 0.36, color=cols, edgecolor=cols)
        b2 = ax.bar(xk + 0.19, pd, 0.36, color='white', edgecolor=cols, hatch='////', lw=0.8)
        top = max(fe + pd)
        bot = min(0.0, min(fe))
        for bars, vals in ((b1, fe), (b2, pd)):
            for r_, v in zip(bars, vals):
                ax.annotate('%.1f' % v, (r_.get_x() + r_.get_width() / 2, max(v, 0.0)),
                            textcoords='offset points', xytext=(0, 2), ha='center', va='bottom',
                            fontsize=7.5)
        ax.axhline(0, color='k', lw=0.6)
        ax.set_ylim(bot - 0.08 * top, 1.22 * top)
        ax.set_xticks(xk)
        ax.set_xticklabels(['E190', 'A320', 'A350-900'], fontsize=8.5)
        ax.set_title('(%s) %s' % ('abcd'[k], LC_SHORT[lc].replace('\n', ' ')), loc='left', fontsize=10)
        ax.grid(True, axis='y')
        ax.tick_params(axis='y', labelsize=8)
        ax.set_xlabel('Aircraft')
        ax.set_ylabel('Stress rise from fill 0.10\nto 0.90, $S_{\\phi}$ [%]')
    hs = [plt.Rectangle((0, 0), 1, 1, fc='0.35', ec='0.35'),
          plt.Rectangle((0, 0), 1, 1, fc='white', ec='0.35', hatch='////')]
    legend_below(fig, hs, ['Stress rise computed by the finite element model',
                           'Pressure rise at the deepest point, $\\Delta p_{\\mathrm{h}}/(p + p_{\\mathrm{h},0.10})$'], 1)
    save(fig, 'fig5_sensitivity', 0.08)


def fig_thermal_panels():
    """One panel per aircraft: wall temperature difference against the fill
    level for the Baseline condition and LC7, with the asymptote q_MLI/h_V."""
    fig, axs = plt.subplots(1, 3, figsize=(W, 2.9))
    for k, ac in enumerate(AC):
        ax = axs[k]
        asym = D[ac]['dT_q_over_hdry_K']
        b = series(ac, 'Baseline', 'dT'); e = series(ac, 'LC7_Emergency', 'dT')
        ax.axhline(asym, color='k', lw=0.8, ls=':')
        ax.plot(FILLS, e, '--D', color=C[ac], ms=3.5, lw=1.0)
        ax.plot(FILLS, b, '-o', color=C[ac], ms=4)
        tfs.label_points(ax, FILLS, b, C[ac], '%.2f', side=['below', 'below', 'below', 'll', 'below'], fs=6.5)
        tfs.label_points(ax, [FILLS[-1]], [e[-1]], C[ac], '%.2f', side='above', fs=6.5)
        lo, hi = min(b), max(e + [asym])
        span = hi - lo
        ax.set_ylim(lo - 0.35 * span, hi + 0.9 * span)
        tfs.note(ax, 'Insulation heat flux %.2f W/m$^2$' % D[ac]['q_Wm2'] + '\n'
                 + 'Limit $q_{\\mathrm{MLI}}/h_{\\mathrm{V}}$ = %.2f K' % asym + '\n'
                 + 'Baseline change 0.10 to 0.90: %+.1f %%' % (100.0 * (b[-1] / b[0] - 1.0)), fs=6.5)
        ax.set_title('(%s) %s' % ('abc'[k], L[ac]), loc='left', fontsize=10)
        fill_axis(ax)
        ax.set_ylabel('Wall temperature\ndifference, $\\Delta T_{\\mathrm{w}}$ [K]')
        ax.grid(True)
    hs = [Line2D([], [], color='0.3', marker='o', ls='-'), Line2D([], [], color='0.3', marker='D', ls='--', ms=3.5),
          Line2D([], [], color='k', lw=0.8, ls=':')]
    legend_below(fig, hs, ['Baseline condition', 'Emergency landing',
                           r'Asymptote $q_{\mathrm{MLI}}/h_{\mathrm{V}}$'], 3)
    save(fig, 'fig5_thermal', 0.1)


def fig_margins_panels():
    """LC7, one column per aircraft. Top: yield reserve factor at 1.5 bar.
    Bottom: minimum axial membrane stress of the far-field barrel at 0.187 bar,
    with the axial tension pR/2t and the allowable compression -sigma_cr/1.5."""
    lc = 'LC7_Emergency'
    fig, axs = plt.subplots(2, 3, figsize=(W, 4.6))
    for k, ac in enumerate(AC):
        top, bot = axs[0, k], axs[1, k]
        rf = series(ac, lc, 'fos')
        sx = [D[ac]['stability']['%s|%.2f' % (lc, f)]['axial_mem_min_MPa'] for f in FILLS]
        allow = -D[ac]['sigma_cr_MPa'] / 1.5
        tens = D[ac]['axial_pR_2t_min_MPa']
        top.axhline(1.0, color='k', lw=0.8, ls=':')
        top.plot(FILLS, rf, '-o', color=C[ac], ms=4)
        tfs.label_points(top, FILLS, rf, C[ac], '%.3f', side='ur', fs=6.5)
        tfs.scale(top, min(rf + [1.0]), max(rf), 0.3, top=0.55, bottom=0.15)
        tfs.note(top, 'At fill 0.90: %.3f' % rf[-1], loc='ur', fs=6.5)
        bot.axhline(0.0, color='k', lw=0.6)
        bot.axhline(allow, color='k', lw=0.8, ls='-.')
        bot.axhline(tens, color='0.45', lw=0.8, ls='--')
        bot.plot(FILLS, sx, '-o', color=C[ac], ms=4)
        tfs.label_points(bot, FILLS, sx, C[ac], '%.2f', side='auto', fs=6.5)
        lo, hi = min(sx + [allow]), max(sx + [tens])
        span = hi - lo
        bot.set_ylim(lo - 0.2 * span, hi + 0.75 * span)
        rfb = abs(allow / sx[-1]) if sx[-1] < 0 else None
        tfs.note(bot, 'Dashed: %+.1f MPa' % tens + '\n' + 'Dash-dot: %.1f MPa' % allow + '\n'
                 + ('Buckling RF at 0.90: %.3f' % rfb if rfb else 'In tension: no buckling'),
                 loc='ur', fs=6.5)
        top.set_title('(%s) %s' % ('abc'[k], L[ac]), loc='left', fontsize=10)
        bot.set_title('(%s) %s' % ('def'[k], L[ac]), loc='left', fontsize=10)
        for ax in (top, bot):
            fill_axis(ax)
            ax.grid(True)
        top.set_ylabel('Yield reserve factor\nat 1.5 bar, $\\mathrm{RF}_{\\mathrm{y}}$ [-]')
        bot.set_ylabel('Minimum axial stress at\n0.187 bar, $\\sigma_{x,\\mathrm{min}}$ [MPa]')
    hs = [Line2D([], [], color='k', lw=0.8, ls=':'), Line2D([], [], color='0.45', lw=0.8, ls='--'),
          Line2D([], [], color='k', lw=0.8, ls='-.')]
    legend_below(fig, hs, ['Requirement, reserve factor = 1 (a to c)',
                           'Axial tension from pressure alone, $pR/2t$ (d to f)',
                           'Allowable compression, $-\\sigma_{\\mathrm{cr}}/1.5$ (d to f)'], 1)
    save(fig, 'fig5_margins', 0.1)


if __name__ == '__main__':
    # Figures used in the report (October 2026 panel convention). The older
    # combined figures (fig_e190_matrix, fig_head_mechanism, fig_stress_vs_head,
    # fig_three_aircraft, fig_sensitivity_vs_tilt, fig_thermal, fig_margins)
    # remain above for reference but are no longer called.
    for f in (fig_sizing_history, fig_e190_matrix_panels,
              fig_sensitivity_bars, fig_thermal_panels, fig_margins_panels, fig_convergence):
        f()
    for lc in ('Baseline', 'LC3_Combined', 'LC7_Emergency'):
        fig_case(lc)
    for lc in LCS:
        fig_stress_vs_pressure_case(lc)
