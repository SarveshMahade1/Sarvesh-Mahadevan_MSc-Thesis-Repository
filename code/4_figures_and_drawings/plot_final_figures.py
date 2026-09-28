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

    python python_scripts/final_numbers.py
    python python_scripts/hand_calcs.py
    python python_scripts/plot_final_figures.py
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
LC_STYLE = {'Baseline': dict(ls='-', marker='o'), 'LC1_Maneuver': dict(ls='--', marker='s'),
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
    ax.set_title(s, loc='left', fontsize=10)


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
    fig, axs = plt.subplots(2, 3, figsize=(W, 4.3), sharex='col', gridspec_kw=dict(height_ratios=[1.35, 1]))
    for k, ac in enumerate(AC):
        rr = [r for r in rows[ac] if r['stage'].startswith('iteration')]
        it = lambda r: int(r['stage'].split(';')[0].split()[1])
        y = [(it(r), r['fos']) for r in rr if r['p_gauge_bar'] > 1.0]
        b = [(it(r), r['buckling_fos']) for r in rr if r['p_gauge_bar'] < 1.0 and math.isfinite(r['buckling_fos'])]
        t = [(it(r), r['t_wall_mm']) for r in rr if r['p_gauge_bar'] > 1.0]
        ax = axs[0, k]
        ax.plot([i for i, _ in y], [v for _, v in y], '-o', color='#1F3864')
        if b:
            ax.plot([i for i, _ in b], [v for _, v in b], '--s', color='#C00000')
        ax.axhline(1.0, color='#1F3864', lw=0.9, ls=':')
        ax.axhline(1.5, color='#C00000', lw=0.9, ls=':')
        ax.set_ylim(0.5, 4.2)
        ax.grid(True)
        title(ax, '(%s) %s' % ('abc'[k], L[ac]))
        ax2 = axs[1, k]
        ax2.plot([i for i, _ in t], [v for _, v in t], '-^', color='#555555')
        ax2.set_xlabel('Iteration')
        ax2.set_xticks(sorted(set(i for i, _ in t)))
        ax2.grid(True)
    axs[0, 0].set_ylabel('Factor of safety [-]')
    axs[1, 0].set_ylabel('Wall $t$ [mm]')
    hs = [Line2D([], [], color='#1F3864', marker='o'), Line2D([], [], color='#C00000', ls='--', marker='s'),
          Line2D([], [], color='#555555', marker='^'), Line2D([], [], color='#1F3864', ls=':'),
          Line2D([], [], color='#C00000', ls=':')]
    legend_below(fig, hs, ['Yield FOS at 1.5 bar', 'Buckling FOS at 0.187 bar', 'Trial wall thickness',
                           'Yield target 1.0', 'Buckling target 1.5'], 3)
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
    axs[0].set_ylabel(r'Far-field $\sigma_{vM}$ [MPa]')
    axs[1].axhline(H['E190']['thermal']['dT_asymptote_K'], color='#555555', lw=1.0, ls=(0, (5, 2)))
    axs[1].set_ylabel(r'$\Delta T$ [K]')
    axs[1].set_ylim(3.35, 3.78)
    for ax, s in zip(axs, ['(a) Stress', '(b) Wetted-to-dry temperature difference']):
        ax.set_xlabel(r'Fill ratio $\phi_f$ [-]')
        ax.set_xticks(FILLS)
        ax.grid(True)
        title(ax, s)
    hs = [Line2D([], [], color=LC_COL[lc], **LC_STYLE[lc]) for lc in LCS] + \
         [Line2D([], [], **REF), Line2D([], [], color='#555555', lw=1.0, ls=(0, (5, 2)))]
    legend_below(fig, hs, [LCL[lc] for lc in LCS] + [r'Allowable 262 MPa (a)', r'Asymptote $q_{\mathrm{MLI}}/h_{\mathrm{dry}}$ (b)'], 3)
    save(fig, 'fig5_e190_matrix', 0.2)


def fig_head_mechanism():
    fig, axs = plt.subplots(1, 2, figsize=(W, 2.9))
    for ac in AC:
        for lc, ls in (('Baseline', '--'), ('LC3_Combined', '-.'), ('LC7_Emergency', '-')):
            dep = [D[ac]['heads']['%s|%.2f' % (lc, f)]['depth_m'] for f in FILLS]
            hp = [100 * D[ac]['heads']['%s|%.2f' % (lc, f)]['head_over_p'] for f in FILLS]
            kw = dict(color=C[ac], ls=ls, marker='o', ms=3.5)
            axs[0].plot(FILLS, dep, **kw)
            axs[1].plot(FILLS, hp, **kw)
    axs[0].set_ylabel(r'Liquid depth along $\hat{\mathbf{g}}$ [m]')
    axs[1].set_ylabel(r'Head $p_h$ / gauge pressure [%]')
    axs[1].set_yscale('log')
    for ax, s in zip(axs, ['(a) Depth of liquid', '(b) Hydrostatic head']):
        ax.set_xlabel(r'Fill ratio $\phi_f$ [-]')
        ax.set_xticks(FILLS)
        ax.grid(True, which='both')
        title(ax, s)
    hs, ls_ = ac_handles()
    hs += [Line2D([], [], color='k', ls=l) for l in ('--', '-.', '-')]
    ls_ += ['Baseline and LC1 (vertical)', r'LC3 Combined (31.0$^\circ$)', r'LC7 Emergency (57.7$^\circ$)']
    legend_below(fig, hs, ls_, 3)
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
    ax.set_xlabel(r'$1 + p_h/p_{\mathrm{gauge}}$ at the deepest wetted point [-]')
    ax.set_ylabel(r'Far-field $\sigma_{vM}\,/\,(pR/t)$ [-]')
    ax.grid(True)
    hs, ls_ = ac_handles()
    hs = [Line2D([], [], color=C[a], marker='o', ls='none') for a in AC]
    hs += [Line2D([], [], color='k', marker=mk[l], ls='none', mfc='none' if l in ('Baseline', 'LC1_Maneuver') else 'k')
           for l in LCS] + [Line2D([], [], **REF)]
    ls_ += [LCL[l] for l in LCS] + [r'Barrel membrane, $0.866\,(1 + p_h/p)$']
    legend_below(fig, hs, ls_, 4)
    save(fig, 'fig5_stress_vs_head', 0.16)


def fig_sensitivity_vs_tilt():
    tilt = {'Baseline': 0.0, 'LC1_Maneuver': 0.0, 'LC3_Combined': 30.96, 'LC7_Emergency': 57.69}
    fig, axs = plt.subplots(1, 2, figsize=(W, 2.9))
    for ac in AC:
        s = D[ac]['fill_sensitivity_pct']
        axs[0].plot([tilt[l] for l in ('Baseline', 'LC3_Combined', 'LC7_Emergency')],
                    [s[l] for l in ('Baseline', 'LC3_Combined', 'LC7_Emergency')], '-o', color=C[ac])
        axs[0].plot([tilt['LC1_Maneuver']], [s['LC1_Maneuver']], 's', mfc='none', color=C[ac], ms=7)
        axs[1].plot([D[ac]['slenderness_L_2R']], [s['LC7_Emergency']], 'D', color=C[ac], ms=7)
    axs[0].set_xlabel(r'Tilt of the resultant from vertical [$^\circ$]')
    axs[0].set_ylabel(r'Change in $\sigma_{vM}$, $\phi_f$ 0.1 to 0.9 [%]')
    axs[1].set_xlabel(r'Barrel slenderness $L/2R$ [-]')
    axs[1].set_ylabel(r'Change under LC7 [%]')
    axs[1].set_xlim(2.4, 7.2)
    axs[1].set_ylim(0, 115)
    for ax, s in zip(axs, ['(a) Against load direction', '(b) Against slenderness, LC7']):
        ax.grid(True)
        title(ax, s)
    hs, ls_ = ac_handles()
    hs += [Line2D([], [], color='k', marker='o'), Line2D([], [], color='k', marker='s', mfc='none', ls='none'),
           Line2D([], [], color='k', marker='D', ls='none')]
    ls_ += ['Baseline, LC3, LC7 (a)', 'LC1 Manoeuvre, 2.5 g (a)', 'LC7 Emergency (b)']
    legend_below(fig, hs, ls_, 3)
    save(fig, 'fig5_sensitivity', 0.2)


def fig_three_aircraft():
    fig, axs = plt.subplots(1, 2, figsize=(W, 2.9), sharey=True)
    for ax, lc in zip(axs, ('Baseline', 'LC7_Emergency')):
        for ac in AC:
            ax.plot(FILLS, [pr(ac, lc, f, 'mises_MPa') / ALLOW for f in FILLS], '-o', color=C[ac])
        ax.axhline(1.0, **REF)
        title(ax, '(%s) %s' % ('a' if lc == 'Baseline' else 'b', LCL[lc]))
        ax.set_xlabel(r'Fill ratio $\phi_f$ [-]')
        ax.set_xticks(FILLS)
        ax.grid(True)
    axs[0].set_ylabel(r'$\sigma_{vM}\,/\,\sigma_{\mathrm{allow}}$ [-]')
    axs[0].set_ylim(0.0, 1.1)
    hs, ls_ = ac_handles()
    legend_below(fig, hs + [Line2D([], [], **REF)], ls_ + ['Allowable, 262 MPa'], 4)
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
    axs[0].set_xlabel(r'Fill ratio $\phi_f$ [-]')
    axs[0].set_ylabel(r'$\Delta T$ [K]')
    axs[0].set_xticks(FILLS)
    axs[1].set_xlabel(r'Insulation flux $q_{\mathrm{MLI}}$ [W/m$^2$]')
    axs[1].set_ylabel(r'$\Delta T$ at $\phi_f = 0.1$ [K]')
    for ax, s in zip(axs, ['(a) Against fill', '(b) Against insulation flux']):
        ax.grid(True)
        title(ax, s)
    hs, ls_ = ac_handles()
    hs += [Line2D([], [], color='k', ls='-', marker='o', ms=4), Line2D([], [], color='k', ls='--', marker='D', ms=4),
           Line2D([], [], color='k', lw=0.9, ls=':')]
    ls_ += ['Baseline', 'LC7 Emergency', r'$\Delta T = q_{\mathrm{MLI}}/h_{\mathrm{dry}}$']
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
        ax.set_xlabel(r'Fill ratio $\phi_f$ [-]')
        ax.set_xticks([0.1, 0.5, 0.9])
        ax.grid(True)
        ax.set_ylim(-100, 90)
    axs[0].set_ylabel('Mid-ring share of reaction [%]')
    hs = [Line2D([], [], color='k', ls='--', marker='o', ms=4), Line2D([], [], color='k', ls='-', marker='o', ms=4),
          Line2D([], [], **REF)]
    legend_below(fig, hs, ['Baseline (1.0 g)', 'LC1 Manoeuvre (2.5 g)', 'Uniform beam on three supports, 62.5 %'], 3)
    save(fig, 'fig5_ring_share', 0.13)


def fig_margins():
    fig, axs = plt.subplots(1, 2, figsize=(W, 3.0))
    for ac in AC:
        axs[0].plot(FILLS, series(ac, 'LC7_Emergency', 'fos'), '-o', color=C[ac])
        axs[1].plot(FILLS, series(ac, 'LC7_Emergency', 'axial_mem_min_MPa', 'stability'), '-D', color=C[ac])
        axs[1].axhline(-D[ac]['sigma_cr_MPa'] / 1.5, color=C[ac], lw=1.0, ls=':')
    axs[0].axhline(1.0, **REF)
    axs[0].set_ylabel('Yield factor of safety [-]')
    axs[0].set_ylim(0.8, 4.0)
    axs[1].axhline(0.0, color='k', lw=0.6)
    axs[1].set_ylabel('Min. axial membrane stress [MPa]')
    for ax, s in zip(axs, ['(a) Strength, 1.5 bar', '(b) Stability, 0.187 bar']):
        ax.set_xlabel(r'Fill ratio $\phi_f$ [-]')
        ax.set_xticks(FILLS)
        ax.grid(True)
        title(ax, s)
    hs, ls_ = ac_handles()
    hs += [Line2D([], [], **REF), Line2D([], [], color='#555555', lw=1.0, ls=':')]
    ls_ += ['Yield target 1.0 (a)', r'Allowable compression $-\sigma_{cr}/1.5$, colour of each aircraft (b)']
    legend_below(fig, hs, ls_, 3)
    save(fig, 'fig5_margins', 0.2)


def fig_convergence():
    fig, axs = plt.subplots(1, 2, figsize=(W, 2.9))
    for ac in AC:
        cd, cm = D[ac]['convergence']['pdesign'], D[ac]['convergence']['pmin']
        n = [c['n_el'] for c in cd]
        axs[0].plot(n, [100 * (c['mises_MPa'] / cd[-1]['mises_MPa'] - 1) for c in cd], '-o', color=C[ac])
        axs[0].plot([n[1]], [100 * (cd[1]['mises_MPa'] / cd[-1]['mises_MPa'] - 1)], 'o', ms=10, mfc='none', color=C[ac])
        if cm[-1]['axial_mem_min_MPa'] < 0:
            axs[1].plot(n, [c['buckling_fos'] for c in cm], '-D', color=C[ac])
            axs[1].plot([n[1]], [cm[1]['buckling_fos']], 'D', ms=10, mfc='none', color=C[ac])
    axs[0].axhline(0, color='k', lw=0.6)
    axs[0].set_ylabel(r'$\sigma_{vM}$ relative to finest [%]')
    axs[1].axhline(1.5, **REF)
    axs[1].set_ylabel('Buckling factor of safety [-]')
    axs[1].set_ylim(1.4, 1.75)
    for ax, s in zip(axs, ['(a) Strength, 1.5 bar', '(b) Stability, 0.187 bar']):
        ax.set_xscale('log')
        ax.set_xlabel('Elements [-]')
        ax.grid(True, which='both')
        title(ax, s)
    hs, ls_ = ac_handles()
    hs += [Line2D([], [], color='k', marker='o', ms=10, mfc='none', ls='none'), Line2D([], [], **REF)]
    ls_ += ['Production mesh', 'Buckling target 1.5 (b)']
    legend_below(fig, hs, ls_, 5)
    save(fig, 'fig5_convergence', 0.13)


def fig_ovality():
    fig, ax = plt.subplots(figsize=(W, 2.9))
    for ac in AC:
        ax.plot(FILLS, series(ac, 'Baseline', 'ovality_mm'), '--o', color=C[ac], ms=4)
        ax.plot(FILLS, series(ac, 'LC7_Emergency', 'ovality_mm'), '-D', color=C[ac], ms=4)
    ax.set_xlabel(r'Fill ratio $\phi_f$ [-]')
    ax.set_ylabel('Section distortion [mm]')
    ax.set_xticks(FILLS)
    ax.grid(True)
    hs, ls_ = ac_handles()
    hs += [Line2D([], [], color='k', ls='--', marker='o', ms=4), Line2D([], [], color='k', ls='-', marker='D', ms=4)]
    ls_ += ['Baseline', 'LC7 Emergency']
    legend_below(fig, hs, ls_, 5)
    save(fig, 'fig5_ovality', 0.13)


if __name__ == '__main__':
    for f in (fig_sizing_history, fig_e190_matrix, fig_head_mechanism, fig_stress_vs_head,
              fig_sensitivity_vs_tilt, fig_three_aircraft, fig_thermal, fig_reactions,
              fig_margins, fig_convergence, fig_ovality):
        f()
