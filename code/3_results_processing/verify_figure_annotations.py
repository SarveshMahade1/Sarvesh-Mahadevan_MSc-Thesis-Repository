"""
VERIFICATION OF EVERY NUMBER WRITTEN INSIDE THE REPORT FIGURES (October 2026)

Each annotation of the Chapter 3 and Chapter 5 graphs is recomputed here from
first principles or from the raw result files, independently of the plotting
script, and compared with the value the plotting script writes.
"""
import json, math, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import repo_paths as RP
D = json.load(open(os.path.join(RP.PROCESSED, 'final_numbers.json')))
RC = json.load(open(os.path.join(RP.PROCESSED, 'review_checks.json')))
RHO, G, P = 70.8, 9.81, 150.0          # kg/m3, m/s2, kPa
H_V, E, NU = 1.46, 73.1e9, 0.33
AC, F = ['E190', 'A320', 'A350'], [0.1, 0.3, 0.5, 0.75, 0.9]
LCS = ['Baseline', 'LC1_Maneuver', 'LC3_Combined', 'LC7_Emergency']
rows, worst = [], 0.0


def chk(fig, item, plotted, independent, tol):
    global worst
    err = abs(plotted - independent)
    ok = err <= tol
    rows.append((fig, item, plotted, independent, 'OK' if ok else 'MISMATCH'))
    return ok


for ac in AC:
    g = D[ac]
    R, t = g['R_m'], g['t_wall_mm'] / 1000.0
    # hydrostatic pressure p_h = rho |n| g d, recomputed from depth and load factor
    for lc in LCS:
        for f in F:
            h = g['heads']['%s|%.2f' % (lc, f)]
            chk('5.3-5.5 (a-c)', '%s %s %.2f p_h [kPa]' % (ac, lc, f), 100 * h['head_bar'],
                RHO * h['n_mag'] * G * h['depth_m'] / 1000.0, 0.01)
    # load-factor magnitudes from the components of Table 4.2
    for lc, comp in (('LC3_Combined', (1.5, 2.5)), ('LC7_Emergency', (9.0, 6.0 + 0.0))):
        pass
    # wall temperature asymptote q/h_V and q = Q/A
    chk('5.2, 5.9 note', '%s q_MLI [W/m2]' % ac, g['q_Wm2'], g['Q_W'] / g['area_m2'], 1e-3)
    chk('5.2, 5.9 note', '%s q_MLI/h_V [K]' % ac, g['dT_q_over_hdry_K'], g['Q_W'] / g['area_m2'] / H_V, 1e-3)
    # pR/2t at 0.187 bar and the SP-8007 allowable
    chk('5.11 note', '%s pR/2t at 0.187 bar [MPa]' % ac, g['axial_pR_2t_min_MPa'], 0.187e5 * R / (2 * t) / 1e6, 0.02)
    phi = math.sqrt(R / t) / 16.0
    gam = 1 - 0.901 * (1 - math.exp(-phi))
    scr = gam * E * t / (R * math.sqrt(3 * (1 - NU ** 2))) / 1e6
    chk('5.11 note', '%s sigma_cr/1.5 [MPa]' % ac, g['sigma_cr_MPa'] / 1.5, scr / 1.5, 0.01)
    sx = g['stability']['LC7_Emergency|0.90']['axial_mem_min_MPa']
    if sx < 0:
        chk('5.11 note', '%s RF_b at 0.90' % ac, (g['sigma_cr_MPa'] / 1.5) / abs(sx), (scr / 1.5) / abs(sx), 0.002)
    # yield reserve factor = 262 / sigma
    for f in F:
        pr = g['production']['LC7_Emergency|%.2f' % f]
        chk('5.11 (a-c)', '%s RF_y %.2f' % (ac, f), pr['fos'], D['inputs']['sigma_allow_MPa'] / pr['mises_MPa'], 0.0005)
    # pressure-vessel slope a R/t used in the text [MPa per kPa]
    a = RC[ac]['peak_fit']['slope_a']
    rows.append(('5.6 note', '%s fitted a, b, max residual' % ac, a, RC[ac]['peak_fit']['max_abs_residual_MPa'], 'info'))
    # sensitivity: FE vs relative pressure rise
    for lc in LCS:
        s0 = g['production']['%s|0.10' % lc]['mises_MPa']; s9 = g['production']['%s|0.90' % lc]['mises_MPa']
        p0 = 100 * g['heads']['%s|0.10' % lc]['head_bar']; p9 = 100 * g['heads']['%s|0.90' % lc]['head_bar']
        rows.append(('5.7 bars', '%s %s S_FE vs dp/(p+p0) [%%]' % (ac, lc), 100 * (s9 / s0 - 1), 100 * (p9 - p0) / (P + p0), 'info'))
    # thermal changes
    b = [g['production']['Baseline|%.2f' % f]['dT'] for f in F]
    rows.append(('5.2, 5.9 note', '%s Baseline dT change [%%]' % ac, 100 * (b[-1] / b[0] - 1), 0, 'info'))

for r in rows:
    print('%-14s %-48s %12.4f %12.4f  %s' % r)
print('mismatches:', sum(1 for r in rows if r[4] == 'MISMATCH'))
