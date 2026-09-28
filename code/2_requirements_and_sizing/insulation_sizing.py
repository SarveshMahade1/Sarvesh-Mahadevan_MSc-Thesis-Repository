# -*- coding: utf-8 -*-
"""
INSULATION ARCHITECTURE: WHAT EACH OPTION WOULD NEED (28 Sep 2026)
===================================================================

Evaluates the three insulation architectures of Chapter 3 against the same
heat-leak budget, so that the selection rests on numbers from one script.

  1. Vacuum + MLI : the resistance network of boiloff_holdtime_sensitivity.py
                    (the same function, imported, so the numbers match Figure
                    3.1 exactly).
  2. Foam only    : the same cylindrical network with the MLI layer replaced
                    by polyurethane foam, k = 0.0163 W/m-K (Bagarello et al.,
                    Table 1, the foam of the benchmark tank).
  3. Vacuum only  : radiation between the inner wall and the jacket, grey
                    parallel surfaces, q = sigma (T_h^4 - T_c^4) eps_eff,
                    eps_eff = 1 / (1/eps1 + 1/eps2 - 1). The script returns the
                    emissivity both surfaces would need to meet the budget.

    python python_scripts/insulation_sizing.py
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
import boiloff_holdtime_sensitivity as bh                     # noqa: E402

SIGMA = 5.670374419e-8
K_FOAM = 0.0163          # W/m-K, PUR-96, Bagarello et al. Table 1
GEO = {'E190': (1.3015, 13.9675, 702.6), 'A320': (1.7080, 10.0120, 919.3),
       'A350': (2.5771, 33.8240, 6583.2)}
T_H, T_C = bh.T_AMBIENT, bh.T_LIQUID


def foam_thickness(R, Lt, Q):
    """Cylindrical conduction through foam plus the same outer film."""
    t = 0.5
    for _ in range(200):
        r1 = R + bh.T_WALL_INNER / 2.0
        r2 = r1 + t
        R_conv = 1.0 / (bh.H_CONV * 2 * math.pi * r2 * Lt)
        R_need = (T_H - T_C) / Q - R_conv
        r2n = r1 * math.exp(R_need * 2 * math.pi * K_FOAM * Lt)
        tn = r2n - r1
        if abs(tn - t) < 1e-9:
            break
        t = tn
    return t


out = {}
for ac, (R, L, Q) in GEO.items():
    h = R / 1.6
    Lt = L + 2 * h
    mli_t, mli_m = bh.required_mli_thickness(R, Lt, Q)
    A = 2 * math.pi * R * L + 2 * math.pi * R * R * (1 + (1 - (1 - (h / R) ** 2)) / math.sqrt(1 - (h / R) ** 2)
                                                     * math.atanh(math.sqrt(1 - (h / R) ** 2)))
    q_need = Q / A
    eps_eff = q_need / (SIGMA * (T_H ** 4 - T_C ** 4))
    eps_each = 2.0 / (1.0 / eps_eff + 1.0)        # equal emissivity on both faces
    out[ac] = {'Q_W': Q, 'area_m2': A, 'q_budget_Wm2': q_need,
               'mli_mm': mli_t * 1e3, 'mli_mass_kg': mli_m,
               'foam_m': foam_thickness(R, Lt, Q), 'foam_over_R': foam_thickness(R, Lt, Q) / R,
               'blackbody_Wm2': SIGMA * (T_H ** 4 - T_C ** 4),
               'eps_eff_needed': eps_eff, 'eps_each_needed': eps_each}

with open(os.path.join(RP.PROCESSED, 'insulation.json'), 'w') as fh:
    json.dump(out, fh, indent=1)
for ac, v in out.items():
    print('%s  q=%.2f W/m2  MLI %.2f mm (%.0f kg)  foam %.3f m (%.2f R)  vacuum-only needs eps %.4f on each face'
          % (ac, v['q_budget_Wm2'], v['mli_mm'], v['mli_mass_kg'], v['foam_m'], v['foam_over_R'], v['eps_each_needed']))
