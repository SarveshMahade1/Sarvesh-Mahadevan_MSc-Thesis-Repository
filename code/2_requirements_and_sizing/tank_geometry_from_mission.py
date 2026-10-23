# -*- coding: utf-8 -*-
"""
TANK GEOMETRY FROM THE MISSION (thesis Section 3.1, Table 3.2)

The first link of the workflow: from the kerosene load of a reference aircraft
to the inner-vessel geometry that the Abaqus model builds.

    m_LH2 = m_JetA * LHV_JetA / LHV_LH2              energy equivalence
    V     = m_LH2 / (rho_LH2 * phi_mission)          volume at the mission fill
    R     = 0.8648 * D_fuselage / 2                  size-class scaling
    h_cap = R / (a/b)                                semi-ellipsoidal caps
    L     = (V - 4/3 pi R^2 h_cap) / (pi R^2)        the barrel length follows

All inputs come from code/design_basis.json. For each aircraft the script prints
the chain, compares the derived geometry with the geometry stored in the design
basis (the thesis values), and writes results/processed/tank_geometry.json.

To add an aircraft: add a record with jet_a_mass_kg and fuselage_diameter_m
(R_mid_m and L_barrel_m null) to design_basis.json and run this script.

    python code/2_requirements_and_sizing/tank_geometry_from_mission.py
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import repo_paths as RP          # noqa: E402
import design_basis as DB        # noqa: E402


def main():
    out = {}
    print('Design basis: %s\n' % DB.PATH)
    for ac in DB.AIRCRAFT:
        r, mg, g = DB.record(ac), DB.mission_geometry(ac), DB.geometry(ac)
        print('%s' % DB.label(ac))
        print('  kerosene replaced        %10.0f kg' % r['jet_a_mass_kg'])
        print('  hydrogen, equal energy   %10.1f kg' % mg['m_lh2_kg'])
        print('  volume at fill %.2f      %10.2f m3' % (DB.MISSION['mission_fill_level'], mg['V_required_m3']))
        print('  radius from fuselage     %10.4f m   (D = %.2f m)' % (mg['R'], r['fuselage_diameter_m']))
        print('  barrel length (derived)  %10.4f m' % mg['L'])
        print('  model geometry           R %.4f m, L %.4f m  (%s)' % (g['R'], g['L'], g['source']))
        dl = (g['L'] - mg['L']) * 1000.0
        print('  derived minus model      L %+.1f mm' % (-dl))
        print('  total length %.2f m, L/2R %.2f, V %.2f m3, area %.1f m2%s\n'
              % (g['L_total'], g['slenderness'], g['V'], g['A'],
                 (', %.0f %% of the fuselage' % (100 * g['share_of_fuselage'])) if 'share_of_fuselage' in g else ''))
        out[ac] = dict(label=DB.label(ac), m_lh2_kg=round(mg['m_lh2_kg'], 1),
                       V_required_m3=round(mg['V_required_m3'], 3), R_derived_m=mg['R'], L_derived_m=mg['L'],
                       R_model_m=g['R'], L_model_m=g['L'], h_cap_m=round(g['h_cap'], 4),
                       L_total_m=round(g['L_total'], 3), V_model_m3=round(g['V'], 3), A_m2=round(g['A'], 2),
                       slenderness=round(g['slenderness'], 3), geometry_source=g['source'])
    if not os.path.isdir(RP.PROCESSED):
        os.makedirs(RP.PROCESSED)
    path = os.path.join(RP.PROCESSED, 'tank_geometry.json')
    with open(path, 'w') as f:
        json.dump(out, f, indent=1)
    print('written %s' % path)


if __name__ == '__main__':
    main()
