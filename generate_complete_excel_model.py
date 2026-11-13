import os
os.environ['PANDAS_USE_PYARROW'] = '0'

import math
import numpy as np
import pandas as pd
import openpyxl
from openpyxl.styles import Font, PatternFill, Alignment, Border, Side
from openpyxl.utils import get_column_letter
from openpyxl.chart import LineChart, Reference

# =====================================================================
# 1. PHYSICAL CONSTANTS & MOLAR MASSES
# =====================================================================
R = 8.314       # J / (mol * K)
V_M = 22.4     # m³ / kmol

M = {
    'C': 12.011, 'Si': 28.086, 'Mn': 54.938, 'Cr': 51.996, 'Fe': 55.845,
    'O': 16.000, 'O2': 31.998, 'CO': 28.010, 'CO2': 44.010,
    'FeO': 71.844, 'SiO2': 60.084, 'MnO': 70.937, 'Cr2O3': 151.990,
    'CaO': 56.077, 'MgO': 40.304, 'Al2O3': 101.960
}

interaction_params = {
    'C':  {'C': 0.14,   'Cr': -0.024,  'Mn': -0.012, 'Si': 0.08,   'O': -0.34,  'Ni': 0.012,  'Mo': -0.0083},
    'Cr': {'C': -0.12,  'Cr': -0.0003, 'Mn': 0.0,    'Si': -0.0043,'O': -0.14,  'Ni': 0.0002, 'Mo': 0.0018},
    'Mn': {'C': -0.0538,'Cr': 0.0039,  'Mn': 0.0,    'Si': 0.0,    'O': -0.083, 'Ni': 0.0,    'Mo': 0.0},
    'Si': {'C': 0.24,   'Cr': 0.015,   'Mn': 0.002,  'Si': 0.37,   'O': 0.0,    'Ni': 0.005,  'Mo': 0.0}
}

K_COEFFS = {
    'C':  (-3795.0,  4.164),
    'Si': (17770.0, -6.122),
    'Mn': (8695.0,  -3.930),
    'Cr': (24025.0, -10.566),
}

CP = {
    'metal': 815.9,   # J / (kg * K)
    'slag':  1196.6,
    'O2':    918.4,
    'CO':    1044.7,
    'CaO':   920.5
}

RHO = {'O2': 1.4277, 'metal': 7370.0, 'slag': 3100.0}
T_IN = 298.15

# Heats of reaction
def delta_H_C(T):
    return (11852.0 - (2.367 * T + 1.708e-4 * (T**2) + 3.835e3 / max(T, 100.0))) * 1000.0

def delta_H_Cr(T):
    return (11519.0 - (1.148 * T + 4.4e-5 * (T**2) + 1.5e4 / max(T, 100.0))) * 1000.0

def delta_H_Si(T):
    return (30658.0 - (2.15 * T + 1.45e-4 * (T**2))) * 1000.0

def delta_H_Mn(T):
    return (7581.0 - (0.845 * T + 7.38e-5 * (T**2) + 6.69e4 / max(T, 100.0))) * 1000.0

# Initial conditions
W_HM = 158000.0
T_bath_0 = 1500.0
metal_pct_0 = {'C': 2.0, 'Si': 0.5, 'Mn': 0.3, 'Cr': 18.0, 'Fe': 79.2}

N_slag_0 = 30.0
T_slag_0 = 1800.0
slag_mol_frac_0 = {
    'FeO': 0.20, 'SiO2': 0.2488, 'CaO': 0.5312,
    'MnO': 0.02, 'Cr2O3': 0.0, 'MgO': 0.0, 'Al2O3': 0.0
}

Q_O2 = 6.871
eta = 0.95
P_total = 1.05
lime_step_kg = 10.0
dt = 2.0
t_end = 600
num_steps = int(t_end / dt)

# Thermodynamic functions
def calculate_metal_activities(metal_comp):
    f, a = {}, {}
    for elem in ['C', 'Si', 'Mn', 'Cr']:
        log_f = sum(interaction_params[elem].get(j, 0.0) * metal_comp.get(j, 0.0) for j in ['C', 'Si', 'Mn', 'Cr'])
        f[elem] = 10.0 ** log_f
        a[elem] = max(f[elem] * metal_comp[elem], 1e-12)
    return f, a

def calculate_slag_activities(T_slag, slag_mol_frac, slag_pct_Cr2O3=0.0):
    N_CaO   = slag_mol_frac.get('CaO', 0.0)
    N_MgO   = slag_mol_frac.get('MgO', 0.0)
    N_SiO2  = slag_mol_frac.get('SiO2', 0.0)
    N_FeO   = slag_mol_frac.get('FeO', 0.0)
    N_MnO   = slag_mol_frac.get('MnO', 0.0)
    N_AlO15 = slag_mol_frac.get('Al2O3', 0.0)
    N_CrO15 = slag_mol_frac.get('Cr2O3', 0.0)
    N_basic = N_CaO + N_MgO

    lg_g_FeO = (
        (4130.0 / T_slag) * N_basic * (N_SiO2 + 0.25 * N_AlO15)
        + (1720.0 / T_slag) * N_MnO * (N_SiO2 + 0.45 * N_CrO15)
        + (1246.0 / T_slag) * N_AlO15 * N_SiO2
        + (42.0 / T_slag) * N_MnO * N_AlO15
        + (692.0 / T_slag) * N_CrO15 * N_SiO2
    )
    lg_g_SiO2 = (
        lg_g_FeO
        - (4130.0 / T_slag) * N_basic
        - (1720.0 / T_slag) * N_MnO
        - (1246.0 / T_slag) * N_AlO15
        - (692.0 / T_slag) * N_CrO15
    )
    lg_g_MnO = (
        lg_g_FeO
        - (1720.0 / T_slag) * (N_SiO2 + 0.45 * N_CrO15)
        - (42.0 / T_slag) * N_AlO15
    )
    lg_g_Cr2O3 = (
        lg_g_FeO
        - (1859.0 / T_slag) * N_basic
        - (774.0 / T_slag) * N_MnO
        - (692.0 / T_slag) * N_SiO2
    )

    gamma = {
        'FeO': 10.0 ** lg_g_FeO,
        'SiO2': 10.0 ** lg_g_SiO2,
        'MnO': 10.0 ** lg_g_MnO,
        'Cr2O3': 10.0 ** lg_g_Cr2O3,
    }

    eps = 1e-15
    a_FeO = max(gamma['FeO'] * N_FeO, eps)
    a_SiO2 = max(gamma['SiO2'] * N_SiO2, eps)
    a_MnO = max(gamma['MnO'] * N_MnO, eps)

    if slag_pct_Cr2O3 >= 5.0:
        a_Cr2O3 = 1.0
    else:
        a_Cr2O3 = max(gamma['Cr2O3'] * N_CrO15, eps)

    return gamma, {'FeO': a_FeO, 'SiO2': a_SiO2, 'MnO': a_MnO, 'Cr2O3': a_Cr2O3}

def calculate_equilibrium_constants(T_bath):
    delta_G_std = {}
    for elem, (A, B) in K_COEFFS.items():
        log_K = (A / T_bath) + B
        delta_G_std[elem] = -2.302585 * R * T_bath * log_K
    return delta_G_std

def calculate_reaction_gibbs(T_bath, a_metal, a_slag, delta_G_std, P_CO):
    eps = 1e-15
    RT = R * T_bath
    Q = {
        'C':  max(P_CO, eps) / (a_metal['C'] * a_slag['FeO']),
        'Si': a_slag['SiO2'] / (a_metal['Si'] * (a_slag['FeO'] ** 2)),
        'Mn': a_slag['MnO'] / (a_metal['Mn'] * a_slag['FeO']),
        'Cr': a_slag['Cr2O3'] / ((a_metal['Cr'] ** 2) * (a_slag['FeO'] ** 3)),
    }
    G = {elem: delta_G_std[elem] + RT * math.log(Q[elem]) for elem in Q}
    return Q, G

def calculate_oxygen_distribution(G, metal_comp):
    w = {}
    w['C']  = abs(G['C']) if (G['C'] < 0 and metal_comp['C'] > 0.001) else 0.0
    w['Si'] = (abs(G['Si']) / 2.0) if (G['Si'] < 0 and metal_comp['Si'] > 0.001) else 0.0
    w['Mn'] = abs(G['Mn']) if (G['Mn'] < 0 and metal_comp['Mn'] > 0.001) else 0.0
    w['Cr'] = (abs(G['Cr']) / 3.0) if (G['Cr'] < 0 and metal_comp['Cr'] > 0.001) else 0.0

    total_w = sum(w.values())
    if total_w <= 0.0:
        return {'C': 0.25, 'Si': 0.25, 'Mn': 0.25, 'Cr': 0.25}
    return {elem: w[elem] / total_w for elem in w}

def calculate_oxidation_rates(Q_O2, X, metal_kg, dt, eta=0.95):
    W_m = sum(metal_kg.values())
    n_dot_O2 = Q_O2 / V_M
    stoich = {'C': 2.0, 'Si': 1.0, 'Mn': 2.0, 'Cr': 4.0 / 3.0}
    rate_pct_s, delta_pct, delta_kg = {}, {}, {}

    for elem in ['C', 'Si', 'Mn', 'Cr']:
        rate = (
            stoich[elem] * eta * n_dot_O2 * X[elem] * (100.0 * M[elem] / W_m)
            if metal_kg[elem] > 0.0 else 0.0
        )
        d_pct = rate * dt
        d_kg = (d_pct / 100.0) * W_m
        d_kg = min(d_kg, max(metal_kg[elem], 0.0))
        d_pct = (d_kg / W_m) * 100.0

        rate_pct_s[elem] = rate
        delta_pct[elem] = d_pct
        delta_kg[elem] = d_kg

    return rate_pct_s, delta_pct, delta_kg

def calculate_dynamic_P_CO(delta_kg_C, Q_O2, dt, eta, P_total=1.05):
    n_CO = delta_kg_C / M['C']
    n_O_unreact = (Q_O2 * (1.0 - eta) * dt) / V_M
    total_gas = n_CO + n_O_unreact
    if total_gas <= 0.0:
        return 0.90
    P_CO = (n_CO / total_gas) * P_total
    return max(min(P_CO, P_total), 0.05)

def calculate_bath_temperature_rise(T_curr, metal_kg_curr, slag_kg_curr, delta_kg, dt, Q_O2=6.871, eta=0.95, lime_add_kg=0.0):
    W_m = sum(metal_kg_curr.values())
    W_s = sum(slag_kg_curr.values())

    r = {el: (delta_kg.get(el, 0.0) / W_m) * (100.0 / dt) for el in delta_kg}

    q_rxn = (
        delta_H_C(T_curr) * r.get('C', 0.0)
        + delta_H_Cr(T_curr) * r.get('Cr', 0.0)
        + delta_H_Si(T_curr) * r.get('Si', 0.0)
        + delta_H_Mn(T_curr) * r.get('Mn', 0.0)
    )

    stoich_slag = (
        (M['Cr2O3'] / (2.0 * M['Cr'])) * r.get('Cr', 0.0)
        + (M['MnO'] / M['Mn']) * r.get('Mn', 0.0)
        + (M['SiO2'] / M['Si']) * r.get('Si', 0.0)
    )
    term_slag_form = -CP['slag'] * T_curr * stoich_slag

    term_metal_loss = CP['metal'] * T_curr * (r.get('C', 0.0) + r.get('Cr', 0.0) + r.get('Mn', 0.0) + r.get('Si', 0.0))
    term_CO_gas = -CP['CO'] * T_curr * (M['CO'] / M['C']) * r.get('C', 0.0)

    m_dot_O2_unreact = Q_O2 * (1.0 - eta) * RHO['O2']
    q_gas = m_dot_O2_unreact * CP['O2'] * (T_curr - T_IN)

    lime_rate_kg_s = lime_add_kg / dt
    q_lime = lime_rate_kg_s * CP['CaO'] * (T_curr - T_IN)

    term_losses = (100.0 / W_m) * (q_gas + q_lime)

    numerator = q_rxn + term_slag_form + term_metal_loss + term_CO_gas - term_losses
    denominator = 100.0 * CP['metal'] + 100.0 * CP['slag'] * (W_s / W_m)

    dT_dt = numerator / denominator
    delta_T = dT_dt * dt
    T_new = T_curr + delta_T

    metal_kg_new = {el: max(metal_kg_curr[el] - delta_kg.get(el, 0.0), 0.0) for el in metal_kg_curr}

    slag_kg_new = dict(slag_kg_curr)
    slag_kg_new['SiO2']  = slag_kg_new.get('SiO2', 0.0) + delta_kg['Si'] * (M['SiO2'] / M['Si'])
    slag_kg_new['MnO']   = slag_kg_new.get('MnO', 0.0) + delta_kg['Mn'] * (M['MnO'] / M['Mn'])
    slag_kg_new['Cr2O3'] = slag_kg_new.get('Cr2O3', 0.0) + delta_kg['Cr'] * (M['Cr2O3'] / (2.0 * M['Cr']))
    slag_kg_new['CaO']   = slag_kg_new.get('CaO', 0.0) + lime_add_kg

    return T_new, delta_T, dT_dt, metal_kg_new, slag_kg_new


# =====================================================================
# SIMULATION EXECUTION
# =====================================================================
print("Running complete 10-minute simulation...")
T_current = T_bath_0
metal_kg = {el: (pct / 100.0) * W_HM for el, pct in metal_pct_0.items()}
slag_kg = {ox: (N_slag_0 * frac) * M[ox] for ox, frac in slag_mol_frac_0.items() if frac > 0}
P_CO_current = 0.90

history = []

# Step 0
W_m0 = sum(metal_kg.values())
W_s0 = sum(slag_kg.values())
f_0, a_m0 = calculate_metal_activities(metal_pct_0)
g_0, a_s0 = calculate_slag_activities(T_current, slag_mol_frac_0, 0.0)
dG_0 = calculate_equilibrium_constants(T_current)
Q_0, G_0 = calculate_reaction_gibbs(T_current, a_m0, a_s0, dG_0, P_CO_current)

history.append({
    'step': 0, 'time_sec': 0.0, 'time_min': 0.0,
    'T_bath_C': T_current - 273.15, 'T_bath_K': T_current,
    'W_metal_kg': W_m0, 'W_slag_kg': W_s0,
    'C_pct': metal_pct_0['C'], 'Si_pct': metal_pct_0['Si'],
    'Mn_pct': metal_pct_0['Mn'], 'Cr_pct': metal_pct_0['Cr'], 'Fe_pct': metal_pct_0['Fe'],
    'm_C_kg': metal_kg['C'], 'm_Si_kg': metal_kg['Si'], 'm_Mn_kg': metal_kg['Mn'], 'm_Cr_kg': metal_kg['Cr'], 'm_Fe_kg': metal_kg['Fe'],
    'CaO_pct': (slag_kg.get('CaO', 0.0) / W_s0) * 100.0,
    'SiO2_pct': (slag_kg.get('SiO2', 0.0) / W_s0) * 100.0,
    'FeO_pct': (slag_kg.get('FeO', 0.0) / W_s0) * 100.0,
    'MnO_pct': (slag_kg.get('MnO', 0.0) / W_s0) * 100.0,
    'Cr2O3_pct': (slag_kg.get('Cr2O3', 0.0) / W_s0) * 100.0,
    'm_CaO_kg': slag_kg.get('CaO', 0.0), 'm_SiO2_kg': slag_kg.get('SiO2', 0.0),
    'm_FeO_kg': slag_kg.get('FeO', 0.0), 'm_MnO_kg': slag_kg.get('MnO', 0.0), 'm_Cr2O3_kg': slag_kg.get('Cr2O3', 0.0),
    'basicity': slag_kg.get('CaO', 0.0) / max(slag_kg.get('SiO2', 1e-4), 1e-4),
    'f_C': f_0['C'], 'f_Si': f_0['Si'], 'f_Mn': f_0['Mn'], 'f_Cr': f_0['Cr'],
    'a_C': a_m0['C'], 'a_Si': a_m0['Si'], 'a_Mn': a_m0['Mn'], 'a_Cr': a_m0['Cr'],
    'a_FeO': a_s0['FeO'], 'a_SiO2': a_s0['SiO2'], 'a_MnO': a_s0['MnO'], 'a_Cr2O3': a_s0['Cr2O3'],
    'G_C': G_0['C'], 'G_Si': G_0['Si'], 'G_Mn': G_0['Mn'], 'G_Cr': G_0['Cr'],
    'X_C': 0.0, 'X_Si': 0.0, 'X_Mn': 0.0, 'X_Cr': 0.0,
    'rate_C_pct_s': 0.0, 'rate_Si_pct_s': 0.0, 'rate_Mn_pct_s': 0.0, 'rate_Cr_pct_s': 0.0,
    'delta_C_kg': 0.0, 'delta_Si_kg': 0.0, 'delta_Mn_kg': 0.0, 'delta_Cr_kg': 0.0,
    'dH_C_J_kg': delta_H_C(T_current), 'dH_Si_J_kg': delta_H_Si(T_current),
    'dH_Mn_J_kg': delta_H_Mn(T_current), 'dH_Cr_J_kg': delta_H_Cr(T_current),
    'P_CO_atm': P_CO_current, 'dT_dt_K_s': 0.0, 'dT_dt_K_min': 0.0
})

step = 0
for time_sec in range(2, t_end + int(dt), int(dt)):
    step += 1
    time_min = time_sec / 60.0

    W_m = sum(metal_kg.values())
    W_s = sum(slag_kg.values())
    curr_metal_pct = {el: (kg / W_m) * 100.0 for el, kg in metal_kg.items()}
    curr_slag_pct  = {ox: (kg / W_s) * 100.0 for ox, kg in slag_kg.items()}

    slag_moles = {ox: kg / M[ox] for ox, kg in slag_kg.items() if ox in M}
    tot_n = sum(slag_moles.values())
    curr_slag_mol_frac = {ox: n / tot_n for ox, n in slag_moles.items()}

    f_curr, a_metal = calculate_metal_activities(curr_metal_pct)
    _, a_slag = calculate_slag_activities(T_current, curr_slag_mol_frac, curr_slag_pct.get('Cr2O3', 0.0))
    dG_std = calculate_equilibrium_constants(T_current)
    Q, G = calculate_reaction_gibbs(T_current, a_metal, a_slag, dG_std, P_CO_current)
    X = calculate_oxygen_distribution(G, curr_metal_pct)

    rate_s, d_pct, d_kg = calculate_oxidation_rates(Q_O2, X, metal_kg, dt, eta)
    P_CO_next = calculate_dynamic_P_CO(d_kg['C'], Q_O2, dt, eta, P_total)

    T_new, delta_T, dT_dt, metal_kg_new, slag_kg_new = calculate_bath_temperature_rise(
        T_curr=T_current, metal_kg_curr=metal_kg, slag_kg_curr=slag_kg,
        delta_kg=d_kg, dt=dt, Q_O2=Q_O2, eta=eta, lime_add_kg=lime_step_kg
    )

    metal_kg = metal_kg_new
    slag_kg = slag_kg_new
    T_current = T_new
    P_CO_current = P_CO_next

    W_m_up = sum(metal_kg.values())
    W_s_up = sum(slag_kg.values())
    new_metal_pct = {el: (kg / W_m_up) * 100.0 for el, kg in metal_kg.items()}
    new_slag_pct  = {ox: (kg / W_s_up) * 100.0 for ox, kg in slag_kg.items()}

    history.append({
        'step': step, 'time_sec': time_sec, 'time_min': time_min,
        'T_bath_C': T_current - 273.15, 'T_bath_K': T_current,
        'W_metal_kg': W_m_up, 'W_slag_kg': W_s_up,
        'C_pct': new_metal_pct['C'], 'Si_pct': new_metal_pct['Si'],
        'Mn_pct': new_metal_pct['Mn'], 'Cr_pct': new_metal_pct['Cr'], 'Fe_pct': new_metal_pct['Fe'],
        'm_C_kg': metal_kg['C'], 'm_Si_kg': metal_kg['Si'], 'm_Mn_kg': metal_kg['Mn'], 'm_Cr_kg': metal_kg['Cr'], 'm_Fe_kg': metal_kg['Fe'],
        'CaO_pct': new_slag_pct.get('CaO', 0.0), 'SiO2_pct': new_slag_pct.get('SiO2', 0.0),
        'FeO_pct': new_slag_pct.get('FeO', 0.0), 'MnO_pct': new_slag_pct.get('MnO', 0.0), 'Cr2O3_pct': new_slag_pct.get('Cr2O3', 0.0),
        'm_CaO_kg': slag_kg.get('CaO', 0.0), 'm_SiO2_kg': slag_kg.get('SiO2', 0.0),
        'm_FeO_kg': slag_kg.get('FeO', 0.0), 'm_MnO_kg': slag_kg.get('MnO', 0.0), 'm_Cr2O3_kg': slag_kg.get('Cr2O3', 0.0),
        'basicity': new_slag_pct.get('CaO', 0.0) / max(new_slag_pct.get('SiO2', 1e-4), 1e-4),
        'f_C': f_curr['C'], 'f_Si': f_curr['Si'], 'f_Mn': f_curr['Mn'], 'f_Cr': f_curr['Cr'],
        'a_C': a_metal['C'], 'a_Si': a_metal['Si'], 'a_Mn': a_metal['Mn'], 'a_Cr': a_metal['Cr'],
        'a_FeO': a_slag['FeO'], 'a_SiO2': a_slag['SiO2'], 'a_MnO': a_slag['MnO'], 'a_Cr2O3': a_slag['Cr2O3'],
        'G_C': G['C'], 'G_Si': G['Si'], 'G_Mn': G['Mn'], 'G_Cr': G['Cr'],
        'X_C': X['C'], 'X_Si': X['Si'], 'X_Mn': X['Mn'], 'X_Cr': X['Cr'],
        'rate_C_pct_s': rate_s['C'], 'rate_Si_pct_s': rate_s['Si'], 'rate_Mn_pct_s': rate_s['Mn'], 'rate_Cr_pct_s': rate_s['Cr'],
        'delta_C_kg': d_kg['C'], 'delta_Si_kg': d_kg['Si'], 'delta_Mn_kg': d_kg['Mn'], 'delta_Cr_kg': d_kg['Cr'],
        'dH_C_J_kg': delta_H_C(T_current), 'dH_Si_J_kg': delta_H_Si(T_current),
        'dH_Mn_J_kg': delta_H_Mn(T_current), 'dH_Cr_J_kg': delta_H_Cr(T_current),
        'P_CO_atm': P_CO_current, 'dT_dt_K_s': dT_dt, 'dT_dt_K_min': dT_dt * 60.0
    })

df_sim = pd.DataFrame(history)
print(f"Simulation completed! Rows generated: {len(df_sim)}")

# =====================================================================
# CREATE REPRODUCED EXCEL WORKBOOK
# =====================================================================
wb = openpyxl.Workbook()

# Sheet 1: Inputs_and_Constants
ws_inputs = wb.active
ws_inputs.title = 'Inputs_and_Constants'

# Sheet 2: Simulation_10min
ws_sim = wb.create_sheet('Simulation_10min')

# Sheet 3: Summary_Milestones
ws_mile = wb.create_sheet('Summary_Milestones')

# Styling palettes (Executive Engineering)
navy_fill = PatternFill(start_color="1F4E79", end_color="1F4E79", fill_type="solid")
accent_fill = PatternFill(start_color="2F5597", end_color="2F5597", fill_type="solid")
sub_header_fill = PatternFill(start_color="D9E1F2", end_color="D9E1F2", fill_type="solid")
gray_fill = PatternFill(start_color="F2F2F2", end_color="F2F2F2", fill_type="solid")
border_thin = Side(style='thin', color="D3D3D3")
cell_border = Border(left=border_thin, right=border_thin, top=border_thin, bottom=border_thin)

header_font = Font(name='Calibri', size=11, bold=True, color="FFFFFF")
sub_font = Font(name='Calibri', size=11, bold=True, color="1F4E79")
regular_font = Font(name='Calibri', size=10)
bold_font = Font(name='Calibri', size=10, bold=True)
title_font = Font(name='Calibri', size=16, bold=True, color="1F4E79")

# ---------------------------------------------------------------------
# BUILD SHEET 1: Inputs_and_Constants
# ---------------------------------------------------------------------
ws_inputs.views.sheetView[0].showGridLines = True
ws_inputs['A1'] = "AOD REFINING MODEL WITH CHROMIUM KINETICS & THERMODYNAMICS"
ws_inputs['A1'].font = title_font
ws_inputs['A2'] = "Reproduced & Expanded from Professor's Mathematical Model (Wei Zhu 2002 Framework)"
ws_inputs['A2'].font = Font(name='Calibri', size=11, italic=True, color="595959")

# Table 1: Operating Inputs
ws_inputs['A4'] = "1. Operating & Initial Metal/Slag Inputs"
ws_inputs['A4'].font = sub_font
inputs_data = [
    ["Parameter", "Symbol", "Value", "Unit", "Notes & References"],
    ["Hot Metal Mass", "W_HM", 158000.0, "kg", "158 tonnes initial molten steel"],
    ["Initial Bath Temperature", "T_bath_0", 1500.0, "K", "1226.85 °C molten hot metal"],
    ["Initial Carbon", "[%C]_0", 2.00, "%", "High carbon molten bath"],
    ["Initial Silicon", "[%Si]_0", 0.50, "%", "Dissolved silicon"],
    ["Initial Manganese", "[%Mn]_0", 0.30, "%", "Dissolved manganese"],
    ["Initial Chromium", "[%Cr]_0", 18.00, "%", "Chromium added (AOD Stainless Steel Grade)"],
    ["Initial Balance Iron", "[%Fe]_0", 79.20, "%", "Balance iron"],
    ["Initial Slag Quantity", "N_slag_0", 30.0, "kmol", "Initial retained slag mass (~1815 kg)"],
    ["Initial Slag Temperature", "T_slag_0", 1800.0, "K", "1526.85 °C"],
    ["Slag Initial FeO", "(%FeO)_0", 20.00, "%", "Slag mole percent"],
    ["Slag Initial SiO2", "(%SiO2)_0", 24.88, "%", "Slag mole percent"],
    ["Slag Initial CaO", "(%CaO)_0", 53.12, "%", "Slag mole percent"],
    ["Slag Initial MnO", "(%MnO)_0", 2.00, "%", "Slag mole percent"],
    ["Slag Initial Cr2O3", "(%Cr2O3)_0", 0.00, "%", "Baseline fresh slag"],
    ["Oxygen Blowing Flow Rate", "Q_O2", 6.871, "m³/s", "Unified oxygen flow (Wei Zhu 2002)"],
    ["Argon Flow Rate", "Q_Ar", 0.0, "m³/s", "Pure O2 blowing stage"],
    ["Oxygen Utilization Efficiency", "eta", 0.95, "-", "Unified oxygen efficiency"],
    ["Vessel Total Pressure", "P_total", 1.05, "atm", "Total converter operating pressure"],
    ["Time Step", "dt", 2.0, "s", "Differential time increment"],
    ["Simulation Duration", "t_end", 600.0, "s", "10.0 minutes total (300 steps)"],
    ["Lime Addition Rate", "W_lime", 10.0, "kg/step", "5 kg/s continuous lime feed for basicity control"]
]

for r_idx, row in enumerate(inputs_data, 5):
    for c_idx, val in enumerate(row, 1):
        cell = ws_inputs.cell(r_idx, c_idx, val)
        cell.border = cell_border
        if r_idx == 5:
            cell.fill = navy_fill
            cell.font = header_font
            cell.alignment = Alignment(horizontal='center', vertical='center')
        else:
            cell.font = regular_font
            if c_idx == 3:
                cell.alignment = Alignment(horizontal='right')
                cell.number_format = '#,##0.00' if isinstance(val, float) else '#,##0'
            elif c_idx in (2, 4):
                cell.alignment = Alignment(horizontal='center')

# Table 2: Wagner Interaction Parameters
ws_inputs['G4'] = "2. Wagner 1st-Order Interaction Parameters (e_i^j)"
ws_inputs['G4'].font = sub_font
wagner_headers = ["Element i \\ j", "C", "Si", "Mn", "Cr"]
ws_inputs.cell(5, 7, wagner_headers[0]).fill = navy_fill
ws_inputs.cell(5, 7).font = header_font
for j_idx, elem_j in enumerate(wagner_headers[1:], 8):
    c = ws_inputs.cell(5, j_idx, elem_j)
    c.fill = navy_fill
    c.font = header_font
    c.alignment = Alignment(horizontal='center')

wagner_rows = [
    ["C",  0.14,    0.08,   -0.012, -0.024],
    ["Si", 0.24,    0.37,    0.002,  0.015],
    ["Mn", -0.0538, 0.0,     0.0,    0.0039],
    ["Cr", -0.12,  -0.0043,  0.0,   -0.0003]
]
for r_idx, row in enumerate(wagner_rows, 6):
    for c_idx, val in enumerate(row, 7):
        cell = ws_inputs.cell(r_idx, c_idx, val)
        cell.border = cell_border
        cell.font = bold_font if c_idx == 7 else regular_font
        if c_idx > 7:
            cell.alignment = Alignment(horizontal='right')
            cell.number_format = '0.0000'

# Table 3: Thermochemical Specific Heats & Heats of Reaction
ws_inputs['G12'] = "3. Specific Heat Capacities (Wei Zhu 2002)"
ws_inputs['G12'].font = sub_font
cp_data = [
    ["Species", "Symbol", "c_p (J/kg.K)", "Notes"],
    ["Liquid Steel", "c_p,m", 815.9, "Molten bath thermal capacity"],
    ["Liquid Slag", "c_p,s", 1196.6, "Molten slag thermal capacity"],
    ["Oxygen Gas", "c_p,O2", 918.4, "Cold gas sensible heating"],
    ["Carbon Monoxide", "c_p,CO", 1044.7, "Off-gas heat removal"],
    ["Lime (CaO)", "c_p,CaO", 920.5, "Solid addition heat absorption"]
]
for r_idx, row in enumerate(cp_data, 13):
    for c_idx, val in enumerate(row, 7):
        cell = ws_inputs.cell(r_idx, c_idx, val)
        cell.border = cell_border
        if r_idx == 13:
            cell.fill = navy_fill
            cell.font = header_font
            cell.alignment = Alignment(horizontal='center')
        else:
            cell.font = regular_font
            if c_idx == 9:
                cell.alignment = Alignment(horizontal='right')
                cell.number_format = '#,##0.0'

# Auto-fit columns in Sheet 1
for col in ws_inputs.columns:
    max_len = max(len(str(cell.value or '')) for cell in col)
    col_letter = get_column_letter(col[0].column)
    ws_inputs.column_dimensions[col_letter].width = max(max_len + 3, 12)

# ---------------------------------------------------------------------
# BUILD SHEET 2: Simulation_10min (Full Telemetry Table)
# ---------------------------------------------------------------------
ws_sim.views.sheetView[0].showGridLines = True

headers_sim = [
    ("Step", "step", "#,##0"),
    ("Time (s)", "time_sec", "#,##0.0"),
    ("Time (min)", "time_min", "0.00"),
    ("Temp (°C)", "T_bath_C", "#,##0.1"),
    ("Temp (K)", "T_bath_K", "#,##0.1"),
    ("Metal Mass (kg)", "W_metal_kg", "#,##0"),
    ("Slag Mass (kg)", "W_slag_kg", "#,##0"),
    ("[%C] (%)", "C_pct", "0.0000"),
    ("[%Si] (%)", "Si_pct", "0.0000"),
    ("[%Mn] (%)", "Mn_pct", "0.0000"),
    ("[%Cr] (%)", "Cr_pct", "0.0000"),
    ("[%Fe] (%)", "Fe_pct", "0.0000"),
    ("Mass C (kg)", "m_C_kg", "#,##0.1"),
    ("Mass Si (kg)", "m_Si_kg", "#,##0.1"),
    ("Mass Mn (kg)", "m_Mn_kg", "#,##0.1"),
    ("Mass Cr (kg)", "m_Cr_kg", "#,##0.1"),
    ("Mass Fe (kg)", "m_Fe_kg", "#,##0.1"),
    ("(%CaO) (%)", "CaO_pct", "0.00"),
    ("(%SiO2) (%)", "SiO2_pct", "0.00"),
    ("(%FeO) (%)", "FeO_pct", "0.00"),
    ("(%MnO) (%)", "MnO_pct", "0.00"),
    ("(%Cr2O3) (%)", "Cr2O3_pct", "0.00"),
    ("Basicity (CaO/SiO2)", "basicity", "0.000"),
    ("O2 Fract X_C", "X_C", "0.0000"),
    ("O2 Fract X_Si", "X_Si", "0.0000"),
    ("O2 Fract X_Mn", "X_Mn", "0.0000"),
    ("O2 Fract X_Cr", "X_Cr", "0.0000"),
    ("Rate C (%/s)", "rate_C_pct_s", "0.00000"),
    ("Rate Si (%/s)", "rate_Si_pct_s", "0.00000"),
    ("Rate Mn (%/s)", "rate_Mn_pct_s", "0.00000"),
    ("Rate Cr (%/s)", "rate_Cr_pct_s", "0.00000"),
    ("Delta C (kg)", "delta_C_kg", "0.00"),
    ("Delta Si (kg)", "delta_Si_kg", "0.00"),
    ("Delta Mn (kg)", "delta_Mn_kg", "0.00"),
    ("Delta Cr (kg)", "delta_Cr_kg", "0.00"),
    ("Delta G C (J/mol)", "G_C", "#,##0"),
    ("Delta G Si (J/mol)", "G_Si", "#,##0"),
    ("Delta G Mn (J/mol)", "G_Mn", "#,##0"),
    ("Delta G Cr (J/mol)", "G_Cr", "#,##0"),
    ("a_C", "a_C", "0.00000"),
    ("a_Si", "a_Si", "0.00000"),
    ("a_Mn", "a_Mn", "0.00000"),
    ("a_Cr", "a_Cr", "0.00000"),
    ("P_CO (atm)", "P_CO_atm", "0.000"),
    ("dT/dt (K/min)", "dT_dt_K_min", "0.00")
]

# Write header row
for c_idx, (col_name, _, _) in enumerate(headers_sim, 1):
    c = ws_sim.cell(1, c_idx, col_name)
    c.fill = navy_fill
    c.font = header_font
    c.alignment = Alignment(horizontal='center', vertical='center', wrap_text=True)
    c.border = cell_border
ws_sim.row_dimensions[1].height = 28

# Write data rows
for r_idx, row_dict in enumerate(history, 2):
    bg_color = gray_fill if r_idx % 2 == 1 else None
    for c_idx, (_, key, num_fmt) in enumerate(headers_sim, 1):
        val = row_dict.get(key, 0.0)
        c = ws_sim.cell(r_idx, c_idx, val)
        c.font = regular_font
        c.border = cell_border
        c.number_format = num_fmt
        c.alignment = Alignment(horizontal='right')
        if bg_color:
            c.fill = bg_color

# Set column widths
for c_idx, (col_name, _, _) in enumerate(headers_sim, 1):
    col_letter = get_column_letter(c_idx)
    ws_sim.column_dimensions[col_letter].width = max(len(col_name) + 3, 11)

# ---------------------------------------------------------------------
# BUILD SHEET 3: Summary_Milestones
# ---------------------------------------------------------------------
ws_mile.views.sheetView[0].showGridLines = True
ws_mile['A1'] = "AOD REFINING SIMULATION - SUMMARY MILESTONES"
ws_mile['A1'].font = title_font
ws_mile['A2'] = "Key process milestones extracted at 0, 1, 2, 5, and 10 minutes of refining"
ws_mile['A2'].font = Font(name='Calibri', size=11, italic=True, color="595959")

milestone_indices = [0, 30, 60, 150, 300]
mile_cols = [
    ("Milestone", lambda r: f"Step {r['step']} ({r['time_min']:.0f} min)", "@"),
    ("Time (min)", lambda r: r['time_min'], "0.0"),
    ("Bath Temp (°C)", lambda r: r['T_bath_C'], "#,##0.1"),
    ("[%C] (wt.%)", lambda r: r['C_pct'], "0.0000"),
    ("[%Si] (wt.%)", lambda r: r['Si_pct'], "0.0000"),
    ("[%Mn] (wt.%)", lambda r: r['Mn_pct'], "0.0000"),
    ("[%Cr] (wt.%)", lambda r: r['Cr_pct'], "0.0000"),
    ("Cr Retention (%)", lambda r: (r['Cr_pct'] / metal_pct_0['Cr']) * 100.0, "0.00%"),
    ("(%CaO) (%)", lambda r: r['CaO_pct'], "0.00"),
    ("(%SiO2) (%)", lambda r: r['SiO2_pct'], "0.00"),
    ("(%FeO) (%)", lambda r: r['FeO_pct'], "0.00"),
    ("(%MnO) (%)", lambda r: r['MnO_pct'], "0.00"),
    ("(%Cr2O3) (%)", lambda r: r['Cr2O3_pct'], "0.00"),
    ("Slag Basicity", lambda r: r['basicity'], "0.000"),
    ("P_CO (atm)", lambda r: r['P_CO_atm'], "0.000"),
    ("Heating Rate (K/min)", lambda r: r['dT_dt_K_min'], "0.00")
]

for c_idx, (name, _, _) in enumerate(mile_cols, 1):
    c = ws_mile.cell(4, c_idx, name)
    c.fill = navy_fill
    c.font = header_font
    c.alignment = Alignment(horizontal='center', vertical='center')
    c.border = cell_border
ws_mile.row_dimensions[4].height = 25

for r_offset, m_idx in enumerate(milestone_indices, 5):
    m_row = history[m_idx]
    for c_idx, (_, fn, fmt) in enumerate(mile_cols, 1):
        val = fn(m_row)
        c = ws_mile.cell(r_offset, c_idx, val)
        c.font = bold_font if c_idx in (1, 3, 4, 7) else regular_font
        c.border = cell_border
        c.number_format = fmt
        c.alignment = Alignment(horizontal='center' if c_idx == 1 else 'right')
    ws_mile.row_dimensions[r_offset].height = 20

for col in ws_mile.columns:
    max_len = max(len(str(cell.value or '')) for cell in col)
    col_letter = get_column_letter(col[0].column)
    ws_mile.column_dimensions[col_letter].width = max(max_len + 3, 14)

# ---------------------------------------------------------------------
# EMBED NATIVE EXCEL CHARTS (Charts Sheet)
# ---------------------------------------------------------------------
ws_charts = wb.create_sheet('Visual_Charts')
ws_charts.views.sheetView[0].showGridLines = True
ws_charts['A1'] = "AOD PROCESS VISUALIZATION PLOTS"
ws_charts['A1'].font = title_font

# Chart 1: Metal Bath Solutes
c1 = LineChart()
c1.title = "1. Metal Bath Composition vs Time (Cr, C, Si, Mn)"
c1.style = 13
c1.y_axis.title = "Concentration (wt.%)"
c1.x_axis.title = "Blowing Time (min)"
c1.width = 18
c1.height = 11

data_c1 = Reference(ws_sim, min_col=8, min_row=1, max_col=11, max_row=len(history)+1) # [%C] to [%Cr]
cats = Reference(ws_sim, min_col=3, min_row=2, max_row=len(history)+1) # Time (min)
c1.add_data(data_c1, titles_from_data=True)
c1.set_categories(cats)
ws_charts.add_chart(c1, 'B3')

# Chart 2: Bath Temperature Rise
c2 = LineChart()
c2.title = "2. Bath Temperature Rise vs Time"
c2.style = 10
c2.y_axis.title = "Bath Temperature (°C)"
c2.x_axis.title = "Blowing Time (min)"
c2.width = 18
c2.height = 11

data_c2 = Reference(ws_sim, min_col=4, min_row=1, max_col=4, max_row=len(history)+1) # Temp (°C)
c2.add_data(data_c2, titles_from_data=True)
c2.set_categories(cats)
ws_charts.add_chart(c2, 'L3')

# Chart 3: Complete Slag Evolution
c3 = LineChart()
c3.title = "3. Complete Slag Composition (%CaO, %SiO2, %FeO, %MnO, %Cr2O3)"
c3.style = 12
c3.y_axis.title = "Slag Oxide (wt.%)"
c3.x_axis.title = "Blowing Time (min)"
c3.width = 18
c3.height = 11

data_c3 = Reference(ws_sim, min_col=18, min_row=1, max_col=22, max_row=len(history)+1) # Slag oxides
c3.add_data(data_c3, titles_from_data=True)
c3.set_categories(cats)
ws_charts.add_chart(c3, 'B25')

# Save workbook
output_path = r"c:\Users\rosha\OneDrive\Desktop\BTP\AOD_Model_Refined_With_Chromium.xlsx"
wb.save(output_path)
print(f"Saved complete reproduced Excel workbook to: {output_path}")

# Also update Wei_Zhu_AOD_Model_Refined.xlsx for seamless continuity
wb.save(r"c:\Users\rosha\OneDrive\Desktop\BTP\Wei_Zhu_AOD_Model_Refined.xlsx")
print("Also updated Wei_Zhu_AOD_Model_Refined.xlsx successfully!")
