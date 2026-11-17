"""
app.py
------
Interactive Web Application for AOD Refining Simulation
Based on Wei Zhu (2002) Model with Chromium Kinetics, Mid-Step Oxygen Redistribution & Heat Balance.

Run locally:
    streamlit run app.py
"""

import os
os.environ['PANDAS_USE_PYARROW'] = '0'

import math
import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# =============================================================================
# PAGE CONFIGURATION
# =============================================================================
st.set_page_config(
    page_title="AOD Refining Simulator | Interactive Model",
    page_icon="🔥",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom Styling
st.markdown("""
<style>
    .main-title {
        font-size: 2.2rem;
        font-weight: 800;
        color: #1E3A8A;
        margin-bottom: 0.2rem;
    }
    .sub-title {
        font-size: 1.05rem;
        color: #4B5563;
        margin-bottom: 1.2rem;
    }
    .stMetric {
        background-color: #F8FAFC;
        border: 1px solid #E2E8F0;
        padding: 10px;
        border-radius: 8px;
    }
</style>
""", unsafe_allow_html=True)

st.markdown('<div class="main-title">🔥 AOD Stainless Steel Refining Simulator</div>', unsafe_allow_html=True)
st.markdown('<div class="sub-title">Interactive Mathematical Model with Chromium Kinetics, Dynamic Slag & Thermochemical Heat Balance</div>', unsafe_allow_html=True)

# =============================================================================
# SIMULATION ENGINE FUNCTION (EXACT FINAL USER MODEL)
# ==========================================================
@st.cache_data(show_spinner=False)
def run_aod_simulation(initial_C, initial_Cr, initial_Si, initial_Mn, initial_Ni, initial_T, oxygen_rate, dt=2, duration=600):
    
    # 1. CONSTANTS & INITIAL INPUTS
    W_m = 158000000.0            # Mass of liquid steel (g) = 158,000 kg
    P_t = 1.05                   # Total pressure in the AOD vessel (atm)
    eta = 0.95                   # Oxygen utilization efficiency (95%)
    R = 8.314                    # Universal gas constant (J / mol*K)
    rho_O2 = 1.4277 * (10**-3)   # Density of oxygen gas (g/cm^3)
    c_p_O2 = 0.918               # Specific heat of oxygen gas (J/g*K)
    c_p_CaO = 0.9205             # Specific heat of calcium oxide (J/g*K)
    c_p_scrap = 0.7113           # Specific heat of scrap (J/g*K)
    c_p_s = 1.1966               # Specific heat of slag (J/g*K)
    c_p_m = 0.8159               # Specific heat of molten steel (J/g*K)
    c_p_CO = 1.0447              # Specific heat of carbon monoxide (J/g*K)
    T_g0 = 298.15                # Reference temperature for gas (K)

    # Atomic Weights (g/mol)
    M_C = 12.011
    M_Cr = 51.996
    M_Si = 28.085
    M_Mn = 54.938
    M_Fe = 55.845
    M_Ni = 58.69 
    
    # Molecular Weights for Slag Oxides (g/mol)
    M_Cr2O3 = 151.99
    M_SiO2 = 60.08
    M_MnO = 70.94
    M_FeO = 71.84
    M_CaO = 56.08
    M_MgO = 40.30
    M_CO = 28.01
    
    # Initial State of Metal Bath
    pct_C = initial_C    
    pct_Cr = initial_Cr  
    pct_Si = initial_Si  
    pct_Mn = initial_Mn  
    pct_Ni = initial_Ni  
    pct_Fe = max(100.0 - (pct_C + pct_Cr + pct_Si + pct_Mn + pct_Ni), 0.0)
    T = initial_T
    T_Slag = 1800.0
         
    # Initial Slag Masses in Grams
    W_CaO = 2000.0 * 1000.0
    W_MgO = 500.0 * 1000.0
    W_SiO2 = 100.0 * 1000.0
    W_MnO = 50.0 * 1000.0
    W_Cr2O3 = 50.0 * 1000.0
    W_FeO = 100.0 * 1000.0
    W_s = W_CaO + W_MgO + W_SiO2 + W_MnO + W_Cr2O3 + W_FeO
    
    Q_O2 = (oxygen_rate * 1e6) / 60.0 
    Q_Ar = 0.0 
    
    results = []
    
    # Record Initial State at t = 0
    results.append({
        'Time': 0,
        'Time (min)': 0.0,
        'Carbon': pct_C,
        'Chromium': pct_Cr,
        'Silicon': pct_Si,
        'Manganese': pct_Mn,
        'Iron': pct_Fe,
        'Temperature': T,
        'Temperature (°C)': T - 273.15,
        'Metal Mass (kg)': W_m / 1000.0,
        'Slag Mass (kg)': W_s / 1000.0,
        'Slag_CaO(%)': (W_CaO / W_s) * 100.0,
        'Slag_MgO(%)': (W_MgO / W_s) * 100.0,
        'Slag_SiO2(%)': (W_SiO2 / W_s) * 100.0,
        'Slag_Cr2O3(%)': (W_Cr2O3 / W_s) * 100.0,
        'Slag_MnO(%)': (W_MnO / W_s) * 100.0,
        'Slag_FeO(%)': (W_FeO / W_s) * 100.0,
        'Basicity': W_CaO / max(W_SiO2, 1e-3),
        'P_CO': 0.75,
        'a_Cr2O3': 1.0
    })
    
    dC_dt = 0.0

    # 2. SIMULATION LOOP
    for t in range(dt, int(duration) + dt, dt):
            
        # Activity Coefficients (f_i)
        log_f_C = (0.14 * pct_C) - (0.024 * pct_Cr) + (0.08 * pct_Si) + (0.012 * pct_Ni)
        f_C = 10 ** log_f_C
            
        log_f_Cr = -(0.12 * pct_C) - (0.0003 * pct_Cr) - (0.0043 * pct_Si) + (0.0002 * pct_Ni)
        f_Cr = 10 ** log_f_Cr
            
        log_f_Si = (0.24 * pct_C) + (0.015 * pct_Cr) + (0.002 * pct_Mn) + (0.37 * pct_Si) + (0.005 * pct_Ni)
        f_Si = 10 ** log_f_Si
            
        log_f_Mn = -(0.0538 * pct_C) + (0.0039 * pct_Cr)
        f_Mn = 10 ** log_f_Mn
            
        # Thermodynamic Activities (a_i)
        a_C  = max(f_C * pct_C, 1e-12)
        a_Cr = max(f_Cr * pct_Cr, 1e-12)
        a_Si = max(f_Si * pct_Si, 1e-12)
        a_Mn = max(f_Mn * pct_Mn, 1e-12)
            
        # Dynamic Slag Thermodynamics
        W_s = W_CaO + W_MgO + W_SiO2 + W_MnO + W_Cr2O3 + W_FeO
            
        pct_Cr2O3_slag = (W_Cr2O3 / W_s) * 100.0
        pct_CaO_slag   = (W_CaO / W_s) * 100.0
        pct_MgO_slag   = (W_MgO / W_s) * 100.0
        pct_SiO2_slag  = (W_SiO2 / W_s) * 100.0
        pct_MnO_slag   = (W_MnO / W_s) * 100.0
        pct_FeO_slag   = (W_FeO / W_s) * 100.0
            
        n_CaO   = W_CaO / M_CaO
        n_MgO   = W_MgO / M_MgO
        n_SiO2  = W_SiO2 / M_SiO2
        n_MnO   = W_MnO / M_MnO
        n_FeO   = W_FeO / M_FeO
        n_Cr2O3 = W_Cr2O3 / M_Cr2O3
        n_Al2O3 = 0.0 
            
        n_total_std = n_CaO + n_MgO + n_SiO2 + n_MnO + n_FeO + n_Cr2O3 + n_Al2O3
        N_SiO2_std  = n_SiO2 / n_total_std
        N_MnO_std   = n_MnO / n_total_std
        N_FeO_std   = n_FeO / n_total_std
        N_Cr2O3_std = n_Cr2O3 / n_total_std
        
        # Moles of Pseudo-Components
        n_CrO15 = 2.0 * n_Cr2O3
        n_AlO15 = 2.0 * n_Al2O3
        n_total_pseudo = n_CaO + n_MgO + n_SiO2 + n_MnO + n_FeO + n_CrO15 + n_AlO15
        
        N_CaO   = n_CaO / n_total_pseudo
        N_MgO   = n_MgO / n_total_pseudo
        N_SiO2  = n_SiO2 / n_total_pseudo
        N_MnO   = n_MnO / n_total_pseudo
        N_CrO15 = n_CrO15 / n_total_pseudo
        N_AlO15 = n_AlO15 / n_total_pseudo
        
        log_gamma_FeO = ( (4130.0/T)*(N_CaO + N_MgO)*(N_SiO2 + 0.25*N_AlO15) + 
                            (1720.0/T)*N_MnO*(N_SiO2 + 0.45*N_CrO15) + 
                            (1246.0/T)*N_AlO15*N_SiO2 + 
                            (42.0/T)*N_MnO*N_AlO15 + 
                            (692.0/T)*N_CrO15*N_SiO2 )
        gamma_FeO = 10 ** log_gamma_FeO
        
        log_gamma_Cr2O3 = log_gamma_FeO - (1859.0/T)*(N_CaO + N_MgO) - (774.0/T)*N_MnO - (692.0/T)*N_SiO2
        gamma_Cr2O3 = 10 ** log_gamma_Cr2O3
        
        log_gamma_SiO2 = log_gamma_FeO - (4130.0/T)*(N_CaO + N_MgO) - (1720.0/T)*N_MnO - (1246.0/T)*N_AlO15 - (692.0/T)*N_CrO15
        gamma_SiO2 = 10 ** log_gamma_SiO2
        
        log_gamma_MnO = log_gamma_FeO - (1720.0/T)*(N_SiO2 + 0.45*N_CrO15) - (42.0/T)*N_AlO15
        gamma_MnO = 10 ** log_gamma_MnO
        
        a_FeO  = max(gamma_FeO * N_FeO_std, 1e-15)
        a_SiO2 = max(gamma_SiO2 * N_SiO2_std, 1e-15)
        a_MnO  = max(gamma_MnO * N_MnO_std, 1e-15)
        
        if pct_Cr2O3_slag >= 5.0:
            a_Cr2O3 = 1.0
        else:
            a_Cr2O3 = max(gamma_Cr2O3 * N_Cr2O3_std, 1e-15)
        
        # Equilibrium Constants
        K_C  = 10 ** ((-3795.0 / T) + 4.164)
        K_Cr = 10 ** ((24025.0 / T) - 10.566)
        K_Si = 10 ** ((17770.0 / T) - 6.122)
        K_Mn = 10 ** ((8695.0 / T) - 3.93)
        
        if t == dt:
            P_CO = 0.75 
        else:
            n_CO_sec = (W_m / (100.0 * M_C)) * abs(dC_dt)
            n_O2_unreacted_sec = (Q_O2 * (1.0 - eta)) / 22400.0
            n_inert_sec = Q_Ar / 22400.0
            n_total_gas_sec = n_CO_sec + n_O2_unreacted_sec + n_inert_sec
            if n_total_gas_sec > 0:
                P_CO = (n_CO_sec / n_total_gas_sec) * P_t
            else:
                P_CO = 0.75
            P_CO = max(min(P_CO, P_t), 0.05)
                
        # Gibbs Free Energy
        RT = R * T
        min_pct = 1e-5
        
        delG0_C  = -RT * math.log(K_C)
        delG0_Cr = -RT * math.log(K_Cr)
        delG0_Si = -RT * math.log(K_Si)
        delG0_Mn = -RT * math.log(K_Mn)
        
        Q_C  = P_CO / (a_C * a_FeO)
        Q_Cr = a_Cr2O3 / ((a_Cr ** 2) * (a_FeO ** 3))
        Q_Si = a_SiO2 / (a_Si * (a_FeO ** 2))
        Q_Mn = a_MnO / (a_Mn * a_FeO)
        
        delG_C  = delG0_C  + (RT * math.log(Q_C))
        delG_Cr = delG0_Cr + (RT * math.log(Q_Cr))
        delG_Si = delG0_Si + (RT * math.log(Q_Si))
        delG_Mn = delG0_Mn + (RT * math.log(Q_Mn))
        
        g_C  = abs(delG_C)        if pct_C  > min_pct else 0.0
        g_Cr = abs(delG_Cr / 3.0) if pct_Cr > min_pct else 0.0
        g_Si = abs(delG_Si / 2.0) if pct_Si > min_pct else 0.0
        g_Mn = abs(delG_Mn)       if pct_Mn > min_pct else 0.0
        
        sum_delG = g_C + g_Cr + g_Si + g_Mn
        if sum_delG > 0:
            x_C  = g_C  / sum_delG
            x_Cr = g_Cr / sum_delG
            x_Si = g_Si / sum_delG
            x_Mn = g_Mn / sum_delG
        else:
            x_C = x_Cr = x_Si = x_Mn = 0.25
        
        # Elemental Oxidation Rates with Mid-Step Redistribution
        mol_O2_eff = (eta * Q_O2) / 22400.0
        
        x_C_max  = max(0.0, (pct_C  - min_pct) * W_m / (200.0 * M_C  * mol_O2_eff * dt))
        x_Cr_max = max(0.0, (pct_Cr - min_pct) * 3.0 * W_m / (400.0 * M_Cr * mol_O2_eff * dt))
        x_Si_max = max(0.0, (pct_Si - min_pct) * W_m / (100.0 * M_Si * mol_O2_eff * dt))
        x_Mn_max = max(0.0, (pct_Mn - min_pct) * W_m / (200.0 * M_Mn * mol_O2_eff * dt))
        
        x_vals = {'C': x_C, 'Cr': x_Cr, 'Si': x_Si, 'Mn': x_Mn}
        g_vals = {'C': g_C, 'Cr': g_Cr, 'Si': g_Si, 'Mn': g_Mn}
        x_max  = {'C': x_C_max, 'Cr': x_Cr_max, 'Si': x_Si_max, 'Mn': x_Mn_max}
        capped = set()
        
        for _ in range(4):
            leftover_x = 0.0
            for el in ['C', 'Cr', 'Si', 'Mn']:
                if el not in capped and x_vals[el] > x_max[el]:
                    leftover_x += (x_vals[el] - x_max[el])
                    x_vals[el] = x_max[el]
                    capped.add(el)
            if leftover_x <= 0.0:
                break
            rem_g = sum(g_vals[el] for el in ['C', 'Cr', 'Si', 'Mn'] if el not in capped)
            if rem_g > 0.0:
                for el in ['C', 'Cr', 'Si', 'Mn']:
                    if el not in capped:
                        x_vals[el] += leftover_x * (g_vals[el] / rem_g)
                        
        x_C, x_Cr, x_Si, x_Mn = x_vals['C'], x_vals['Cr'], x_vals['Si'], x_vals['Mn']
        
        dC_dt  = - ( (200.0 * M_C * mol_O2_eff) / W_m ) * x_C
        dCr_dt = - ( (400.0 * M_Cr * mol_O2_eff) / (3.0 * W_m) ) * x_Cr
        dSi_dt = - ( (100.0 * M_Si * mol_O2_eff) / W_m ) * x_Si
        dMn_dt = - ( (200.0 * M_Mn * mol_O2_eff) / W_m ) * x_Mn

        # Enthalpies & Heat Balance
        T_g = T
        delH_C  = 11852.0 - (2.367 * T_g + 1.708e-4 * (T_g ** 2) + 3.835e3 / T_g)
        delH_Cr = 11519.0 - (1.148 * T + 4.4e-5 * (T ** 2) + 1.5e4 / T)
        delH_Si = 30658.0 - (2.15 * T + 1.45e-4 * (T ** 2))
        delH_Mn = 7581.0 - (0.845 * T + 7.38e-5 * (T ** 2) + 6.69e4 / T)
        
        term1 = c_p_s * T * (
            (M_Cr2O3 / (2.0 * M_Cr)) * dCr_dt +
            (M_MnO / M_Mn) * dMn_dt +
            (M_SiO2 / M_Si) * dSi_dt
        )
        term2 = -c_p_m * T * (dC_dt + dCr_dt + dMn_dt + dSi_dt)
        
        q_loss = 0.0
        q_5 = 0.0
        term3 = -(100.0 / W_m) * (
            Q_O2 * rho_O2 * c_p_O2 * ((1.0 - eta) * T_g - T_g0) + q_loss + q_5
        )
        term4 = c_p_CO * T_g * (M_CO / M_C) * dC_dt
        term5 = -(
            delH_C * dC_dt +
            delH_Cr * dCr_dt +
            delH_Mn * dMn_dt +
            delH_Si * dSi_dt
        )
        
        denom = (100.0 * c_p_m) + (100.0 * c_p_s * (W_s / W_m))
        dT_dt = (term1 + term2 + term3 + term4 + term5) / denom
        
        # State Update
        dW_Cr2O3 = -(W_m * dt / 100.0) * dCr_dt * (M_Cr2O3 / (2.0 * M_Cr))
        dW_SiO2  = -(W_m * dt / 100.0) * dSi_dt * (M_SiO2 / M_Si)
        dW_MnO   = -(W_m * dt / 100.0) * dMn_dt * (M_MnO / M_Mn)
        
        W_Cr2O3 += dW_Cr2O3
        W_SiO2  += dW_SiO2
        W_MnO   += dW_MnO
        W_s = W_CaO + W_MgO + W_SiO2 + W_Cr2O3 + W_MnO + W_FeO
        
        W_m = W_m * (1.0 + (dC_dt + dCr_dt + dMn_dt + dSi_dt) * (dt / 100.0))
        
        pct_C  += dC_dt * dt
        pct_Cr += dCr_dt * dt
        pct_Si += dSi_dt * dt
        pct_Mn += dMn_dt * dt
        pct_Fe = max(100.0 - (pct_C + pct_Cr + pct_Si + pct_Mn + pct_Ni), 0.0)
        
        T += dT_dt * dt
        
        results.append({
            'Time': t,
            'Time (min)': t / 60.0,
            'Carbon': pct_C,
            'Chromium': pct_Cr,
            'Silicon': pct_Si,
            'Manganese': pct_Mn,
            'Iron': pct_Fe,
            'Temperature': T,
            'Temperature (°C)': T - 273.15,
            'Metal Mass (kg)': W_m / 1000.0,
            'Slag Mass (kg)': W_s / 1000.0,
            'Slag_CaO(%)': pct_CaO_slag,
            'Slag_MgO(%)': pct_MgO_slag,
            'Slag_SiO2(%)': pct_SiO2_slag,
            'Slag_Cr2O3(%)': pct_Cr2O3_slag,
            'Slag_MnO(%)': pct_MnO_slag,
            'Slag_FeO(%)': pct_FeO_slag,
            'Basicity': W_CaO / max(W_SiO2, 1e-3),
            'P_CO': P_CO,
            'a_Cr2O3': a_Cr2O3
        })

    return pd.DataFrame(results)

# =============================================================================
# SIDEBAR CONTROLS (USER INPUTS)
# =============================================================================
with st.sidebar:
    st.header("⚙️ Simulation Inputs")
    
    with st.expander("1. Hot Metal Composition (wt.%)", expanded=True):
        col1, col2 = st.columns(2)
        with col1:
            inp_C = st.number_input("[%C]", min_value=0.05, max_value=5.0, value=2.0, step=0.1)
            inp_Si = st.number_input("[%Si]", min_value=0.01, max_value=2.0, value=0.5, step=0.05)
            inp_Ni = st.number_input("[%Ni]", min_value=0.0, max_value=25.0, value=8.0, step=0.5)
        with col2:
            inp_Cr = st.number_input("[%Cr]", min_value=0.0, max_value=35.0, value=18.0, step=0.5)
            inp_Mn = st.number_input("[%Mn]", min_value=0.01, max_value=3.0, value=0.8, step=0.05)
            
        balance_Fe = max(100.0 - (inp_C + inp_Cr + inp_Si + inp_Mn + inp_Ni), 0.0)
        st.info(f"Balance [%Fe]: **{balance_Fe:.2f}%**")

    with st.expander("2. Initial Bath Temperature", expanded=True):
        temp_unit = st.radio("Input Unit", ["Kelvin (K)", "Celsius (°C)"], horizontal=True)
        if temp_unit == "Kelvin (K)":
            inp_T = st.number_input("Bath Temperature (K)", min_value=1400.0, max_value=1900.0, value=1500.0, step=25.0)
        else:
            temp_c = st.number_input("Bath Temperature (°C)", min_value=1126.85, max_value=1626.85, value=1226.85, step=25.0)
            inp_T = temp_c + 273.15
        st.caption(f"Current initial: **{inp_T:.1f} K** ({inp_T-273.15:.1f} °C)")

    with st.expander("3. Gas Blowing & Timing", expanded=True):
        inp_o2_rate = st.number_input("Oxygen Rate (Nm³/min)", min_value=50.0, max_value=1000.0, value=450.0, step=25.0)
        col_t1, col_t2 = st.columns(2)
        with col_t1:
            inp_duration = st.number_input("Duration (s)", min_value=60, max_value=1800, value=600, step=60)
        with col_t2:
            inp_dt = st.number_input("Time Step dt (s)", min_value=1, max_value=5, value=2, step=1)
        st.caption(f"Total Steps: **{int(inp_duration / inp_dt)} discrete steps**")

    run_sim = st.button("🚀 Run AOD Simulation", type="primary", use_container_width=True)

# =============================================================================
# RUN SIMULATION & DISPLAY
# =============================================================================
with st.spinner("Computing AOD Thermochemical Kinetics & Heat Balance..."):
    df_sim = run_aod_simulation(
        initial_C=inp_C,
        initial_Cr=inp_Cr,
        initial_Si=inp_Si,
        initial_Mn=inp_Mn,
        initial_Ni=inp_Ni,
        initial_T=inp_T,
        oxygen_rate=inp_o2_rate,
        dt=int(inp_dt),
        duration=int(inp_duration)
    )

# Summary KPI Cards
row_0 = df_sim.iloc[0]
row_end = df_sim.iloc[-1]

c1, c2, c3, c4, c5 = st.columns(5)
with c1:
    delta_T = row_end['Temperature (°C)'] - row_0['Temperature (°C)']
    st.metric("Bath Temperature", f"{row_end['Temperature (°C)']:.1f} °C", f"+{delta_T:.1f} °C")
with c2:
    delta_C = row_end['Carbon'] - row_0['Carbon']
    st.metric("Final Carbon", f"{row_end['Carbon']:.3f} %", f"{delta_C:.3f} %")
with c3:
    cr_ret = (row_end['Chromium'] / row_0['Chromium']) * 100.0
    st.metric("Final Chromium", f"{row_end['Chromium']:.2f} %", f"{cr_ret:.1f}% retained")
with c4:
    st.metric("Final Slag Mass", f"{row_end['Slag Mass (kg)']:,.0f} kg", f"+{row_end['Slag Mass (kg)'] - row_0['Slag Mass (kg)']:,.0f} kg")
with c5:
    st.metric("Final Basicity", f"{row_end['Basicity']:.2f}", f"Initial {row_0['Basicity']:.1f}")

st.markdown("---")

# =============================================================================
# THE 4 USER PLOTS (INTERACTIVE PLOTLY & MATPLOTLIB)
# =============================================================================
st.subheader("📊 Simulation Trajectory (The 4 Key Engineering Plots)")

tab_plotly, tab_mpl = st.tabs(["Interactive Charts (Plotly)", "Publication 4-Panel Plot (Matplotlib)"])

with tab_plotly:
    col_p1, col_p2 = st.columns(2)
    
    # 1. Minor elements (C, Si, Mn)
    with col_p1:
        fig1 = go.Figure()
        fig1.add_trace(go.Scatter(x=df_sim['Time'], y=df_sim['Carbon'], mode='lines', name='Carbon (C %)', line=dict(color='#1f77b4', width=2.5)))
        fig1.add_trace(go.Scatter(x=df_sim['Time'], y=df_sim['Manganese'], mode='lines', name='Manganese (Mn %)', line=dict(color='#d62728', width=2.5)))
        fig1.add_trace(go.Scatter(x=df_sim['Time'], y=df_sim['Silicon'], mode='lines', name='Silicon (Si %)', line=dict(color='#2ca02c', width=2.5)))
        fig1.update_layout(title="<b>1. Metal Bath Composition: C, Si, Mn vs. Time</b>", xaxis_title="Time (s)", yaxis_title="Weight Percentage (%)", template="plotly_white", hovermode="x unified", height=400)
        st.plotly_chart(fig1, use_container_width=True)

    # 2. Major elements (Cr & Fe Dual Axis)
    with col_p2:
        fig2 = make_subplots(specs=[[{"secondary_y": True}]])
        fig2.add_trace(go.Scatter(x=df_sim['Time'], y=df_sim['Chromium'], mode='lines', name='Chromium (Cr %)', line=dict(color='#ff7f0e', width=2.5)), secondary_y=False)
        fig2.add_trace(go.Scatter(x=df_sim['Time'], y=df_sim['Iron'], mode='lines', name='Iron (Fe %)', line=dict(color='#7f7f7f', width=2.0, dash='dash')), secondary_y=True)
        fig2.update_layout(title="<b>2. Metal Bath Composition: Cr & Fe vs. Time</b>", xaxis_title="Time (s)", template="plotly_white", hovermode="x unified", height=400)
        fig2.update_yaxes(title_text="Chromium Weight %", secondary_y=False)
        fig2.update_yaxes(title_text="Iron Weight %", secondary_y=True)
        st.plotly_chart(fig2, use_container_width=True)

    col_p3, col_p4 = st.columns(2)
    
    # 3. Slag Composition
    with col_p3:
        fig3 = go.Figure()
        fig3.add_trace(go.Scatter(x=df_sim['Time'], y=df_sim['Slag_CaO(%)'], mode='lines', name='CaO (%)', line=dict(color='#1f77b4', width=2.0)))
        fig3.add_trace(go.Scatter(x=df_sim['Time'], y=df_sim['Slag_Cr2O3(%)'], mode='lines', name='Cr2O3 (%)', line=dict(color='#ff7f0e', width=2.5)))
        fig3.add_trace(go.Scatter(x=df_sim['Time'], y=df_sim['Slag_SiO2(%)'], mode='lines', name='SiO2 (%)', line=dict(color='#2ca02c', width=2.0)))
        fig3.add_trace(go.Scatter(x=df_sim['Time'], y=df_sim['Slag_MnO(%)'], mode='lines', name='MnO (%)', line=dict(color='#d62728', width=2.0)))
        fig3.add_trace(go.Scatter(x=df_sim['Time'], y=df_sim['Slag_MgO(%)'], mode='lines', name='MgO (%)', line=dict(color='#9467bd', width=1.8, dash='dashdot')))
        fig3.add_trace(go.Scatter(x=df_sim['Time'], y=df_sim['Slag_FeO(%)'], mode='lines', name='FeO (%)', line=dict(color='#8c564b', width=1.8, dash='dot')))
        fig3.update_layout(title="<b>3. Slag Composition Evolution vs. Time</b>", xaxis_title="Time (s)", yaxis_title="Slag Weight Percentage (%)", template="plotly_white", hovermode="x unified", height=400)
        st.plotly_chart(fig3, use_container_width=True)

    # 4. Bath Temperature
    with col_p4:
        fig4 = make_subplots(specs=[[{"secondary_y": True}]])
        fig4.add_trace(go.Scatter(x=df_sim['Time'], y=df_sim['Temperature'], mode='lines', name='Temperature (K)', line=dict(color='#e31a1c', width=2.5)), secondary_y=False)
        fig4.add_trace(go.Scatter(x=df_sim['Time'], y=df_sim['Temperature (°C)'], mode='lines', name='Temperature (°C)', line=dict(color='#b10026', width=0)), secondary_y=True)
        fig4.update_layout(title="<b>4. Bath Temperature vs. Time (K & °C)</b>", xaxis_title="Time (s)", template="plotly_white", hovermode="x unified", height=400)
        fig4.update_yaxes(title_text="Temperature (K)", secondary_y=False)
        fig4.update_yaxes(title_text="Temperature (°C)", secondary_y=True)
        st.plotly_chart(fig4, use_container_width=True)

with tab_mpl:
    # Exact Matplotlib 4-Panel Plot
    fig_mpl, axes = plt.subplots(2, 2, figsize=(15, 10))
    fig_mpl.suptitle('AOD Process Simulation: Metal Composition, Slag Composition & Bath Temperature vs. Time', fontsize=14, fontweight='bold')

    # 1. C, Si, Mn
    ax1 = axes[0, 0]
    ax1.plot(df_sim['Time'], df_sim['Carbon'], label='Carbon (C %)', color='#1f77b4', linewidth=2.2)
    ax1.plot(df_sim['Time'], df_sim['Manganese'], label='Manganese (Mn %)', color='#d62728', linewidth=2.2)
    ax1.plot(df_sim['Time'], df_sim['Silicon'], label='Silicon (Si %)', color='#2ca02c', linewidth=2.2)
    ax1.set_title('Metal Bath Composition: C, Si, Mn vs. Time', fontsize=11, fontweight='bold')
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('Weight Percentage (%)')
    ax1.legend(loc='upper right')
    ax1.grid(True, linestyle='--', alpha=0.6)

    # 2. Cr & Fe
    ax2 = axes[0, 1]
    ax2_fe = ax2.twinx()
    l_cr = ax2.plot(df_sim['Time'], df_sim['Chromium'], label='Chromium (Cr %)', color='#ff7f0e', linewidth=2.2)
    l_fe = ax2_fe.plot(df_sim['Time'], df_sim['Iron'], label='Iron (Fe % - Balance)', color='#7f7f7f', linewidth=2.0, linestyle='--')
    ax2.set_title('Metal Bath Composition: Cr & Fe vs. Time', fontsize=11, fontweight='bold')
    ax2.set_xlabel('Time (s)')
    ax2.set_ylabel('Chromium Weight %', color='#ff7f0e')
    ax2_fe.set_ylabel('Iron Weight %', color='#7f7f7f')
    lines = l_cr + l_fe
    labels = [l.get_label() for l in lines]
    ax2.legend(lines, labels, loc='center right')
    ax2.grid(True, linestyle='--', alpha=0.6)

    # 3. Slag
    ax3 = axes[1, 0]
    ax3.plot(df_sim['Time'], df_sim['Slag_CaO(%)'], label='CaO (%)', color='#1f77b4', linewidth=2.0)
    ax3.plot(df_sim['Time'], df_sim['Slag_Cr2O3(%)'], label='Cr2O3 (%)', color='#ff7f0e', linewidth=2.2)
    ax3.plot(df_sim['Time'], df_sim['Slag_SiO2(%)'], label='SiO2 (%)', color='#2ca02c', linewidth=2.0)
    ax3.plot(df_sim['Time'], df_sim['Slag_MnO(%)'], label='MnO (%)', color='#d62728', linewidth=2.0)
    ax3.plot(df_sim['Time'], df_sim['Slag_MgO(%)'], label='MgO (%)', color='#9467bd', linewidth=1.8, linestyle='-.')
    ax3.plot(df_sim['Time'], df_sim['Slag_FeO(%)'], label='FeO (%)', color='#8c564b', linewidth=1.8, linestyle=':')
    ax3.set_title('Slag Composition Evolution vs. Time', fontsize=11, fontweight='bold')
    ax3.set_xlabel('Time (s)')
    ax3.set_ylabel('Slag Weight Percentage (%)')
    ax3.legend(loc='center right', ncol=2, fontsize=9)
    ax3.grid(True, linestyle='--', alpha=0.6)

    # 4. Temperature
    ax4 = axes[1, 1]
    ax4.plot(df_sim['Time'], df_sim['Temperature'], label='Bath Temperature (K)', color='#e31a1c', linewidth=2.5)
    ax4_c = ax4.twinx()
    ax4_c.plot(df_sim['Time'], df_sim['Temperature'] - 273.15, color='#e31a1c', alpha=0.0)
    ax4.set_title('Bath Temperature vs. Time', fontsize=11, fontweight='bold')
    ax4.set_xlabel('Time (s)')
    ax4.set_ylabel('Temperature (K)')
    ax4_c.set_ylabel('Temperature (°C)')
    ax4.legend(loc='lower right')
    ax4.grid(True, linestyle='--', alpha=0.6)

    plt.tight_layout(rect=[0, 0.02, 1, 0.95])
    st.pyplot(fig_mpl)

# =============================================================================
# DATA TABLE & EXPORT
# =============================================================================
st.markdown("---")
st.subheader("📋 Step-by-Step Simulation Telemetry Table")

col_d1, col_d2 = st.columns([4, 1])
with col_d1:
    st.caption(f"Showing full step telemetry ({len(df_sim)} rows recorded).")
with col_d2:
    csv_bytes = df_sim.to_csv(index=False).encode('utf-8')
    st.download_button(
        label="📥 Download CSV Table",
        data=csv_bytes,
        file_name=f"AOD_Simulation_C{inp_C:.1f}_Cr{inp_Cr:.1f}.csv",
        mime="text/csv",
        use_container_width=True
    )

st.dataframe(
    df_sim.style.format({
        'Time': '{:.0f}',
        'Time (min)': '{:.2f}',
        'Carbon': '{:.4f}',
        'Chromium': '{:.4f}',
        'Silicon': '{:.4f}',
        'Manganese': '{:.4f}',
        'Iron': '{:.4f}',
        'Temperature': '{:.1f}',
        'Temperature (°C)': '{:.1f}',
        'Metal Mass (kg)': '{:,.1f}',
        'Slag Mass (kg)': '{:,.1f}',
        'Slag_CaO(%)': '{:.2f}',
        'Slag_MgO(%)': '{:.2f}',
        'Slag_SiO2(%)': '{:.2f}',
        'Slag_Cr2O3(%)': '{:.2f}',
        'Slag_MnO(%)': '{:.2f}',
        'Slag_FeO(%)': '{:.2f}',
        'Basicity': '{:.3f}',
        'P_CO': '{:.3f}',
        'a_Cr2O3': '{:.3f}'
    }),
    height=400,
    use_container_width=True
)
