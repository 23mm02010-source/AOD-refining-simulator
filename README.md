# AOD Stainless Steel Refining Simulator

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://share.streamlit.io)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

An interactive, physics-based mathematical simulation of the **Argon Oxygen Decarburization (AOD)** refining process for stainless steelmaking. 

This model reproduces and expands upon the seminal work of **Wei Zhu (2002)** (*"Mathematical Model for AOD Refining Process"*), incorporating **multicomponent solute thermodynamics**, **Turkdogan regular solution slag activities**, **mid-step oxygen redistribution**, and a **closed-form thermochemical heat balance** to track bath temperature and elemental trajectories in real time.

---

## Project Overview

During AOD stainless steelmaking, hot metal containing dissolved Carbon ($\text{C}$), Silicon ($\text{Si}$), Manganese ($\text{Mn}$), and valuable Chromium ($\text{Cr}$) is refined using high-velocity oxygen gas. Because the process is autogenous (it heats itself without external heating arcs), predicting the competition between **decarburization** and **chromium oxidation** requires simultaneously solving:

1. **Multicomponent Solute Thermodynamics:** Henrian activities via Wagner interaction parameters ($e_i^j$).
2. **Slag Thermodynamics:** Raoultian activities of $\text{FeO}, \text{SiO}_2, \text{MnO}, \text{Cr}_2\text{O}_3$ via Turkdogan's pseudo-component model with a 5% chromite saturation cap.
3. **Chemical Driving Forces:** Instantaneous Gibbs free energies ($\Delta G_i = \Delta G_i^\circ + RT \ln Q_i$).
4. **Kinetics & Oxygen Redistribution:** Oxygen allocation based on driving force modulus ($X_i$) with adaptive mid-step reallocation when an element finishes oxidation.
5. **Coupled Heat Balance:** Real-time bath temperature rise rate ($\frac{dT}{dt}$) from dynamic heats of reaction ($\Delta H_i(T)$) and slag/metal sensible heat capacities.

---

## Mathematical Formulation

### 1. Metal Bath Activity Coefficients (Wagner 1st-Order)
$$\log_{10} f_i = \sum_{j} e_i^j \cdot [\% j]$$
$$a_i = f_i \cdot [\% i]$$

### 2. Instantaneous Gibbs Free Energy & Oxygen Distribution
$$\Delta G_i = \Delta G_i^\circ + R \cdot T \cdot \ln(Q_i)$$
$$X_i = \frac{w_i}{\sum w_k} \quad \text{where} \quad w_i = \frac{|\Delta G_i|}{n_{\text{valence}}} \quad (\text{for } \Delta G_i < 0 \text{ and } [\%i] > 0)$$

### 3. Elemental Oxidation Rates (% / s)
$$\frac{d[\% C]}{dt} = - \left( \frac{200 \cdot M_{\text{C}} \cdot \dot{n}_{\text{O}_2, \text{eff}}}{W_m} \right) \cdot X_{\text{C}}$$
$$\frac{d[\% Cr]}{dt} = - \left( \frac{400 \cdot M_{\text{Cr}} \cdot \dot{n}_{\text{O}_2, \text{eff}}}{3 \cdot W_m} \right) \cdot X_{\text{Cr}}$$
$$\frac{d[\% Si]}{dt} = - \left( \frac{100 \cdot M_{\text{Si}} \cdot \dot{n}_{\text{O}_2, \text{eff}}}{W_m} \right) \cdot X_{\text{Si}}$$
$$\frac{d[\% Mn]}{dt} = - \left( \frac{200 \cdot M_{\text{Mn}} \cdot \dot{n}_{\text{O}_2, \text{eff}}}{W_m} \right) \cdot X_{\text{Mn}}$$

### 4. Bath Temperature Rising Rate (Wei Zhu Eq. [22])
$$\frac{dT}{dt} = \frac{\sum r_i \Delta H_i(T) - c_{p,s} T \sum \left(\frac{M_{\text{oxide}}}{n M_i}\right) r_i + c_{p,m} T \sum r_i - c_{p,\text{CO}} T \left(\frac{M_{\text{CO}}}{M_{\text{C}}}\right) r_{\text{C}} - \frac{100}{W_m} q_{\text{losses}}}{100 \cdot c_{p,m} + 100 \cdot c_{p,s} \left(\frac{W_s}{W_m}\right)}$$

---

## Key Features

* **Real-Time Interactive Web Dashboard:** Built with Streamlit for clean, responsive browser controls.
* **4-Panel Publication Plots:**
  * **Panel 1:** Metal Bath Minor Elements ($[\%C], [\%Si], [\%Mn]$ vs. time).
  * **Panel 2:** Metal Bath Major Elements ($[\%Cr]$ and Balance $[\%Fe]$ with dual Y-axes).
  * **Panel 3:** Complete Slag Composition Evolution ($\% \text{CaO}, \% \text{Cr}_2\text{O}_3, \% \text{SiO}_2, \% \text{MnO}, \% \text{MgO}, \% \text{FeO}$).
  * **Panel 4:** Molten Bath Temperature Rise ($T$ in $\text{K}$ and $^\circ\text{C}$).
* **Dynamic $P_{\text{CO}}$ Calculation:** Accounts for carbon monoxide evolution and dilution.
* **Full Data Telemetry & CSV Export:** Download complete step-by-step telemetry for reporting.

---

## Local Installation & Setup

### Prerequisites
* Python 3.9, 3.10, 3.11, or 3.12

### 1. Clone the Repository
```bash
git clone https://github.com/your-username/aod-refining-simulator.git
cd aod-refining-simulator
```

### 2. Install Dependencies
```bash
pip install -r requirements.txt
```

### 3. Launch the Web Application
```bash
streamlit run app.py
```
Your default web browser will automatically open:
```
Local URL: http://localhost:8501
```

---

## Deploy to Streamlit Community Cloud (Free Public Access)

To share this simulator with your professor or evaluators with a 24/7 public link:

1. Push this folder to your GitHub account (containing `app.py`, `requirements.txt`, and `README.md`).
2. Visit **[share.streamlit.io](https://share.streamlit.io)** and sign in with GitHub.
3. Click **"New app"** and select:
   * **Repository:** `your-username/aod-refining-simulator`
   * **Branch:** `main`
   * **Main file path:** `app.py`
4. Click **"Deploy!"**
5. Streamlit will generate a public URL (e.g. `https://aod-refining-simulator.streamlit.app`) accessible on any computer or mobile device.

---

## Repository File Structure

```
├── app.py                  # Main interactive Streamlit web application
├── requirements.txt        # Python package dependencies
├── README.md               # Project documentation & mathematical manual
└── AOD_Model_postmidsem.ipynb # Jupyter notebook for mathematical derivation
```

---

## References

1. **Wei, J., & Zhu, Y. (2002).** *Mathematical model for AOD refining process.* ISIJ International, 42(1), 14–22.
2. **Turkdogan, E. T. (1996).** *Fundamentals of Steelmaking.* The Institute of Materials, London.
3. **Sigworth, G. K., & Elliott, J. F. (1974).** *The Thermodynamics of Liquid Dilute Iron Alloys.* Metal Science, 8(1), 298–310.

