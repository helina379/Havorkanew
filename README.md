## Contents

- `app.py` — Streamlit interface: patient/simulation inputs, plotting, CSV export
- `model.py` — core ODE model, integrated with RK4
- `constants.py` — model parameters, initial conditions, and the default meal/insulin schedule
- `requirements.txt`

## Setup

```bash
git clone https://github.com/helina379/Havorkanew.git
cd Havorkanew
pip install -r requirements.txt
streamlit run app.py
```

This launches a browser interface where you can set patient body weight and basal insulin rate, configure meals and insulin boluses (time, duration, amount), and run the simulation. Output is a glucose-vs-time plot, with the trace also downloadable as CSV.

## Model Overview

- **Meal absorption** — two-compartment carbohydrate absorption feeding the glucose appearance rate
- **Glucose subsystem** — plasma and interstitial glucose, renal glucose loss above threshold, saturable non-insulin-mediated uptake
- **EGP** — driven by a glycogen store, which is in turn driven by glucagon
- **Glucagon dynamics**
- **Insulin absorption** — classic two-compartment subcutaneous delay model
- **Insulin action** — three action variables affecting glucose distribution, disposal, and EGP suppression
- **Plasma insulin**



