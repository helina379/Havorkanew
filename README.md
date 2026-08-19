# Hovorka Glucose-Insulin Simulator

A Python/Streamlit implementation of the Hovorka glucose-insulin model, extended with an endogenous glucose production (EGP) term driven by glucagon-mediated glycogenolysis. Developed as part of an undergraduate research project (URF) at Thapar Institute of Engineering and Technology, under the supervision of Dr. Sangeeta Kamboj and Dr. Sahaj Saxena, with an IEEE publication as the end goal.

This began as a port of an existing Arduino/C reference implementation. While porting it, I found several places where that implementation didn't actually match the paper it was based on, so I reworked the affected parameters against the paper's own equations and tables. Details are below and in `constants.py`.

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

Each timestep integrates, via RK4:

- **Meal absorption** — two-compartment carbohydrate absorption feeding the glucose appearance rate
- **Glucose subsystem** — plasma and interstitial glucose, renal glucose loss above threshold, saturable non-insulin-mediated uptake
- **EGP** — driven by a glycogen store, which is in turn driven by glucagon
- **Glucagon dynamics**
- **Insulin absorption** — classic two-compartment subcutaneous delay model
- **Insulin action** — three action variables affecting glucose distribution, disposal, and EGP suppression
- **Plasma insulin**

## Corrections to the Original Reference Implementation

The Arduino/C source this was ported from diverged from the underlying paper in a few places. These were identified and corrected in `constants.py` (see the module docstring for full derivations):

1. **`ke` (insulin elimination rate)** — corrected to 0.138/min per the paper (was 0.02). This also happens to be required for the paper's own insulin steady-state equation to hold.
2. **Insulin absorption ODE** — reverted to the paper's classic two-compartment `tau_s`-based delay model. The original code used a different rate-constant formulation not referenced anywhere in the paper.
3. **Renal threshold** — the paper specifies two distinct thresholds (renal: 162 mg/dl, hypoglycemic: 60 mg/dl) that the original code had collapsed into one, causing continuous renal glucose loss even at normal glucose levels.
4. **Interstitial glucose initial condition** — set to 70 mg/dl per the paper, distinct from plasma glucose's initial 90 mg/dl (the original code set both equal).
5. **`EGP(0)`** — corrected to the paper's stated value rather than a derived approximation.
6. **`G6P(0)`** — now computed directly from the paper's formula rather than hardcoded.
7. **`Sit`/`Sid`/`Sie`** — verified correct at their original values by cross-checking derived `kb1`/`kb2`/`kb3` against the paper's table.

**Open questions, flagged but not silently resolved:**
- `sigma` — the paper states this in a different unit system than the code uses; the conversion needs confirming before changing it.
- `Srhb` — likely a transcription error in the paper's parameter table; using the analytically derived value in the meantime.

The default meal/insulin schedule in `constants.py` uses internal minute-offsets rather than literal clock times, since this is what reproduces the shape and timing of the reference plot the model was validated against.

## Context

This simulator supports a broader URF project alongside CarbSnap, a food-image carbohydrate estimation app. A MATLAB port of this model was also produced for collaborators working outside Python.
