# Hovorka + EGP Glucose-Insulin Simulator ("Proposed Model")

A Python/Streamlit implementation of the Hovorka glucose-insulin model extended with a glycogenolysis-driven endogenous glucose production (EGP) term. Built for an undergraduate research project (URF) at Thapar Institute of Engineering and Technology, targeting an IEEE publication.

This started as a port of an Arduino/C reference implementation, then was rewritten and validated directly against the source paper's equations and parameter tables — several real discrepancies between the original Arduino code and the paper were found and fixed along the way (see **Model Notes & Fixes** below).

## Overview

The model simulates blood glucose response to meals (carbohydrate intake) and exogenous insulin infusion, using a system of ODEs solved with 4th-order Runge-Kutta integration. On top of the classic Hovorka structure, it adds a glucagon-driven glycogenolysis mechanism affecting endogenous glucose production.

Run it interactively via the Streamlit app: set patient body weight and basal insulin, configure up to 8 meals and 8 insulin boluses with timing/duration/amount, and get a glucose-vs-time plot plus a downloadable CSV trace.

## Repository Structure

```
.
├── app.py            # Streamlit UI — patient/simulation inputs, plotting, CSV export
├── model.py           # Body class: RK4 ODE integration of the glucose-insulin-glucagon system
├── constants.py        # Model parameters, initial conditions, and default meal/insulin schedule
└── requirements.txt      # streamlit, matplotlib
```

## Getting Started

```bash
git clone https://github.com/helina379/Havorkanew.git
cd Havorkanew
pip install -r requirements.txt
streamlit run app.py
```

This opens a browser dashboard where you can adjust patient parameters, meal/insulin schedules, and simulation settings, then run the simulation and download the resulting glucose trace as CSV.

## Model Structure

`model.py`'s `Body` class integrates, at each timestep:

- **Meal absorption** — two-compartment carbohydrate absorption (`Dm1`, `Dm2`) feeding glucose appearance rate `Ug`
- **Glucose subsystem** — plasma (`G`) and interstitial (`Gt`) glucose, renal glucose loss above a renal threshold, and saturable non-insulin-mediated glucose uptake (`F01`)
- **Endogenous glucose production (EGP)** — driven by glycogen store `G6p`, itself driven by a glucagon-dependent glycogenolysis term (`Ggg`) using a tanh-based time-decay factor
- **Glucagon (`c`, i.e. H(t) in the paper)** — secretion/suppression dynamics feeding into EGP
- **Insulin absorption** — classic two-compartment subcutaneous delay model (`S1`, `S2` → `tau_s`)
- **Insulin action** — three action variables (`x1`, `x2`, `x3`) on glucose distribution, disposal, and EGP suppression
- **Plasma insulin (`I`)** — driven by absorbed insulin `Ui` and elimination rate `ke`

All subsystems use RK4 integration with a shared, configurable step size.

## Model Notes & Fixes

The original Arduino/C source this was ported from diverged from the reference paper in several places. These were identified and corrected in `constants.py` (see the module docstring there for full derivations):

1. **`ke` (insulin elimination rate)** — corrected to 0.138/min per the paper (was 0.02), which is also required for the paper's insulin steady-state equation to actually hold.
2. **Insulin absorption ODE** — reverted to the paper's classic two-compartment `tau_s`-based delay model; the Arduino code used a different, unreferenced rate-constant formulation that didn't match its own paper-derived initial conditions.
3. **Renal threshold** — split into two distinct paper values that the Arduino code had collapsed into one: `Gth = 162 mg/dl` (renal glucose loss) and `Gth1 = 60 mg/dl` (hypoglycemic threshold, glucagon secretion only).
4. **Interstitial glucose initial condition** — `Gt(0) = 70 mg/dl`, distinct from `G(0) = 90 mg/dl`, per the paper (the Arduino code incorrectly set them equal).
5. **`EGP(0)`** — set to the paper's stated `EGPb = 1.23 mg/dl/min` rather than the Arduino code's derived value.
6. **`G6P(0)`** — computed from the paper's own `EGPb`/`kG6P`/`g6po` values rather than a hardcoded constant.
7. **`Sit`/`Sid`/`Sie`** — confirmed correct at their classic, unscaled values by cross-checking derived `kb1`/`kb2`/`kb3` against the paper's directly stated values.

**Still open / unverified** (flagged in code, not silently changed):
- `sigma` — the paper's stated value uses a different unit system than the code's; left as-is pending clarification.
- `SRbH` (`Srhb`) — a possible transcription artifact in the paper's parameter table; using the analytically-derived value instead.

The default meal/insulin schedule in `constants.py` uses the reference C code's internal minute-offsets (not literal clock times converted to minutes) so that simulated output matches the shape and peak timing of the original reference plot.

## Research Context

Developed as part of undergraduate research supervised by Dr. Sangeeta Kamboj and Dr. Sahaj Saxena (TIET Patiala), aimed at an IEEE publication. A MATLAB port of this simulator was also produced for collaborators. This work also underpins related efforts such as CarbSnap, a food-image carbohydrate estimation tool.

## License

Add a license (e.g. MIT) if this is intended to be shared/reused publicly.
