"""
Constants for the Hovorka + glycogenolysis-driven EGP ("Proposed") model.

Rewritten against the actual paper text/tables (New_Endogenous_Glucose_
Production_Model PDFs), not just the Arduino/Egp6_model.txt code, after
cross-checking parameter values. Compared to the earlier Arduino-literal
port, this fixes several real discrepancies between the code and the paper:

  1. ke: paper Table 4 states ke = 0.138 /min ("fractional elimination rate
     of insulin from blood"). The Arduino code used 0.02. 0.138 is also the
     value that makes eq (7)'s I(0) = u(0)/(0.01656*BW) an actual steady
     state (0.138 * VI-coefficient 0.12 = 0.01656), which only holds together
     with fix #2 below.

  2. Insulin absorption ODE: paper eq (6) s1(0)=s2(0)=tau_s*u(0) is the
     steady-state of the CLASSIC 2-compartment delay model:
         ds1/dt = u - s1/tau_s
         ds2/dt = (s1 - s2)/tau_s
         uI = s2/tau_s
     The Arduino code instead used a different absorption ODE (rate
     constants k21, kd, ka) that doesn't appear anywhere in the paper's
     parameter tables. That mismatch is why I(0) never actually held: the
     code's absorption dynamics don't reduce to steady state at the paper's
     initial values. Reverted to the classic tau_s form.

  3. Renal threshold: the paper defines TWO distinct thresholds that the
     Arduino code collapsed into one `Gth=60`:
       - Gth  = 162 mg/dl  (Table 2, "Renal threshold of glucose", used for
         glomerular filtration / renal glucose loss -- matches classic
         Hovorka's 9 mmol/L * 18 = 162 exactly)
       - Gth1 = 60  mg/dl  (Table 4, "Hypoglycemic threshold", used only in
         the glucagon secretion (Srh) piecewise formula)
     Using 60 mg/dl as the renal threshold caused constant renal glucose
     loss even at/below normal fasting glucose, which is not physiological
     and not what the paper specifies.

  4. Interstitial glucose initial condition: paper Table 2 gives
     G1(0) = 70 mg/dl, distinct from G(0) = 90 mg/dl. The Arduino code set
     Gt(0) = G(0) = 90 (per an explicit code comment), which is contrary to
     the paper -- the two compartments are meant to start apart, producing a
     normal, bounded initial transient (not an error).

  5. EGP(0) = EGPb = 1.23 mg/dl/min per Table 3, rather than the code's
     0.0161*BW (=1.127 for 70kg) -- close but not the paper's stated value.

  6. G6P(0): paper eq (13) G6P(0) = EGPb/kG6P + g6po, computed here from the
     paper's own EGPb/kG6P/g6po values rather than the hardcoded 41.897 in
     the Arduino source (numerically close, ~41.68, but derived properly).

  7. Sit/Sid/Sie are the CLASSIC, unscaled values (51.2/10000, 8.2/10000,
     520/10000). These are confirmed correct because kb1=Sit*ka1,
     kb2=Sid*ka2, kb3=Sie*ka3 computed from them match the paper's directly
     stated kb1/kb2/kb3 (Table 4) to 4+ significant figures. My earlier
     rescaling of Sit/Sid/Sie (to patch the ke=0.02 instability) contradicted
     the paper and has been fully reverted -- fixes #1 and #2 above resolve
     the instability without needing to touch Sit/Sid/Sie at all.

STILL OPEN / UNVERIFIED (flagged, not silently changed):
  - sigma: paper Table 4 lists sigma = 1.714410e-11 in units
    "ng/L/min/mg/dL/pmol" -- a very different value and a very different unit
    system (pmol/ng) than the code's sigma=1e-9. Left as the code's original
    value here since plugging in the paper's number without knowing the
    intended unit conversion could silently break the glucagon term instead
    of fixing it. Needs a direct question to your supervisor about the unit
    system used for I(t)/H(t) in that specific equation.
  - SRbH (Srhb): paper Table 4 lists SRbH = 0.0007 mg/dl/min, but the
    derivation in the paper text states SRbH = SRSHb = n*Hb, which computes
    to 5.8e-8, not 0.0007 (0.0007 exactly equals kp2's value in the row
    directly above it in Table 4 -- likely a table transcription artifact).
    Using the analytically-derived n*Hb here; worth confirming visually
    against the actual PDF page.
"""


class Constants:
    def __init__(self, BW: float = 70.0, u_basal: float = 12.9127):
        self.BW = BW

        # ---- meal absorption ----
        self.tau_d = 40.0
        self.Ag = 0.8
        self.Mwg = 180.0

        # ---- glucose subsystem ----
        self.F01 = 0.00097 * BW      # Fii in the paper
        self.Vg = 0.16 * BW
        self.k12 = 0.066
        self.G = 90.0                 # mg/dL, paper Table 2 G(0)
        self.Gb = 90.0                 # mg/dL
        self.Gt = 70.0                 # mg/dL, paper Table 2 G1(0) -- NOT equal to G(0)

        # ---- insulin subsystem (classic tau_s-based absorption, paper eq 6) ----
        self.Vi = 0.12 * BW
        self.tau_s = 55.0
        self.u_basal = u_basal
        u = u_basal

        self.S1 = self.tau_s * u          # paper eq 6
        self.S2 = self.tau_s * u          # paper eq 6

        # ---- insulin action rate constants ----
        self.ka1 = 0.006
        self.ka2 = 0.06
        self.ka3 = 0.03
        self.ke = 0.138                    # paper Table 4 (was 0.02 in the Arduino code)

        self.Sit = 51.2 / 10000            # classic, unscaled -- confirmed via kb1 cross-check
        self.Sid = 8.2 / 10000
        self.Sie = 520 / 10000
        self.kb1 = self.Sit * self.ka1
        self.kb2 = self.Sid * self.ka2
        self.kb3 = self.Sie * self.ka3

        self.I = u / (0.01656 * BW)        # paper eq 7 -- now a true steady state (see fixes 1,2)
        self.x1 = 0.30898 * u / BW         # paper eq 8
        self.x2 = 0.04951 * u / BW          # paper eq 9
        self.x3 = 3.2206 * u / BW          # paper eq 10

        # ---- renal / F01 saturation ----
        self.kp2 = 0.0007
        self.ke1 = 0.007
        self.Gth = 162.0     # renal threshold (Table 2) -- was incorrectly 60 in the Arduino code
        self.Gth1 = 60.0     # hypoglycemic threshold (Table 4), used only in the Srh formula

        # ---- glycogenolysis / glucagon-driven EGP extension ----
        self.K6gp = 0.034
        self.Sc = 297.0
        self.cth = 8 / 1_000_000     # Hth, glucagon threshold
        self.tD = 59.9
        self.z = 23.24
        self.Ggng1b = 0.495
        self.Ggg1b = 0.7425

        EGPb = 1.23
        g6po = 5.50
        self.EGP = EGPb                          # paper Table 3 (was 0.0161*BW in the Arduino code)
        self.G6p = EGPb / self.K6gp + g6po        # paper eq 13

        # ---- glucagon (H(t) in the paper; called `c` in the Arduino code) ----
        self.n = 0.01
        self.cb = 58 / 10_000_000   # Hb
        self.rho = 0.86
        self.sigma = 0.01 / 10_000_000   # UNVERIFIED -- see module docstring
        self.S_glyc = 0.98 / 10_000_000  # delta in the paper, matches code
        self.c = self.cb                 # H(0) = Hb
        self.Srhs = 0.0
        self.Srhd = 0.0
        self.Srhb = 0.0                  # matches Arduino source's `double Srhb=0;` (gets
                                          # overwritten by n*Hb every step anyway, per the
                                          # faithfully-preserved missing-braces bug)
        self.Srh = 0.0

        # ---- meal absorption states ----
        self.Dm1 = 0.0
        self.Dm2 = 0.0
        self.Ug = 0.0


# ---- Default schedule ----
# Table 5/6 give real clock times (8AM/1PM/4PM/7PM/11PM) and the actual carb
# amounts / insulin rates. But the reference "Proposed" plot you were shown
# was generated using the reference C code's internal offsets (210, 390, 600,
# 765, 900 minutes from simulation start), not those clock times converted to
# minutes-since-midnight. Using clock times shifts the whole curve ~4-8 hours
# later and no longer matches the reference plot's peak locations. So: keep
# Table 5/6's carb amounts and insulin rates, but use the original code's
# timing offsets so the shape/timing matches the reference plot.
DEFAULT_MEALS = [
    # time_min, duration_min, amount_g
    (210, 10, 124.17),
    (390, 20, 34.5),
    (600, 10, 34.5),
    (765, 20, 76.95),
    (900, 10, 34.5),
]

DEFAULT_INSULIN = [
    # time_min, duration_min, rate_mU_per_min
    (210, 5, 700.0),
    (390, 5, 250.0),
    (600, 5, 200.0),
    (765, 5, 490.0),
    (900, 5, 210.0),
]

DEFAULT_TOTAL_TIME = 1500   # minutes, matches reference plot x-axis
DEFAULT_STEP = 0.1          # minutes, matches Egp6_model.txt
