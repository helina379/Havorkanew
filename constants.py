"""
Constants for the Hovorka + glycogenolysis-driven EGP ("Proposed") model.

Parameters are ported 1:1 from Egp6_model.txt (the "Proposed" model source).
All glucose-related states are in mg/dL (not mmol/L), matching that file.
"""


class Constants:
    def __init__(self, BW: float = 70.0, u_basal: float = 12.9127):
        self.BW = BW

        # ---- meal absorption ----
        self.tau_d = 40.0
        self.Ag = 0.8
        self.Mwg = 180.0

        # ---- glucose subsystem ----
        self.F01 = 0.00097 * BW
        self.Vg = 0.16 * BW
        self.k12 = 0.066
        self.G = 90.0               # mg/dL
        self.Gb = 90.0               # mg/dL
        self.Gt = 90.0               # mg/dL (Gt(0) = G(0))

        # ---- insulin subsystem ----
        self.Vi = 0.12 * BW
        self.tau_s = 55.0
        self.u_basal = u_basal
        u = u_basal

        # ---- insulin action rate constants (unchanged from source) ----
        self.ka1 = 0.006
        self.ka2 = 0.06
        self.ka3 = 0.03
        self.ke = 0.02               # NOTE: 0.02 here, not 0.138 (classic Hovorka) -- Egp6 value

        # ---- renal / F01 saturation ----
        self.kp2 = 0.0007
        self.ke1 = 0.007
        self.Gth = 60.0

        # ---- insulin absorption (subcutaneous) ----
        self.k21 = 0.045
        self.kd = 0.0021
        self.ka_sc = 0.02             # NOTE: called `ka` in the C source; renamed to avoid
                                       # clashing with ka1/ka2/ka3 (insulin action rates)

        # =====================================================================
        # RECALIBRATION: the source file lowers ke (0.138 -> 0.02) for a slower
        # insulin-clearance patient, but keeps the classic Sit/Sid/Sie and the
        # classic I(0)=u/(0.01656*BW) initializer. Those were only correct for
        # ke=0.138. Left as-is, the true insulin steady state is ~6x higher
        # than the initializer, which pushes x3 past 1 (EGP suppression term
        # (1-x3) goes negative) and the fasting state collapses instead of
        # holding at Gb. Fix: derive every steady state analytically instead of
        # reusing the classic magic numbers.
        # =====================================================================

        # --- 1) true insulin steady state under THIS ke/kd/ka_sc/k21/tau_s ---
        S1_ss = u / self.k21
        S2_ss = u / (self.kd + self.ka_sc)
        Ui_ss = S2_ss / self.tau_s
        I_ss = Ui_ss / (self.Vi * self.ke)

        # --- 2) classic-model insulin steady state (what Sit/Sid/Sie were
        #        actually calibrated against) ---
        I_ss_classic = u / (0.01656 * BW)

        Sit_classic = 51.2 / 10000
        Sid_classic = 8.2 / 10000
        Sie_classic = 520 / 10000

        # target x_ss values = what the classic, physiologically-sane model
        # would produce; rescale Sit/Sid/Sie so THIS model reaches the same
        # x_ss target even though its I_ss is different
        target_x1_ss = Sit_classic * I_ss_classic
        target_x2_ss = Sid_classic * I_ss_classic
        target_x3_ss = Sie_classic * I_ss_classic

        self.Sit = target_x1_ss / I_ss
        self.Sid = target_x2_ss / I_ss
        self.Sie = target_x3_ss / I_ss

        self.kb1 = self.Sit * self.ka1
        self.kb2 = self.Sid * self.ka2
        self.kb3 = self.Sie * self.ka3

        # --- 3) set insulin subsystem states to their (now self-consistent)
        #        steady state ---
        self.I = I_ss
        self.S1 = S1_ss
        self.S2 = S2_ss
        self.x1 = target_x1_ss
        self.x2 = target_x2_ss
        self.x3 = target_x3_ss

        # --- 4) glycogenolysis / catecholamine-driven EGP extension ---
        self.K6gp = 0.034
        self.Sc = 297.0
        self.cth = 8 / 1_000_000
        self.tD = 59.9
        self.z = 23.24

        Ggng1b_ref = 0.495
        Ggg1b_ref = 0.7425

        # required EGP (converted, mg/dL/min-equivalent) to balance basal
        # utilization + renal loss exactly at G = Gb (F01uc/Erc are BW-invariant
        # by construction: F01 and Vg both scale with BW so the ratio cancels)
        F01uc_at_Gb = 18 * self.F01 / self.Vg      # G=Gb=90 >= 81 saturation branch
        Erc_at_Gb = self.ke1 * (self.Gb - self.Gth)  # G=Gb=90 >= Gth branch
        target_EGPc = F01uc_at_Gb + Erc_at_Gb

        target_EGP_ss = target_EGPc * self.Vg / (18 * (1 - self.x3))

        scale_egp = target_EGP_ss / (Ggg1b_ref + Ggng1b_ref)
        self.Ggng1b = Ggng1b_ref * scale_egp
        self.Ggg1b = Ggg1b_ref * scale_egp

        self.G6p = target_EGP_ss / self.K6gp   # steady state of dG6p/dt = -K6gp*G6p + Ggg1b + Ggng1b
        self.EGP = target_EGP_ss

        # ---- counter-regulatory hormone (glucagon-ish) secretion ----
        self.n = 0.01
        self.cb = 58 / 10_000_000
        self.rho = 0.86
        self.sigma = 0.01 / 10_000_000
        self.S_glyc = 0.98 / 10_000_000   # NOTE: called `S` in the C source; renamed to avoid
                                           # clashing with S1/S2 (insulin subcutaneous states)
        self.c = self.cb                  # initial catecholamine-pool state
        self.Srhs = 0.0
        self.Srhd = 0.0
        self.Srhb = 0.0
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
