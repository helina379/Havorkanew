"""
Proposed model: Hovorka glucose/insulin kinetics extended with a
glycogenolysis / catecholamine-driven EGP subsystem.

This is a line-by-line port of the `rungeKutta4Solver_custom()` function in
Egp6_model.txt, generalized to accept arbitrary meal and insulin schedules.

IMPORTANT - faithful bug reproduction:
The original C source has an unbraced if/else around Srh:

    if (G>=Gb)
    Srhs=rho*(Srhs-Srhb);
    else
    Srhb=n*cb;
    Srhs=rho*(Srhs - max(sigma*(Gth-G)/(I+1)+Srhb ,0));

In C, only `Srhb=n*cb;` belongs to the else-branch. The final `Srhs=...` line
runs UNCONDITIONALLY every step, so the `if(G>=Gb)` assignment to Srhs is
always immediately overwritten and never actually has any effect. That
behavior is intentionally preserved below (per your call) so this matches
the reference "Proposed" plot exactly. It's not fixed to be "correct" if/else
behavior.
"""

import math
from constants import Constants


def _max0(x):
    return x if x > 0 else 0.0


class Body:
    def __init__(self, constants: Constants):
        self.c = constants

    def _meal_and_insulin_at(self, t, meals, insulin):
        """meals/insulin: list of (start_min, duration_min, rate) tuples."""
        d_cho = 0.0
        for start, dur, rate in meals:
            if start <= t <= start + dur:
                d_cho = rate
                break

        u = self.c.u_basal
        for start, dur, rate in insulin:
            if start <= t <= start + dur:
                u = rate
                break

        return d_cho, u

    def _step(self, t, step, d_cho, u):
        c = self.c

        D_meal = 1000 * d_cho / c.Mwg

        # ---- meal (two-compartment) absorption, RK4 ----
        rkm11 = step * (c.Ag * D_meal - c.Dm1 / c.tau_d)
        rkm21 = step * (c.Dm1 / c.tau_d - c.Dm2 / c.tau_d)

        rkm12 = step * (c.Ag * D_meal - (c.Dm1 + rkm11 / 2) / c.tau_d)
        rkm22 = step * ((c.Dm1 + rkm11 / 2) / c.tau_d - (c.Dm2 + rkm21 / 2) / c.tau_d)

        rkm13 = step * (c.Ag * D_meal - (c.Dm1 + rkm12 / 2) / c.tau_d)
        rkm23 = step * ((c.Dm1 + rkm12 / 2) / c.tau_d - (c.Dm2 + rkm22 / 2) / c.tau_d)

        rkm14 = step * (c.Ag * D_meal - (c.Dm1 + rkm13) / c.tau_d)
        rkm24 = step * ((c.Dm1 + rkm13) / c.tau_d - (c.Dm2 + rkm23) / c.tau_d)

        Ug = c.Dm2 / c.tau_d   # uses PRE-update Dm2, matches source ordering
        c.Dm1 = c.Dm1 + (rkm11 + 2 * rkm12 + 2 * rkm13 + rkm14) / 6
        c.Dm2 = c.Dm2 + (rkm21 + 2 * rkm22 + 2 * rkm23 + rkm24) / 6

        # ---- glucose subsystem terms (held fixed across the 4 RK stages,
        #      exactly as in the source) ----
        Ugc = 18 * Ug / c.Vg
        EGPc = 18 * (c.EGP / c.Vg) * (1 - c.x3) - c.kp2 * (c.G - c.Gb)

        Erc = 0.0
        if c.G >= c.Gth:
            Erc = c.ke1 * (c.G - c.Gth)

        F01uc = 18 * c.F01 * c.G / (c.Vg * 81)
        if c.G >= 81:
            F01uc = 18 * c.F01 / c.Vg

        G0, Gt0, Gb, x1, x2, k12 = c.G, c.Gt, c.Gb, c.x1, c.x2, c.k12

        rkg11 = step * (Ugc - F01uc - Erc + k12 * (Gt0 - Gb) - x1 * (G0 - Gb) + EGPc)
        rkg21 = step * (x1 * (G0 - Gb) - (k12 + x2) * (Gt0 - Gb))

        rkg12 = step * (Ugc - F01uc - Erc + k12 * ((Gt0 + rkg21 / 2) - Gb) - x1 * ((G0 + rkg11 / 2) - Gb) + EGPc)
        rkg22 = step * (x1 * ((G0 + rkg11 / 2) - Gb) - (k12 + x2) * ((Gt0 + rkg21 / 2) - Gb))

        rkg13 = step * (Ugc - F01uc - Erc + k12 * ((Gt0 + rkg22 / 2) - Gb) - x1 * ((G0 + rkg12 / 2) - Gb) + EGPc)
        rkg23 = step * (x1 * ((G0 + rkg12 / 2) - Gb) - (k12 + x2) * ((Gt0 + rkg22 / 2) - Gb))

        rkg14 = step * (Ugc - F01uc - Erc + k12 * ((Gt0 + rkg23) - Gb) - x1 * ((G0 + rkg13) - Gb) + EGPc)
        rkg24 = step * (x1 * ((G0 + rkg13) - Gb) - (k12 + x2) * ((Gt0 + rkg23) - Gb))

        # ---- glycogenolysis drive ----
        E = (1 - math.tanh((t - c.tD) / c.z)) / 2
        Ggg = c.Ggg1b + c.Sc * _max0(c.c - c.cth) * E

        dGt_scaled = (rkg21 + 2 * rkg22 + 2 * rkg23 + rkg24) / 6 * step
        if dGt_scaled >= 0:
            new_EGP = c.K6gp * c.G6p - c.x3 * dGt_scaled - c.kp2 * (c.G - c.Gb)
        else:
            new_EGP = c.K6gp * c.G6p - c.kp2 * (c.G - c.Gb)

        # ---- Srh (faithful bug-for-bug: see module docstring) ----
        if c.G >= c.Gb:
            c.Srhs = c.rho * (c.Srhs - c.Srhb)
        else:
            c.Srhb = c.n * c.cb
        c.Srhs = c.rho * (c.Srhs - _max0(c.sigma * (c.Gth1 - c.G) / (c.I + 1) + c.Srhb))

        c.Srhd = c.S_glyc * _max0(-dGt_scaled)
        c.Srh = c.Srhs + c.Srhd

        rkc1 = step * (-c.n * c.c + c.Srh)
        rkc2 = step * (-c.n * (c.c + rkc1 / 2) + c.Srh)
        rkc3 = step * (-c.n * (c.c + rkc2 / 2) + c.Srh)
        rkc4 = step * (-c.n * (c.c + rkc3) + c.Srh)
        c.c = c.c + (rkc1 + 2 * rkc2 + 2 * rkc3 + rkc4) / 6

        c.G = c.G + (rkg11 + 2 * rkg12 + 2 * rkg13 + rkg14) / 6
        c.Gt = c.Gt + (rkg21 + 2 * rkg22 + 2 * rkg23 + rkg24) / 6

        rkg61 = step * (-c.K6gp * c.G6p + Ggg + c.Ggng1b)
        rkg62 = step * (-c.K6gp * (c.G6p + rkg61 / 2) + Ggg + c.Ggng1b)
        rkg63 = step * (-c.K6gp * (c.G6p + rkg62 / 2) + Ggg + c.Ggng1b)
        rkg64 = step * (-c.K6gp * (c.G6p + rkg63) + Ggg + c.Ggng1b)
        c.G6p = c.G6p + (rkg61 + 2 * rkg62 + 2 * rkg63 + rkg64) / 6

        c.EGP = new_EGP

        # ---- insulin (subcutaneous absorption), RK4 ----
        # classic tau_s-based 2-compartment delay model (paper eq 6), NOT the
        # k21/kd/ka model that was in the Arduino source -- see module
        # docstring in constants.py for why.
        S1, S2, tau_s = c.S1, c.S2, c.tau_s

        rki11 = step * (u - S1 / tau_s)
        rki21 = step * ((S1 - S2) / tau_s)

        rki12 = step * (u - (S1 + rki11 / 2) / tau_s)
        rki22 = step * (((S1 + rki11 / 2) - (S2 + rki21 / 2)) / tau_s)

        rki13 = step * (u - (S1 + rki12 / 2) / tau_s)
        rki23 = step * (((S1 + rki12 / 2) - (S2 + rki22 / 2)) / tau_s)

        rki14 = step * (u - (S1 + rki13) / tau_s)
        rki24 = step * (((S1 + rki13) - (S2 + rki23)) / tau_s)

        Ui = c.S2 / c.tau_s
        c.S1 = c.S1 + (rki11 + 2 * rki12 + 2 * rki13 + rki14) / 6
        c.S2 = c.S2 + (rki21 + 2 * rki22 + 2 * rki23 + rki24) / 6

        # ---- insulin action, RK4 ----
        I0 = c.I
        rkx11 = step * (-c.ka1 * c.x1 + c.kb1 * I0)
        rkx21 = step * (-c.ka2 * c.x2 + c.kb2 * I0)
        rkx31 = step * (-c.ka3 * c.x3 + c.kb3 * I0)

        rkx12 = step * (-c.ka1 * (c.x1 + rkx11 / 2) + c.kb1 * I0)
        rkx22 = step * (-c.ka2 * (c.x2 + rkx21 / 2) + c.kb2 * I0)
        rkx32 = step * (-c.ka3 * (c.x3 + rkx31 / 2) + c.kb3 * I0)

        rkx13 = step * (-c.ka1 * (c.x1 + rkx12 / 2) + c.kb1 * I0)
        rkx23 = step * (-c.ka2 * (c.x2 + rkx22 / 2) + c.kb2 * I0)
        rkx33 = step * (-c.ka3 * (c.x3 + rkx32 / 2) + c.kb3 * I0)

        rkx14 = step * (-c.ka1 * (c.x1 + rkx13) + c.kb1 * I0)
        rkx24 = step * (-c.ka2 * (c.x2 + rkx23) + c.kb2 * I0)
        rkx34 = step * (-c.ka3 * (c.x3 + rkx33) + c.kb3 * I0)

        c.x1 = c.x1 + (rkx11 + 2 * rkx12 + 2 * rkx13 + rkx14) / 6
        c.x2 = c.x2 + (rkx21 + 2 * rkx22 + 2 * rkx23 + rkx24) / 6
        c.x3 = c.x3 + (rkx31 + 2 * rkx32 + 2 * rkx33 + rkx34) / 6

        # ---- plasma insulin, RK4 ----
        rkia1 = step * (Ui / c.Vi - c.ke * I0)
        rkia2 = step * (Ui / c.Vi - c.ke * (I0 + rkia1 / 2))
        rkia3 = step * (Ui / c.Vi - c.ke * (I0 + rkia2 / 2))
        rkia4 = step * (Ui / c.Vi - c.ke * (I0 + rkia3))
        c.I = c.I + (rkia1 + 2 * rkia2 + 2 * rkia3 + rkia4) / 6

    def simulate(self, meals, insulin, total_time, step):
        """
        meals / insulin: list of (start_min, duration_min, rate) tuples.
        Returns (t_points, G_points) both as plain lists, G in mg/dL.
        """
        n_steps = int(total_time / step)
        t_points = []
        G_points = []
        t = 0.0

        for _ in range(n_steps):
            d_cho, u = self._meal_and_insulin_at(t, meals, insulin)
            t_points.append(t)
            G_points.append(self.c.G)
            self._step(t, step, d_cho, u)
            t = t + step

        return t_points, G_points
