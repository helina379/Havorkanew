import matplotlib.pyplot as plt
import streamlit as st

from constants import Constants, DEFAULT_MEALS, DEFAULT_INSULIN, DEFAULT_TOTAL_TIME, DEFAULT_STEP
from model import Body


def main():
    st.set_page_config(page_title="Proposed Glucose-Insulin Model", layout="centered")
    st.markdown("<h1 style='text-align: center;'>Proposed Model (Hovorka + EGP)</h1>", unsafe_allow_html=True)
    st.markdown("<br>", unsafe_allow_html=True)

    with st.sidebar:
        st.header("Patient")
        BW = st.number_input("Body weight (kg)", min_value=20.0, value=70.0, step=1.0)
        u_basal = st.number_input("Basal insulin infusion (mU/min)", min_value=0.0, value=12.9127, format="%.4f")

        st.header("Simulation")
        total_time = st.number_input("Total simulation time (min)", min_value=60, value=DEFAULT_TOTAL_TIME, step=10)
        step = st.number_input("Solver step (min)", min_value=0.01, value=DEFAULT_STEP, step=0.01, format="%.2f")

    st.markdown("### Meals")
    st.caption("Time = minutes since simulation start (00:00). Amount = total grams of carbohydrate for that meal.")
    num_meals = st.number_input("Number of meals", min_value=0, max_value=8, value=len(DEFAULT_MEALS))

    meals = []
    for i in range(num_meals):
        default = DEFAULT_MEALS[i] if i < len(DEFAULT_MEALS) else (0, 10, 0.0)
        cols = st.columns(3)
        t = cols[0].number_input(f"Meal {i+1} time (min)", min_value=0, value=int(default[0]), key=f"meal_t_{i}")
        d = cols[1].number_input(f"Meal {i+1} duration (min)", min_value=1, value=int(default[1]), key=f"meal_d_{i}")
        g = cols[2].number_input(f"Meal {i+1} carbs (g)", min_value=0.0, value=float(default[2]), key=f"meal_g_{i}")
        rate = g / d if d > 0 else 0.0
        meals.append((t, d, rate))

    st.markdown("### Insulin boluses")
    st.caption("Time = minutes since simulation start. Rate = exogenous insulin infusion rate during the bolus window.")
    num_insulin = st.number_input("Number of insulin boluses", min_value=0, max_value=8, value=len(DEFAULT_INSULIN))

    insulin = []
    for i in range(num_insulin):
        default = DEFAULT_INSULIN[i] if i < len(DEFAULT_INSULIN) else (0, 5, 0.0)
        cols = st.columns(3)
        t = cols[0].number_input(f"Bolus {i+1} time (min)", min_value=0, value=int(default[0]), key=f"ins_t_{i}")
        d = cols[1].number_input(f"Bolus {i+1} duration (min)", min_value=1, value=int(default[1]), key=f"ins_d_{i}")
        rate = cols[2].number_input(f"Bolus {i+1} rate (mU/min)", min_value=0.0, value=float(default[2]), key=f"ins_r_{i}")
        insulin.append((t, d, rate))

    st.markdown("<br>", unsafe_allow_html=True)

    if st.button("Run simulation", type="primary"):
        constants = Constants(BW=BW, u_basal=u_basal)
        body = Body(constants)
        t_points, G_points = body.simulate(meals, insulin, total_time=total_time, step=step)

        fig, ax = plt.subplots(figsize=(9, 5))
        ax.plot(t_points, G_points, label="Proposed", color="#1f77b4", linewidth=2)
        ax.set_title("Glucose Profile")
        ax.set_xlabel("Time (min)")
        ax.set_ylabel("Glucose (mg/dl)")
        ax.legend()
        ax.grid(alpha=0.3)
        st.pyplot(fig)

        st.download_button(
            "Download glucose trace (CSV)",
            data="time_min,glucose_mgdl\n" + "\n".join(f"{t:.2f},{g:.4f}" for t, g in zip(t_points, G_points)),
            file_name="proposed_glucose_trace.csv",
            mime="text/csv",
        )


if __name__ == "__main__":
    main()
