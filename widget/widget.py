# To-do: be very sure we know what rows vs columns mean
import numpy as np
import pandas as pd
import streamlit as st

import ngm


def app():
    st.title("3-Group NGM Calculator")

    # Group names
    group_names = ["Core", "Kids", "General Population"]

    # Sidebar for inputs
    st.sidebar.header("Model Inputs")

    # Population size
    st.sidebar.subheader("Population Sizes")
    default_values = np.array([0.2, 0.2, 0.595]) * 1000000
    N = np.array(
        [
        st.sidebar.number_input(f"Population ({group})", value=int(default_values[i]), min_value=0)
            for i, group in enumerate(group_names)
        ]
    )

    # Vaccine doses
    st.sidebar.subheader("Vaccine Allocation Strategies")
    starting_vax = st.sidebar.selectbox(
        "Vaccine allocation strategy",
        [
            "All core (low)", "All kids (low)", "Even (low)",
            "All core (high)", "All kids (high)", "Even (high)",
        ]
    )
    ndoses_default = 100000
    if starting_vax.split(r" (")[1][:-1] == "high":
        ndoses_default = 200000

    allocation_default = [0, 0, 0]
    if starting_vax.split(r" (")[0] == "All core":
        allocation_default = np.array([100.0, 0, 0])
    elif starting_vax.split(r" (")[0] == "All kids":
        allocation_default = np.array([0, 100.0, 0])
    elif starting_vax.split(r" (")[0] == "Even":
        allocation_default = 100.0 * N / N.sum()

    st.sidebar.subheader("Vaccine Allocation Customization")
    ndoses = st.sidebar.number_input("Total Number of Doses", value=ndoses_default, min_value=0, max_value=10000000, step=1)

    allocation = []
    remaining = 100.0
    for i, group in enumerate(group_names[:-1]):
        allocation.append(
            st.sidebar.number_input(
                f"Percent of Vaccine Doses going to {group}",
                value=allocation_default[i],
                min_value=0.0,
                max_value=100.0 - sum(allocation[:i]),
                step=1.0
            )
        )
        remaining -= allocation[-1]
    st.sidebar.write(f"(Allocating {remaining:.2f}% to {group_names[-1]})")
    allocation.append(remaining)

    V = np.floor(np.array(allocation) / 100.0 * ndoses).astype("int")

    # Contact matrix
    st.sidebar.subheader("High and low contact rates")
    # Define lo and hi using Streamlit inputs
    lo = st.sidebar.number_input("Low value", value=1)
    hi = st.sidebar.number_input("High value", value=10)

    # Create the contact matrix K
    beta = np.array(
        [
            [hi, lo, lo],  # core
            [lo, lo, lo],  # kids
            [lo, lo, lo],  # general
        ]
    )

    with st.sidebar.expander("Advanced Settings"):
        st.sidebar.subheader("Vaccine efficacy")
        VE = st.sidebar.slider("Vaccine Efficacy", 0.0, 1.0, value=0.7, step=0.01)

    # Perform the NGM calculation
    result = ngm.simulate(
        n=N, n_vax=V, beta=beta, p_severe=np.array([0.02, 0.06, 0.02]), ve=VE
    )

    # Display the adjusted contact matrix
    st.subheader("NGM with vaccination")
    st.write("This matrix reflects the impact of vaccine efficacy and numbers of susceptible individuals:")

    R = pd.DataFrame(
        result["R"],
        columns=group_names,
        index=group_names,
    )

    st.dataframe(R)

    # Display results
    st.subheader("Results")
    st.write(f"R-effective: {result['Re']:.2f}")
    st.write("Distribution of Infections by Group:")
    st.json(
        {group: round(inf, 2) for group, inf in zip(group_names, result["infections"])}
    )
    total_severe = sum(result["severe_infections"])
    st.write(f"Total number of severe infections: {total_severe:.2f}")


if __name__ == "__main__":
    app()
