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

    ndoses = 50000
    st.sidebar.subheader("Vaccine Doses (Pre-filled Scenarios)")
    starting_vax = st.sidebar.selectbox("Vaccine allocation strategy", ["No vaccination", "All core", "All kids", "Even"])
    # No vax is default
    allocation = [0, 0, 0]
    if starting_vax == "All core":
        allocation = [1, 0, 0]
    elif starting_vax == "All kids":
        allocation = [0, 1, 0]
    elif starting_vax == "Even":
        allocation = [1/3, 1/3, 1/3]
    allocation = np.floor(np.array(allocation) * ndoses).astype("int")

    st.sidebar.subheader("Vaccine Doses (Customization)")
    V = np.array(
        [
            st.sidebar.number_input(
                f"Vaccine Doses ({group})",
                value=allocation[i],
                min_value=0,
                max_value=N[i],
            )
            for i, group in enumerate(group_names)
        ]
    )

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
