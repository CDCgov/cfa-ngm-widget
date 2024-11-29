# To-do: be very sure we know what rows vs columns mean
import numpy as np
import pandas as pd
import streamlit as st

import ngm


def app():
    st.title("4-Group NGM Calculator")

    # Group names
    group_names = ["Core", "Kids", "Travelers", "General Population"]

    # Sidebar for inputs
    st.sidebar.header("Model Inputs")

    # Population size
    st.sidebar.subheader("Population Sizes")
    default_values = [100, 100, 50, 1000]
    N = np.array(
        [
            st.sidebar.number_input(f"Population ({group})", value=default_values[i], min_value=0)
            for i, group in enumerate(group_names)
        ]
    )

    # Vaccine doses
    st.sidebar.subheader("Vaccine Doses")
    V = np.array(
        [
            st.sidebar.number_input(
                f"Vaccine Doses ({group})",
                value=0,
                min_value=0,
                max_value=N[i],
            )
            for i, group in enumerate(group_names)
        ]
    )

    # Contact matrix
    st.sidebar.subheader("R0 (low) and R0 (high), entries to NGM (K)")
    # Define lo and hi using Streamlit inputs
    lo = st.sidebar.number_input("Low value", value=0.5)
    hi = st.sidebar.number_input("High value", value=3)

    # Create the contact matrix K
    K = np.array(
        [
            [hi, lo, hi, lo],  # core
            [lo, lo, lo, lo],  # kids
            [hi, lo, lo, lo],  # travelers
            [lo, lo, lo, lo],  # general
        ]
    )

    with st.sidebar.expander("Advanced Settings"):
        st.sidebar.subheader("Vaccine efficacy")
        VE = st.sidebar.slider("Vaccine Efficacy", 0.0, 1.0, value=0.7, step=0.01)

    # Perform the NGM calculation
    result = ngm.simulate(
        n=N, n_vax=V, K=K, p_severe=np.array([0.02, 0.06, 0.02, 0.02]), ve=VE
    )

    # Display the adjusted contact matrix
    st.subheader("NGM with vaccination")
    st.write("This matrix reflects the impact of vaccine efficacy and numbers of susceptible individuals:")

    K_adjusted_df = pd.DataFrame(
        result["reduced_K"],
        columns=group_names,
        index=group_names,
    )

    st.dataframe(K_adjusted_df)

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
