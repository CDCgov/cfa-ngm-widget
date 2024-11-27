# To-do: be very sure we know what rows vs columns mean
import numpy as np
import pandas as pd
import streamlit as st

from ngm import ngm_sir


def app():
    st.title("4-Group SIR Model NGM Calculator")

    # Group names
    group_names = ["Core", "Kids", "Travelers", "General Population"]

    # Sidebar for inputs
    st.sidebar.header("Model Inputs")

    # Population size
    st.sidebar.subheader("Population Sizes")
    N = np.array(
        [
            st.sidebar.number_input(
                f"Population ({group})", value=100, min_value=0
            )
            for group in group_names
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
        VE = st.sidebar.slider(
            "Vaccine Efficacy", 0.0, 1.0, value=0.7, step=0.01
        )

    # Perform the NGM calculation
    if st.sidebar.button("Calculate"):
        result = ngm_sir(n=N, doses=V, r=K, v_e=VE)

        # Display the adjusted contact matrix
        st.subheader("NGM with vaccination")
        st.write(
            "This matrix reflects the impact of vaccine efficacy and susceptibility:"
        )

        K_adjusted_df = pd.DataFrame(
            result["ngm_adjusted"],
            columns=group_names,
            index=group_names,
        )

        st.dataframe(K_adjusted_df)

        # Display results
        st.subheader("Results")
        st.write(f"R-effective: {result['dominant_eigenvalue']:.2f}")
        st.write("Distribution of Infections by Group:")
        st.json(
            {
                group: round(inf, 2)
                for group, inf in zip(
                    group_names, result["dominant_eigenvector"]
                )
            }
        )


if __name__ == "__main__":
    app()
