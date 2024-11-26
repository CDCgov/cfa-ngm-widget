# To-do: be very sure we know what rows vs columns mean
import numpy as np
import pandas as pd
import streamlit as st


def ngm_sir(N, V, K, VE, p_severe):
    """
    Function to calculate next-generation matrix (NGM) for a 4-group SIR model

    Args:
        N (array): Population sizes for each group
        V (array): Number of vaccine doses administered to each group
        K (2D array): Contact matrix (4x4); (flexible?)
        VE (float): Vaccine efficacy, all or nothing
        p_severe (array): Group spcfeic probability of severe infection

    Returns:
        dict: Contains R-effective, total infections, severe infections, adjusted K matrix
    """
    # check sizes of inputs? do we need this
    assert np.all(N >= V), "Vaccinated cannot exceed population size"
    assert len(N) == len(V) == len(p_severe), "Input dimensions must match"

    # Calculate susceptibles for NGM
    S = N - VE * V
    K_adjusted = K * S / N

    # Eigenvalue computation for NGM
    eigenvalues, eigenvectors = np.linalg.eig(K_adjusted)
    r_effective = max(eigenvalues)

    # Calculate infections and severe infections
    dominant_vector = eigenvectors[:, np.argmax(eigenvalues)]
    infections = dominant_vector * N
    severe_infections = infections * p_severe

    return {
        "R_effective": r_effective,
        "Infections": infections.real,
        "Severe_Infections": severe_infections.real,
        "K_adjusted": K_adjusted,
    }


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

    # Vaccine efficacy
    VE = st.sidebar.slider("Vaccine Efficacy", 0.0, 1.0, value=0.7, step=0.01)

    # Probability of severe infections
    st.sidebar.subheader("Severe Infection Probabilities")
    p_severe = np.array(
        [
            st.sidebar.slider(
                f"Probability of Severe Infection ({group})",
                0.0,
                1.0,
                value=0.03,
                step=0.01,
            )
            for group in group_names
        ]
    )

    # Contact matrix
    st.sidebar.subheader("NGM (4x4)")
    K = np.array(
        [
            [
                st.sidebar.number_input(
                    f"K[{group_names[i]} → {group_names[j]}]",
                    value=1 if i != j else 3,
                )
                for j in range(4)
            ]
            for i in range(4)
        ]
    )

    # Perform the NGM calculation
    if st.sidebar.button("Calculate"):
        result = ngm_sir(N, V, K, VE, p_severe)

        # Display the adjusted contact matrix
        st.subheader("NGM with vaccination")
        st.write(
            "This matrix reflects the impact of vaccine efficacy and susceptibility:"
        )

        K_adjusted_df = pd.DataFrame(
            result["K_adjusted"],
            columns=group_names,
            index=group_names,
        )

        st.dataframe(K_adjusted_df)

        # Display results
        st.subheader("Results")
        st.write(f"R-effective: {result['R_effective']:.2f}")
        st.write("Infections by Group:")
        st.json(
            {
                group: round(inf, 2)
                for group, inf in zip(group_names, result["Infections"])
            }
        )
        st.write("Severe Infections by Group:")
        st.json(
            {
                group: round(sev, 2)
                for group, sev in zip(group_names, result["Severe_Infections"])
            }
        )


if __name__ == "__main__":
    app()
