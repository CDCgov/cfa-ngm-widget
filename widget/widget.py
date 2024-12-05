# To-do: be very sure we know what rows vs columns mean
import numpy as np
import pandas as pd
import streamlit as st

import ngm


def app():
    st.title("3-Group NGM Calculator")

    # Information we should be getting from scratch/config.yaml
    p_severe = np.array([0.02, 0.06, 0.02])

    # Group names
    group_names = ["Core", "Kids", "General Population"]
    n_groups = len(group_names)
    # Sidebar for inputs
    st.sidebar.header("Model Inputs")

    # Population size
    st.sidebar.subheader("Population Sizes")
    default_values = np.array([0.05, 0.45, 0.5]) * int(1e7)
    N = np.array(
        [
        st.sidebar.number_input(f"Population ({group})", value=int(default_values[i]), min_value=0)
            for i, group in enumerate(group_names)
        ]
    )

    # Vaccine doses
    st.sidebar.subheader("Vaccine Allocation: Doses")
    ndoses_default = int(1e6)
    ndoses = st.sidebar.number_input("Total Number of Doses", value=ndoses_default, min_value=0, max_value=sum(N))

    st.sidebar.subheader("Vaccine Allocation: Strategies")
    strategy = st.sidebar.selectbox(
        "Vaccine allocation strategy",
        [
            "Core", "Kids", "Even",
        ]
    )

    button_to_core = {"Core" : 0, "Kids" : 1, "Even" : "even"}

    allocation_default = [0] * n_groups
    if ndoses > 0:
        allocation_default = ngm.distribute_vaccines(V=ndoses, N_i=N, strategy=button_to_core[strategy])
        allocation_default = 100 * allocation_default / allocation_default.sum()

    st.sidebar.subheader("Vaccine Allocation: Customization")
    allocation = []
    remaining = 100.0
    for i, group in enumerate(group_names[:-1]):
        this_max = 100.0 - sum(allocation[:i])
        if (this_max / 100.0) * ndoses > N[i]:
            this_max = N[i] / ndoses * 100

        allocation.append(
            st.sidebar.number_input(
                f"Percent of Vaccine Doses going to {group}",
                value=allocation_default[i],
                min_value=0.0,
                max_value=this_max,
                step=1.0
            )
        )
        remaining -= allocation[-1]
    st.sidebar.write(f"(Allocating {remaining:.2f}% to {group_names[-1]})")
    allocation.append(remaining)

    V = np.floor(np.array(allocation) / 100.0 * ndoses).astype("int")

    # Contact matrix
    st.sidebar.subheader("Next Generation Matrix")
    st.sidebar.write("For a single new infection of category `from`, specify how many infections of category `to` it will make.")

    from_to = [((i, group_names[i]), (j, group_names[j]),) for i in range(n_groups) for j in range(n_groups)]
    r_default = np.array([[3.0, 0.0, 0.2], [0.10, 1.0, 0.5], [0.25, 1.0, 1.5]])

    r_novax = np.zeros((n_groups, n_groups,))
    for ft in from_to:
        row = ft[1][0]
        col = ft[0][0]
        r_novax[row, col] = st.sidebar.number_input(
            f"From {ft[0][1]} to {ft[1][1]}",
            value=r_default[row, col],
            min_value=0.0,
            max_value=100.0
        )

    with st.sidebar.expander("Advanced Settings"):
        st.sidebar.subheader("Vaccine efficacy")
        VE = st.sidebar.slider("Vaccine Efficacy", 0.0, 1.0, value=0.7, step=0.01)

    # Perform the NGM calculation
    result = ngm.simulate(
        R_novax=r_novax, n=N, n_vax=V, p_severe=p_severe, ve=VE
    )

    # Display the adjusted contact matrix
    st.subheader("Results with vaccination:")

    R = pd.DataFrame(
        np.round(result["R"], 2),
        columns=group_names,
        index=group_names,
    )

    st.write("Next Generation Matrix:")
    st.dataframe(R)

    # Display results
    st.write(f"R-effective: {result['Re']:.2f}")
    st.write("Proportion of infections in each group:")
    st.dataframe(
        pd.DataFrame(
            [np.round(result["infections"], 2)],
            columns=group_names,
        ).style.hide(axis="index")
    )
    st.write(f"Infection fatality ratio: {(p_severe * result['infections']).sum():.3f}")

    st.subheader("Counterfactual (no vaccination):")
    st.write("Next Generation Matrix:")
    st.dataframe(
        pd.DataFrame(
            np.round(r_novax, 2),
            columns=group_names,
            index=group_names,
        )
    )
    novax_eigen = ngm.dominant_eigen(r_novax)
    st.write(f"R0: {novax_eigen.value:.2f}")
    st.write("Proportion of infections in each group:")
    st.dataframe(
        pd.DataFrame(
            [np.round(novax_eigen.vector, 2)],
            columns=group_names,
        ).style.hide(axis="index")
    )
    st.write(f"Infection fatality ratio: {(p_severe * novax_eigen.vector).sum():.3f}")


if __name__ == "__main__":
    app()
