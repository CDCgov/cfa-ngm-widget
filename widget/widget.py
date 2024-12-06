# To-do: be very sure we know what rows vs columns mean
import numpy as np
import pandas as pd
import streamlit as st

import ngm
from scratch.simulate import simulate_scenario

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

        st.sidebar.subheader("Generations of spread")
        G = st.sidebar.slider("Generations", 1, 10, value=10, step=1)


    params = {
        "n_total": N.sum(),
        "pop_props": N/N.sum(),
        "R_novax": r_novax,
        "p_severe": p_severe,
        "n_vax": V,
        "ve": VE,
        "G": G,
    }


    # Run the simulation with vaccination
    result = simulate_scenario(params, return_polars=False)

    st.subheader("Results with vaccination:")

    st.write(f"R-effective: {float(result['Re']):.2f}")
    st.write("Proportion of infections in each group:")
    infections = [1.0, result["infections_core"], result["infections_children"], result["infections_adults"]]
    fatalities_after_G_generations =  [result["deaths_after_G_generations"], result["deaths_after_G_generations_core"], result["deaths_after_G_generations_children"], result["deaths_after_G_generations_adults"]]

    st.dataframe(
        pd.DataFrame(
            [np.round(infections, 2),
             np.round(fatalities_after_G_generations, 1)],
            columns=["Total", *group_names],
        ).style.hide(axis="index")
    )
    st.write(f"Infection fatality ratio: {result['ifr']:.3f}")


    # Run counterfactural scenario
    params_no_vax = params.copy()
    params_no_vax["n_vax_total"] = 0
    result_no_vax = simulate_scenario(params_no_vax, return_polars=False)


    st.subheader("Counterfactual (no vaccination):")
    st.write(f"R-effective: {float(result_no_vax['Re']):.2f}")
    st.write("Proportion of infections in each group:")
    infections = [1.0, result_no_vax["infections_core"], result_no_vax["infections_children"], result_no_vax["infections_adults"]]
    fatalities_after_G_generations =  [result_no_vax["deaths_after_G_generations"], result_no_vax["deaths_after_G_generations_core"], result_no_vax["deaths_after_G_generations_children"], result_no_vax["deaths_after_G_generations_adults"]]

    st.dataframe(
        pd.DataFrame(
            [np.round(infections, 2),
             np.round(fatalities_after_G_generations, 1)],
            columns=["Total", *group_names],
        ).style.hide(axis="index")
    )
    st.write(f"Infection fatality ratio: {result_no_vax['ifr']:.3f}")


if __name__ == "__main__":
    app()
