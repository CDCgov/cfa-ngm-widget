# To-do: be very sure we know what rows vs columns mean
import numpy as np
import polars as pl
import streamlit as st

import ngm
from scratch.simulate import simulate_scenario

def extract_vector(prefix: str, df: pl.DataFrame, index_name: str, sigdigs, groups=["core", "children", "adults"]):
    assert df.shape[0] == 1
    cols = [prefix + grp for grp in groups]
    vec = (
        df
        .with_columns(
            total=pl.sum_horizontal(cols),
        )
        .select(
            pl.col(col).round_sig_figs(sigdigs) for col in ["total", *cols]
        )
        .with_columns(
            summary=pl.lit(index_name),
        )
        .select(["summary", "total", *cols])
        .rename(lambda cname: cname.replace(prefix, "") if prefix in cname else cname)
    )
    return vec

def summarize_scenario(params, sigdigs, display=["infections_", "deaths_per_prior_infection_", "deaths_after_G_generations_"], display_names=["Percent of infections", "Deaths per prior infection", "Deaths after G generations"]):
    # Run the simulation with vaccination
    result = simulate_scenario(params, distributions_as_percents=True)

    st.subheader(params["scenario_title"])

    st.dataframe(
        pl.concat([
            extract_vector(disp, result, disp_name, sigdigs) for disp,disp_name in zip(display, display_names)
        ])
    )

    st.write(f"R-effective: {result['Re'].round_sig_figs(sigdigs)[0]}")

    st.write(f"Infection fatality ratio: {result['ifr'].round_sig_figs(sigdigs)[0]}")


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

    M_novax = np.zeros((n_groups, n_groups,))
    for ft in from_to:
        row = ft[1][0]
        col = ft[0][0]
        M_novax[row, col] = st.sidebar.number_input(
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

        st.sidebar.subheader("Misc")
        sigdigs = st.sidebar.slider("Displayed significant figures", 1, 10, value=3, step=1)

    scenario = {
        "scenario_title": "Results with vaccination",
        "group_names": group_names,
        "n_total": N.sum(),
        "pop_props": N/N.sum(),
        "M_novax": M_novax,
        "p_severe": p_severe,
        "n_vax": V,
        "ve": VE,
        "G": G,
    }

    counterfactual = scenario.copy()
    counterfactual["scenario_title"] = "Counterfactual (no vaccination)"
    counterfactual["n_vax"] = 0 * V

    scenarios = [
        scenario,
        counterfactual,
    ]

    for s in scenarios:
        summarize_scenario(s, sigdigs)


if __name__ == "__main__":
    app()
