# To-do: be very sure we know what rows vs columns mean
import numpy as np
import pandas as pd
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
    p_vax = params["n_vax"] / (params["n_total"] * params["pop_props"])

    # Run the simulation with vaccination
    result = simulate_scenario(params, distributions_as_percents=True)

    st.subheader(params["scenario_title"])

    st.write("Percent of each group vaccinated:")
    st.dataframe(
        (
            pl.DataFrame({
                grp : [prob * 100]
                for grp,prob in zip(params["group_names"], p_vax)
            })
            .select(
                pl.col(col).round_sig_figs(sigdigs) for col in params["group_names"]
            )

        )
    )

    st.write("Summaries of Infections:")
    st.dataframe(
        pl.concat([
            extract_vector(disp, result, disp_name, sigdigs) for disp,disp_name in zip(display, display_names)
        ])
    )

    st.write("Next Generation Matrix:")
    m_vax = ngm.vaccinate_M(params["M_novax"], p_vax, params["ve"])
    ngm_df = (
        pl.DataFrame({
            f"from {grp}": m_vax[:,i]
            for i,grp in enumerate(params["group_names"])
        })
        .with_columns(pl.Series("", [f"to {grp}" for grp in params["group_names"]]))
        .select(["", *[f"from {grp}" for grp in params["group_names"]]])
    )
    st.dataframe(
        ngm_df
    )

    st.write(f"R-effective: {result['Re'].round_sig_figs(sigdigs)[0]}")

    st.write(f"Infection fatality ratio: {result['ifr'].round_sig_figs(sigdigs)[0]}")


def app():

    # set defaults ------------------------------------------------------------
    params_default = pd.DataFrame(
        {
            "Group name": ["Core", "Children", "Non-core adult"],
            "Pop. size": np.array([0.05, 0.45, 0.5]) * int(1e7),
            "No. vaccines": [2.5e5, 2.5e5, 5e5],
            "Prob. severe": [0.02, 0.06, 0.02],
        }
    )

    n_groups_default = params_default.shape[0]

    M_default = pd.DataFrame(
        {"to": params_default["Group name"]}
        | {f"from {x}": [0.0] * n_groups_default for x in params_default["Group name"]}
    )
    M_default.set_index("to")
    M_default.iloc[:, 1:] = np.array(
        [[3.0, 0.0, 0.2], [0.1, 1.0, 0.5], [0.25, 1.0, 1.5]]
    )

    # get parameters ------------------------------------------------------------
    st.title("NGM Calculator")

    st.sidebar.header("Model Inputs", help="If you can't see the full matrices without scrolling, drag the sidebar to make it wider.")

    st.sidebar.subheader("Population Information", help="Define the:\n - Group names\n - Numbers of people in each group\n - Number of vaccines allocated to each group\n - Probability that an infection will produce the severe outcome of interest (e.g. death) in each group")
    params = st.sidebar.data_editor(params_default)

    st.sidebar.subheader("Next Generation Matrix", help="For a single new infection of category `from`, specify how many infections of category `to` it will make.")
    M_df = st.sidebar.data_editor(M_default, disabled=["to"], hide_index=True)
    M_novax = M_df.iloc[:, 1:].to_numpy()

    VE = st.sidebar.slider("Vaccine Efficacy", 0.0, 1.0, value=0.7, step=0.01, help="Vaccines are assumed to be all or nothing, and the probability that someone is successfully immunized is this value.")

    G = st.sidebar.slider("Generations", 1, 10, value=10, step=1, help="Outcomes after this many generations are summarized.")

    sigdigs = st.sidebar.slider("Displayed significant figures", 1, 10, value=3, step=1, help="Values are reported only to this many significant figures.")

    # make and run scenarios ------------------------------------------------------------
    group_names = params["Group name"]
    N = params["Pop. size"].to_numpy()
    V = params["No. vaccines"].to_numpy()
    p_severe = params["Prob. severe"].to_numpy()

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
    counterfactual["p_vax"] = 0.0 * V

    scenarios = [
        scenario,
        counterfactual,
    ]

    # present results ------------------------------------------------------------
    st.write("The...")
    for s in scenarios:
        summarize_scenario(s, sigdigs)


if __name__ == "__main__":
    app()
