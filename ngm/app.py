import altair as alt
import numpy as np
import polars as pl
import streamlit as st
import streamlit.delta_generator

import ngm


def simulate_scenario(params, distributions_as_percents=False):
    assert sum(params["pop_props"]) == 1.0

    mult = 1.0
    if distributions_as_percents:
        mult = 100.0

    # population sizes
    N_i = params["n_total"] * np.array(params["pop_props"])

    M_novax = np.array(params["M_novax"])
    p_severe = np.array(params["p_severe"])

    if "n_vax" in params:
        n_vax = params["n_vax"]
    else:
        n_vax = ngm.distribute_vaccines(
            params["n_vax_total"], N_i, strategy=params["vax_strategy"]
        )

    result = ngm.run_ngm(M_novax=M_novax, n=N_i, n_vax=n_vax, ve=params["ve"])

    Re = result["Re"]
    ifr = np.dot(result["infection_distribution"], p_severe)
    fatalities_per_prior_infection = ngm.severity(
        eigenvalue=Re,
        eigenvector=result["infection_distribution"],
        p_severe=p_severe,
        G=1,
    )
    fatalities_after_G_generations = ngm.severity(
        eigenvalue=Re,
        eigenvector=result["infection_distribution"],
        p_severe=p_severe,
        G=params["G"],
    )

    infection_distribution_dict = {
        f"infections_{group}": result["infection_distribution"][i] * mult
        for i, group in enumerate(params["group_names"])
    }

    deaths_per_prior_infection_dict = {
        f"deaths_per_prior_infection_{group}": fatalities_per_prior_infection[i]
        for i, group in enumerate(params["group_names"])
    }

    deaths_after_G_generations_dict = {
        f"deaths_after_G_generations_{group}": fatalities_after_G_generations[i]
        for i, group in enumerate(params["group_names"])
    }

    # Combine all dictionaries into results_dict
    results_dict = {
        "Re": Re,
        "ifr": ifr,
        "deaths_per_prior_infection": fatalities_per_prior_infection.sum(),
        "deaths_after_G_generations": fatalities_after_G_generations.sum(),
        **infection_distribution_dict,
        **deaths_per_prior_infection_dict,
        **deaths_after_G_generations_dict,
    }
    return pl.DataFrame(results_dict)


def extract_vector(
    prefix: str,
    df: pl.DataFrame,
    index_name: str,
    sigdigs,
    groups=["core", "children", "adults"],
):
    assert df.shape[0] == 1
    cols = [prefix + grp for grp in groups]
    vec = (
        df.with_columns(
            total=pl.sum_horizontal(cols),
        )
        .select(pl.col(col).round_sig_figs(sigdigs) for col in ["total", *cols])
        .with_columns(
            summary=pl.lit(index_name),
        )
        .select(["summary", "total", *cols])
        .rename(lambda cname: cname.replace(prefix, "") if prefix in cname else cname)
    )
    return vec


def summarize_scenario(
    c: streamlit.delta_generator.DeltaGenerator,
    params: dict,
    sigdigs,
    groups,
    display=[
        "infections_",
        "deaths_per_prior_infection_",
        "deaths_after_G_generations_",
    ],
    display_names=[
        "Percent of infections",
        "Severe infections per prior infection",
        "Severe infections after G generations",
    ],
):
    p_vax = params["n_vax"] / (params["n_total"] * params["pop_props"])

    # Run the simulation with vaccination
    result = simulate_scenario(params, distributions_as_percents=True)

    c.header(f"*{params['scenario_title']}*")

    prop_vax_help = f"Based on allocated doses, what percent of each group is vaccinated? In the counter factual scenario, we assume no vaccines are administered. Vaccination of 100% does not guarantee complete immunity if VE is less than 1. VE is {params['ve']}"
    c.subheader("% of each group vaccinated:", help=prop_vax_help)
    c.dataframe(
        (
            pl.DataFrame(
                {grp: [prob * 100] for grp, prob in zip(params["group_names"], p_vax)}
            ).select(
                pl.col(col).round_sig_figs(sigdigs) for col in params["group_names"]
            )
        )
    )
    c.subheader("Summaries of Infections:")

    res = pl.concat(
        [
            extract_vector(disp, result, disp_name, sigdigs, groups=groups)
            for disp, disp_name in zip(display, display_names)
        ]
    )
    c.dataframe(res)

    ngm_help = "This is the Next Generation Matrix accounting for the specified administration of vaccines in this scenario."
    c.subheader("Next Generation Matrix given vaccine scenario:")
    m_vax = ngm.vaccinate_M(params["M_novax"], p_vax, params["ve"])
    ngm_df = (
        pl.DataFrame(
            {f"from {grp}": m_vax[:, i] for i, grp in enumerate(params["group_names"])}
        )
        .with_columns(pl.Series("", [f"to {grp}" for grp in params["group_names"]]))
        .select(["", *[f"from {grp}" for grp in params["group_names"]]])
    )
    c.dataframe(ngm_df)
    c.write(ngm_help)

    re_help = "The effective reproductive number accounting for the specified administration of vaccines in this scenario."
    c.subheader(f"R-effective: {result['Re'].round_sig_figs(sigdigs)[0]}", help=re_help)

    ifr_help = 'The probability that a random infection will result in the severe outcome of interest, e.g. death, accounting for the specified administration of vaccines in this scenario. Here "random" means drawing uniformly across all infections, so the probability that one draws an infection in any class is given by the distribution specified in the summary table above.'
    c.subheader(
        f"Severe infection ratio: {result['ifr'].round_sig_figs(sigdigs)[0]}",
        help=ifr_help,
    )

    c.subheader(
        "Cumulative infections after G generations of infection",
        help="This plot shows how many infections (in total across groups) there will be, both severe and otherwise, cumulatively, up to and including G generations of infection. The first generation is the generation produced by the index case, so G = 1 includes the index infection (generation 0) and one generation of spread",
    )

    percent_infections = np.array(res.select(list(params["group_names"]))[0] / 100)

    growth_df = (
        pl.from_numpy(
            ngm.exp_growth_model_severity(
                result["Re"],
                percent_infections,
                params["p_severe"],
                params["G"],
            ),
            schema=["Generation", "All Infections", "Severe Infections"],
        )
        .with_columns(
            (pl.col("All Infections") - pl.col("Severe Infections")).alias(
                "Non-Severe Infections"
            )
        )
        .drop("All Infections")
        .unpivot(index="Generation", variable_name="Infection Type", value_name="Count")
    )

    # Bar plot
    chart = (
        alt.Chart(growth_df)
        .mark_bar()
        .encode(x="Generation:O", y="Count:Q", color="Infection Type:N")
        .properties(title="")
    )

    c.altair_chart(chart, use_container_width=True)


def app():
    st.info(
        "This interactive application is a prototype designed for software testing and educational purposes."
    )

    st.title("Vaccine Allocation Widget")
    st.write(
        "Uses a Next Generation Matrix (NGM) approach to approximate the dynamics of disease spread around the disease-free equilibrium."
    )

    params_default = pl.DataFrame(
        {
            "Group name": ["Core", "Children", "Adults"],
            "Pop. size": np.array([0.05, 0.45, 0.5]) * int(1e7),
            "No. vaccines": [2.5e5, 2.5e5, 5e5],
            "Prob. severe": [0.02, 0.06, 0.02],
        }
    )

    with st.sidebar:
        st.header(
            "Model Inputs",
            help="If you can't see the full matrices without scrolling, drag the sidebar to make it wider.",
        )

        st.subheader(
            "Population Information",
            help="Edit entries in the following matrix to define the:\n - Group names\n - Numbers of people in each group\n - Number of vaccines allocated to each group\n - Probability that an infection will produce the severe outcome of interest (e.g. death) in each group",
        )
        params = st.sidebar.data_editor(params_default)

        VE = st.slider(
            "Vaccine Efficacy",
            0.0,
            1.0,
            value=0.74,
            step=0.01,
            help="Protection due to vaccination is assumed to be all or nothing with this efficacy.",
        )

        m_def_np = np.array([[3.0, 0.0, 0.2], [0.10, 1.0, 0.5], [0.25, 1.0, 1.5]])
        M_default = (
            pl.DataFrame(
                {
                    f"from {grp}": m_def_np[:, i]
                    for i, grp in enumerate(params["Group name"])
                }
            )
            .with_columns(pl.Series("", [f"to {grp}" for grp in params["Group name"]]))
            .select(["", *[f"from {grp}" for grp in params["Group name"]]])
        )

        st.subheader(
            "Next Generation Matrix",
            help="For a single new infection of category `from`, specify how many infections it will generate of category `to` by editing the corresponding entry in the matrix.",
        )
        M_df = st.data_editor(M_default, disabled=["to"], hide_index=True)
        M_novax = M_df.drop("").to_numpy()

        with st.expander("Advanced Options"):
            G = st.slider(
                "Generations",
                1,
                10,
                value=10,
                step=1,
                help="Outcomes after this many generations are summarized.",
            )
            sigdigs = st.slider(
                "Displayed significant figures",
                1,
                4,
                value=3,
                step=1,
                help="Values are reported only to this many significant figures.",
            )

    # # make and run scenarios ------------------------------------------------------------
    group_names = params["Group name"]
    N = params["Pop. size"].to_numpy()
    V = params["No. vaccines"].to_numpy()
    p_severe = params["Prob. severe"].to_numpy()

    scenario = {
        "scenario_title": "Scenario: Vaccination",
        "group_names": group_names,
        "n_total": N.sum(),
        "pop_props": N / N.sum(),
        "M_novax": M_novax,
        "p_severe": p_severe,
        "n_vax": V,
        "ve": VE,
        "G": G,
    }

    counterfactual = scenario.copy()
    counterfactual["scenario_title"] = "Scenario: Counterfactual (no vaccination)"
    counterfactual["n_vax"] = 0 * V
    counterfactual["p_vax"] = 0.0 * V

    scenarios = [
        scenario,
        counterfactual,
    ]

    # text outside scenarios
    summary_help = (
        "Each scenario below gives a summary of infections, including:\n"
        "- Percent of infections: The percent of all infections which are in the given group.\n"
        "- Severe infections per prior infection: If there is one infection, how many severe infections in each group will there be in the next generation of infections?\n"
        "- Severe infections after G generations: Starting with one index infection, how many severe infections will there have been, cumulatively, in each group after G generations of infection? Note that the index infection is marginalized over the the distribution on infections from the table above.\n"
    )
    st.write(summary_help)

    # present results ------------------------------------------------------------
    c = st.container()
    for s in scenarios:
        summarize_scenario(c=c, params=s, sigdigs=sigdigs, groups=params["Group name"])


if __name__ == "__main__":
    app()
