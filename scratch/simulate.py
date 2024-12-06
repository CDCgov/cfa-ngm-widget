import numpy as np
import ngm as ngm
import polars as pl
import polars.selectors as cs
import griddler
import griddler.griddle

strategy_names = {"even": "even", "0": "core first", "1": "children first"}

parameter_sets = griddler.griddle.read("scratch/config.yaml")


def simulate(params):
    assert sum(params["pop_props"]) == 1.0

    # population sizes
    N_i = params["n_total"] * np.array(params["pop_props"])

    R_novax = np.array(params["R_novax"])
    p_severe = np.array(params["p_severe"])

    n_vax = ngm.distribute_vaccines(
        params["n_vax_total"], N_i, strategy=params["vax_strategy"]
    )

    result = ngm.simulate(
        R_novax=R_novax, n=N_i, n_vax=n_vax, p_severe=p_severe, ve=params["ve"]
    )

    Re = result["Re"]
    ifr = np.dot(result["infection_distribution"], p_severe)
    fatalities_per_prior_infection = Re * result["infection_distribution"] * p_severe

    # there are two ways to get the Re*IFR
    assert np.isclose(
        Re * result["severe_infection_ratio"], fatalities_per_prior_infection.sum()
    )

    return pl.DataFrame(
        {
            "Re": Re,
            "ifr": ifr,
            "deaths_per_prior_infection": fatalities_per_prior_infection.sum(),
            "deaths_per_prior_infection_core": fatalities_per_prior_infection[0],
            "deaths_per_prior_infection_children": fatalities_per_prior_infection[1],
            "deaths_per_prior_infection_adults": fatalities_per_prior_infection[2],
        }
    )


results_all = griddler.run_squash(simulate, parameter_sets).with_columns(
    pl.col("vax_strategy").replace_strict(strategy_names)
)

results = (
    results_all.with_columns(cs.float().round(3))
    .select(
        [
            "n_vax_total",
            "vax_strategy",
            "Re",
            "ifr",
            "deaths_per_prior_infection",
            "deaths_per_prior_infection_core",
            "deaths_per_prior_infection_children",
            "deaths_per_prior_infection_adults",
        ]
    )
    .sort(["n_vax_total", "vax_strategy"])
)

with pl.Config(tbl_rows=-1):
    print(results)

# save results
results.write_csv("scratch/results.csv")
