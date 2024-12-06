import numpy as np
import ngm as ngm
import polars as pl
import polars.selectors as cs
import griddler
import griddler.griddle

strategy_names = {"even": "even", "0": "core first", "1": "children first"}

parameter_sets = griddler.griddle.read("scratch/config.yaml")


def simulate_scenario(params, return_polars = True):
    assert sum(params["pop_props"]) == 1.0

    # population sizes
    N_i = params["n_total"] * np.array(params["pop_props"])

    R_novax = np.array(params["R_novax"])
    p_severe = np.array(params["p_severe"])

    if "n_vax" in params:
        n_vax = params["n_vax"]
    else:
        n_vax = ngm.distribute_vaccines(
            params["n_vax_total"], N_i, strategy=params["vax_strategy"]
        )

    result = ngm.run_ngm(
        M_novax=R_novax, n=N_i, n_vax=n_vax, ve=params["ve"]
    )

    Re = result["Re"]
    ifr = np.dot(result["infection_distribution"], p_severe)
    fatalities_per_prior_infection = ngm.severity(eigenvalue = Re, eigenvector = result["infection_distribution"],
                                               p_severe = p_severe, G = 1)
    fatalities_after_G_generations = ngm.severity(eigenvalue = Re, eigenvector = result["infection_distribution"],
                                                  p_severe = p_severe, G = 10)



    results_dict = {
            "Re": Re,
            "infections_core": result["infection_distribution"][0],
            "infections_children": result["infection_distribution"][1],
            "infections_adults": result["infection_distribution"][2],
            "ifr": ifr,
            "deaths_per_prior_infection": fatalities_per_prior_infection.sum(),
            "deaths_per_prior_infection_core": fatalities_per_prior_infection[0],
            "deaths_per_prior_infection_children": fatalities_per_prior_infection[1],
            "deaths_per_prior_infection_adults": fatalities_per_prior_infection[2],
            "deaths_after_G_generations": fatalities_after_G_generations.sum(),
            "deaths_after_G_generations_core": fatalities_after_G_generations[0],
            "deaths_after_G_generations_children": fatalities_after_G_generations[1],
            "deaths_after_G_generations_adults": fatalities_after_G_generations[2],

        }

    if return_polars:
        return pl.DataFrame(results_dict)
    else:
        return(results_dict)


if __name__ == "__main__":

    results_all = griddler.run_squash(simulate_scenario, parameter_sets).with_columns(
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
