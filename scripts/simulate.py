import numpy as np
import ngm as ngm
import polars as pl
import polars.selectors as cs
import griddler
import griddler.griddle

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

    result = ngm.run_ngm(
        M_novax=M_novax, n=N_i, n_vax=n_vax, ve=params["ve"]
    )

    Re = result["Re"]
    ifr = np.dot(result["infection_distribution"], p_severe)
    fatalities_per_prior_infection = ngm.severity(eigenvalue = Re, eigenvector = result["infection_distribution"],
                                               p_severe = p_severe, G = 1)
    fatalities_after_G_generations = ngm.severity(eigenvalue = Re, eigenvector = result["infection_distribution"],
                                                  p_severe = p_severe, G = params["G"])

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


if __name__ == "__main__":
    parameter_sets = griddler.griddle.read("scripts/config.yaml")

    strategy_names = {"even": "even", "0": "core first", "1": "group 1 first", "2": "group 2 first", "0_1": "core and group 1 first"}
    results_all = griddler.run_squash(simulate_scenario, parameter_sets).with_columns(
        pl.col("vax_strategy").replace_strict(strategy_names)
    )

    scen = results_all.select("scenario").row(0)[0]

    cols_to_select = [
        "n_vax_total",
        "vax_strategy",
        "Re",
        "ifr"]

    results = (
        results_all.with_columns(cs.float().round(3))
        .select(cs.by_name(cols_to_select) | cs.starts_with("deaths_per_prior", "infections_"))
        .sort(["n_vax_total", "vax_strategy"])
    )

    with pl.Config(tbl_rows=-1):
        print(results)

    # save results
    results.write_csv(f"scripts/results_{scen}.csv")
