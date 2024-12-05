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
    ifr = result["severe_infection_ratio"]

    return pl.DataFrame({"Re": Re, "ifr": ifr, "ifr_times_Re": ifr * Re})


results_all = griddler.run_squash(simulate, parameter_sets).with_columns(
    pl.col("vax_strategy").replace_strict(strategy_names)
)

results = (
    results_all.select(["n_vax_total", "vax_strategy", "Re", "ifr", "ifr_times_Re"])
    .with_columns(cs.float().round(3))
    .sort(["n_vax_total", "vax_strategy"])
)

print(results_all)

with pl.Config(tbl_rows=-1):
    print(results)

# save results
results.write_csv("scratch/results.csv")
