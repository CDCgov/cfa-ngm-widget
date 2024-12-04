import numpy as np
import ngm as ngm
import polars as pl
import polars.selectors as cs
import griddler
import griddler.griddle

strategy_names = {"even": "even", "0": "core first", "1": "children first"}

parameter_sets = griddler.griddle.read("config.yaml")

R = np.array([[]])


def make_beta(lo, hi):
    return np.array([[hi, lo, lo], [lo, lo, lo], [lo, lo, lo]])


def simulate(params):
    assert sum(params["pop_props"]) == 1.0
    # population sizes
    N_i = params["n_total"] * np.array(params["pop_props"])

    beta = make_beta(lo=params["beta_low"], hi=params["beta_high"])
    p_severe = np.array(params["p_severe"])

    n_vax = ngm.distribute_vaccines(
        params["n_vax_total"], N_i, strategy=params["vax_strategy"]
    )

    result = ngm.simulate(
        n=N_i, n_vax=n_vax, beta=beta, p_severe=p_severe, ve=params["ve"]
    )

    Re = result["Re"]
    ifr = np.dot(result["infections"], result["severe_infections"])

    return pl.DataFrame({"Re": Re, "ifr": ifr, "ifr_times_Re": ifr * Re})


results = griddler.run_squash(simulate, parameter_sets).with_columns(
    pl.col("vax_strategy").replace_strict(strategy_names)
)


print(results.to_pandas().to_markdown())

# save results
results.select(
    ["n_vax_total", "vax_strategy", "Re", "ifr", "Re_times_ifr"]
).with_columns(cs.float().round(3)).write_csv("results.csv")
