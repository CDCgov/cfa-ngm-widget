import griddler
import griddler.griddle
import polars as pl
import polars.selectors as cs

from ngm.widget import simulate_scenario

if __name__ == "__main__":
    parameter_sets = griddler.griddle.read("scripts/config.yaml")

    strategy_names = {
        "even": "even",
        "0": "core first",
        "1": "group 1 first",
        "2": "group 2 first",
        "0_1": "core and group 1 first",
    }
    results_all = griddler.run_squash(simulate_scenario, parameter_sets).with_columns(
        pl.col("vax_strategy").replace_strict(strategy_names)
    )

    scen = results_all.select("scenario").row(0)[0]

    cols_to_select = ["n_vax_total", "vax_strategy", "Re", "ifr"]

    results = (
        results_all.with_columns(cs.float().round(3))
        .select(
            cs.by_name(cols_to_select)
            | cs.starts_with("deaths_per_prior", "infections_")
        )
        .sort(["n_vax_total", "vax_strategy"])
    )

    with pl.Config(tbl_rows=-1):
        print(results)

    # save results
    results.write_csv(f"scripts/results_{scen}.csv")
