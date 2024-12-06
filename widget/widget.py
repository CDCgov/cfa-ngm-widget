# To-do: be very sure we know what rows vs columns mean
import numpy as np
import pandas as pd
import streamlit as st

import ngm


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

    K_default = pd.DataFrame(
        {"to": params_default["Group name"]}
        | {f"from {x}": [0.0] * n_groups_default for x in params_default["Group name"]}
    )
    K_default.set_index("to")
    K_default.iloc[:, 1:] = np.array(
        [[3.0, 0.0, 0.2], [0.1, 1.0, 0.5], [0.25, 1.0, 1.5]]
    )

    st.title("NGM Calculator")

    st.header("Model inputs")

    params = st.data_editor(params_default)

    K_df = st.data_editor(K_default, disabled=["to"], hide_index=True)
    K = K_df.iloc[:, 1:].to_numpy()

    VE = st.slider("Vaccine Efficacy", 0.0, 1.0, value=0.7, step=0.01)

    # Perform the NGM calculation
    result = ngm.simulate(
        R_novax=K,
        n=params["Pop. size"].to_numpy(),
        n_vax=params["No. vaccines"].to_numpy(),
        p_severe=params["Prob. severe"].to_numpy(),
        ve=VE,
    )

    st.header("Results with vaccination")

    st.write("Next Generation Matrix:")
    st.dataframe(result["R"])

    # Display results
    st.write(f"R-effective: {result['Re']:.2f}")
    st.write("Proportion of infections in each group:")
    st.dataframe(
        pd.DataFrame(
            [np.round(result["infections"], 2)],
            columns=params["Group name"],
        ).style.hide(axis="index")
    )
    st.write(
        f"Infection fatality ratio: {(params['Prob. severe'].to_numpy() * result['infections']).sum():.3f}"
    )

    st.subheader("Counterfactual (no vaccination):")
    st.write("Next Generation Matrix:")

    novax_eigen = ngm.dominant_eigen(K)
    st.write(f"R0: {novax_eigen.value:.2f}")
    st.write("Proportion of infections in each group:")
    st.dataframe(
        pd.DataFrame(
            [np.round(novax_eigen.vector, 2)],
            columns=params["Group name"],
        ).style.hide(axis="index")
    )
    st.write(
        f"Infection fatality ratio: {(params['Prob. severe'].to_numpy() * novax_eigen.vector).sum():.3f}"
    )


if __name__ == "__main__":
    app()
