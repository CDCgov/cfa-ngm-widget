import numpy as np
import pandas as pd
import streamlit as st
import altair as alt

def exp_growth_model(G, inf_distribution, p_severe, R_e):
    total_severe_infections_list = []
    total_infections_list = []

    total_severe_infections = 0
    total_infections = 0
    for g in range(1, G + 1):
        total_infections += pow(R_e, g)
        total_severe_infections += sum(pow(R_e, g) * inf_distribution * p_severe)

        total_infections_list.append(total_infections)
        total_severe_infections_list.append(round(total_severe_infections))

    non_severe_infections_list = np.array(total_infections_list) - np.array(total_severe_infections_list)
    return total_infections_list, total_severe_infections_list, non_severe_infections_list

def app():
    st.title("Infection Simulation")

    # User input
    G = st.slider("Number of Generations (G)", min_value=1, max_value=20, value=10)
    inf_distribution = st.text_input("Infection Distribution (comma-separated)", "0.25,0.25,0.5")
    inf_distribution = np.array([float(x) for x in inf_distribution.split(",")])
    p_severe = st.text_input("Probability of Severe Outcome (comma-separated)", "0.02,0.06,0.02")
    p_severe = np.array([float(x) for x in p_severe.split(",")])
    R_e = st.number_input("Reproduction Number (R_e)", min_value=0.0, value=3.0)

    # Make data objects
    total_infections_list, total_severe_infections_list, non_severe_infections_list = exp_growth_model(G, inf_distribution, p_severe, R_e)

    # wrangle data
    data = pd.DataFrame({
        'Generation': np.arange(1, G + 1),
        'Severe Infections': total_severe_infections_list,
        'Non-Severe Infections': non_severe_infections_list
    })

    data_melted = data.melt(id_vars='Generation', var_name='Infection Type', value_name='Count')

    # Bar plot
    chart = alt.Chart(data_melted).mark_bar().encode(
        x='Generation:O',
        y='Count:Q',
        color='Infection Type:N'
    ).properties(
        title='Total Non-Severe and Severe Infections by Generation'
    )

    st.altair_chart(chart, use_container_width=True)

    # Display result table
    st.subheader("Infections by Generation")
    st.table(data)


if __name__ == "__main__":
    app()
