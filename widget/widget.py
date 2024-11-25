import numpy as np
import streamlit as st

# 2 groups as test
# To-do: be very sure we know what rows vs columns mean
r_mat = np.array([[1.0, 2.0], [2.0, 1.0]])

st.sidebar.markdown("## Parameters")
high_r0 = st.sidebar.slider("High R0", 0.1, 10.0, 1.0)  # min, max, default
low_r0 = st.sidebar.slider("Low R0", 0.1, 10.0, 1.0)  # min, max, default
