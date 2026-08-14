# -*- coding: utf-8 -*-
"""Small helper script to generate sample_input.xlsx, a demo workbook
with 2 sheets (2 rainfall series) using the exact column names the
app expects: 'Year' and 'Rainfall (mm)'. Run this once with:
    python make_sample_input.py
"""
import numpy as np
import pandas as pd

np.random.seed(42)

def make_series(start_year, n, shape, scale, base):
    years = list(range(start_year, start_year + n))
    values = np.round(np.random.gamma(shape=shape, scale=scale, size=n) + base, 1)
    return pd.DataFrame({"Year": years, "Rainfall (mm)": values})

station_a = make_series(1995, 28, shape=6, scale=28, base=60)
station_b = make_series(2000, 22, shape=5, scale=35, base=45)

with pd.ExcelWriter("sample_input.xlsx", engine="xlsxwriter") as writer:
    station_a.to_excel(writer, sheet_name="Station A", index=False)
    station_b.to_excel(writer, sheet_name="Station B", index=False)

print("sample_input.xlsx created.")
