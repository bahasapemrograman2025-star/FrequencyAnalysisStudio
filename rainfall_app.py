# -*- coding: utf-8 -*-
"""
=======================================================================
 RAINFALL / ANNUAL DISCHARGE FREQUENCY ANALYSIS
 6 Distributions: Normal, LogNormal, 3-Param LogNormal, Gumbel,
 Pearson III, Log Pearson III, GEV (PWM / Hosking)

 v4.0 - "Aurora Studio" Edition
   - Dual data source: Import from Excel  OR  Manual entry in-app
   - Unified "Series Manager" with live KPI dashboard cards
   - Redesigned, colour-coded workflow sidebar with numbered steps
     and a live workspace snapshot (series loaded / ready)
   - Interactive Studio (Sidebar Tab 4) with Station Filter,
     an at-a-glance statistics strip, colour-coded Data Quality
     cards with an overview banner, and a log-scale Return Period
     chart with per-distribution toggles and Show/Hide All controls.

 v3.0 - "Data Studio" Edition
   - Dual data source: Import from Excel  OR  Manual entry in-app
   - Unified "Series Manager"
   - Interactive Studio (Sidebar Tab 4) with Station Filter,
     Detailed Data Quality Cards, and Interactive Chart.
=======================================================================
"""
import sys
import os
import re
import io
import csv
import traceback
import queue
import threading
import subprocess
import numpy as np
import pandas as pd
from scipy import stats, special

import tkinter as tk
from tkinter import ttk, filedialog, messagebox, scrolledtext
import tkinter.font as tkfont

# --- PUSTAKA TAMBAHAN UNTUK UI/UX MODERN ---
import customtkinter as ctk
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.ticker
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import mplcursors

# When built with PyInstaller in --windowed mode (no console window)
class _NullWriter:
    def write(self, *args, **kwargs): pass
    def flush(self): pass

if sys.stdout is None:
    sys.stdout = _NullWriter()
if sys.stderr is None:
    sys.stderr = _NullWriter()

# =======================================================================
#                            BASIC STATISTICS
# =======================================================================

def basic_moments(x):
    x = np.asarray(x, dtype=float)
    n = x.size
    mean = x.mean()
    std = x.std(ddof=1)
    cv = std / mean if mean != 0 else np.nan
    cs = (n / ((n - 1) * (n - 2))) * np.sum((x - mean) ** 3) / std ** 3 if n > 2 else np.nan
    ck = (n ** 2 / ((n - 1) * (n - 2) * (n - 3))) * np.sum((x - mean) ** 4) / std ** 4 if n > 3 else np.nan
    return dict(n=n, mean=mean, std=std, cv=cv, cs=cs, ck=ck)

def lmoments(x):
    x = np.sort(np.asarray(x, dtype=float))
    n = x.size
    i = np.arange(1, n + 1, dtype=float)
    b0 = x.mean()
    b1 = np.sum((i - 1) / (n - 1) * x) / n if n > 1 else np.nan
    b2 = (np.sum((i - 1) * (i - 2) / ((n - 1) * (n - 2)) * x) / n if n > 2 else np.nan)
    b3 = (np.sum((i - 1) * (i - 2) * (i - 3) / ((n - 1) * (n - 2) * (n - 3)) * x) / n if n > 3 else np.nan)

    L1 = b0
    L2 = 2 * b1 - b0 if n > 1 else np.nan
    L3 = 6 * b2 - 6 * b1 + b0 if n > 2 else np.nan
    L4 = 20 * b3 - 30 * b2 + 12 * b1 - b0 if n > 3 else np.nan

    t3 = L3 / L2 if (n > 2 and L2 not in (0, np.nan)) else np.nan
    t4 = L4 / L2 if (n > 3 and L2 not in (0, np.nan)) else np.nan
    return dict(L1=L1, L2=L2, L3=L3, L4=L4, t3=t3, t4=t4)

# =======================================================================
#                DATA QUALITY TESTS
# =======================================================================

def data_quality_tests(x, years=None):
    x = np.asarray(x, dtype=float)
    n = len(x)
    no = np.arange(1, n + 1)
    if years is None:
        years = no
    years = np.asarray(years)

    # 1. OUTLIER TEST
    log_x = np.log10(x)
    mean_log = np.mean(log_x)
    std_log = np.std(log_x, ddof=1)
    t_crit_out = stats.t.ppf(1 - 0.1 / (2 * n), n - 2)
    kn = ((n - 1) / np.sqrt(n)) * np.sqrt(t_crit_out ** 2 / (n - 2 + t_crit_out ** 2))
    yh = mean_log + kn * std_log
    yl = mean_log - kn * std_log
    xh = 10 ** yh
    xl = 10 ** yl
    outliers_high = x[x > xh]
    outliers_low = x[x < xl]
    outlier_status = "No Outliers" if (len(outliers_high) == 0 and len(outliers_low) == 0) else "Outliers Present"

    df_outlier_detail = pd.DataFrame({"No": no, "Year": years, "Rainfall (mm)": x, "Log Rainfall": log_x})
    outlier_summary = {
        "Mean": mean_log, "Std Dev": std_log, "Kn (Table)": kn,
        "Upper Limit Yh": yh, "Lower Limit Yl": yl,
        "Upper Limit Xh": xh, "Lower Limit Xl": xl,
        "Conclusion": outlier_status,
    }

    # 2. TREND TEST
    t = np.arange(1, n + 1, dtype=float)
    rank_t = stats.rankdata(t)
    rank_x = stats.rankdata(x)
    d = rank_t - rank_x
    d2 = d ** 2
    rho, p_trend = stats.spearmanr(t, x)
    t_trend = rho * np.sqrt((n - 2) / (1 - rho ** 2)) if rho != 1.0 else np.inf
    t_trend_crit = stats.t.ppf(0.975, n - 2)
    trend_status = "Trend Detected (Fail)" if abs(t_trend) > t_trend_crit else "No Trend (Pass)"

    df_trend_detail = pd.DataFrame({
        "No": no, "Year": years, "Rainfall (mm)": x,
        "Rank Time (t)": rank_t, "Rank Rainfall (X)": rank_x, "d = Rt-Rx": d, "d^2": d2,
    })
    trend_summary = {
        "Sum d^2": np.sum(d2), "Rho (Spearman)": rho,
        "t Calculated": t_trend, "t Critical": t_trend_crit, "Conclusion": trend_status,
    }

    # 3. VARIANCE HOMOGENEITY
    n1 = n // 2
    idx1 = np.arange(0, n1)
    idx2 = np.arange(n1, n)
    x1, x2 = x[idx1], x[idx2]
    mean1, mean2 = np.mean(x1), np.mean(x2)
    v1, v2 = np.var(x1, ddof=1), np.var(x2, ddof=1)
    dev1 = x1 - mean1
    dev2 = x2 - mean2
    f_calc = max(v1, v2) / min(v1, v2) if min(v1, v2) > 0 else np.inf
    f_crit = stats.f.ppf(0.975, len(x1) - 1, len(x2) - 1)
    var_status = "Variance Not Homogeneous (Fail)" if f_calc > f_crit else "Variance Homogeneous (Pass)"

    df_group1 = pd.DataFrame({"No": no[idx1], "Year": years[idx1], "Rainfall (mm)": x1, "(Xi-X1bar)": dev1, "(Xi-X1bar)^2": dev1 ** 2})
    df_group2 = pd.DataFrame({"No": no[idx2], "Year": years[idx2], "Rainfall (mm)": x2, "(Xi-X2bar)": dev2, "(Xi-X2bar)^2": dev2 ** 2})
    var_summary = {
        "n1": len(x1), "n2": len(x2), "Mean Group 1": mean1, "Mean Group 2": mean2,
        "Variance Group 1 (S1^2)": v1, "Variance Group 2 (S2^2)": v2,
        "F Calculated": f_calc, "F Critical": f_crit, "Conclusion": var_status,
    }

    # 4. MEAN HOMOGENEITY
    sp2 = ((len(x1) - 1) * v1 + (len(x2) - 1) * v2) / (len(x1) + len(x2) - 2)
    sp = np.sqrt(sp2)
    t_mean_calc = abs(mean1 - mean2) / np.sqrt(sp2 * (1 / len(x1) + 1 / len(x2)))
    t_mean_crit = stats.t.ppf(0.975, len(x1) + len(x2) - 2)
    mean_status = "Mean Not Homogeneous (Fail)" if t_mean_calc > t_mean_crit else "Mean Homogeneous (Pass)"

    mean_summary = {
        "n1": len(x1), "n2": len(x2), "Mean Group 1": mean1, "Mean Group 2": mean2,
        "Sp (Pooled)": sp, "Sp^2": sp2,
        "T Calculated": t_mean_calc, "T Critical": t_mean_crit, "Conclusion": mean_status,
    }

    # 5. INDEPENDENCE TEST
    mean_x = np.mean(x)
    dev = x - mean_x
    num_terms = dev[:-1] * dev[1:]
    den = np.sum(dev ** 2)
    num = np.sum(num_terms)
    r1 = num / den if den > 0 else 0
    ll = (-1 - 1.96 * np.sqrt(n - 2)) / (n - 1)
    ul = (-1 + 1.96 * np.sqrt(n - 2)) / (n - 1)
    indep_status = "Not Independent (Fail)" if not (ll <= r1 <= ul) else "Independent (Pass)"

    df_indep_detail = pd.DataFrame({
        "No": no[:-1], "Year (i)": years[:-1], "Xi": x[:-1], "Xi+1": x[1:],
        "(Xi-Xbar)": dev[:-1], "(Xi+1-Xbar)": dev[1:], "(Xi-Xbar)(Xi+1-Xbar)": num_terms, "(Xi-Xbar)^2": dev[:-1] ** 2,
    })
    indep_summary = {
        "Mean (Xbar)": mean_x, "Sum Numerator": num, "Sum Denominator": den,
        "r1 (Correlation)": r1, "Lower Limit (Ll)": ll, "Upper Limit (Ul)": ul,
        "Conclusion": indep_status,
    }

    return {
        "outlier": {"kn": kn, "xh": xh, "xl": xl, "high": outliers_high, "low": outliers_low,
                    "detail": df_outlier_detail, "summary": outlier_summary},
        "trend": {"rho": rho, "t_calc": t_trend, "t_crit": t_trend_crit, "status": trend_status,
                  "detail": df_trend_detail, "summary": trend_summary},
        "var": {"f_calc": f_calc, "f_crit": f_crit, "status": var_status,
                "detail1": df_group1, "detail2": df_group2, "summary": var_summary},
        "mean": {"t_calc": t_mean_calc, "t_crit": t_mean_crit, "status": mean_status,
                 "detail1": df_group1, "detail2": df_group2, "summary": mean_summary},
        "indep": {"r1": r1, "ll": ll, "ul": ul, "status": indep_status,
                  "detail": df_indep_detail, "summary": indep_summary},
    }

def gev_pwm_ks_tables(x, gev_fit):
    x_sort_desc = np.sort(x)[::-1]
    n = len(x)
    b1_list, b2_list, b3_list, emp_list, teo_list, sim_list = [], [], [], [], [], []

    for i in range(n):
        val = x_sort_desc[i]
        rank_asc = n - i
        term_b1 = val * (rank_asc - 1) / (n - 1) if n > 1 else 0
        term_b2 = val * (rank_asc - 1) * (rank_asc - 2) / ((n - 1) * (n - 2)) if n > 2 else 0
        term_b3 = val * (rank_asc - 1) * (rank_asc - 2) * (rank_asc - 3) / ((n - 1) * (n - 2) * (n - 3)) if n > 3 else 0
        b1_list.append(term_b1); b2_list.append(term_b2); b3_list.append(term_b3)

        emp = ((i + 1) / (n + 1)) * 100
        teo = (1.0 - gev_fit.cdf(val)) * 100
        sim = abs(emp - teo)
        emp_list.append(emp); teo_list.append(teo); sim_list.append(sim)

    df_pwm = pd.DataFrame({"X (Sorted)": x_sort_desc, "b1": b1_list, "b2": b2_list, "b3": b3_list,
                           "Empirical (%)": emp_list, "Theoretical (%)": teo_list, "Deviation (%)": sim_list})
    pct = np.arange(99, 0, -1)
    px = np.arange(1, 100) / 100.0
    fx = [float(gev_fit.ppf(1 - p)) for p in px]
    df_ks_gev = pd.DataFrame({"%": pct, "P(x)": px, "F(x)": fx})
    return df_pwm, df_ks_gev


# =======================================================================
#                      DISTRIBUTIONS (6 METHODS)
# =======================================================================
EULER_GAMMA = 0.5772156649015329

class FittedDist:
    def __init__(self, name, nparams, cdf_func, ppf_func, pdf_func=None, params=None, ok=True, note=""):
        self.name = name; self.nparams = nparams; self._cdf = cdf_func; self._ppf = ppf_func
        self._pdf = pdf_func; self.params = params or {}; self.ok = ok; self.note = note

    def cdf(self, x):
        try: return np.clip(self._cdf(np.asarray(x, dtype=float)), 1e-12, 1 - 1e-12)
        except Exception: return np.full(np.shape(x), np.nan)

    def ppf(self, p):
        try: return self._ppf(np.asarray(p, dtype=float))
        except Exception: return np.full(np.shape(p), np.nan)

    def pdf(self, x):
        if self._pdf is None: return None
        try: return self._pdf(np.asarray(x, dtype=float))
        except Exception: return None

def _failed(name, note="data does not meet requirements / fit failed"):
    return FittedDist(name, nparams=np.nan, cdf_func=None, ppf_func=None, ok=False, note=note)

def _from_scipy(name, dist, params, nparams):
    frozen = dist(**params) if isinstance(params, dict) else dist(*params)
    return FittedDist(name, nparams, frozen.cdf, frozen.ppf, frozen.pdf, params=params)

def fit_normal_mom(x):
    return _from_scipy("Normal", stats.norm, dict(loc=np.mean(x), scale=np.std(x, ddof=1)), 2)

def fit_lognormal_mom(x):
    x = np.asarray(x, dtype=float)
    if np.any(x <= 0): return _failed("LogNormal")
    y = np.log(x)
    return _from_scipy("LogNormal", stats.lognorm, dict(s=y.std(ddof=1), loc=0, scale=np.exp(y.mean())), 2)

def fit_lognormal3_mom(x, name="3-Parameter Log Normal"):
    x = np.asarray(x, dtype=float)
    n = x.size; mean, std = x.mean(), x.std(ddof=1)
    cs = (n / ((n - 1) * (n - 2))) * np.sum((x - mean) ** 3) / std ** 3 if n > 2 else np.nan
    if not np.isfinite(cs) or cs == 0: return _failed(name)

    cs_abs = abs(cs); half = cs_abs / 2.0; disc = np.sqrt(half ** 2 + 1.0)
    cv = np.cbrt(half + disc) + np.cbrt(half - disc)
    if not np.isfinite(cv) or cv <= 0: return _failed(name)

    sigma2 = np.log(cv ** 2 + 1.0); sigma = np.sqrt(sigma2)
    skala = std / cv
    if skala <= 0 or not np.isfinite(skala): return _failed(name)
    mu = np.log(skala) - sigma2 / 2.0

    xo = mean - skala if cs >= 0 else mean + skala
    sign = 1.0 if cs >= 0 else -1.0
    base = stats.lognorm(s=sigma, loc=0, scale=np.exp(mu))

    def cdf(v): return base.cdf(np.asarray(v, dtype=float) - xo) if sign > 0 else base.sf(xo - np.asarray(v, dtype=float))
    def ppf(p): return xo + base.ppf(np.asarray(p, dtype=float)) if sign > 0 else xo - base.ppf(1.0 - np.asarray(p, dtype=float))
    def pdf(v): return base.pdf(np.asarray(v, dtype=float) - xo) if sign > 0 else base.pdf(xo - np.asarray(v, dtype=float))
    return FittedDist(name, 3, cdf, ppf, pdf, params=dict(xo=xo, mu=mu, sigma=sigma, sign=sign))

def fit_gumbel_max_mom(x):
    m, s = np.mean(x), np.std(x, ddof=1)
    alpha = s * np.sqrt(6) / np.pi
    xi = m - EULER_GAMMA * alpha
    return _from_scipy("Gumbel", stats.gumbel_r, dict(loc=xi, scale=alpha), 2)

def fit_pearson3_mom(x, name="Pearson III"):
    x = np.asarray(x, dtype=float)
    n = x.size; mean, std = x.mean(), x.std(ddof=1)
    cs = (n / ((n - 1) * (n - 2))) * np.sum((x - mean) ** 3) / std ** 3 if n > 2 else np.nan
    if not np.isfinite(cs) or cs == 0: return _failed(name)
    alpha = 4.0 / cs ** 2; beta = 0.5 * std * cs; xi = mean - alpha * beta
    if beta > 0:
        frozen = stats.gamma(a=alpha, loc=xi, scale=beta)
        return FittedDist(name, 3, frozen.cdf, frozen.ppf, frozen.pdf, params=dict(alpha=alpha, beta=beta, xi=xi))
    else:
        g = stats.gamma(a=alpha, scale=abs(beta))
        return FittedDist(name, 3, lambda v: 1.0 - g.cdf(xi - v), lambda p: xi - g.ppf(1.0 - p), lambda v: g.pdf(xi - v), params=dict(alpha=alpha, beta=beta, xi=xi))

def fit_logpearson3_mom(x):
    x = np.asarray(x, dtype=float)
    if np.any(x <= 0): return _failed("Log Pearson III")
    fy = fit_pearson3_mom(np.log(x), name="Log Pearson III")
    if not fy.ok: return fy
    return FittedDist("Log Pearson III", 3, lambda v: fy.cdf(np.log(v)), lambda p: np.exp(fy.ppf(p)), params=fy.params)

def fit_gev_pwm(x, lm, name="GEV"):
    t3 = lm["t3"]
    if not np.isfinite(t3): return _failed(name)
    c_h = 2.0 / (3.0 + t3) - np.log(2) / np.log(3)
    k = 7.8590 * c_h + 2.9554 * c_h ** 2

    if abs(k) < 1e-8:
        alpha = lm["L2"] / np.log(2)
        xi = lm["L1"] - EULER_GAMMA * alpha
        return _from_scipy(name, stats.gumbel_r, dict(loc=xi, scale=alpha), 2)

    alpha = k * lm["L2"] / ((1 - 2 ** (-k)) * special.gamma(1 + k))
    xi = lm["L1"] - (alpha / k) * (1 - special.gamma(1 + k))
    if alpha <= 0: return _failed(name)

    def ppf(p):
        with np.errstate(invalid="ignore", divide="ignore"): return xi + (alpha / k) * (1.0 - (-np.log(np.asarray(p, dtype=float))) ** k)
    def cdf(v):
        with np.errstate(invalid="ignore", divide="ignore"):
            arg = 1.0 - k * (np.asarray(v, dtype=float) - xi) / alpha
            return np.exp(-(np.where(arg > 0, arg, np.nan) ** (1.0 / k)))
    def pdf(v):
        with np.errstate(invalid="ignore", divide="ignore"):
            arg = 1.0 - k * (np.asarray(v, dtype=float) - xi) / alpha
            arg = np.where(arg > 0, arg, np.nan)
            return (1.0 / alpha) * (arg ** (1.0 / k - 1.0)) * np.exp(-(arg ** (1.0 / k)))
    return FittedDist(name, 3, cdf, ppf, pdf, params=dict(alpha=alpha, k=k, xi=xi))

def build_all_fits(x, lm):
    return {
        "Normal": fit_normal_mom(x), "LogNormal": fit_lognormal_mom(x),
        "3-Parameter Log Normal": fit_lognormal3_mom(x), "Gumbel": fit_gumbel_max_mom(x),
        "Pearson III": fit_pearson3_mom(x), "Log Pearson III": fit_logpearson3_mom(x),
        "GEV": fit_gev_pwm(x, lm, name="GEV")
    }

# =======================================================================
#             GOODNESS OF FIT TESTS & QUANTILES
# =======================================================================
ALPHAS = [0.01, 0.05, 0.10]
RATIO_LOW, RATIO_HIGH = 1.7, 3.2

def n_classes(n): return max(int(round(1 + 3.322 * np.log10(n))), 4)

def chi_square_test(x, fitted, k=None):
    x = np.asarray(x, dtype=float); n = x.size; k = k or n_classes(n)
    if not fitted.ok: return dict(statistic=np.nan, dof=np.nan, pvalue=np.nan, k=k)
    edges = np.concatenate([[-np.inf], np.sort(fitted.ppf(np.linspace(0, 1, k + 1)[1:-1])), [np.inf]])
    observed, _ = np.histogram(x, bins=edges)
    expected = n / k; dof = k - 1 - fitted.nparams
    if dof <= 0 or not np.isfinite(fitted.nparams): return dict(statistic=np.nan, dof=np.nan, pvalue=np.nan, k=k)
    stat = np.sum((observed - expected) ** 2 / expected)
    return dict(statistic=stat, dof=dof, pvalue=stats.chi2.sf(stat, dof), k=k)

def ks_test(x, fitted):
    if not fitted.ok: return dict(dmax=np.nan, pvalue=np.nan)
    try:
        res = stats.kstest(np.asarray(x, dtype=float), fitted.cdf)
        return dict(dmax=res.statistic, pvalue=res.pvalue)
    except Exception: return dict(dmax=np.nan, pvalue=np.nan)

def decision(pvalue, alpha):
    return "" if not np.isfinite(pvalue) else ("ACCEPT" if pvalue >= alpha else "REJECT")

def quantiles_for_periods(fitted, periods, x=None):
    if not fitted.ok: return {T: np.nan for T in periods}
    lo, hi = -np.inf, np.inf
    if x is not None:
        mean, std = np.mean(x), np.std(x, ddof=1)
        lo, hi = mean - 20 * std, mean + 200 * std
    out = {}
    for T in periods:
        try:
            val = float(fitted.ppf(np.array(1.0 - 1.0 / T)))
            out[T] = val if np.isfinite(val) and lo <= val <= hi else np.nan
        except Exception: out[T] = np.nan
    return out

def empirical_plotting_position(x):
    x_sorted = np.sort(np.asarray(x, dtype=float))[::-1]
    return (np.arange(1, x_sorted.size + 1) / (x_sorted.size + 1.0)) ** -1, x_sorted

def build_return_period_curve_data(x, fits_dict, periods, n_curve_points=100):
    x = np.asarray(x, dtype=float)
    T_emp, x_emp = empirical_plotting_position(x)
    Tline = np.linspace(1.01, max(periods), n_curve_points)
    p_line = 1.0 - 1.0 / Tline
    mean, std = np.mean(x), np.std(x, ddof=1)
    ylo, yhi = mean - 20 * std, mean + 200 * std
    curves = {}
    for name, f in fits_dict.items():
        if f.ok:
            try:
                y = f.ppf(p_line)
                curves[name] = np.where((y >= ylo) & (y <= yhi), y, np.nan)
            except Exception: curves[name] = np.full_like(Tline, np.nan)
        else: curves[name] = np.full_like(Tline, np.nan)
    return T_emp, x_emp, Tline, curves

# =======================================================================
#                          MAIN PROGRAM / CONSTANTS
# =======================================================================
COL_TAHUN = "Year"
COL_DATA = "Rainfall (mm)"
RETURN_PERIODS = [1.01, 1.1, 2, 5, 10, 20, 25, 50, 100, 200, 500, 1000]
ROUND_DECIMALS = 4
CREDENTIAL_TEXT = "Created and developed by Ali Nursamsi Dahlan in 2026. Please report any errors immediately."

# =======================================================================
#                              DATA I/O
# =======================================================================
def read_input_dataframe(path, sheet, col_tahun=COL_TAHUN, col_data=COL_DATA):
    df = pd.read_excel(path, sheet_name=sheet)
    df.columns = [str(c).strip() for c in df.columns]
    if col_tahun not in df.columns or col_data not in df.columns:
        raise ValueError(f"Column '{col_tahun}' and/or '{col_data}' not found in sheet.")
    df = df[[col_tahun, col_data]].dropna()
    df[col_data] = pd.to_numeric(df[col_data], errors="coerce")
    df = df.dropna()
    df[col_tahun] = df[col_tahun].apply(lambda v: int(v) if float(v).is_integer() else v)
    return df.reset_index(drop=True)

def read_excel_all_sheets(path, col_tahun=COL_TAHUN, col_data=COL_DATA):
    xls = pd.ExcelFile(path)
    series, errors = {}, {}
    for sheet in xls.sheet_names:
        try:
            df = read_input_dataframe(path, sheet, col_tahun, col_data)
            if len(df) == 0:
                errors[sheet] = "No valid numeric rows found."
                continue
            series[sheet] = df
        except Exception as e: errors[sheet] = str(e)
    return series, errors

def empty_series_df():
    return pd.DataFrame({COL_TAHUN: pd.Series(dtype="Int64"), COL_DATA: pd.Series(dtype="float")})

def validate_series_df(df):
    n = len(df)
    if n == 0: return False, "No data rows yet"
    if n < 5: return False, f"Only {n} rows (minimum 5 recommended)"
    return True, "Ready"

def parse_bulk_text(raw_text):
    rows, bad_lines = [], []
    for raw_line in raw_text.splitlines():
        line = raw_line.strip()
        if not line: continue
        parts = [p for p in re.split(r"[,\t;]+|\s+", line) if p != ""]
        if len(parts) < 2:
            bad_lines.append(raw_line); continue
        try:
            year_val = float(parts[0].replace(",", "."))
            data_val = float(parts[1].replace(",", "."))
            rows.append((int(year_val) if float(year_val).is_integer() else year_val, data_val))
        except Exception: bad_lines.append(raw_line)
    return rows, bad_lines

# =======================================================================
#                        MAIN ANALYSIS PIPELINE
# =======================================================================
def run_analysis(df, col_tahun=COL_TAHUN, col_data=COL_DATA):
    x = df[col_data].to_numpy(dtype=float)
    years = df[col_tahun].to_numpy()
    n = x.size
    if n < 5: raise ValueError("Insufficient data (minimum suggested >= 5 values).")

    stats_dasar = basic_moments(x)
    lm = lmoments(x)
    s = pd.Series(x)
    pandas_stats = dict(n=len(x), min=s.min(), max=s.max(), mean=s.mean(), std=s.std(ddof=1), ck=s.kurt(), cs=s.skew())

    fits = build_all_fits(x, lm)
    k = n_classes(n)
    dq_tests = data_quality_tests(x, years=years)
    df_pwm, df_ks_gev = gev_pwm_ks_tables(x, fits["GEV"])

    chi_rows, ks_rows, quant_rows = [], [], []
    for name, f in fits.items():
        chi = chi_square_test(x, f, k=k)
        chi_rows.append(dict(Distribution=name, **{f"a={int(a*100)}%": decision(chi["pvalue"], a) for a in ALPHAS},
                             **{"Attained alpha": chi["pvalue"], "Chi-Square Statistic": chi["statistic"], "dof": chi["dof"]}))
        ks = ks_test(x, f)
        ks_rows.append(dict(Distribution=name, **{f"a={int(a*100)}%": decision(ks["pvalue"], a) for a in ALPHAS},
                            **{"Attained alpha": ks["pvalue"], "Dmax": ks["dmax"]}))
        q = quantiles_for_periods(f, RETURN_PERIODS, x=x)
        quant_rows.append(dict(Distribution=name, **{f"T={T}": q[T] for T in RETURN_PERIODS}))

    chi_df, ks_df, quant_df = pd.DataFrame(chi_rows), pd.DataFrame(ks_rows), pd.DataFrame(quant_rows)
    ks_df["Rank"] = ks_df["Dmax"].rank(method="min", na_option="bottom").astype("Int64")
    ks_df = ks_df.sort_values("Distribution", key=lambda s: s.map({n: i for i, n in enumerate(fits)}))
    
    num_cols = [c for c in quant_df.columns if c != "Distribution"]
    quant_df[num_cols] = quant_df[num_cols].round(ROUND_DECIMALS)

    r100, r2 = pd.to_numeric(quant_df["T=100"], errors="coerce"), pd.to_numeric(quant_df["T=2"], errors="coerce")
    with np.errstate(divide="ignore", invalid="ignore"): ratio = r100 / r2
    ratio = ratio.where(r2 != 0, np.nan)
    quant_df["R100/R2"] = ratio.round(ROUND_DECIMALS)
    quant_df["R100/R2 OK?"] = np.where(quant_df["R100/R2"].between(RATIO_LOW, RATIO_HIGH), "Yes", "No")
    chi_df["Attained alpha"] = chi_df["Attained alpha"].round(6)
    ks_df["Attained alpha"] = ks_df["Attained alpha"].round(6)
    ks_df["Dmax"] = ks_df["Dmax"].round(5)

    try: ks_crit_05 = stats.kstwo.ppf(0.95, n)
    except AttributeError: ks_crit_05 = stats.ksone.ppf(0.975, n)
    ks_df["Pass KS 5%"] = np.where(ks_df["Dmax"].notna() & (ks_df["Dmax"] < ks_crit_05), "Pass", "Fail")

    return dict(x=x, stats_dasar=stats_dasar, pandas_stats=pandas_stats, lmoments=lm, fits=fits,
                chi_df=chi_df, ks_df=ks_df, quant_df=quant_df, k_classes=k,
                dq_tests=dq_tests, df_pwm=df_pwm, df_ks_gev=df_ks_gev)

def _stamp_all_sheets(wb):
    fmt_credential = wb.add_format({"font_size": 9, "italic": True, "bold": True, "font_color": "#595959"})
    for ws in wb.worksheets(): ws.write(0, 0, CREDENTIAL_TEXT, fmt_credential)

def _write_df_table(ws, df, start_row, start_col, fmt_hdr, fmt_num, fmt_text, n_decimals=4):
    for c, col_name in enumerate(df.columns): ws.write(start_row, start_col + c, col_name, fmt_hdr)
    for r in range(len(df)):
        for c, col_name in enumerate(df.columns):
            val = df.iat[r, c]
            if isinstance(val, (int, np.integer)): ws.write_number(start_row + r + 1, start_col + c, int(val), fmt_text)
            elif isinstance(val, (float, np.floating)): ws.write_number(start_row + r + 1, start_col + c, float(val), fmt_num)
            else: ws.write(start_row + r + 1, start_col + c, str(val), fmt_text)
    return start_row + len(df) + 2

def _write_summary_kv(ws, summary_dict, start_row, start_col, fmt_label, fmt_value, fmt_value_text):
    r = start_row
    for label, val in summary_dict.items():
        ws.write(r, start_col, label, fmt_label)
        if isinstance(val, (int, np.integer)): ws.write_number(r, start_col + 1, int(val), fmt_value)
        elif isinstance(val, (float, np.floating)): ws.write_number(r, start_col + 1, float(val), fmt_value)
        else: ws.write(r, start_col + 1, str(val), fmt_value_text)
        r += 1
    return r + 1

def export_excel(df_input, result, out_path, col_tahun=COL_TAHUN, col_data=COL_DATA):
    with pd.ExcelWriter(out_path, engine="xlsxwriter") as writer:
        wb = writer.book
        fmt_title = wb.add_format({"bold": True, "font_size": 13})
        fmt_title2 = wb.add_format({"bold": True, "font_size": 11, "bg_color": "#F2F2F2"})
        fmt_hdr = wb.add_format({"bold": True, "bg_color": "#D9E1F2", "border": 1, "align": "center", "valign": "vcenter"})
        fmt_accept = wb.add_format({"bg_color": "#C6EFCE", "font_color": "#006100", "align": "center", "border": 1})
        fmt_reject = wb.add_format({"bg_color": "#FFC7CE", "font_color": "#9C0006", "align": "center", "border": 1})
        fmt_num = wb.add_format({"num_format": "0.0000", "border": 1, "align": "center"})
        fmt_num_new = wb.add_format({"num_format": "0.000", "border": 1, "align": "center"})
        fmt_border = wb.add_format({"border": 1})
        fmt_border_left = wb.add_format({"border": 1, "align": "left"})
        fmt_border_center = wb.add_format({"border": 1, "align": "center"})
        fmt_label = wb.add_format({"border": 1, "align": "left", "bold": True})
        fmt_value_text = wb.add_format({"border": 1, "align": "left"})

        # --- Data Input SHEET ---
        df_input.to_excel(writer, sheet_name="Input Data", index=False, startrow=1)
        ws = writer.sheets["Input Data"]
        ws.set_column(0, 1, 16)
        for c, name in enumerate(df_input.columns): ws.write(1, c, name, fmt_hdr)

        # --- Data Quality Test SHEET ---
        ws = wb.add_worksheet("Data Quality Test")
        writer.sheets["Data Quality Test"] = ws
        dq = result["dq_tests"]

        row = 1
        ws.write(row, 0, "RAINFALL DATA QUALITY TEST (DATA CONSISTENCY TEST)", fmt_title); row += 2
        ws.write(row, 0, "1. RAINFALL DATA OUTLIER TEST (Grubbs / WRC, \u03b1 = 10%)", fmt_title2); row += 1
        row = _write_df_table(ws, dq["outlier"]["detail"], row, 0, fmt_hdr, fmt_num_new, fmt_border_center)
        row = _write_summary_kv(ws, dq["outlier"]["summary"], row, 0, fmt_label, fmt_num_new, fmt_value_text)
        ws.write(row - 2, 1, dq["outlier"]["summary"]["Conclusion"], fmt_accept if "No Outliers" in dq["outlier"]["summary"]["Conclusion"] else fmt_reject); row += 1

        ws.write(row, 0, "2. RAINFALL DATA NO-TREND TEST (Spearman's Rho)", fmt_title2); row += 1
        row = _write_df_table(ws, dq["trend"]["detail"], row, 0, fmt_hdr, fmt_num_new, fmt_border_center)
        row = _write_summary_kv(ws, dq["trend"]["summary"], row, 0, fmt_label, fmt_num_new, fmt_value_text)
        ws.write(row - 2, 1, dq["trend"]["summary"]["Conclusion"], fmt_accept if "No Trend" in dq["trend"]["summary"]["Conclusion"] else fmt_reject); row += 1

        ws.write(row, 0, "3. RAINFALL DATA VARIANCE STATIONARITY TEST (F-Test)", fmt_title2); row += 1
        ws.write(row, 0, "Group 1 (First Half)", fmt_border_left); row += 1
        row = _write_df_table(ws, dq["var"]["detail1"], row, 0, fmt_hdr, fmt_num_new, fmt_border_center)
        ws.write(row, 0, "Group 2 (Second Half)", fmt_border_left); row += 1
        row = _write_df_table(ws, dq["var"]["detail2"], row, 0, fmt_hdr, fmt_num_new, fmt_border_center)
        row = _write_summary_kv(ws, dq["var"]["summary"], row, 0, fmt_label, fmt_num_new, fmt_value_text)
        ws.write(row - 2, 1, dq["var"]["summary"]["Conclusion"], fmt_accept if "Homogeneous" in dq["var"]["summary"]["Conclusion"] and "Not" not in dq["var"]["summary"]["Conclusion"] else fmt_reject); row += 1

        ws.write(row, 0, "4. RAINFALL DATA MEAN STATIONARITY TEST (T-Test)", fmt_title2); row += 1
        ws.write(row, 0, "Group 1 (First Half)", fmt_border_left); row += 1
        row = _write_df_table(ws, dq["mean"]["detail1"], row, 0, fmt_hdr, fmt_num_new, fmt_border_center)
        ws.write(row, 0, "Group 2 (Second Half)", fmt_border_left); row += 1
        row = _write_df_table(ws, dq["mean"]["detail2"], row, 0, fmt_hdr, fmt_num_new, fmt_border_center)
        row = _write_summary_kv(ws, dq["mean"]["summary"], row, 0, fmt_label, fmt_num_new, fmt_value_text)
        ws.write(row - 2, 1, dq["mean"]["summary"]["Conclusion"], fmt_accept if "Homogeneous" in dq["mean"]["summary"]["Conclusion"] and "Not" not in dq["mean"]["summary"]["Conclusion"] else fmt_reject); row += 1

        ws.write(row, 0, "5. RAINFALL DATA INDEPENDENCE TEST (Lag-1 Serial Correlation)", fmt_title2); row += 1
        row = _write_df_table(ws, dq["indep"]["detail"], row, 0, fmt_hdr, fmt_num_new, fmt_border_center)
        row = _write_summary_kv(ws, dq["indep"]["summary"], row, 0, fmt_label, fmt_num_new, fmt_value_text)
        ws.write(row - 2, 1, dq["indep"]["summary"]["Conclusion"], fmt_accept if ("Independent" in dq["indep"]["summary"]["Conclusion"] and "Not" not in dq["indep"]["summary"]["Conclusion"]) else fmt_reject); row += 1

        row += 1; ws.write(row, 0, "RAINFALL DATA QUALITY TEST - QUICK RECAP", fmt_title); row += 2
        ws.write(row, 0, "1. RAINFALL DATA OUTLIER TEST (Grubbs / WRC)", fmt_title2); row += 1
        for label, val in [("Upper Threshold (Xh)", dq["outlier"]["summary"]["Upper Limit Xh"]), ("Lower Threshold (Xl)", dq["outlier"]["summary"]["Lower Limit Xl"]), ("Kn Value", dq["outlier"]["summary"]["Kn (Table)"])]:
            ws.write(row, 0, label, fmt_border_left); ws.write_number(row, 1, float(val), fmt_num_new); row += 1
        ws.write(row, 0, "Outlier Status", fmt_label)
        status = "Safe (No Outliers)" if "No Outliers" in dq["outlier"]["summary"]["Conclusion"] else "Outliers Detected"
        ws.write(row, 1, status, fmt_accept if "Safe" in status else fmt_reject); row += 2

        ws.write(row, 0, "2. RAINFALL DATA NO-TREND TEST (Spearman's Rho)", fmt_title2); row += 1
        for label, val in [("Correlation Value (Rho)", dq["trend"]["summary"]["Rho (Spearman)"]), ("Calculated T", dq["trend"]["summary"]["t Calculated"]), ("Critical T", dq["trend"]["summary"]["t Critical"])]:
            ws.write(row, 0, label, fmt_border_left); ws.write_number(row, 1, float(val), fmt_num_new); row += 1
        ws.write(row, 0, "Trend Status", fmt_label)
        ws.write(row, 1, dq["trend"]["summary"]["Conclusion"], fmt_accept if "No Trend" in dq["trend"]["summary"]["Conclusion"] else fmt_reject); row += 2

        ws.write(row, 0, "3. RAINFALL DATA VARIANCE / F-TEST", fmt_title2); row += 1
        for label, val in [("Calculated F", dq["var"]["summary"]["F Calculated"]), ("Critical F", dq["var"]["summary"]["F Critical"])]:
            ws.write(row, 0, label, fmt_border_left); ws.write_number(row, 1, float(val), fmt_num_new); row += 1
        ws.write(row, 0, "Variance Status", fmt_label)
        var_short = "Homogeneous Variance (Pass)" if ("Homogeneous" in dq["var"]["summary"]["Conclusion"] and "Not" not in dq["var"]["summary"]["Conclusion"]) else "Non-Homogeneous Variance (Fail)"
        ws.write(row, 1, var_short, fmt_accept if "Pass" in var_short else fmt_reject); row += 2

        ws.write(row, 0, "4. RAINFALL DATA MEAN / T-TEST", fmt_title2); row += 1
        for label, val in [("Calculated T", dq["mean"]["summary"]["T Calculated"]), ("Critical T", dq["mean"]["summary"]["T Critical"])]:
            ws.write(row, 0, label, fmt_border_left); ws.write_number(row, 1, float(val), fmt_num_new); row += 1
        ws.write(row, 0, "Mean Status", fmt_label)
        mean_short = "Homogeneous Mean (Pass)" if ("Homogeneous" in dq["mean"]["summary"]["Conclusion"] and "Not" not in dq["mean"]["summary"]["Conclusion"]) else "Non-Homogeneous Mean (Fail)"
        ws.write(row, 1, mean_short, fmt_accept if "Pass" in mean_short else fmt_reject); row += 2

        ws.write(row, 0, "5. RAINFALL DATA INDEPENDENCE TEST (Lag-1 Serial Correlation)", fmt_title2); row += 1
        for label, val in [("Correlation (r1)", dq["indep"]["summary"]["r1 (Correlation)"]), ("Lower Limit", dq["indep"]["summary"]["Lower Limit (Ll)"]), ("Upper Limit", dq["indep"]["summary"]["Upper Limit (Ul)"])]:
            ws.write(row, 0, label, fmt_border_left); ws.write_number(row, 1, float(val), fmt_num_new); row += 1
        ws.write(row, 0, "Independence Status", fmt_label)
        ws.write(row, 1, dq["indep"]["summary"]["Conclusion"], fmt_accept if ("Independent" in dq["indep"]["summary"]["Conclusion"] and "Not" not in dq["indep"]["summary"]["Conclusion"]) else fmt_reject); row += 1

        ws.set_column(0, 0, 26); ws.set_column(1, 9, 16)

        # --- Basic Statistics SHEET ---
        sd, lm, ps = result["stats_dasar"], result["lmoments"], result["pandas_stats"]
        ws = wb.add_worksheet("Basic Statistics")
        writer.sheets["Basic Statistics"] = ws

        ws.write(1, 0, "Data Statistics", fmt_title)
        for i, (label, val) in enumerate([("Data Name", "Rainfall (mm/day)"), ("Data Count", ps["n"]), ("Minimum", ps["min"]), ("Maximum", ps["max"]), ("Mean", ps["mean"]), ("Standard Deviation", ps["std"]), ("Kurtosis", ps["ck"]), ("Skewness", ps["cs"])], start=2):
            ws.write(i, 0, label, fmt_border_left)
            if isinstance(val, str): ws.write(i, 1, val, fmt_border_left)
            else: ws.write_number(i, 1, float(val), fmt_num_new)

        ws.write(12, 0, "BASIC STATISTICS (Method of Moments)", fmt_title)
        for i, (label, val) in enumerate([("Data Count (n)", sd["n"]), ("Mean", sd["mean"]), ("Standard Deviation (Std Dev)", sd["std"]), ("Coef. of Variation (Cv)", sd["cv"]), ("Coef. of Skewness (Cs)", sd["cs"]), ("Coef. of Kurtosis (Ck)", sd["ck"])], start=14):
            ws.write(i, 0, label, fmt_border); ws.write_number(i, 1, float(val), fmt_num)

        ws.write(22, 0, "L-MOMENTS", fmt_title)
        for i, (label, val) in enumerate([("L1", lm["L1"]), ("L2", lm["L2"]), ("L3", lm["L3"]), ("L4", lm["L4"]), ("t3 (L-Skewness)", lm["t3"]), ("t4 (L-Kurtosis)", lm["t4"])], start=24):
            ws.write(i, 0, label, fmt_border); ws.write_number(i, 1, float(val), fmt_num)
        ws.set_column(0, 0, 26); ws.set_column(1, 1, 20)

        req_dist = ["Normal", "LogNormal", "3-Parameter Log Normal", "Gumbel", "Pearson III", "Log Pearson III", "GEV"]
        disp_dist = ["Normal", "Log Normal", "3-Parameter Log Normal", "Gumbel", "Pearson III", "Log Pearson III", "GEV"]

        # --- Chi-Square Test SHEET ---
        chi_df = result["chi_df"].drop(columns=["Chi-Square Statistic", "dof"]).rename(columns={"Attained alpha": "Attained a"})
        chi_full = result["chi_df"]
        chi_df.to_excel(writer, sheet_name="Chi-Square Test", index=False, startrow=3)
        ws = writer.sheets["Chi-Square Test"]

        ws.write(1, 0, f"ORIGINAL CHI-SQUARE TEST (number of classes k = {result['k_classes']})", fmt_title)
        for c, name in enumerate(chi_df.columns): ws.write(3, c, name, fmt_hdr)
        for r in range(len(chi_df)):
            for c in [1, 2, 3]:
                val = chi_df.iat[r, c]
                ws.write(r + 4, c, val, fmt_accept if val == "ACCEPT" else (fmt_reject if val == "REJECT" else fmt_border))
            ws.write_number(r + 4, 4, 0 if pd.isna(chi_df.iat[r, 4]) else float(chi_df.iat[r, 4]), fmt_num) if not pd.isna(chi_df.iat[r, 4]) else ws.write_blank(r + 4, 4, None, fmt_border)
            pearson_val = chi_full["Chi-Square Statistic"].iat[r]
            ws.write_blank(r + 4, 5, None, fmt_border) if pd.isna(pearson_val) else ws.write_number(r + 4, 5, float(pearson_val), fmt_num)
            ws.write(r + 4, 0, chi_df.iat[r, 0], fmt_border)
        ws.write(3, 5, "Chi-Square Statistic", fmt_hdr)

        row_offset = len(chi_df) + 7
        for alpha in ALPHAS:
            ws.write(row_offset, 0, f"Goodness of fit test (Chi-Square) at \u03b1 = {int(round(alpha * 100))}%", fmt_title)
            ws.write(row_offset + 1, 0, "Distribution", fmt_hdr)
            for c, dname in enumerate(disp_dist, start=1): ws.write(row_offset + 1, c, dname, fmt_hdr)
            ws.write(row_offset + 2, 0, "Maximum Delta", fmt_border_left)
            ws.write(row_offset + 3, 0, "Critical Delta", fmt_border_left)
            ws.write(row_offset + 4, 0, "Test Result", fmt_border_left)

            for c, r_name in enumerate(req_dist, start=1):
                row2 = chi_full[chi_full["Distribution"] == r_name].iloc[0]
                stat, dof = row2["Chi-Square Statistic"], row2["dof"]
                if pd.isna(stat) or pd.isna(dof) or dof <= 0:
                    ws.write(row_offset + 2, c, "-", fmt_border_center); ws.write(row_offset + 3, c, "-", fmt_border_center); ws.write(row_offset + 4, c, "Fail", fmt_reject)
                else:
                    crit = stats.chi2.ppf(1 - alpha, dof)
                    hasil = "Pass" if stat < crit else "Fail"
                    ws.write_number(row_offset + 2, c, stat, fmt_num_new); ws.write_number(row_offset + 3, c, crit, fmt_num_new); ws.write(row_offset + 4, c, hasil, fmt_accept if hasil == "Pass" else fmt_reject)
            row_offset += 7
        ws.set_column(0, 0, 30); ws.set_column(1, 7, 16)

        # --- Kolmogorov-Smirnov Test SHEET ---
        ks_df = result["ks_df"].rename(columns={"Attained alpha": "Attained a"})
        ks_df.to_excel(writer, sheet_name="Kolmogorov-Smirnov Test", index=False, startrow=3)
        ws = writer.sheets["Kolmogorov-Smirnov Test"]

        ws.write(1, 0, "ORIGINAL KOLMOGOROV-SMIRNOV TEST", fmt_title)
        for c, name in enumerate(ks_df.columns): ws.write(3, c, name, fmt_hdr)
        for r in range(len(ks_df)):
            ws.write(r + 4, 0, ks_df.iat[r, 0], fmt_border)
            for c in [1, 2, 3]:
                val = ks_df.iat[r, c]
                ws.write(r + 4, c, val, fmt_accept if val == "ACCEPT" else (fmt_reject if val == "REJECT" else fmt_border))
            for c in [4, 5]:
                val = ks_df.iat[r, c]
                ws.write_blank(r + 4, c, None, fmt_border) if pd.isna(val) else ws.write_number(r + 4, c, float(val), fmt_num)
            rank_val = ks_df.iat[r, 6]
            ws.write(r + 4, 6, "" if pd.isna(rank_val) else int(rank_val), fmt_border)

        row_offset = len(ks_df) + 7; n_data = ps["n"]
        for alpha in ALPHAS:
            ws.write(row_offset, 0, f"Goodness of fit test (Smirnov-Kolmogorov) at \u03b1 = {int(round(alpha * 100))}%", fmt_title)
            ws.write(row_offset + 1, 0, "Distribution", fmt_hdr)
            for c, dname in enumerate(disp_dist, start=1): ws.write(row_offset + 1, c, dname, fmt_hdr)
            ws.write(row_offset + 2, 0, "Maximum Delta", fmt_border_left)
            ws.write(row_offset + 3, 0, "Critical Delta", fmt_border_left)
            ws.write(row_offset + 4, 0, "Test Result", fmt_border_left)
            try: ks_crit = stats.kstwo.ppf(1 - alpha, n_data)
            except AttributeError: ks_crit = stats.ksone.ppf(1 - alpha / 2, n_data)

            for c, r_name in enumerate(req_dist, start=1):
                dmax = ks_df[ks_df["Distribution"] == r_name].iloc[0]["Dmax"]
                if pd.isna(dmax):
                    ws.write(row_offset + 2, c, "-", fmt_border_center); ws.write(row_offset + 3, c, ks_crit, fmt_num_new); ws.write(row_offset + 4, c, "Fail", fmt_reject)
                else:
                    hasil = "Pass" if dmax < ks_crit else "Fail"
                    ws.write_number(row_offset + 2, c, dmax, fmt_num_new); ws.write_number(row_offset + 3, c, ks_crit, fmt_num_new); ws.write(row_offset + 4, c, hasil, fmt_accept if hasil == "Pass" else fmt_reject)
            row_offset += 7
        ws.set_column(0, 0, 30); ws.set_column(1, 7, 16)

        # --- GEV - PWM & KS SHEET ---
        ws = wb.add_worksheet("GEV - PWM & KS")
        writer.sheets["GEV - PWM & KS"] = ws
        df_pwm, df_ks_gev = result["df_pwm"], result["df_ks_gev"]

        ws.write(1, 0, "Probability in % (GEV Calculation using PWM Method)", fmt_title)
        for c, col_name in enumerate(df_pwm.columns): ws.write(2, c, col_name, fmt_hdr)
        for r in range(len(df_pwm)):
            ws.write(r + 3, 0, df_pwm.iat[r, 0], fmt_border_center)
            for c in range(1, len(df_pwm.columns)): ws.write_number(r + 3, c, df_pwm.iat[r, c], fmt_num)
        ws.set_column(0, 6, 15)
        row_ks_offset = len(df_pwm) + 7
        ws.write(row_ks_offset, 0, "Kolmogorov-Smirnov GEV", fmt_title)
        for c, col_name in enumerate(df_ks_gev.columns): ws.write(row_ks_offset + 2, c, col_name, fmt_hdr)
        for r in range(len(df_ks_gev)):
            ws.write_number(row_ks_offset + r + 3, 0, df_ks_gev.iat[r, 0], fmt_border_center)
            ws.write_number(row_ks_offset + r + 3, 1, df_ks_gev.iat[r, 1], fmt_border_center)
            ws.write_number(row_ks_offset + r + 3, 2, df_ks_gev.iat[r, 2], fmt_num)

        # --- Return Period SHEET ---
        quant_df = result["quant_df"]
        n_periods = len(RETURN_PERIODS)
        quant_df.to_excel(writer, sheet_name="Return Period", index=False, startrow=3)
        ws = writer.sheets["Return Period"]
        ws.write(1, 0, "RETURN PERIOD VALUES (Return Period, T years)", fmt_title)
        ws.write(2, 0, f"Note: R100/R2 is considered reasonable within the range of {RATIO_LOW}\u2013{RATIO_HIGH}", fmt_border_left)
        for c, name in enumerate(quant_df.columns): ws.write(3, c, name, fmt_hdr)
        for r in range(len(quant_df)):
            for c, col_name in enumerate(quant_df.columns):
                val = quant_df.iat[r, c]
                if c == 0: ws.write(r + 4, c, val, fmt_border)
                elif col_name == "R100/R2 OK?": ws.write(r + 4, c, val, fmt_accept if val == "Yes" else fmt_reject)
                elif pd.isna(val): ws.write_blank(r + 4, c, None, fmt_border)
                else: ws.write_number(r + 4, c, float(val), fmt_num)
        ws.set_column(0, 0, 25); ws.set_column(1, len(quant_df.columns) - 1, 12)

        helper_row = len(quant_df) + 6
        ws.write(helper_row, 0, "T (chart helper)", fmt_border)
        for c, T in enumerate(RETURN_PERIODS, start=1): ws.write_number(helper_row, c, float(T))
        chart = wb.add_chart({"type": "scatter", "subtype": "straight_with_markers"})
        for r in range(len(quant_df)):
            chart.add_series({"name": [ws.get_name(), r + 3, 0], "categories": [ws.get_name(), helper_row, 1, helper_row, n_periods], "values": [ws.get_name(), r + 3, 1, r + 3, n_periods], "marker": {"type": "circle", "size": 4}, "line": {"width": 1}})
        chart.set_title({"name": "Return Period Curve - All Methods"}); chart.set_x_axis({"name": "Return Period T (years)"}); chart.set_y_axis({"name": "Value"}); chart.set_size({"width": 900, "height": 500}); chart.set_legend({"font": {"size": 8}})
        ws.insert_chart(3, len(quant_df.columns) + 1, chart)

        # --- Graph Data (hidden helper sheet) ---
        ws_data = wb.add_worksheet("Graph Data")
        writer.sheets["Graph Data"] = ws_data
        T_emp, x_emp, Tline, curves = build_return_period_curve_data(result["x"], result["fits"], RETURN_PERIODS, n_curve_points=100)
        n_emp, n_curve = len(T_emp), len(Tline)
        ws_data.write(0, 0, "T (Empirical)", fmt_hdr); ws_data.write(0, 1, "X (Empirical)", fmt_hdr)
        for i in range(n_emp): ws_data.write_number(i + 1, 0, float(T_emp[i])); ws_data.write_number(i + 1, 1, float(x_emp[i]))
        curve_col = {}; col = 3; names_ordered = list(result["fits"].keys())
        for name in names_ordered:
            ws_data.write(0, col, f"T_{name}"); ws_data.write(0, col + 1, f"X_{name}")
            y = curves[name]
            for i in range(n_curve):
                ws_data.write_number(i + 1, col, float(Tline[i]))
                if np.isfinite(y[i]): ws_data.write_number(i + 1, col + 1, float(y[i]))
                else: ws_data.write_blank(i + 1, col + 1, None)
            curve_col[name] = col; col += 3
        ws_data.hide()

        # --- Graph SHEET (native Excel charts) ---
        ws = wb.add_worksheet("Graph")
        writer.sheets["Graph"] = ws
        ws.write(1, 0, "Frequency Curve - All Methods (native Excel charts, linear T axis)", fmt_title)
        ncols_grid, chart_width, chart_height, row_start, row_step, col_step = 3, 480, 300, 3, 16, 9
        data_sheet_name = ws_data.get_name()
        for idx, name in enumerate(names_ordered):
            r_idx, c_idx = divmod(idx, ncols_grid); col = curve_col[name]
            chart = wb.add_chart({"type": "scatter", "subtype": "straight_with_markers"})
            chart.add_series({"name": "Empirical Data", "categories": [data_sheet_name, 1, 0, n_emp, 0], "values": [data_sheet_name, 1, 1, n_emp, 1], "marker": {"type": "circle", "size": 5, "fill": {"color": "black"}, "border": {"color": "black"}}, "line": {"none": True}})
            chart.add_series({"name": f"{name} Fit", "categories": [data_sheet_name, 1, col, n_curve, col], "values": [data_sheet_name, 1, col + 1, n_curve, col + 1], "marker": {"type": "none"}, "line": {"color": "red", "width": 1.5}})
            chart.set_title({"name": name}); chart.set_x_axis({"name": "Return Period T (years)"}); chart.set_y_axis({"name": "Value"})
            chart.set_size({"width": chart_width, "height": chart_height}); chart.set_legend({"font": {"size": 7}})
            ws.insert_chart(row_start + r_idx * row_step, c_idx * col_step, chart)
        _stamp_all_sheets(wb)

def pilih_distribusi_terbaik_ks(result):
    ks_df_all = result["ks_df"].dropna(subset=["Rank"]).sort_values("Rank")
    quant_df = result["quant_df"]
    if len(ks_df_all) == 0: return None, None, "No", None, "Fail"
    def jarak(ratio): return RATIO_LOW - ratio if ratio < RATIO_LOW else (ratio - RATIO_HIGH if ratio > RATIO_HIGH else 0.0)
    def kumpulkan(ks_subset):
        candidates = []
        for _, row in ks_subset.iterrows():
            match = quant_df[quant_df["Distribution"] == row["Distribution"]]
            if len(match) > 0 and pd.notna(match.iloc[0].get("R100/R2", np.nan)):
                candidates.append((row["Rank"], row["Distribution"], match.iloc[0], match.iloc[0].get("R100/R2")))
        return candidates
    def lulus_ks5(name):
        row = ks_df_all[ks_df_all["Distribution"] == name]
        return "Fail" if len(row) == 0 else row.iloc[0]["Pass KS 5%"]
    
    candidates = kumpulkan(ks_df_all[ks_df_all["Pass KS 5%"] == "Pass"]) or kumpulkan(ks_df_all)
    if not candidates:
        best_name = ks_df_all.iloc[0]["Distribution"]
        match = quant_df[quant_df["Distribution"] == best_name]
        return best_name, match.iloc[0] if len(match) else None, "No", ks_df_all.iloc[0]["Rank"], lulus_ks5(best_name)

    for rank_val, name, qrow, ratio in candidates:
        if RATIO_LOW <= ratio <= RATIO_HIGH: return name, qrow, "Yes", rank_val, lulus_ks5(name)
    candidates_sorted = sorted(candidates, key=lambda c: (jarak(c[3]), c[0]))
    return candidates_sorted[0][1], candidates_sorted[0][2], "Approaching", candidates_sorted[0][0], lulus_ks5(candidates_sorted[0][1])

def export_rekap_excel(rekap_rows, out_path, dq_rekap_rows=None):
    rekap_df = pd.DataFrame(rekap_rows)
    with pd.ExcelWriter(out_path, engine="xlsxwriter") as writer:
        rekap_df.to_excel(writer, sheet_name="Best Return Period Summary", index=False, startrow=3)
        wb, ws = writer.book, writer.sheets["Best Return Period Summary"]
        fmt_title = wb.add_format({"bold": True, "font_size": 13})
        fmt_hdr = wb.add_format({"bold": True, "bg_color": "#D9E1F2", "border": 1, "align": "center", "valign": "vcenter"})
        fmt_border = wb.add_format({"border": 1})
        fmt_num = wb.add_format({"num_format": "0.0000", "border": 1, "align": "center"})
        fmt_num3 = wb.add_format({"num_format": "0.000", "border": 1, "align": "center"})
        fmt_accept = wb.add_format({"bg_color": "#C6EFCE", "font_color": "#006100", "align": "center", "border": 1})
        fmt_reject = wb.add_format({"bg_color": "#FFC7CE", "font_color": "#9C0006", "align": "center", "border": 1})
        fmt_warning = wb.add_format({"bg_color": "#FFEB9C", "font_color": "#9C6500", "align": "center", "border": 1})
        ws.write(1, 0, "BEST DISTRIBUTION SUMMARY (KS TEST RANKING = 1) & RETURN PERIOD VALUES", fmt_title)
        for c, name in enumerate(rekap_df.columns): ws.write(3, c, name, fmt_hdr)
        for r in range(len(rekap_df)):
            for c, col in enumerate(rekap_df.columns):
                val = rekap_df.iat[r, c]
                if col in ("Series", "Selected Distribution", "Selected KS Rank"): ws.write(r + 4, c, val, fmt_border)
                elif col == "Pass KS 5%?": ws.write(r + 4, c, val, fmt_accept if val == "Pass" else fmt_reject)
                elif col == "R100/R2 In Range?": ws.write(r + 4, c, val, fmt_accept if val == "Yes" else (fmt_warning if val == "Approaching" else fmt_reject))
                elif pd.isna(val): ws.write_blank(r + 4, c, None, fmt_border)
                else: ws.write_number(r + 4, c, float(val), fmt_num)
        ws.set_column(0, 1, 28); ws.set_column(2, len(rekap_df.columns) - 1, 12)

        if dq_rekap_rows:
            ws2 = wb.add_worksheet("Data Quality Test Summary")
            writer.sheets["Data Quality Test Summary"] = ws2
            fmt_hdr_out = wb.add_format({"bold": True, "bg_color": "#FCE4D6", "border": 1, "align": "center", "valign": "vcenter", "text_wrap": True})
            fmt_hdr_trn = wb.add_format({"bold": True, "bg_color": "#E2EFDA", "border": 1, "align": "center", "valign": "vcenter", "text_wrap": True})
            fmt_hdr_var = wb.add_format({"bold": True, "bg_color": "#DDEBF7", "border": 1, "align": "center", "valign": "vcenter", "text_wrap": True})
            fmt_hdr_mea = wb.add_format({"bold": True, "bg_color": "#EDD9F5", "border": 1, "align": "center", "valign": "vcenter", "text_wrap": True})
            fmt_hdr_ind = wb.add_format({"bold": True, "bg_color": "#FFF2CC", "border": 1, "align": "center", "valign": "vcenter", "text_wrap": True})
            COL_DEFS = [("No", "No", fmt_hdr, False, None), ("Series / Station", "Series", fmt_hdr, False, None),
                        ("Upper Threshold\nXh", "out_Xh", fmt_hdr_out, False, None), ("Lower Threshold\nXl", "out_Xl", fmt_hdr_out, False, None), ("Kn Value", "out_Kn", fmt_hdr_out, False, None), ("Outlier Status", "out_Status", fmt_hdr_out, True, "No Outliers"),
                        ("Rho\n(Spearman)", "trn_Rho", fmt_hdr_trn, False, None), ("T Calculated\n(Trend)", "trn_T_calc", fmt_hdr_trn, False, None), ("T Critical\n(Trend)", "trn_T_crit", fmt_hdr_trn, False, None), ("Trend Status", "trn_Status", fmt_hdr_trn, True, "No Trend"),
                        ("F Calculated", "var_F_calc", fmt_hdr_var, False, None), ("F Critical", "var_F_crit", fmt_hdr_var, False, None), ("Variance Status", "var_Status", fmt_hdr_var, True, "Homogeneous"),
                        ("T Calculated\n(Mean)", "mea_T_calc", fmt_hdr_mea, False, None), ("T Critical\n(Mean)", "mea_T_crit", fmt_hdr_mea, False, None), ("Mean Status", "mea_Status", fmt_hdr_mea, True, "Homogeneous"),
                        ("r1\n(Lag-1)", "ind_r1", fmt_hdr_ind, False, None), ("Lower Limit\n(Ll)", "ind_Ll", fmt_hdr_ind, False, None), ("Upper Limit\n(Ul)", "ind_Ul", fmt_hdr_ind, False, None), ("Independence Status", "ind_Status", fmt_hdr_ind, True, "Independent")]
            ws2.write(1, 0, "DATA QUALITY TEST SUMMARY - ALL SERIES", fmt_title)
            ws2.write(2, 0, "One row per series/station. Green = Pass, Red = Fail.", fmt_border)
            groups = [(0, 1, ""), (2, 5, "OUTLIER TEST (Grubbs/WRC)"), (6, 9, "NO-TREND TEST (Spearman's Rho)"), (10, 12, "VARIANCE TEST (F-Test)"), (13, 15, "MEAN TEST (T-Test)"), (16, 19, "INDEPENDENCE TEST (Lag-1)")]
            group_fmt = {"OUTLIER TEST (Grubbs/WRC)": wb.add_format({"bold": True, "bg_color": "#FCE4D6", "border": 1, "align": "center", "valign": "vcenter"}), "NO-TREND TEST (Spearman's Rho)": wb.add_format({"bold": True, "bg_color": "#E2EFDA", "border": 1, "align": "center", "valign": "vcenter"}), "VARIANCE TEST (F-Test)": wb.add_format({"bold": True, "bg_color": "#DDEBF7", "border": 1, "align": "center", "valign": "vcenter"}), "MEAN TEST (T-Test)": wb.add_format({"bold": True, "bg_color": "#EDD9F5", "border": 1, "align": "center", "valign": "vcenter"}), "INDEPENDENCE TEST (Lag-1)": wb.add_format({"bold": True, "bg_color": "#FFF2CC", "border": 1, "align": "center", "valign": "vcenter"})}
            for (c_start, c_end, label) in groups:
                if label == "":
                    for c in range(c_start, c_end + 1): ws2.write(4, c, "", wb.add_format({"border": 1}))
                elif c_start == c_end: ws2.write(4, c_start, label, group_fmt.get(label, fmt_hdr))
                else: ws2.merge_range(4, c_start, 4, c_end, label, group_fmt.get(label, fmt_hdr))
            for c, (hdr, key, hfmt, is_status, pass_kw) in enumerate(COL_DEFS): ws2.write(5, c, hdr, hfmt)
            for r, drow in enumerate(dq_rekap_rows):
                for c, (hdr, key, hfmt, is_status, pass_kw) in enumerate(COL_DEFS):
                    val = drow.get(key, "")
                    if key == "No": ws2.write_number(r + 6, c, r + 1, fmt_num3)
                    elif is_status: ws2.write(r + 6, c, val, fmt_accept if isinstance(val, str) and pass_kw in val and "Not" not in val and "Non" not in val else fmt_reject)
                    elif isinstance(val, (int, float, np.integer, np.floating)) and not isinstance(val, bool): ws2.write_number(r + 6, c, float(val), fmt_num3)
                    else: ws2.write(r + 6, c, str(val) if val is not None else "", fmt_border)
            ws2.set_column(0, 0, 5); ws2.set_column(1, 1, 28); ws2.set_column(2, 3, 14); ws2.set_column(4, 4, 10); ws2.set_column(5, 5, 22); ws2.set_column(6, 8, 13); ws2.set_column(9, 9, 22); ws2.set_column(10, 11, 13); ws2.set_column(12, 12, 26); ws2.set_column(13, 14, 13); ws2.set_column(15, 15, 26); ws2.set_column(16, 18, 13); ws2.set_column(19, 19, 24); ws2.set_row(4, 20); ws2.set_row(5, 35)
        _stamp_all_sheets(wb)

def _safe_sheet_filename(name):
    name = re.sub(r'[\\/*?:"<>|]', "_", str(name).strip())
    return name if name else "Series"

def run_batch_from_series(series_dict, out_dir, log, progress):
    os.makedirs(out_dir, exist_ok=True)
    names = list(series_dict.keys()); total = len(names)
    success, failed, rekap_rows, dq_rekap_rows = [], [], [], []
    for idx, name in enumerate(names, start=1):
        try:
            df = series_dict[name]
            if not validate_series_df(df)[0]: raise ValueError("Not ready")
            result = run_analysis(df)
            out_path = os.path.join(out_dir, f"{_safe_sheet_filename(name)}.xlsx")
            export_excel(df, result, out_path)
            best_name, quant_row, in_range, rank_val, lulus_ks5_status = pilih_distribusi_terbaik_ks(result)
            entry = {"Series": name, "Selected Distribution": best_name or "-", "Selected KS Rank": int(rank_val) if pd.notna(rank_val) else "-", "Pass KS 5%?": lulus_ks5_status or "-", "R100/R2": quant_row.get("R100/R2", np.nan) if quant_row is not None else np.nan, "R100/R2 In Range?": in_range or "-"}
            for T in RETURN_PERIODS: entry[f"T={T}"] = quant_row[f"T={T}"] if quant_row is not None else np.nan
            rekap_rows.append(entry)
            dq = result["dq_tests"]
            dq_rekap_rows.append({"Series": name, "out_Xh": dq["outlier"]["summary"]["Upper Limit Xh"], "out_Xl": dq["outlier"]["summary"]["Lower Limit Xl"], "out_Kn": dq["outlier"]["summary"]["Kn (Table)"], "out_Status": "Safe (No Outliers)" if "No Outliers" in dq["outlier"]["summary"]["Conclusion"] else "Outliers Detected", "trn_Rho": dq["trend"]["summary"]["Rho (Spearman)"], "trn_T_calc": dq["trend"]["summary"]["t Calculated"], "trn_T_crit": dq["trend"]["summary"]["t Critical"], "trn_Status": dq["trend"]["summary"]["Conclusion"], "var_F_calc": dq["var"]["summary"]["F Calculated"], "var_F_crit": dq["var"]["summary"]["F Critical"], "var_Status": "Homogeneous Variance (Pass)" if ("Homogeneous" in dq["var"]["summary"]["Conclusion"] and "Not" not in dq["var"]["summary"]["Conclusion"]) else "Non-Homogeneous Variance (Fail)", "mea_T_calc": dq["mean"]["summary"]["T Calculated"], "mea_T_crit": dq["mean"]["summary"]["T Critical"], "mea_Status": "Homogeneous Mean (Pass)" if ("Homogeneous" in dq["mean"]["summary"]["Conclusion"] and "Not" not in dq["mean"]["summary"]["Conclusion"]) else "Non-Homogeneous Mean (Fail)", "ind_r1": dq["indep"]["summary"]["r1 (Correlation)"], "ind_Ll": dq["indep"]["summary"]["Lower Limit (Ll)"], "ind_Ul": dq["indep"]["summary"]["Upper Limit (Ul)"], "ind_Status": dq["indep"]["summary"]["Conclusion"]})
            success.append(name)
        except Exception as e: failed.append((name, str(e)))
        finally: progress(idx, total)
    if rekap_rows: export_rekap_excel(rekap_rows, os.path.join(out_dir, "Best_Return_Period_Summary.xlsx"), dq_rekap_rows=dq_rekap_rows)
    return success, failed

def _open_in_explorer(path):
    try:
        if sys.platform.startswith("win"): os.startfile(path)
        elif sys.platform == "darwin": subprocess.Popen(["open", path])
        else: subprocess.Popen(["xdg-open", path])
    except Exception: pass

# =======================================================================
#                             GUI  -  THEME
# =======================================================================
# -----------------------------------------------------------------------
#  DESIGN SYSTEM  -  "Aurora Studio" palette
#  A richer, more contemporary colour system: deep indigo/navy sidebar,
#  teal + amber + violet accents, soft tinted surfaces for status states,
#  and a dedicated, colour-blind-friendly palette for charts.
# -----------------------------------------------------------------------
COLOR_BG            = "#F3F6FC"   # app background (very light blue-grey)
COLOR_BG_ALT        = "#E9EFFB"   # secondary background tint
COLOR_SIDEBAR       = "#0B1F3F"   # deep navy sidebar (top of gradient)
COLOR_SIDEBAR_2     = "#132C55"   # sidebar gradient bottom tone
COLOR_SIDEBAR_ACTIVE= "#1D4ED8"   # active nav highlight (vivid blue)
COLOR_SIDEBAR_HOVER = "#16345F"   # nav hover tone
COLOR_PRIMARY       = "#1B4F91"   # primary action blue
COLOR_PRIMARY_DARK  = "#123A6B"
COLOR_SECONDARY     = "#7C3AED"   # violet secondary accent
COLOR_SECONDARY_DARK= "#6D28D9"
COLOR_ACCENT        = "#0EA5A5"   # teal accent
COLOR_ACCENT_DARK   = "#0B8484"
COLOR_AMBER         = "#F59E0B"   # amber highlight accent
COLOR_AMBER_DARK    = "#B45309"
COLOR_CARD          = "#FFFFFF"
COLOR_CARD_ALT      = "#F8FAFC"
COLOR_TEXT          = "#0F172A"
COLOR_MUTED         = "#64748B"
COLOR_MUTED_LIGHT   = "#94A3B8"
COLOR_SUCCESS       = "#15803D"
COLOR_SUCCESS_BG    = "#DCFCE7"
COLOR_SUCCESS_BORDER= "#86EFAC"
COLOR_DANGER        = "#B91C1C"
COLOR_DANGER_BG     = "#FEE2E2"
COLOR_DANGER_BORDER = "#FCA5A5"
COLOR_WARN          = "#B45309"
COLOR_WARN_BG       = "#FEF3C7"
COLOR_WARN_BORDER   = "#FCD34D"
COLOR_INFO          = "#0369A1"
COLOR_INFO_BG       = "#EFF8FF"
COLOR_INFO_BORDER   = "#7DD3FC"
COLOR_BORDER        = "#E1E7F2"
COLOR_BORDER_STRONG = "#C7D2E5"
FONT_FAMILY = "Segoe UI" if sys.platform.startswith("win") else "Helvetica"

# Distinct, print-friendly palette shared by every chart/legend in the app
CHART_PALETTE = ["#2563EB", "#F59E0B", "#10B981", "#EF4444", "#8B5CF6",
                  "#EC4899", "#0EA5A5", "#F97316", "#84CC16", "#6366F1"]

# Small helper: blend two hex colours (used for subtle gradient fills)
def _blend_hex(c1, c2, t):
    c1 = c1.lstrip("#"); c2 = c2.lstrip("#")
    r1, g1, b1 = int(c1[0:2], 16), int(c1[2:4], 16), int(c1[4:6], 16)
    r2, g2, b2 = int(c2[0:2], 16), int(c2[2:4], 16), int(c2[4:6], 16)
    r = int(r1 + (r2 - r1) * t); g = int(g1 + (g2 - g1) * t); b = int(b1 + (b2 - b1) * t)
    return f"#{r:02x}{g:02x}{b:02x}"

class Card(tk.Frame):
    """A soft, elevated content card with an optional accent-coloured top
    strip, an icon-badge, a title and a subtitle. Gently highlights its
    border on hover to give the whole app a more tactile, alive feel."""
    def __init__(self, parent, title=None, subtitle=None, icon=None, accent=COLOR_PRIMARY, **kwargs):
        outer_kwargs = {k: v for k, v in kwargs.items() if k not in ("padx", "pady")}
        super().__init__(parent, bg=COLOR_BORDER, bd=0, **outer_kwargs)
        self._accent = accent
        self._strip = tk.Frame(self, bg=accent, height=3)
        self._strip.pack(fill="x", side="top")
        self._inner = tk.Frame(self, bg=COLOR_CARD)
        self._inner.pack(fill="both", expand=True, padx=1, pady=(0, 1))
        self._pad = tk.Frame(self._inner, bg=COLOR_CARD)
        self._pad.pack(fill="both", expand=True, padx=18, pady=16)
        if title:
            head = tk.Frame(self._pad, bg=COLOR_CARD); head.pack(anchor="w", fill="x")
            if icon:
                badge = tk.Label(head, text=icon, bg=accent, fg="white", font=(FONT_FAMILY, 11, "bold"), width=3, height=1)
                badge.pack(side="left", padx=(0, 8))
            tk.Label(head, text=title, bg=COLOR_CARD, fg=COLOR_TEXT, font=(FONT_FAMILY, 12, "bold")).pack(side="left", anchor="w")
        if subtitle: tk.Label(self._pad, text=subtitle, bg=COLOR_CARD, fg=COLOR_MUTED, font=(FONT_FAMILY, 9)).pack(anchor="w", pady=(2, 6))
        self.bind("<Enter>", self._on_enter); self.bind("<Leave>", self._on_leave)
    def _on_enter(self, e): self.configure(bg=self._accent)
    def _on_leave(self, e): self.configure(bg=COLOR_BORDER)
    @property
    def body(self): return self._pad


def make_kpi_card(parent, icon, value, label, accent=COLOR_PRIMARY, tint=None):
    """A compact 'dashboard' stat tile: big number, small caption, coloured
    accent bar and icon badge. Used to surface useful at-a-glance metrics
    (counts, readiness, statistics) throughout the app."""
    tint = tint or COLOR_CARD
    outer = tk.Frame(parent, bg=COLOR_BORDER)
    top = tk.Frame(outer, bg=accent, height=3); top.pack(fill="x")
    body = tk.Frame(outer, bg=tint); body.pack(fill="both", expand=True, padx=1, pady=(0, 1))
    row = tk.Frame(body, bg=tint); row.pack(fill="x", padx=14, pady=(12, 2))
    tk.Label(row, text=icon, bg=accent, fg="white", font=(FONT_FAMILY, 12, "bold"), width=3).pack(side="left")
    val_lbl = tk.Label(row, text=str(value), bg=tint, fg=COLOR_TEXT, font=(FONT_FAMILY, 22, "bold"))
    val_lbl.pack(side="left", padx=(10, 0))
    lbl = tk.Label(body, text=label, bg=tint, fg=COLOR_MUTED, font=(FONT_FAMILY, 9))
    lbl.pack(anchor="w", padx=14, pady=(0, 12))
    outer.value_label = val_lbl
    return outer


def make_gradient_canvas(parent, width, height, color_top, color_bottom, **kwargs):
    """Draws a smooth vertical gradient onto a Canvas - used to give the
    sidebar a modern, layered look instead of a flat block of colour."""
    cv = tk.Canvas(parent, width=width, height=height, highlightthickness=0, bd=0, **kwargs)
    steps = max(int(height), 2)
    for i in range(steps):
        t = i / (steps - 1)
        cv.create_line(0, i, width, i, fill=_blend_hex(color_top, color_bottom, t))
    return cv

def make_pill_button(parent, text, command, bg=COLOR_PRIMARY, fg="white", hover=None, font_size=10, padx=14, pady=7, state="normal"):
    hover = hover or bg
    btn = tk.Label(parent, text=text, bg=bg, fg=fg, font=(FONT_FAMILY, font_size, "bold"), padx=padx, pady=pady, cursor="hand2")
    btn._enabled = (state == "normal"); btn._bg = bg; btn._hover = hover
    btn.bind("<Button-1>", lambda e: command() if btn._enabled else None)
    btn.bind("<Enter>", lambda e: btn.configure(bg=btn._hover) if btn._enabled else None)
    btn.bind("<Leave>", lambda e: btn.configure(bg=btn._bg) if btn._enabled else None)
    def set_state(st):
        btn._enabled = (st == "normal")
        btn.configure(bg=btn._bg if btn._enabled else "#B7C3D6", fg="white" if btn._enabled else "#EEF2F9", cursor="hand2" if btn._enabled else "arrow")
    btn.set_state = set_state
    if state != "normal": set_state(state)
    return btn

class ManualSeriesEditor(tk.Toplevel):
    def __init__(self, master, series_name, df, on_save):
        super().__init__(master)
        self.title(f"Manual Data Entry - {series_name}"); self.geometry("620x620"); self.minsize(520, 480)
        self.configure(bg=COLOR_BG); self.transient(master); self.grab_set()
        self.series_name, self.on_save, self.df = series_name, on_save, df.copy().reset_index(drop=True)
        self._build_ui(); self._refresh_table(); self._refresh_stats()

    def _build_ui(self):
        header = tk.Frame(self, bg=COLOR_PRIMARY); header.pack(fill="x")
        tk.Label(header, text=f"\u270E  {self.series_name}", bg=COLOR_PRIMARY, fg="white", font=(FONT_FAMILY, 13, "bold")).pack(anchor="w", padx=16, pady=10)
        toolbar = tk.Frame(self, bg=COLOR_BG); toolbar.pack(fill="x", padx=14, pady=(10, 4))
        make_pill_button(toolbar, "+ Add Row", self._add_row, bg=COLOR_ACCENT, hover=COLOR_ACCENT_DARK, font_size=9).pack(side="left", padx=(0, 6))
        make_pill_button(toolbar, "Edit Selected", self._edit_row, bg="#475569", hover="#334155", font_size=9).pack(side="left", padx=6)
        make_pill_button(toolbar, "Delete Selected", self._delete_row, bg=COLOR_DANGER, hover="#991B1B", font_size=9).pack(side="left", padx=6)
        make_pill_button(toolbar, "Sort by Year", self._sort_rows, bg="#475569", hover="#334155", font_size=9).pack(side="left", padx=6)
        make_pill_button(toolbar, "\u2398 Paste / Bulk Add", self._open_bulk_paste, bg="#7C3AED", hover="#6D28D9", font_size=9).pack(side="left", padx=6)
        table_card = Card(self, title="Data Table (Year / Rainfall)", subtitle="Double-click a row to edit it.", icon="\U0001F4CB", accent=COLOR_PRIMARY)
        table_card.pack(fill="both", expand=True, padx=14, pady=6)
        self.tree = ttk.Treeview(table_card.body, columns=("no", "year", "value"), show="headings", height=14)
        for c, t, w in [("no", "No", 50), ("year", "Year", 110), ("value", COL_DATA, 140)]:
            self.tree.heading(c, text=t); self.tree.column(c, width=w, anchor="center")
        vsb = ttk.Scrollbar(table_card.body, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=vsb.set); self.tree.pack(side="left", fill="both", expand=True); vsb.pack(side="right", fill="y")
        self.tree.bind("<Double-1>", lambda e: self._edit_row())
        stats_card = Card(self, title="Live Preview Statistics", icon="\U0001F4C8", accent=COLOR_ACCENT); stats_card.pack(fill="x", padx=14, pady=(0, 6))
        self.lbl_stats = tk.Label(stats_card.body, text="-", bg=COLOR_CARD, fg=COLOR_MUTED, font=(FONT_FAMILY, 9), justify="left"); self.lbl_stats.pack(anchor="w")
        footer = tk.Frame(self, bg=COLOR_BG); footer.pack(fill="x", padx=14, pady=(0, 14))
        make_pill_button(footer, "Save & Close", self._save_and_close, bg=COLOR_PRIMARY, hover="#123A6B").pack(side="right")
        make_pill_button(footer, "Cancel", self.destroy, bg="#94A3B8", hover="#64748B").pack(side="right", padx=(0, 8))

    def _refresh_table(self):
        for i in self.tree.get_children(): self.tree.delete(i)
        for i, row in self.df.iterrows(): self.tree.insert("", "end", iid=str(i), values=(i + 1, row[COL_TAHUN], row[COL_DATA]))
    def _refresh_stats(self):
        n = len(self.df)
        if n == 0: self.lbl_stats.configure(text="No data yet - use '+ Add Row' or 'Paste / Bulk Add' to begin."); return
        ok, msg = validate_series_df(self.df)
        txt = f"N = {n}"
        if n >= 2:
            x = self.df[COL_DATA].to_numpy(dtype=float)
            txt += f"   |   Mean = {x.mean():.3f}   |   Std Dev = {x.std(ddof=1):.3f}"
        if n >= 3:
            m = basic_moments(self.df[COL_DATA].to_numpy(dtype=float))
            txt += f"   |   Cv = {m['cv']:.3f}   |   Cs = {m['cs']:.3f}"
        self.lbl_stats.configure(text=txt + f"\nStatus: {msg}")
    def _row_dialog(self, year=None, value=None):
        dlg = tk.Toplevel(self); dlg.title("Row"); dlg.configure(bg=COLOR_BG); dlg.geometry("300x180"); dlg.transient(self); dlg.grab_set()
        tk.Label(dlg, text="Year", bg=COLOR_BG, fg=COLOR_TEXT, font=(FONT_FAMILY, 10)).pack(pady=(16, 2)); e_year = tk.Entry(dlg, font=(FONT_FAMILY, 11), justify="center"); e_year.pack(padx=30, fill="x")
        if year is not None: e_year.insert(0, str(year))
        tk.Label(dlg, text=COL_DATA, bg=COLOR_BG, fg=COLOR_TEXT, font=(FONT_FAMILY, 10)).pack(pady=(12, 2)); e_val = tk.Entry(dlg, font=(FONT_FAMILY, 11), justify="center"); e_val.pack(padx=30, fill="x")
        if value is not None: e_val.insert(0, str(value))
        res = {}
        def on_ok():
            try:
                y = float(e_year.get().strip().replace(",", ".")); v = float(e_val.get().strip().replace(",", "."))
                res["year"] = int(y) if float(y).is_integer() else y; res["value"] = v; dlg.destroy()
            except Exception: messagebox.showerror("Invalid Input", "Please enter valid numbers for Year and Value.", parent=dlg)
        btns = tk.Frame(dlg, bg=COLOR_BG); btns.pack(pady=16)
        make_pill_button(btns, "OK", on_ok, bg=COLOR_PRIMARY).pack(side="left", padx=6); make_pill_button(btns, "Cancel", dlg.destroy, bg="#94A3B8").pack(side="left", padx=6)
        e_year.focus_set(); dlg.wait_window(); return res.get("year"), res.get("value")
    def _add_row(self):
        y, v = self._row_dialog()
        if y is not None: self.df.loc[len(self.df)] = [y, v]; self._refresh_table(); self._refresh_stats()
    def _selected_index(self):
        sel = self.tree.selection()
        if not sel: messagebox.showinfo("No Selection", "Please select a row first.", parent=self)
        return int(sel[0]) if sel else None
    def _edit_row(self):
        idx = self._selected_index()
        if idx is not None:
            y, v = self._row_dialog(self.df.at[idx, COL_TAHUN], self.df.at[idx, COL_DATA])
            if y is not None: self.df.at[idx, COL_TAHUN], self.df.at[idx, COL_DATA] = y, v; self._refresh_table(); self._refresh_stats()
    def _delete_row(self):
        idx = self._selected_index()
        if idx is not None: self.df = self.df.drop(index=idx).reset_index(drop=True); self._refresh_table(); self._refresh_stats()
    def _sort_rows(self):
        try: self.df = self.df.sort_values(COL_TAHUN).reset_index(drop=True); self._refresh_table()
        except Exception: pass
    def _open_bulk_paste(self):
        dlg = tk.Toplevel(self); dlg.title("Paste / Bulk Add Data"); dlg.configure(bg=COLOR_BG); dlg.geometry("460x420"); dlg.transient(self); dlg.grab_set()
        tk.Label(dlg, text="Paste rows below (one per line).", bg=COLOR_BG, fg=COLOR_TEXT, font=(FONT_FAMILY, 10, "bold")).pack(anchor="w", padx=16, pady=(14, 0))
        tk.Label(dlg, text="Accepted formats:  2010,120.5   /   2010\\t120.5   /   2010 120.5", bg=COLOR_BG, fg=COLOR_MUTED, font=(FONT_FAMILY, 9)).pack(anchor="w", padx=16, pady=(0, 8))
        txt = scrolledtext.ScrolledText(dlg, height=14, font=("Consolas", 10)); txt.pack(fill="both", expand=True, padx=16)
        def on_add():
            rows, bad = parse_bulk_text(txt.get("1.0", "end"))
            if not rows: messagebox.showwarning("Nothing to Add", "No valid Year / Value pairs were found.", parent=dlg); return
            for y, v in rows: self.df.loc[len(self.df)] = [y, v]
            self._refresh_table(); self._refresh_stats()
            msg = f"Added {len(rows)} row(s)."
            if bad: msg += f"\n{len(bad)} line(s) could not be parsed and were skipped."
            messagebox.showinfo("Bulk Add Complete", msg, parent=dlg); dlg.destroy()
        btns = tk.Frame(dlg, bg=COLOR_BG); btns.pack(pady=12)
        make_pill_button(btns, "Add These Rows", on_add, bg=COLOR_ACCENT, hover=COLOR_ACCENT_DARK).pack(side="left", padx=6)
        make_pill_button(btns, "Cancel", dlg.destroy, bg="#94A3B8").pack(side="left", padx=6)
    def _save_and_close(self):
        self.on_save(self.series_name, self.df); self.destroy()

# =======================================================================
#                          MAIN APPLICATION
# =======================================================================
APP_TITLE = "Rainfall Frequency Analysis"
APP_VERSION = "4.0"

class RainfallApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title(f"{APP_TITLE}  v{APP_VERSION}")
        self.geometry("1180x760")
        self.minsize(980, 640)
        self.configure(bg=COLOR_BG)
        
        ctk.set_appearance_mode("Light")
        ctk.set_default_color_theme("blue")

        self.series_store = {}
        self.series_source = {}
        self.out_dir = None
        self.worker_thread = None
        self.msg_queue = queue.Queue()
        self.current_interactive_canvas = None # Track canvas to avoid memory leaks

        self._build_style()
        self._build_layout()
        self.protocol("WM_DELETE_WINDOW", self._on_close)
        self.after(100, self._poll_queue)
        self._show_page("data")

    def _build_style(self):
        style = ttk.Style(self)
        try: style.theme_use("clam")
        except Exception: pass
        style.configure("Treeview", rowheight=28, font=(FONT_FAMILY, 9), background="white", fieldbackground="white", bordercolor=COLOR_BORDER, borderwidth=0)
        style.configure("Treeview.Heading", font=(FONT_FAMILY, 9, "bold"), background="#DCE6F5", foreground=COLOR_TEXT, relief="flat", padding=6)
        style.map("Treeview.Heading", background=[("active", "#C7D9F2")])
        style.map("Treeview", background=[("selected", COLOR_ACCENT)], foreground=[("selected", "white")])
        style.configure("TProgressbar", background=COLOR_ACCENT, troughcolor="#DCE3EE", bordercolor=COLOR_BG, lightcolor=COLOR_ACCENT, darkcolor=COLOR_ACCENT)
        style.configure("Aurora.Horizontal.TProgressbar", thickness=14, background=COLOR_ACCENT, troughcolor="#E2E8F0", bordercolor="#E2E8F0", lightcolor=COLOR_ACCENT, darkcolor=COLOR_SECONDARY)
        style.configure("TCombobox", fieldbackground="white", background="white", foreground=COLOR_TEXT, arrowcolor=COLOR_PRIMARY, padding=6)
        style.map("TCombobox", fieldbackground=[("readonly", "white")])
        style.configure("Vertical.TScrollbar", background="#C7D2E5", troughcolor=COLOR_BG, arrowsize=12)

    def _build_layout(self):
        self.sidebar = tk.Frame(self, bg=COLOR_SIDEBAR, width=250)
        self.sidebar.pack(side="left", fill="y"); self.sidebar.pack_propagate(False)

        # Soft vertical gradient strip along the very left edge for depth
        make_gradient_canvas(self.sidebar, 250, 760, COLOR_SIDEBAR, COLOR_SIDEBAR_2).place(x=0, y=0, relheight=1)

        brand = tk.Frame(self.sidebar, bg=COLOR_SIDEBAR); brand.place(x=0, y=0, relwidth=1)
        brand_pad = tk.Frame(brand, bg=COLOR_SIDEBAR); brand_pad.pack(fill="x", pady=(24, 22), padx=20)
        logo_row = tk.Frame(brand_pad, bg=COLOR_SIDEBAR); logo_row.pack(anchor="w", fill="x")
        logo_badge = tk.Label(logo_row, text="\U0001F327", bg=COLOR_ACCENT, fg="white", font=(FONT_FAMILY, 18), width=2, height=1)
        logo_badge.pack(side="left")
        title_col = tk.Frame(logo_row, bg=COLOR_SIDEBAR); title_col.pack(side="left", padx=(10, 0))
        tk.Label(title_col, text="Rainfall Frequency", bg=COLOR_SIDEBAR, fg="white", font=(FONT_FAMILY, 13, "bold")).pack(anchor="w")
        tk.Label(title_col, text="Analysis Studio", bg=COLOR_SIDEBAR, fg="#8FA8D6", font=(FONT_FAMILY, 9)).pack(anchor="w")
        tk.Frame(brand_pad, bg="#1E3A66", height=1).pack(fill="x", pady=(16, 0))

        # --- Workflow step navigation: numbered circular badges + connecting rail ---
        nav_wrap = tk.Frame(self.sidebar, bg=COLOR_SIDEBAR); nav_wrap.place(x=0, y=118, relwidth=1)
        self.nav_buttons = {}
        self.nav_badges = {}
        nav_items = [
            ("data", "1", "\U0001F4C2", "Data Source"),
            ("manager", "2", "\U0001F5C2", "Series Manager"),
            ("process", "3", "\u25B6", "Process & Results"),
            ("interactive", "4", "\u2728", "Interactive Studio"),
        ]
        for idx, (key, num, icon, label) in enumerate(nav_items):
            row_wrap = tk.Frame(nav_wrap, bg=COLOR_SIDEBAR, cursor="hand2")
            row_wrap.pack(fill="x")
            row = tk.Frame(row_wrap, bg=COLOR_SIDEBAR, padx=20, pady=11)
            row.pack(fill="x")
            badge = tk.Label(row, text=num, bg="#1E3A66", fg="#B9CBEA", font=(FONT_FAMILY, 10, "bold"), width=2, height=1)
            badge.pack(side="left")
            lbl = tk.Label(row, text=f"{icon}  {label}", bg=COLOR_SIDEBAR, fg="#C7D6EF", font=(FONT_FAMILY, 11), anchor="w")
            lbl.pack(side="left", padx=(10, 0), fill="x", expand=True)
            for w in (row_wrap, row, badge, lbl):
                w.bind("<Button-1>", lambda e, k=key: self._show_page(k))
                w.bind("<Enter>", lambda e, k=key: self._nav_hover(k, True))
                w.bind("<Leave>", lambda e, k=key: self._nav_hover(k, False))
            self.nav_buttons[key] = (row_wrap, row, lbl)
            self.nav_badges[key] = badge
            if idx < len(nav_items) - 1:
                rail = tk.Frame(nav_wrap, bg=COLOR_SIDEBAR, height=14); rail.pack(fill="x")
                tk.Frame(rail, bg="#1E3A66", width=2).place(x=29, y=0, relheight=1)

        # Workspace snapshot mini-card near the bottom of the sidebar
        self.sidebar_stat = tk.Frame(self.sidebar, bg="#0F274E", highlightbackground="#1E3A66", highlightthickness=1)
        self.sidebar_stat.place(x=20, rely=1.0, y=-92, relwidth=1, width=-40)
        tk.Label(self.sidebar_stat, text="WORKSPACE", bg="#0F274E", fg="#6E89BE", font=(FONT_FAMILY, 8, "bold")).pack(anchor="w", padx=12, pady=(10, 0))
        self.lbl_sidebar_series = tk.Label(self.sidebar_stat, text="0 series loaded", bg="#0F274E", fg="white", font=(FONT_FAMILY, 11, "bold")); self.lbl_sidebar_series.pack(anchor="w", padx=12)
        self.lbl_sidebar_ready = tk.Label(self.sidebar_stat, text="0 ready to process", bg="#0F274E", fg=COLOR_ACCENT, font=(FONT_FAMILY, 9)); self.lbl_sidebar_ready.pack(anchor="w", padx=12, pady=(0, 10))

        about_btn = tk.Label(self.sidebar, text="\u2139  About this App", bg=COLOR_SIDEBAR, fg="#8FA8D6", font=(FONT_FAMILY, 9), anchor="w", padx=20, pady=10, cursor="hand2")
        about_btn.place(x=0, rely=1.0, y=-34, relwidth=1)
        about_btn.bind("<Button-1>", lambda e: self._show_about())
        about_btn.bind("<Enter>", lambda e: about_btn.configure(fg="white"))
        about_btn.bind("<Leave>", lambda e: about_btn.configure(fg="#8FA8D6"))

        main = tk.Frame(self, bg=COLOR_BG)
        main.pack(side="left", fill="both", expand=True)
        self.header = tk.Frame(main, bg=COLOR_BG); self.header.pack(fill="x", padx=28, pady=(22, 4))
        head_top = tk.Frame(self.header, bg=COLOR_BG); head_top.pack(fill="x")
        title_side = tk.Frame(head_top, bg=COLOR_BG); title_side.pack(side="left", fill="x", expand=True)
        self.lbl_page_title = tk.Label(title_side, text="", bg=COLOR_BG, fg=COLOR_TEXT, font=(FONT_FAMILY, 19, "bold")); self.lbl_page_title.pack(anchor="w")
        self.lbl_page_subtitle = tk.Label(title_side, text="", bg=COLOR_BG, fg=COLOR_MUTED, font=(FONT_FAMILY, 10)); self.lbl_page_subtitle.pack(anchor="w", pady=(2, 0))
        self.lbl_header_badge = tk.Label(head_top, text="\u2728 v4.0 Aurora Studio", bg=COLOR_INFO_BG, fg=COLOR_INFO, font=(FONT_FAMILY, 9, "bold"), padx=10, pady=5)
        self.lbl_header_badge.pack(side="right", anchor="n")
        accent_bar = tk.Frame(self.header, bg=COLOR_ACCENT, height=3); accent_bar.pack(fill="x", pady=(14, 0))
        self.page_container = tk.Frame(main, bg=COLOR_BG); self.page_container.pack(fill="both", expand=True, padx=28, pady=(10, 18))

        self.pages = {}
        self.pages["data"] = self._build_page_data(self.page_container)
        self.pages["manager"] = self._build_page_manager(self.page_container)
        self.pages["process"] = self._build_page_process(self.page_container)
        self.pages["interactive"] = self._build_page_interactive(self.page_container)
        self.current_page = None

    def _nav_hover(self, key, entering):
        if self.current_page == key: return
        row_wrap, row, lbl = self.nav_buttons[key]
        color = COLOR_SIDEBAR_HOVER if entering else COLOR_SIDEBAR
        for w in (row_wrap, row): w.configure(bg=color)
        lbl.configure(bg=color)

    def _update_sidebar_stat(self):
        total = len(self.series_store)
        ready = sum(1 for df in self.series_store.values() if validate_series_df(df)[0])
        self.lbl_sidebar_series.configure(text=f"{total} series loaded")
        self.lbl_sidebar_ready.configure(text=f"\u2713 {ready} ready to process" if ready else "No series ready yet")

    def _show_page(self, key):
        titles = {
            "data": ("Data Source", "Import an Excel workbook and/or create series manually."),
            "manager": ("Series Manager", "Every series currently loaded, from any source."),
            "process": ("Process & Results", "Choose an output folder and run the full analysis."),
            "interactive": ("Interactive Studio", "Live preview of Data Quality, Goodness of Fit, and the Return Period curve, filterable by station.")
        }
        for k, (row_wrap, row, lbl) in self.nav_buttons.items():
            active = (k == key)
            bg = COLOR_SIDEBAR_ACTIVE if active else COLOR_SIDEBAR
            for w in (row_wrap, row): w.configure(bg=bg)
            lbl.configure(bg=bg, fg="white" if active else "#C7D6EF", font=(FONT_FAMILY, 11, "bold" if active else "normal"))
            self.nav_badges[k].configure(bg="white" if active else "#1E3A66", fg=COLOR_SIDEBAR_ACTIVE if active else "#B9CBEA")
        for k, frame in self.pages.items(): frame.pack_forget()
        self.pages[key].pack(fill="both", expand=True)
        self.lbl_page_title.configure(text=titles[key][0])
        self.lbl_page_subtitle.configure(text=titles[key][1])
        self.current_page = key
        self._update_sidebar_stat()

        if key == "manager": self._refresh_manager_table()
        if key == "process": self._refresh_process_summary()
        if key == "interactive": self._refresh_interactive_page()

    # --- PAGE 1: DATA ---
    def _build_page_data(self, parent):
        page = tk.Frame(parent, bg=COLOR_BG)
        cols_frame = tk.Frame(page, bg=COLOR_BG); cols_frame.pack(fill="both", expand=True); cols_frame.columnconfigure(0, weight=1); cols_frame.columnconfigure(1, weight=1); cols_frame.rowconfigure(0, weight=1)
        left = Card(cols_frame, title="Import from Excel", subtitle="Each sheet becomes one data series.", icon="\U0001F4E5", accent=COLOR_PRIMARY); left.grid(row=0, column=0, sticky="nsew", padx=(0, 10))
        make_pill_button(left.body, "Browse Excel File...", self._on_import_excel, bg=COLOR_PRIMARY, hover="#123A6B").pack(anchor="w", pady=(4, 10))
        self.lbl_import_file = tk.Label(left.body, text="No file imported yet.", bg=COLOR_CARD, fg=COLOR_MUTED, font=(FONT_FAMILY, 9)); self.lbl_import_file.pack(anchor="w", pady=(0, 8))
        preview_frame = tk.Frame(left.body, bg=COLOR_CARD); preview_frame.pack(fill="both", expand=True)
        self.import_tree = ttk.Treeview(preview_frame, columns=("sheet", "n", "mean", "std", "status"), show="headings", height=10)
        for c, t, w in [("sheet", "Sheet", 150), ("n", "N", 50), ("mean", "Mean", 80), ("std", "Std Dev", 80), ("status", "Status", 150)]: self.import_tree.heading(c, text=t); self.import_tree.column(c, width=w, anchor="center" if c != "sheet" else "w")
        vsb = ttk.Scrollbar(preview_frame, orient="vertical", command=self.import_tree.yview)
        self.import_tree.configure(yscrollcommand=vsb.set); self.import_tree.pack(side="left", fill="both", expand=True); vsb.pack(side="right", fill="y")

        right = Card(cols_frame, title="Manual Entry", subtitle="Type or paste your own rainfall series directly.", icon="\u270E", accent=COLOR_SECONDARY); right.grid(row=0, column=1, sticky="nsew", padx=(10, 0))
        row_btns = tk.Frame(right.body, bg=COLOR_CARD); row_btns.pack(fill="x", pady=(4, 10))
        make_pill_button(row_btns, "+ New Manual Series", self._on_new_manual_series, bg=COLOR_ACCENT, hover=COLOR_ACCENT_DARK).pack(side="left")
        tk.Label(right.body, text="Manual series created in this session:", bg=COLOR_CARD, fg=COLOR_MUTED, font=(FONT_FAMILY, 9)).pack(anchor="w", pady=(4, 4))
        manual_frame = tk.Frame(right.body, bg=COLOR_CARD); manual_frame.pack(fill="both", expand=True)
        self.manual_tree = ttk.Treeview(manual_frame, columns=("name", "n", "status"), show="headings", height=10)
        for c, t, w in [("name", "Series Name", 160), ("n", "N", 50), ("status", "Status", 150)]: self.manual_tree.heading(c, text=t); self.manual_tree.column(c, width=w, anchor="center" if c != "name" else "w")
        vsb2 = ttk.Scrollbar(manual_frame, orient="vertical", command=self.manual_tree.yview)
        self.manual_tree.configure(yscrollcommand=vsb2.set); self.manual_tree.pack(side="left", fill="both", expand=True); vsb2.pack(side="right", fill="y")
        self.manual_tree.bind("<Double-1>", lambda e: self._on_edit_manual_series())
        manual_btns2 = tk.Frame(right.body, bg=COLOR_CARD); manual_btns2.pack(fill="x", pady=(8, 0))
        make_pill_button(manual_btns2, "Edit Selected", self._on_edit_manual_series, bg="#475569", hover="#334155", font_size=9).pack(side="left", padx=(0, 6))
        make_pill_button(manual_btns2, "Delete Selected", self._on_delete_manual_series, bg=COLOR_DANGER, hover="#991B1B", font_size=9).pack(side="left")
        return page

    def _on_import_excel(self):
        path = filedialog.askopenfilename(title="Select Input Excel File", filetypes=[("Excel files", "*.xlsx *.xls"), ("All files", "*.*")])
        if not path: return
        self.lbl_import_file.configure(text=f"Reading: {os.path.basename(path)} ..."); self.update_idletasks()
        try: series, errors = read_excel_all_sheets(path)
        except Exception as e: messagebox.showerror("Import Failed", f"Could not read this Excel file:\n\n{e}"); self.lbl_import_file.configure(text="Import failed."); return
        for i in self.import_tree.get_children(): self.import_tree.delete(i)
        added = 0
        for name, df in series.items():
            u_name = self._unique_series_name(name)
            self.series_store[u_name] = df; self.series_source[u_name] = f"Excel: {os.path.basename(path)}"
            msg = validate_series_df(df)[1]; m = basic_moments(df[COL_DATA].to_numpy(dtype=float)) if len(df) >= 2 else None
            self.import_tree.insert("", "end", values=(u_name, len(df), f"{m['mean']:.2f}" if m else "-", f"{m['std']:.2f}" if m else "-", msg)); added += 1
        for name, err in errors.items(): self.import_tree.insert("", "end", values=(name, "-", "-", "-", f"Skipped: {err}"))
        self.lbl_import_file.configure(text=f"Imported {added} series from '{os.path.basename(path)}'" + (f" ({len(errors)} skipped)." if errors else "."))
        if self.current_page == "manager": self._refresh_manager_table()
    def _unique_series_name(self, base_name):
        name = str(base_name).strip() or "Series"
        if name not in self.series_store: return name
        i = 2
        while f"{name} ({i})" in self.series_store: i += 1
        return f"{name} ({i})"
    def _on_new_manual_series(self):
        dlg = tk.Toplevel(self); dlg.title("New Manual Series"); dlg.configure(bg=COLOR_BG); dlg.geometry("340x160"); dlg.transient(self); dlg.grab_set()
        tk.Label(dlg, text="Series / Station Name", bg=COLOR_BG, fg=COLOR_TEXT, font=(FONT_FAMILY, 10, "bold")).pack(pady=(20, 4))
        e_name = tk.Entry(dlg, font=(FONT_FAMILY, 11), justify="center"); e_name.pack(padx=30, fill="x"); e_name.focus_set()
        def on_ok():
            n = e_name.get().strip()
            if not n: messagebox.showerror("Invalid Name", "Please enter a series name.", parent=dlg); return
            un = self._unique_series_name(n); dlg.destroy(); self.series_store[un] = empty_series_df(); self.series_source[un] = "Manual"
            self._refresh_manual_tree(); self._open_manual_editor(un)
        btns = tk.Frame(dlg, bg=COLOR_BG); btns.pack(pady=18)
        make_pill_button(btns, "Create & Open Editor", on_ok, bg=COLOR_ACCENT, hover=COLOR_ACCENT_DARK).pack(side="left", padx=6)
        make_pill_button(btns, "Cancel", dlg.destroy, bg="#94A3B8").pack(side="left", padx=6)
        dlg.bind("<Return>", lambda e: on_ok())
    def _open_manual_editor(self, name):
        def on_save(sn, ndf):
            self.series_store[sn] = ndf.reset_index(drop=True); self._refresh_manual_tree()
            if self.current_page == "manager": self._refresh_manager_table()
            if self.current_page == "process": self._refresh_process_summary()
        ManualSeriesEditor(self, name, self.series_store.get(name, empty_series_df()), on_save)
    def _refresh_manual_tree(self):
        for i in self.manual_tree.get_children(): self.manual_tree.delete(i)
        for n, src in self.series_source.items():
            if src == "Manual": self.manual_tree.insert("", "end", values=(n, len(self.series_store[n]), validate_series_df(self.series_store[n])[1]))
    def _selected_manual_name(self):
        sel = self.manual_tree.selection()
        if not sel: messagebox.showinfo("No Selection", "Please select a manual series first."); return None
        return self.manual_tree.item(sel[0])["values"][0]
    def _on_edit_manual_series(self):
        n = self._selected_manual_name()
        if n: self._open_manual_editor(str(n))
    def _on_delete_manual_series(self):
        n = self._selected_manual_name()
        if n and messagebox.askyesno("Confirm Delete", f"Remove manual series '{n}'?"):
            self.series_store.pop(n, None); self.series_source.pop(n, None); self._refresh_manual_tree()
            if self.current_page == "manager": self._refresh_manager_table()

    # --- PAGE 2: MANAGER ---
    def _build_page_manager(self, parent):
        page = tk.Frame(parent, bg=COLOR_BG)
        summary = tk.Frame(page, bg=COLOR_BG); summary.pack(fill="x", pady=(0, 12))
        kpi_specs = [
            ("\U0001F4E6", COLOR_PRIMARY, COLOR_CARD, "Total series loaded"),
            ("\u2705", COLOR_SUCCESS, COLOR_SUCCESS_BG, "Ready to process"),
            ("\u26A0", COLOR_AMBER, COLOR_WARN_BG, "Not ready (needs data)"),
        ]
        self.manager_summary_labels = []
        for i, (icon, accent, tint, caption) in enumerate(kpi_specs):
            kpi = make_kpi_card(summary, icon, "0", caption, accent=accent, tint=tint)
            kpi.pack(side="left", fill="x", expand=True, padx=(0 if i == 0 else 8, 0))
            self.manager_summary_labels.append((kpi.value_label, None))
        card = Card(page, title="All Loaded Series", subtitle="Combined list. Select a row to rename, edit, or remove it.", icon="\U0001F5C2", accent=COLOR_PRIMARY); card.pack(fill="both", expand=True)
        toolbar = tk.Frame(card.body, bg=COLOR_CARD); toolbar.pack(fill="x", pady=(0, 8))
        make_pill_button(toolbar, "\u21BB Refresh", self._refresh_manager_table, bg="#475569", hover="#334155", font_size=9).pack(side="left", padx=(0, 6))
        make_pill_button(toolbar, "Edit / View", self._manager_edit_selected, bg=COLOR_PRIMARY, hover="#123A6B", font_size=9).pack(side="left", padx=6)
        make_pill_button(toolbar, "Rename", self._manager_rename_selected, bg="#475569", hover="#334155", font_size=9).pack(side="left", padx=6)
        make_pill_button(toolbar, "Remove", self._manager_remove_selected, bg=COLOR_DANGER, hover="#991B1B", font_size=9).pack(side="left", padx=6)
        make_pill_button(toolbar, "Remove All", self._manager_remove_all, bg=COLOR_DANGER, hover="#991B1B", font_size=9).pack(side="left", padx=6)
        table_frame = tk.Frame(card.body, bg=COLOR_CARD); table_frame.pack(fill="both", expand=True)
        cols = ("name", "source", "n", "mean", "std", "cv", "cs", "status")
        self.manager_tree = ttk.Treeview(table_frame, columns=cols, show="headings", height=14)
        for c, t, w in [("name", "Series / Station", 190), ("source", "Source", 170), ("n", "N", 50), ("mean", "Mean", 80), ("std", "Std Dev", 80), ("cv", "Cv", 60), ("cs", "Cs", 60), ("status", "Status", 170)]:
            self.manager_tree.heading(c, text=t); self.manager_tree.column(c, width=w, anchor="w" if c in ("name", "source") else "center")
        vsb = ttk.Scrollbar(table_frame, orient="vertical", command=self.manager_tree.yview)
        self.manager_tree.configure(yscrollcommand=vsb.set); self.manager_tree.pack(side="left", fill="both", expand=True); vsb.pack(side="right", fill="y")
        self.manager_tree.bind("<Double-1>", lambda e: self._manager_edit_selected())
        return page

    def _refresh_manager_table(self):
        for i in self.manager_tree.get_children(): self.manager_tree.delete(i)
        ready_count = 0
        for name, df in self.series_store.items():
            src = self.series_source.get(name, "-")
            ok, msg = validate_series_df(df)
            if ok: ready_count += 1
            if len(df) >= 2: x = df[COL_DATA].to_numpy(dtype=float); mean_txt, std_txt = f"{x.mean():.2f}", f"{x.std(ddof=1):.2f}"
            else: mean_txt = std_txt = "-"
            if len(df) >= 3: m = basic_moments(df[COL_DATA].to_numpy(dtype=float)); cv_txt, cs_txt = f"{m['cv']:.2f}", f"{m['cs']:.2f}"
            else: cv_txt = cs_txt = "-"
            self.manager_tree.insert("", "end", values=(name, src, len(df), mean_txt, std_txt, cv_txt, cs_txt, msg))
        total = len(self.series_store)
        self.manager_summary_labels[0][0].configure(text=str(total)); self.manager_summary_labels[1][0].configure(text=str(ready_count)); self.manager_summary_labels[2][0].configure(text=str(total - ready_count))

    def _selected_manager_name(self):
        sel = self.manager_tree.selection()
        if not sel: messagebox.showinfo("No Selection", "Please select a series first."); return None
        return self.manager_tree.item(sel[0])["values"][0]
    def _manager_edit_selected(self):
        n = self._selected_manager_name()
        if n: self._open_manual_editor(str(n))
    def _manager_rename_selected(self):
        name = self._selected_manager_name()
        if not name: return
        new_name = tk.simpledialog.askstring("Rename Series", "New name:", initialvalue=str(name), parent=self) if hasattr(tk, "simpledialog") else None
        if not new_name or not new_name.strip() or new_name.strip() == str(name): return
        new_name = new_name.strip()
        if new_name in self.series_store: messagebox.showerror("Name In Use", "A series with that name already exists."); return
        self.series_store[new_name] = self.series_store.pop(str(name)); self.series_source[new_name] = self.series_source.pop(str(name), "-")
        self._refresh_manager_table(); self._refresh_manual_tree()
    def _manager_remove_selected(self):
        n = self._selected_manager_name()
        if n and messagebox.askyesno("Confirm Remove", f"Remove series '{n}' from the workspace?"):
            self.series_store.pop(str(n), None); self.series_source.pop(str(n), None)
            self._refresh_manager_table(); self._refresh_manual_tree()
    def _manager_remove_all(self):
        if self.series_store and messagebox.askyesno("Confirm Remove All", "Remove ALL loaded series from the workspace?"):
            self.series_store.clear(); self.series_source.clear()
            self._refresh_manager_table(); self._refresh_manual_tree()

    # --- PAGE 3: PROCESS ---
    def _build_page_process(self, parent):
        page = tk.Frame(parent, bg=COLOR_BG)
        kpi_row = tk.Frame(page, bg=COLOR_BG); kpi_row.pack(fill="x", pady=(0, 12))
        kpi_specs = [
            ("\U0001F4E6", COLOR_PRIMARY, COLOR_CARD, "Series loaded"),
            ("\u2705", COLOR_SUCCESS, COLOR_SUCCESS_BG, "Ready to process"),
            ("\U0001F4C1", COLOR_SECONDARY, COLOR_CARD, "Output folder"),
        ]
        self.process_kpi_labels = []
        for i, (icon, accent, tint, caption) in enumerate(kpi_specs):
            kpi = make_kpi_card(kpi_row, icon, "0", caption, accent=accent, tint=tint)
            kpi.pack(side="left", fill="x", expand=True, padx=(0 if i == 0 else 8, 0))
            self.process_kpi_labels.append(kpi.value_label)

        top = tk.Frame(page, bg=COLOR_BG); top.pack(fill="x")
        card_out = Card(top, title="Output Folder", subtitle="Where the Excel reports will be saved.", icon="\U0001F4C1", accent=COLOR_SECONDARY); card_out.pack(side="left", fill="both", expand=True, padx=(0, 10))
        row = tk.Frame(card_out.body, bg=COLOR_CARD); row.pack(fill="x")
        make_pill_button(row, "Choose Folder...", self._on_choose_outdir, bg=COLOR_PRIMARY, hover="#123A6B", font_size=9).pack(side="left")
        self.lbl_outdir = tk.Label(card_out.body, text="(not selected)", bg=COLOR_CARD, fg=COLOR_MUTED, font=(FONT_FAMILY, 9), wraplength=380, justify="left"); self.lbl_outdir.pack(anchor="w", pady=(8, 0))
        card_run = Card(top, title="Run Analysis", subtitle="Processes every series currently loaded.", icon="\u25B6", accent=COLOR_ACCENT); card_run.pack(side="left", fill="both", expand=True, padx=(10, 0))
        self.btn_process = make_pill_button(card_run.body, "\u25B6  Process All Series", self._on_process, bg=COLOR_ACCENT, hover=COLOR_ACCENT_DARK, state="disabled"); self.btn_process.pack(anchor="w")
        self.lbl_process_summary = tk.Label(card_run.body, text="0 series ready.", bg=COLOR_CARD, fg=COLOR_MUTED, font=(FONT_FAMILY, 9)); self.lbl_process_summary.pack(anchor="w", pady=(8, 0))
        progress_card = Card(page, title="Progress", icon="\u23F1", accent=COLOR_AMBER); progress_card.pack(fill="x", pady=(10, 10))
        self.progress = ttk.Progressbar(progress_card.body, mode="determinate", style="Aurora.Horizontal.TProgressbar"); self.progress.pack(fill="x", pady=(0, 6))
        self.lbl_status = tk.Label(progress_card.body, text="Waiting to start.", bg=COLOR_CARD, fg=COLOR_MUTED, font=(FONT_FAMILY, 9), anchor="w"); self.lbl_status.pack(anchor="w")
        log_card = Card(page, title="Processing Log", subtitle="Live console output while your series are being analysed.", icon="\U0001F4DC", accent="#334155"); log_card.pack(fill="both", expand=True)
        self.txt_log = scrolledtext.ScrolledText(log_card.body, height=14, state="disabled", font=("Consolas", 9), wrap="word", bg="#0f172a", fg="#5EEAD4", insertbackground="#e2e8f0", relief="flat"); self.txt_log.pack(fill="both", expand=True)
        bottom = tk.Frame(page, bg=COLOR_BG); bottom.pack(fill="x", pady=(10, 0))
        self.btn_open_result = make_pill_button(bottom, "\U0001F5C2  Open Result Folder", self._on_open_result, bg="#475569", hover="#334155", state="disabled"); self.btn_open_result.pack(side="left")
        return page

    def _refresh_process_summary(self):
        total = len(self.series_store)
        ready = sum(1 for df in self.series_store.values() if validate_series_df(df)[0])
        self.lbl_process_summary.configure(text=f"{ready} of {total} loaded series are ready to process." if total else "No series loaded yet - go to 'Data Source' first.")
        if hasattr(self, "process_kpi_labels"):
            self.process_kpi_labels[0].configure(text=str(total))
            self.process_kpi_labels[1].configure(text=str(ready))
            self.process_kpi_labels[2].configure(text="Set" if self.out_dir else "Not set")
        self._update_process_button()
        self._update_sidebar_stat()
    def _on_choose_outdir(self):
        d = filedialog.askdirectory(title="Select Output Folder")
        if d: self.out_dir = d; self.lbl_outdir.configure(text=d); self._update_process_button()
    def _update_process_button(self):
        self.btn_process.set_state("normal" if self.out_dir and any(validate_series_df(df)[0] for df in self.series_store.values()) else "disabled")
    def _append_log_direct(self, msg):
        self.txt_log.configure(state="normal"); self.txt_log.insert("end", str(msg) + "\n"); self.txt_log.see("end"); self.txt_log.configure(state="disabled")
    def _clear_log(self):
        self.txt_log.configure(state="normal"); self.txt_log.delete("1.0", "end"); self.txt_log.configure(state="disabled")
    def _log(self, msg): self.msg_queue.put(("log", msg))
    def _set_progress(self, done, total): self.msg_queue.put(("progress", (done, total)))
    def _poll_queue(self):
        try:
            while True:
                kind, payload = self.msg_queue.get_nowait()
                if kind == "log": self._append_log_direct(payload)
                elif kind == "progress":
                    done, total = payload
                    self.progress["maximum"] = max(total, 1); self.progress["value"] = done
                    self.lbl_status.configure(text=f"Processing... ({done}/{total} series)")
                elif kind == "done":
                    success, failed, out_dir = payload
                    self._on_process_finished(success, failed, out_dir)
        except queue.Empty: pass
        self.after(100, self._poll_queue)

    def _on_process(self):
        if not self.out_dir: return
        ready_series = {name: df for name, df in self.series_store.items() if validate_series_df(df)[0]}
        if not ready_series: messagebox.showwarning(APP_TITLE, "No series are ready to process (each series needs at least 5 rows)."); return
        if self.worker_thread and self.worker_thread.is_alive(): return
        self.btn_process.set_state("disabled"); self.btn_open_result.set_state("disabled")
        self.progress["value"] = 0; self.lbl_status.configure(text="Processing..."); self._clear_log()
        self._append_log_direct(">>> Starting analysis for all ready series...\n")
        def worker():
            try:
                success, failed = run_batch_from_series(ready_series, self.out_dir, self._log, self._set_progress)
                self.msg_queue.put(("done", (success, failed, self.out_dir)))
            except Exception:
                self._log(traceback.format_exc())
                self.msg_queue.put(("done", ([], [("FATAL", "See the Processing Log panel for details")], self.out_dir)))
        self.worker_thread = threading.Thread(target=worker, daemon=True); self.worker_thread.start()

    def _on_process_finished(self, success, failed, out_dir):
        self.btn_process.set_state("normal"); self.btn_open_result.set_state("normal")
        if failed and not success:
            self.lbl_status.configure(text="Finished with errors.")
            messagebox.showerror(APP_TITLE, f"Processing finished with errors.\n\nSucceeded: {len(success)}\nFailed: {len(failed)}\n\nSee the 'Processing Log' panel for details.")
        elif failed:
            self.lbl_status.configure(text="Finished (some series failed).")
            messagebox.showwarning(APP_TITLE, f"Processing finished.\n\nSucceeded: {len(success)}\nFailed: {len(failed)}\n\nSee the 'Processing Log' panel for details.")
        else:
            self.lbl_status.configure(text="Finished. All series processed successfully.")
            messagebox.showinfo(APP_TITLE, f"Processing complete!\n\nSucceeded: {len(success)} series\n\nResults saved to:\n{out_dir}")

    def _on_open_result(self):
        if self.out_dir: _open_in_explorer(self.out_dir)

    # --- PAGE 4: INTERACTIVE STUDIO ---
    def _build_page_interactive(self, parent):
        page = tk.Frame(parent, bg=COLOR_BG)

        # Station selector, styled as a compact toolbar card
        ctrl_card = tk.Frame(page, bg=COLOR_BORDER)
        ctrl_card.pack(fill="x", pady=(0, 12))
        strip = tk.Frame(ctrl_card, bg=COLOR_SECONDARY, height=3); strip.pack(fill="x")
        ctrl_inner = tk.Frame(ctrl_card, bg=COLOR_CARD); ctrl_inner.pack(fill="x", padx=1, pady=(0, 1))
        ctrl_frame = tk.Frame(ctrl_inner, bg=COLOR_CARD); ctrl_frame.pack(fill="x", padx=16, pady=12)

        icon_badge = tk.Label(ctrl_frame, text="\U0001F4CD", bg=COLOR_SECONDARY, fg="white", font=(FONT_FAMILY, 11, "bold"), width=3)
        icon_badge.pack(side="left")
        tk.Label(ctrl_frame, text="Select Station Data:", bg=COLOR_CARD, fg=COLOR_TEXT, font=(FONT_FAMILY, 11, "bold")).pack(side="left", padx=(10, 0))
        self.combo_station = ttk.Combobox(ctrl_frame, state="readonly", width=40, font=(FONT_FAMILY, 10))
        self.combo_station.pack(side="left", padx=10)
        self.combo_station.bind("<<ComboboxSelected>>", self._on_interactive_station_changed)

        # KPI strip: quick at-a-glance statistics for the selected station
        self.interactive_kpi_frame = tk.Frame(page, bg=COLOR_BG)
        self.interactive_kpi_frame.pack(fill="x", pady=(0, 12))

        # Area Konten Tab (CTkFrame)
        self.interactive_content = ctk.CTkFrame(page, fg_color="transparent")
        self.interactive_content.pack(fill="both", expand=True)
        
        return page

    def _refresh_interactive_kpis(self, result):
        for w in self.interactive_kpi_frame.winfo_children(): w.destroy()
        ps = result["pandas_stats"]
        best_name, *_ = pilih_distribusi_terbaik_ks(result)
        specs = [
            ("\U0001F522", f"{ps['n']:.0f}", "Data points (N)", COLOR_PRIMARY, COLOR_CARD),
            ("\U0001F4CF", f"{ps['mean']:.2f}", "Mean (mm)", COLOR_ACCENT, COLOR_CARD),
            ("\U0001F4C9", f"{ps['std']:.2f}", "Std. Deviation", COLOR_SECONDARY, COLOR_CARD),
            ("\u2696", f"{ps['cs']:.3f}", "Skewness (Cs)", COLOR_AMBER, COLOR_CARD),
            ("\U0001F3C6", best_name or "-", "Recommended Fit", COLOR_INFO, COLOR_INFO_BG),
        ]
        for i, (icon, val, cap, accent, tint) in enumerate(specs):
            kpi = make_kpi_card(self.interactive_kpi_frame, icon, val, cap, accent=accent, tint=tint)
            kpi.pack(side="left", fill="x", expand=True, padx=(0 if i == 0 else 6, 0))
            kpi.value_label.configure(font=(FONT_FAMILY, 16 if isinstance(val, str) and len(str(val)) > 6 else 20, "bold"))

    def _refresh_interactive_page(self):
        ready_series = [name for name, df in self.series_store.items() if validate_series_df(df)[0]]
        self.combo_station["values"] = ready_series
        if ready_series:
            if self.combo_station.get() not in ready_series:
                self.combo_station.current(0)
            self._on_interactive_station_changed(None)
        else:
            self.combo_station.set("")
            for widget in self.interactive_content.winfo_children(): widget.destroy()
            for widget in self.interactive_kpi_frame.winfo_children(): widget.destroy()
            if self.current_interactive_canvas:
                self.current_interactive_canvas = None
            empty = tk.Frame(self.interactive_content, bg=COLOR_BG); empty.pack(fill="both", expand=True)
            tk.Label(empty, text="\U0001F327", bg=COLOR_BG, font=(FONT_FAMILY, 34)).pack(pady=(60, 6))
            tk.Label(empty, text="No station data ready for analysis yet.", bg=COLOR_BG, fg=COLOR_TEXT, font=(FONT_FAMILY, 13, "bold")).pack()
            tk.Label(empty, text="Import or type at least 5 rows for a station on the Data Source page first.", bg=COLOR_BG, fg=COLOR_MUTED, font=(FONT_FAMILY, 9)).pack(pady=(2, 0))

    def _on_interactive_station_changed(self, event):
        name = self.combo_station.get()
        if not name or name not in self.series_store: return
        
        # Clear existing tabs and charts to free memory
        for widget in self.interactive_content.winfo_children():
            widget.destroy()
        if self.current_interactive_canvas:
            self.current_interactive_canvas = None

        df = self.series_store[name]
        try:
            result = run_analysis(df)
        except Exception as e:
            tk.Label(self.interactive_content, text=f"Error analyzing data:\n{e}", bg=COLOR_BG, fg="red").pack(pady=40)
            return

        self._refresh_interactive_kpis(result)

        # Main Tabview for Interactive
        self.interactive_tabview = ctk.CTkTabview(self.interactive_content)
        self.interactive_tabview.pack(fill="both", expand=True)
        
        tab_dq = self.interactive_tabview.add("📊 Data Quality")
        tab_gof = self.interactive_tabview.add("🏆 GOF & Selection")
        tab_freq = self.interactive_tabview.add("📈 Return Period")
        
        self._build_interactive_dq(tab_dq, result)
        self._build_interactive_gof(tab_gof, result)
        self._build_interactive_freq(tab_freq, name, result)

    def _make_dq_card(self, parent, icon, title, details, status, pass_keywords):
        is_pass = any(kw in status for kw in pass_keywords) and "Not" not in status and "Non" not in status

        # Desain warna dinamis (Hijau Pastel untuk Pass, Merah Pastel untuk Fail)
        bg_color = "#F0FDF4" if is_pass else "#FEF2F2"
        border_color = "#BBF7D0" if is_pass else "#FECACA"
        status_color = "#15803D" if is_pass else "#B91C1C"
        status_bg = "#DCFCE7" if is_pass else "#FEE2E2"
        result_icon = "✅" if is_pass else "❌"

        card = ctk.CTkFrame(parent, fg_color=bg_color, corner_radius=12, border_width=1, border_color=border_color)
        card.pack(fill="x", pady=6)

        icon_col = ctk.CTkFrame(card, fg_color=status_color, corner_radius=10, width=42, height=42)
        icon_col.pack(side="left", padx=(14, 10), pady=14)
        icon_col.pack_propagate(False)
        ctk.CTkLabel(icon_col, text=icon, font=("Segoe UI", 16), text_color="white").pack(expand=True)

        left_frame = ctk.CTkFrame(card, fg_color="transparent")
        left_frame.pack(side="left", fill="both", expand=True, padx=(0, 15), pady=12)

        ctk.CTkLabel(left_frame, text=title, font=("Segoe UI", 14, "bold"), text_color="#0F172A").pack(anchor="w")
        ctk.CTkLabel(left_frame, text=details, font=("Consolas", 12), text_color="#475569", justify="left").pack(anchor="w", pady=(4,0))

        pill = ctk.CTkFrame(card, fg_color=status_bg, corner_radius=20)
        pill.pack(side="right", padx=20)
        ctk.CTkLabel(pill, text=f" {'PASS' if is_pass else 'FAIL'} {result_icon} ", font=("Segoe UI", 13, "bold"), text_color=status_color).pack(padx=6, pady=4)
        return is_pass

    def _build_interactive_dq(self, parent, result):
        scroll = ctk.CTkScrollableFrame(parent, fg_color="transparent")
        scroll.pack(fill="both", expand=True, padx=10, pady=10)
        dq = result["dq_tests"]

        def _passed(status, keywords):
            return any(kw in status for kw in keywords) and "Not" not in status and "Non" not in status

        tests = [
            ("outlier", "No Outliers"), ("trend", "No Trend"), ("var", "Homogeneous"),
            ("mean", "Homogeneous"), ("indep", "Independent"),
        ]
        n_pass = sum(1 for key, kw in tests if _passed(dq[key]["summary"]["Conclusion"], [kw]))
        summary_color = COLOR_SUCCESS if n_pass == len(tests) else (COLOR_AMBER if n_pass >= len(tests) - 1 else COLOR_DANGER)
        summary = ctk.CTkFrame(scroll, fg_color=summary_color, corner_radius=10)
        summary.pack(fill="x", pady=(0, 10))
        ctk.CTkLabel(summary, text=f"\U0001F4CB  Data Quality Overview: {n_pass} of {len(tests)} tests passed",
                     font=("Segoe UI", 13, "bold"), text_color="white").pack(padx=16, pady=10, anchor="w")

        self._make_dq_card(scroll, "\U0001F50E", "1. Outlier Test (Grubbs / WRC)",
                           f"Lower Limit (Xl) : {dq['outlier']['summary']['Lower Limit Xl']:.3f}\n"
                           f"Upper Limit (Xh) : {dq['outlier']['summary']['Upper Limit Xh']:.3f}\n"
                           f"Table Kn Value   : {dq['outlier']['summary']['Kn (Table)']:.3f}",
                           dq["outlier"]["summary"]["Conclusion"], ["No Outliers"])

        self._make_dq_card(scroll, "\U0001F4C8", "2. Trend Test (Spearman's Rho)",
                           f"Rho (Spearman)   : {dq['trend']['summary']['Rho (Spearman)']:.3f}\n"
                           f"T Calculated     : {dq['trend']['summary']['t Calculated']:.3f}\n"
                           f"T Critical       : {dq['trend']['summary']['t Critical']:.3f}",
                           dq["trend"]["summary"]["Conclusion"], ["No Trend"])

        self._make_dq_card(scroll, "\u2696", "3. Variance Homogeneity Test (F-Test)",
                           f"F Calculated     : {dq['var']['summary']['F Calculated']:.3f}\n"
                           f"F Critical       : {dq['var']['summary']['F Critical']:.3f}",
                           dq["var"]["summary"]["Conclusion"], ["Homogeneous"])

        self._make_dq_card(scroll, "\U0001F4CA", "4. Mean Homogeneity Test (T-Test)",
                           f"T Calculated     : {dq['mean']['summary']['T Calculated']:.3f}\n"
                           f"T Critical       : {dq['mean']['summary']['T Critical']:.3f}",
                           dq["mean"]["summary"]["Conclusion"], ["Homogeneous"])

        self._make_dq_card(scroll, "\U0001F517", "5. Independence Test (Lag-1 Serial Correlation)",
                           f"Correlation (r1) : {dq['indep']['summary']['r1 (Correlation)']:.3f}\n"
                           f"Lower Limit (Ll) : {dq['indep']['summary']['Lower Limit (Ll)']:.3f}\n"
                           f"Upper Limit (Ul) : {dq['indep']['summary']['Upper Limit (Ul)']:.3f}",
                           dq["indep"]["summary"]["Conclusion"], ["Independent"])

    def _build_interactive_gof(self, parent, result):
        scroll = ctk.CTkScrollableFrame(parent, fg_color="transparent")
        scroll.pack(fill="both", expand=True, padx=10, pady=10)
        
        # 1. Kartu Highlight Distribusi Terbaik (Warna Biru Modern)
        best_name, quant_row, in_range, rank_val, ks5_pass = pilih_distribusi_terbaik_ks(result)
        
        best_card = ctk.CTkFrame(scroll, fg_color="#F0F9FF", corner_radius=12, border_width=2, border_color="#7DD3FC")
        best_card.pack(fill="x", pady=(0, 20))
        
        header_frame = ctk.CTkFrame(best_card, fg_color="transparent")
        header_frame.pack(fill="x", padx=20, pady=(15, 5))
        ctk.CTkLabel(header_frame, text="🏆 RECOMMENDED DISTRIBUTION", font=("Segoe UI", 16, "bold"), text_color="#0369A1").pack(side="left")
        
        content_frame = ctk.CTkFrame(best_card, fg_color="transparent")
        content_frame.pack(fill="x", padx=20, pady=(0, 15))
        
        r100_ratio = quant_row.get('R100/R2', np.nan) if quant_row is not None else np.nan
        ratio_status = "✅ Valid Range" if in_range == "Yes" else ("⚠️ Approaching" if in_range == "Approaching" else "❌ Invalid Range")
        
        ctk.CTkLabel(content_frame, text=f"Selected Method : {best_name or 'None'}", font=("Segoe UI", 14, "bold"), text_color="#0F172A").grid(row=0, column=0, sticky="w", pady=2)
        ctk.CTkLabel(content_frame, text=f"KS Test Rank    : {rank_val or '-'}", font=("Consolas", 13), text_color="#334155").grid(row=1, column=0, sticky="w", pady=2)
        ctk.CTkLabel(content_frame, text=f"R100/R2 Ratio   : {r100_ratio:.3f} ({ratio_status})", font=("Consolas", 13), text_color="#334155").grid(row=2, column=0, sticky="w", pady=2)

        # 2. Detail GOF Seluruh Distribusi
        ctk.CTkLabel(scroll, text="📊 Detailed Goodness of Fit (α = 5%)", font=("Segoe UI", 16, "bold"), text_color="#1E293B").pack(anchor="w", pady=(10, 10))
        
        chi_df = result["chi_df"]
        ks_df = result["ks_df"]
        n_data = result["pandas_stats"]["n"]
        
        # Hitung Nilai Kritis KS secara universal
        try: ks_crit_05 = stats.kstwo.ppf(0.95, n_data)
        except AttributeError: ks_crit_05 = stats.ksone.ppf(0.975, n_data)
        
        for _, row in ks_df.iterrows():
            dist = row["Distribution"]
            chi_row = chi_df[chi_df["Distribution"] == dist].iloc[0] if len(chi_df[chi_df["Distribution"] == dist]) else None
            
            card = ctk.CTkFrame(scroll, fg_color="white", corner_radius=10, border_width=1, border_color="#CBD5E1")
            card.pack(fill="x", pady=6)
            
            # Judul Distribusi
            ctk.CTkLabel(card, text=f"📌 {dist}", font=("Segoe UI", 15, "bold"), text_color="#0F172A").pack(anchor="w", padx=15, pady=(10, 5))
            
            # Membagi Layar Jadi 2 Kolom (Kiri: KS, Kanan: Chi-Square)
            grid_frame = ctk.CTkFrame(card, fg_color="transparent")
            grid_frame.pack(fill="x", padx=15, pady=(0, 10))
            grid_frame.columnconfigure(0, weight=1)
            grid_frame.columnconfigure(1, weight=1)
            
            # ----- DATA KOLMOGOROV-SMIRNOV -----
            ks_dmax = row['Dmax']
            ks_pass = row.get('Pass KS 5%', 'Fail')
            ks_color = "#15803D" if ks_pass == "Pass" else "#B91C1C"
            ks_icon = "✅" if ks_pass == "Pass" else "❌"
            
            ks_frame = ctk.CTkFrame(grid_frame, fg_color="#F8FAFC", corner_radius=6, border_width=1, border_color="#E2E8F0")
            ks_frame.grid(row=0, column=0, sticky="nsew", padx=(0, 5))
            ctk.CTkLabel(ks_frame, text="Kolmogorov-Smirnov", font=("Segoe UI", 12, "bold"), text_color="#475569").pack(anchor="w", padx=10, pady=(5,0))
            
            # Pisahkan teks baris atas dan baris Result agar sejajar sempurna
            ks_text_top = f"Dmax Calc : {ks_dmax:.4f}\nDmax Crit : {ks_crit_05:.4f}"
            ctk.CTkLabel(ks_frame, text=ks_text_top, font=("Consolas", 11), text_color="#334155", justify="left").pack(anchor="w", padx=10, pady=(2, 0))
            
            ks_res_frame = ctk.CTkFrame(ks_frame, fg_color="transparent")
            ks_res_frame.pack(anchor="w", padx=10, pady=(0, 5))
            ctk.CTkLabel(ks_res_frame, text="Result    : ", font=("Consolas", 11), text_color="#334155").pack(side="left")
            ctk.CTkLabel(ks_res_frame, text=f"{ks_pass} {ks_icon}", font=("Consolas", 11, "bold"), text_color=ks_color).pack(side="left")

            # ----- DATA CHI-SQUARE -----
            chi_stat = chi_row['Chi-Square Statistic'] if chi_row is not None else np.nan
            dof = chi_row['dof'] if chi_row is not None else np.nan
            
            if pd.isna(chi_stat) or pd.isna(dof) or dof <= 0:
                chi_crit = np.nan
                chi_pass = "Fail"
            else:
                chi_crit = stats.chi2.ppf(0.95, dof)
                chi_pass = "Pass" if chi_stat < chi_crit else "Fail"
                
            chi_color = "#15803D" if chi_pass == "Pass" else "#B91C1C"
            chi_icon = "✅" if chi_pass == "Pass" else "❌"
            
            chi_frame = ctk.CTkFrame(grid_frame, fg_color="#F8FAFC", corner_radius=6, border_width=1, border_color="#E2E8F0")
            chi_frame.grid(row=0, column=1, sticky="nsew", padx=(5, 0))
            ctk.CTkLabel(chi_frame, text=f"Chi-Square (dof={dof})", font=("Segoe UI", 12, "bold"), text_color="#475569").pack(anchor="w", padx=10, pady=(5,0))
            
            if pd.isna(chi_stat):
                chi_text_top = f"Chi² Calc : N/A\nChi² Crit : N/A"
            else:
                chi_text_top = f"Chi² Calc : {chi_stat:.4f}\nChi² Crit : {chi_crit:.4f}"
                
            ctk.CTkLabel(chi_frame, text=chi_text_top, font=("Consolas", 11), text_color="#334155", justify="left").pack(anchor="w", padx=10, pady=(2, 0))
            
            chi_res_frame = ctk.CTkFrame(chi_frame, fg_color="transparent")
            chi_res_frame.pack(anchor="w", padx=10, pady=(0, 5))
            ctk.CTkLabel(chi_res_frame, text="Result    : ", font=("Consolas", 11), text_color="#334155").pack(side="left")
            ctk.CTkLabel(chi_res_frame, text=f"{chi_pass} {chi_icon}", font=("Consolas", 11, "bold"), text_color=chi_color).pack(side="left")

    def _build_interactive_freq(self, parent, name, result):
        sidebar = ctk.CTkFrame(parent, width=230, fg_color="white", corner_radius=10, border_width=1, border_color="#D9E1EC")
        sidebar.pack(side="left", fill="y", padx=(10, 5), pady=10)
        sidebar.pack_propagate(False)

        main_graph = ctk.CTkFrame(parent, fg_color="transparent")
        main_graph.pack(side="right", fill="both", expand=True, padx=(5, 10), pady=10)

        header_bar = ctk.CTkFrame(sidebar, fg_color=COLOR_SECONDARY, corner_radius=10)
        header_bar.pack(fill="x", padx=12, pady=(16, 4))
        ctk.CTkLabel(header_bar, text="\U0001F3AF  Filter Distributions", font=("Segoe UI", 13, "bold"), text_color="white").pack(padx=12, pady=10, anchor="w")

        T_emp, x_emp, Tline, curves = build_return_period_curve_data(
            result["x"], result["fits"], RETURN_PERIODS, n_curve_points=100
        )

        fig = Figure(figsize=(7, 4.2), dpi=100)
        fig.patch.set_facecolor(COLOR_BG)
        ax = fig.add_subplot(111)
        ax.set_facecolor('white')

        ax.set_title("Return Period Curve - Fitted Distributions vs. Empirical Data", color=COLOR_TEXT, pad=14, fontweight="bold", fontsize=11)
        ax.set_xlabel("Return Period T (Years, log scale)", color=COLOR_TEXT, fontsize=9)
        ax.set_ylabel("Rainfall (mm)", color=COLOR_TEXT, fontsize=9)
        ax.tick_params(colors=COLOR_TEXT, labelsize=8)
        ax.set_xscale("log")
        ax.set_xticks(RETURN_PERIODS)
        ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.grid(True, which="major", linestyle="--", alpha=0.55, color="#CBD5E1")
        ax.grid(True, which="minor", linestyle=":", alpha=0.25, color="#CBD5E1")
        for spine in ax.spines.values(): spine.set_color("#CBD5E1")

        scatter = ax.scatter(T_emp, x_emp, color=COLOR_TEXT, edgecolor="white", linewidth=0.6,
                              label="Empirical Data", zorder=5, s=55, marker="o")
        colors = CHART_PALETTE
        self.inter_lines = {}

        canvas = FigureCanvasTkAgg(fig, master=main_graph)

        toggles_frame = ctk.CTkFrame(sidebar, fg_color="transparent")
        toggles_frame.pack(fill="both", expand=True, padx=6, pady=(6, 6))
        self.inter_vars = {}

        for (dist_name, y), color in zip(curves.items(), colors):
            line, = ax.plot(Tline, y, label=dist_name, color=color, linewidth=2.2, alpha=0.95)
            self.inter_lines[dist_name] = line

            row = ctk.CTkFrame(toggles_frame, fg_color="transparent")
            row.pack(fill="x", padx=6, pady=3)
            swatch = ctk.CTkFrame(row, fg_color=color, width=14, height=14, corner_radius=4)
            swatch.pack(side="left", padx=(2, 6)); swatch.pack_propagate(False)
            var = ctk.BooleanVar(value=True)
            chk = ctk.CTkCheckBox(
                row, text=dist_name, variable=var,
                command=lambda n=dist_name, v=var: self._toggle_line(n, v.get(), canvas),
                fg_color=color, hover_color=color, text_color=COLOR_TEXT, font=("Segoe UI", 11)
            )
            chk.pack(side="left")
            self.inter_vars[dist_name] = var

        btn_row = ctk.CTkFrame(sidebar, fg_color="transparent"); btn_row.pack(fill="x", padx=12, pady=(0, 14))
        ctk.CTkButton(btn_row, text="Show All", width=90, height=26, fg_color=COLOR_ACCENT, hover_color=COLOR_ACCENT_DARK,
                      command=lambda: self._toggle_all_lines(True, canvas)).pack(side="left", padx=(0, 6))
        ctk.CTkButton(btn_row, text="Hide All", width=90, height=26, fg_color="#94A3B8", hover_color="#64748B",
                      command=lambda: self._toggle_all_lines(False, canvas)).pack(side="left")

        legend = ax.legend(facecolor='white', edgecolor='#D9E1EC', labelcolor=COLOR_TEXT, fontsize=8, loc="upper left")
        legend.get_frame().set_alpha(0.92)

        cursor = mplcursors.cursor(ax.lines + [scatter], hover=True)
        @cursor.connect("add")
        def on_add(sel):
            sel.annotation.set_text(f"Return Period: {sel.target[0]:.1f} yrs\nRainfall: {sel.target[1]:.1f} mm")
            sel.annotation.get_bbox_patch().set(boxstyle="round,pad=0.6", fc="white", ec="#94A3B8", alpha=0.95)
            sel.annotation.set_color(COLOR_TEXT)

        fig.tight_layout()
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)
        self.current_interactive_canvas = canvas 
        
        toolbar = NavigationToolbar2Tk(canvas, main_graph)
        toolbar.update()
        toolbar.pack(side="bottom", fill="x")

    def _toggle_all_lines(self, visible, canvas):
        for name, var in getattr(self, "inter_vars", {}).items():
            var.set(visible)
            self._toggle_line(name, visible, canvas)

    def _toggle_line(self, name, is_visible, canvas):
        if name in self.inter_lines:
            self.inter_lines[name].set_visible(is_visible)
            canvas.draw()

    def _show_about(self):
        messagebox.showinfo(APP_TITLE, f"{APP_TITLE} v{APP_VERSION} - \"Aurora Studio\" Edition\n\nRainfall frequency analysis using 7 distributions:\nNormal, Log Normal, 3-Parameter Log Normal, Gumbel,\nPearson III, Log Pearson III, and GEV.\n\nData can be imported from a multi-sheet Excel workbook, typed\nmanually inside the app, or both combined.\n\nWhat's new in v4.0:\n - Redesigned interface with a colour-coded workflow sidebar\n - At-a-glance KPI dashboards on every page\n - Log-scale Return Period chart with Show/Hide All controls\n - Data Quality overview banner and clearer pass/fail badges\n\n{CREDENTIAL_TEXT}")

    def _on_close(self):
        if self.worker_thread and self.worker_thread.is_alive():
            if not messagebox.askyesno(APP_TITLE, "A process is still running. Are you sure you want to close the application?\n(The unfinished process will be force-stopped.)"):
                return
        self.destroy()

def main():
    import tkinter.simpledialog  # noqa: F401
    app = RainfallApp()
    app.mainloop()

if __name__ == "__main__":
    main()
