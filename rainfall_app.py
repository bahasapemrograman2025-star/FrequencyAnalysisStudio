# -*- coding: utf-8 -*-
"""
=======================================================================
 RAINFALL / ANNUAL DISCHARGE FREQUENCY ANALYSIS
 6 Distributions: Normal, LogNormal, 3-Param LogNormal, Gumbel,
 Pearson III, Log Pearson III, GEV (PWM / Hosking)

 v3.0 - "Data Studio" Edition
   - Dual data source: Import from Excel  OR  Manual entry in-app
     (multiple series, spreadsheet-like editor, exactly like working
     with multiple sheets in Excel - but without ever leaving the app)
   - Unified "Series Manager" so imported and manually typed series
     can be mixed, edited, renamed and processed together
   - Redesigned UI/UX: card-based layout, color-coded status, a
     dashboard-style workspace with clear step-by-step navigation
   - Same analytical engine & native Excel export/report as before
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

# When built with PyInstaller in --windowed mode (no console window),
# sys.stdout / sys.stderr are None. Some helper functions call print()
# as a harmless debug echo, so make that a no-op instead of crashing.
class _NullWriter:
    def write(self, *args, **kwargs):
        pass

    def flush(self):
        pass

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
#                DATA QUALITY TESTS (WITH FULL DETAIL TABLES)
# =======================================================================

def data_quality_tests(x, years=None):
    x = np.asarray(x, dtype=float)
    n = len(x)
    no = np.arange(1, n + 1)
    if years is None:
        years = no
    years = np.asarray(years)

    # ---------------- 1. OUTLIER TEST (Grubbs / WRC, alpha = 10%) ----------------
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

    df_outlier_detail = pd.DataFrame({
        "No": no, "Year": years, "Rainfall (mm)": x, "Log Rainfall": log_x,
    })
    outlier_summary = {
        "Mean": mean_log, "Std Dev": std_log, "Kn (Table)": kn,
        "Upper Limit Yh": yh, "Lower Limit Yl": yl,
        "Upper Limit Xh": xh, "Lower Limit Xl": xl,
        "Conclusion": outlier_status,
    }

    # ---------------- 2. TREND TEST (Spearman's Rho) ----------------
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
        "Rank Time (t)": rank_t, "Rank Rainfall (X)": rank_x,
        "d = Rt-Rx": d, "d^2": d2,
    })
    trend_summary = {
        "Sum d^2": np.sum(d2), "Rho (Spearman)": rho,
        "t Calculated": t_trend, "t Critical": t_trend_crit, "Conclusion": trend_status,
    }

    # ---------------- 3. VARIANCE HOMOGENEITY TEST (F-Test) ----------------
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

    df_group1 = pd.DataFrame({
        "No": no[idx1], "Year": years[idx1], "Rainfall (mm)": x1,
        "(Xi-X1bar)": dev1, "(Xi-X1bar)^2": dev1 ** 2,
    })
    df_group2 = pd.DataFrame({
        "No": no[idx2], "Year": years[idx2], "Rainfall (mm)": x2,
        "(Xi-X2bar)": dev2, "(Xi-X2bar)^2": dev2 ** 2,
    })
    var_summary = {
        "n1": len(x1), "n2": len(x2), "Mean Group 1": mean1, "Mean Group 2": mean2,
        "Variance Group 1 (S1^2)": v1, "Variance Group 2 (S2^2)": v2,
        "F Calculated": f_calc, "F Critical": f_crit, "Conclusion": var_status,
    }

    # ---------------- 4. MEAN HOMOGENEITY TEST (T-Test) ----------------
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

    # ---------------- 5. INDEPENDENCE TEST (Lag-1 Serial Correlation) ----------------
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
        "(Xi-Xbar)": dev[:-1], "(Xi+1-Xbar)": dev[1:],
        "(Xi-Xbar)(Xi+1-Xbar)": num_terms, "(Xi-Xbar)^2": dev[:-1] ** 2,
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

    b1_list, b2_list, b3_list = [], [], []
    emp_list, teo_list, sim_list = [], [], []

    for i in range(n):
        val = x_sort_desc[i]
        rank_asc = n - i

        term_b1 = val * (rank_asc - 1) / (n - 1) if n > 1 else 0
        term_b2 = val * (rank_asc - 1) * (rank_asc - 2) / ((n - 1) * (n - 2)) if n > 2 else 0
        term_b3 = val * (rank_asc - 1) * (rank_asc - 2) * (rank_asc - 3) / ((n - 1) * (n - 2) * (n - 3)) if n > 3 else 0

        b1_list.append(term_b1)
        b2_list.append(term_b2)
        b3_list.append(term_b3)

        emp = ((i + 1) / (n + 1)) * 100
        teo = (1.0 - gev_fit.cdf(val)) * 100
        sim = abs(emp - teo)

        emp_list.append(emp)
        teo_list.append(teo)
        sim_list.append(sim)

    df_pwm = pd.DataFrame({
        "X (Sorted)": x_sort_desc, "b1": b1_list, "b2": b2_list, "b3": b3_list,
        "Empirical (%)": emp_list, "Theoretical (%)": teo_list, "Deviation (%)": sim_list
    })

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
        self.name = name
        self.nparams = nparams
        self._cdf = cdf_func
        self._ppf = ppf_func
        self._pdf = pdf_func
        self.params = params or {}
        self.ok = ok
        self.note = note

    def cdf(self, x):
        try:
            v = self._cdf(np.asarray(x, dtype=float))
            return np.clip(v, 1e-12, 1 - 1e-12)
        except Exception:
            return np.full(np.shape(x), np.nan)

    def ppf(self, p):
        try:
            return self._ppf(np.asarray(p, dtype=float))
        except Exception:
            return np.full(np.shape(p), np.nan)

    def pdf(self, x):
        if self._pdf is None:
            return None
        try:
            return self._pdf(np.asarray(x, dtype=float))
        except Exception:
            return None


def _failed(name, note="data does not meet requirements / fit failed"):
    return FittedDist(name, nparams=np.nan, cdf_func=None, ppf_func=None, ok=False, note=note)


def _from_scipy(name, dist, params, nparams):
    frozen = dist(**params) if isinstance(params, dict) else dist(*params)
    return FittedDist(name, nparams, frozen.cdf, frozen.ppf, frozen.pdf, params=params)


def fit_normal_mom(x):
    m, s = np.mean(x), np.std(x, ddof=1)
    return _from_scipy("Normal", stats.norm, dict(loc=m, scale=s), 2)


def fit_lognormal_mom(x):
    x = np.asarray(x, dtype=float)
    if np.any(x <= 0):
        return _failed("LogNormal")
    y = np.log(x)
    mu, sigma = y.mean(), y.std(ddof=1)
    return _from_scipy("LogNormal", stats.lognorm, dict(s=sigma, loc=0, scale=np.exp(mu)), 2)


def fit_lognormal3_mom(x, name="3-Parameter Log Normal"):
    x = np.asarray(x, dtype=float)
    n = x.size
    mean, std = x.mean(), x.std(ddof=1)
    cs = (n / ((n - 1) * (n - 2))) * np.sum((x - mean) ** 3) / std ** 3 if n > 2 else np.nan
    if not np.isfinite(cs) or cs == 0:
        return _failed(name)

    cs_abs = abs(cs)
    half = cs_abs / 2.0
    disc = np.sqrt(half ** 2 + 1.0)
    cv = np.cbrt(half + disc) + np.cbrt(half - disc)
    if not np.isfinite(cv) or cv <= 0:
        return _failed(name)

    sigma2 = np.log(cv ** 2 + 1.0)
    sigma = np.sqrt(sigma2)
    skala = std / cv
    if skala <= 0 or not np.isfinite(skala):
        return _failed(name)
    mu = np.log(skala) - sigma2 / 2.0

    if cs >= 0:
        xo = mean - skala
        sign = 1.0
    else:
        xo = mean + skala
        sign = -1.0

    base = stats.lognorm(s=sigma, loc=0, scale=np.exp(mu))

    def cdf(v):
        v = np.asarray(v, dtype=float)
        return base.cdf(v - xo) if sign > 0 else base.sf(xo - v)

    def ppf(p):
        p = np.asarray(p, dtype=float)
        return xo + base.ppf(p) if sign > 0 else xo - base.ppf(1.0 - p)

    def pdf(v):
        v = np.asarray(v, dtype=float)
        return base.pdf(v - xo) if sign > 0 else base.pdf(xo - v)

    return FittedDist(name, 3, cdf, ppf, pdf, params=dict(xo=xo, mu=mu, sigma=sigma, sign=sign))


def fit_gumbel_max_mom(x):
    m, s = np.mean(x), np.std(x, ddof=1)
    alpha = s * np.sqrt(6) / np.pi
    xi = m - EULER_GAMMA * alpha
    return _from_scipy("Gumbel", stats.gumbel_r, dict(loc=xi, scale=alpha), 2)


def _pearson3_params_from_moments(mean, std, cs):
    if not np.isfinite(cs) or cs == 0:
        return None
    alpha = 4.0 / cs ** 2
    beta = 0.5 * std * cs
    xi = mean - alpha * beta
    return alpha, beta, xi


def _pearson3_dist(alpha, beta, xi):
    if beta > 0:
        frozen = stats.gamma(a=alpha, loc=xi, scale=beta)
        return frozen.cdf, frozen.ppf, frozen.pdf
    else:
        b = abs(beta)
        g = stats.gamma(a=alpha, scale=b)

        def cdf(v):
            return 1.0 - g.cdf(xi - v)

        def ppf(p):
            return xi - g.ppf(1.0 - p)

        def pdf(v):
            return g.pdf(xi - v)

        return cdf, ppf, pdf


def fit_pearson3_mom(x, name="Pearson III"):
    x = np.asarray(x, dtype=float)
    n = x.size
    mean, std = x.mean(), x.std(ddof=1)
    cs = (n / ((n - 1) * (n - 2))) * np.sum((x - mean) ** 3) / std ** 3 if n > 2 else np.nan
    p = _pearson3_params_from_moments(mean, std, cs)
    if p is None:
        return _failed(name)
    alpha, beta, xi = p
    cdf, ppf, pdf = _pearson3_dist(alpha, beta, xi)
    return FittedDist(name, 3, cdf, ppf, pdf, params=dict(alpha=alpha, beta=beta, xi=xi))


def fit_logpearson3_mom(x):
    x = np.asarray(x, dtype=float)
    if np.any(x <= 0):
        return _failed("Log Pearson III")
    y = np.log(x)
    fy = fit_pearson3_mom(y, name="Log Pearson III")
    if not fy.ok:
        return fy

    def cdf(v):
        return fy.cdf(np.log(v))

    def ppf(p):
        return np.exp(fy.ppf(p))

    return FittedDist("Log Pearson III", 3, cdf, ppf, params=fy.params)


def fit_gev_pwm(x, lm, name="GEV"):
    t3 = lm["t3"]
    if not np.isfinite(t3):
        return _failed(name)

    c_h = 2.0 / (3.0 + t3) - np.log(2) / np.log(3)
    k = 7.8590 * c_h + 2.9554 * c_h ** 2

    if abs(k) < 1e-8:
        alpha = lm["L2"] / np.log(2)
        xi = lm["L1"] - EULER_GAMMA * alpha
        return _from_scipy(name, stats.gumbel_r, dict(loc=xi, scale=alpha), 2)

    alpha = k * lm["L2"] / ((1 - 2 ** (-k)) * special.gamma(1 + k))
    xi = lm["L1"] - (alpha / k) * (1 - special.gamma(1 + k))
    if alpha <= 0:
        return _failed(name)

    def ppf(p):
        p = np.asarray(p, dtype=float)
        with np.errstate(invalid="ignore", divide="ignore"):
            return xi + (alpha / k) * (1.0 - (-np.log(p)) ** k)

    def cdf(v):
        v = np.asarray(v, dtype=float)
        with np.errstate(invalid="ignore", divide="ignore"):
            arg = 1.0 - k * (v - xi) / alpha
            arg = np.where(arg > 0, arg, np.nan)
            return np.exp(-(arg ** (1.0 / k)))

    def pdf(v):
        v = np.asarray(v, dtype=float)
        with np.errstate(invalid="ignore", divide="ignore"):
            arg = 1.0 - k * (v - xi) / alpha
            arg = np.where(arg > 0, arg, np.nan)
            return (1.0 / alpha) * (arg ** (1.0 / k - 1.0)) * np.exp(-(arg ** (1.0 / k)))

    return FittedDist(name, 3, cdf, ppf, pdf, params=dict(alpha=alpha, k=k, xi=xi))


def build_all_fits(x, lm):
    fits = {}
    fits["Normal"] = fit_normal_mom(x)
    fits["LogNormal"] = fit_lognormal_mom(x)
    fits["3-Parameter Log Normal"] = fit_lognormal3_mom(x)
    fits["Gumbel"] = fit_gumbel_max_mom(x)
    fits["Pearson III"] = fit_pearson3_mom(x)
    fits["Log Pearson III"] = fit_logpearson3_mom(x)
    fits["GEV"] = fit_gev_pwm(x, lm, name="GEV")
    return fits


# =======================================================================
#             GOODNESS OF FIT TESTS & QUANTILES
# =======================================================================

ALPHAS = [0.01, 0.05, 0.10]
RATIO_LOW, RATIO_HIGH = 1.7, 3.2


def n_classes(n):
    k = int(round(1 + 3.322 * np.log10(n)))
    return max(k, 4)


def chi_square_test(x, fitted, k=None):
    x = np.asarray(x, dtype=float)
    n = x.size
    if k is None:
        k = n_classes(n)
    if not fitted.ok:
        return dict(statistic=np.nan, dof=np.nan, pvalue=np.nan, k=k)

    probs = np.linspace(0, 1, k + 1)
    edges = fitted.ppf(probs[1:-1])
    edges = np.concatenate([[-np.inf], np.sort(edges), [np.inf]])
    observed, _ = np.histogram(x, bins=edges)
    expected = n / k

    nparams = fitted.nparams
    dof = k - 1 - nparams
    if dof <= 0 or not np.isfinite(nparams):
        return dict(statistic=np.nan, dof=np.nan, pvalue=np.nan, k=k)

    stat = np.sum((observed - expected) ** 2 / expected)
    pvalue = stats.chi2.sf(stat, dof)
    return dict(statistic=stat, dof=dof, pvalue=pvalue, k=k)


def ks_test(x, fitted):
    if not fitted.ok:
        return dict(dmax=np.nan, pvalue=np.nan)
    x = np.asarray(x, dtype=float)
    try:
        res = stats.kstest(x, fitted.cdf)
        return dict(dmax=res.statistic, pvalue=res.pvalue)
    except Exception:
        return dict(dmax=np.nan, pvalue=np.nan)


def decision(pvalue, alpha):
    if not np.isfinite(pvalue):
        return ""
    return "ACCEPT" if pvalue >= alpha else "REJECT"


def quantiles_for_periods(fitted, periods, x=None):
    if not fitted.ok:
        return {T: np.nan for T in periods}
    lo, hi = -np.inf, np.inf
    if x is not None:
        x = np.asarray(x, dtype=float)
        mean, std = x.mean(), x.std(ddof=1)
        lo, hi = mean - 20 * std, mean + 200 * std
    out = {}
    for T in periods:
        p = 1.0 - 1.0 / T
        try:
            val = float(fitted.ppf(np.array(p)))
            if not np.isfinite(val) or val < lo or val > hi:
                val = np.nan
        except Exception:
            val = np.nan
        out[T] = val
    return out


def empirical_plotting_position(x):
    x_sorted = np.sort(np.asarray(x, dtype=float))[::-1]
    n = x_sorted.size
    m = np.arange(1, n + 1)
    P = m / (n + 1.0)
    T = 1.0 / P
    return T, x_sorted


def build_return_period_curve_data(x, fits_dict, periods, n_curve_points=100):
    x = np.asarray(x, dtype=float)
    T_emp, x_emp = empirical_plotting_position(x)

    max_T = max(periods)
    Tline = np.linspace(1.01, max_T, n_curve_points)
    p_line = 1.0 - 1.0 / Tline

    mean, std = np.mean(x), np.std(x, ddof=1)
    ylo, yhi = mean - 20 * std, mean + 200 * std

    curves = {}
    for name, f in fits_dict.items():
        if f.ok:
            try:
                y = f.ppf(p_line)
                y = np.where((y >= ylo) & (y <= yhi), y, np.nan)
            except Exception:
                y = np.full_like(Tline, np.nan)
        else:
            y = np.full_like(Tline, np.nan)
        curves[name] = y

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
    """Read one sheet of an Excel workbook into a clean (Year, Rainfall) DataFrame."""
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
    """Read every sheet of an Excel workbook into a dict {sheet_name: DataFrame}.
    Sheets that cannot be parsed are skipped and reported separately."""
    xls = pd.ExcelFile(path)
    series = {}
    errors = {}
    for sheet in xls.sheet_names:
        try:
            df = read_input_dataframe(path, sheet, col_tahun, col_data)
            if len(df) == 0:
                errors[sheet] = "No valid numeric rows found."
                continue
            series[sheet] = df
        except Exception as e:
            errors[sheet] = str(e)
    return series, errors


def empty_series_df():
    return pd.DataFrame({COL_TAHUN: pd.Series(dtype="Int64"), COL_DATA: pd.Series(dtype="float")})


def validate_series_df(df):
    """Return (is_ready, message) describing whether a series is ready to be analyzed."""
    n = len(df)
    if n == 0:
        return False, "No data rows yet"
    if n < 5:
        return False, f"Only {n} rows (minimum 5 recommended)"
    return True, "Ready"


def parse_bulk_text(raw_text):
    """Parse pasted / typed bulk text into a list of (year, value) pairs.
    Accepts comma, tab, semicolon or whitespace separated values, one pair
    per line. Lines that cannot be parsed are skipped and returned as
    'bad_lines' so the caller can warn the user."""
    rows = []
    bad_lines = []
    for raw_line in raw_text.splitlines():
        line = raw_line.strip()
        if not line:
            continue
        # normalize separators to a single comma
        parts = re.split(r"[,\t;]+|\s+", line)
        parts = [p for p in parts if p != ""]
        if len(parts) < 2:
            bad_lines.append(raw_line)
            continue
        year_raw, val_raw = parts[0], parts[1]
        try:
            year_raw_clean = year_raw.replace(",", ".")
            val_raw_clean = val_raw.replace(",", ".")
            year_val = float(year_raw_clean)
            data_val = float(val_raw_clean)
            year_val = int(year_val) if float(year_val).is_integer() else year_val
            rows.append((year_val, data_val))
        except Exception:
            bad_lines.append(raw_line)
    return rows, bad_lines


# =======================================================================
#                        MAIN ANALYSIS PIPELINE
# =======================================================================

def run_analysis(df, col_tahun=COL_TAHUN, col_data=COL_DATA):
    x = df[col_data].to_numpy(dtype=float)
    years = df[col_tahun].to_numpy()
    n = x.size
    if n < 5:
        raise ValueError("Insufficient data (minimum suggested >= 5 values).")

    stats_dasar = basic_moments(x)
    lm = lmoments(x)

    s = pd.Series(x)
    pandas_stats = dict(
        n=len(x), min=s.min(), max=s.max(), mean=s.mean(),
        std=s.std(ddof=1), ck=s.kurt(), cs=s.skew()
    )

    fits = build_all_fits(x, lm)
    k = n_classes(n)

    dq_tests = data_quality_tests(x, years=years)
    df_pwm, df_ks_gev = gev_pwm_ks_tables(x, fits["GEV"])

    chi_rows, ks_rows, quant_rows = [], [], []

    for name, f in fits.items():
        chi = chi_square_test(x, f, k=k)
        chi_rows.append(dict(
            Distribution=name,
            **{f"a={int(a*100)}%": decision(chi["pvalue"], a) for a in ALPHAS},
            **{"Attained alpha": chi["pvalue"], "Chi-Square Statistic": chi["statistic"], "dof": chi["dof"]}
        ))

        ks = ks_test(x, f)
        ks_rows.append(dict(
            Distribution=name,
            **{f"a={int(a*100)}%": decision(ks["pvalue"], a) for a in ALPHAS},
            **{"Attained alpha": ks["pvalue"], "Dmax": ks["dmax"]}
        ))

        q = quantiles_for_periods(f, RETURN_PERIODS, x=x)
        quant_rows.append(dict(Distribution=name, **{f"T={T}": q[T] for T in RETURN_PERIODS}))

    chi_df = pd.DataFrame(chi_rows)
    ks_df = pd.DataFrame(ks_rows)
    ks_df["Rank"] = ks_df["Dmax"].rank(method="min", na_option="bottom").astype("Int64")
    ks_df = ks_df.sort_values("Distribution", key=lambda s: s.map({n: i for i, n in enumerate(fits)}))

    quant_df = pd.DataFrame(quant_rows)
    num_cols = [c for c in quant_df.columns if c != "Distribution"]
    quant_df[num_cols] = quant_df[num_cols].round(ROUND_DECIMALS)

    r100 = pd.to_numeric(quant_df["T=100"], errors="coerce")
    r2 = pd.to_numeric(quant_df["T=2"], errors="coerce")
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = r100 / r2
    ratio = ratio.where(r2 != 0, np.nan)
    quant_df["R100/R2"] = ratio.round(ROUND_DECIMALS)
    quant_df["R100/R2 OK?"] = np.where(
        quant_df["R100/R2"].between(RATIO_LOW, RATIO_HIGH), "Yes", "No"
    )
    chi_df["Attained alpha"] = chi_df["Attained alpha"].round(6)
    ks_df["Attained alpha"] = ks_df["Attained alpha"].round(6)
    ks_df["Dmax"] = ks_df["Dmax"].round(5)

    try:
        ks_crit_05 = stats.kstwo.ppf(0.95, n)
    except AttributeError:
        ks_crit_05 = stats.ksone.ppf(0.975, n)
    ks_df["Pass KS 5%"] = np.where(
        ks_df["Dmax"].notna() & (ks_df["Dmax"] < ks_crit_05), "Pass", "Fail"
    )

    return dict(x=x, stats_dasar=stats_dasar, pandas_stats=pandas_stats, lmoments=lm, fits=fits,
                chi_df=chi_df, ks_df=ks_df, quant_df=quant_df, k_classes=k,
                dq_tests=dq_tests, df_pwm=df_pwm, df_ks_gev=df_ks_gev)


def _stamp_all_sheets(wb):
    fmt_credential = wb.add_format({"font_size": 9, "italic": True, "bold": True, "font_color": "#595959"})
    for ws in wb.worksheets():
        ws.write(0, 0, CREDENTIAL_TEXT, fmt_credential)


def _write_df_table(ws, df, start_row, start_col, fmt_hdr, fmt_num, fmt_text, n_decimals=4):
    for c, col_name in enumerate(df.columns):
        ws.write(start_row, start_col + c, col_name, fmt_hdr)
    for r in range(len(df)):
        for c, col_name in enumerate(df.columns):
            val = df.iat[r, c]
            if isinstance(val, (int, np.integer)):
                ws.write_number(start_row + r + 1, start_col + c, int(val), fmt_text)
            elif isinstance(val, (float, np.floating)):
                ws.write_number(start_row + r + 1, start_col + c, float(val), fmt_num)
            else:
                ws.write(start_row + r + 1, start_col + c, str(val), fmt_text)
    return start_row + len(df) + 2


def _write_summary_kv(ws, summary_dict, start_row, start_col, fmt_label, fmt_value, fmt_value_text):
    r = start_row
    for label, val in summary_dict.items():
        ws.write(r, start_col, label, fmt_label)
        if isinstance(val, (int, np.integer)):
            ws.write_number(r, start_col + 1, int(val), fmt_value)
        elif isinstance(val, (float, np.floating)):
            ws.write_number(r, start_col + 1, float(val), fmt_value)
        else:
            ws.write(r, start_col + 1, str(val), fmt_value_text)
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
        for c, name in enumerate(df_input.columns):
            ws.write(1, c, name, fmt_hdr)

        # --- Data Quality Test SHEET ---
        ws = wb.add_worksheet("Data Quality Test")
        writer.sheets["Data Quality Test"] = ws
        dq = result["dq_tests"]

        row = 1
        ws.write(row, 0, "RAINFALL DATA QUALITY TEST (DATA CONSISTENCY TEST)", fmt_title)
        row += 2

        ws.write(row, 0, "1. RAINFALL DATA OUTLIER TEST (Grubbs / WRC, \u03b1 = 10%)", fmt_title2)
        row += 1
        row = _write_df_table(ws, dq["outlier"]["detail"], row, 0, fmt_hdr, fmt_num_new, fmt_border_center)
        row = _write_summary_kv(ws, dq["outlier"]["summary"], row, 0, fmt_label, fmt_num_new, fmt_value_text)
        ws.write(row - 2, 1, dq["outlier"]["summary"]["Conclusion"],
                 fmt_accept if "No Outliers" in dq["outlier"]["summary"]["Conclusion"] else fmt_reject)
        row += 1

        ws.write(row, 0, "2. RAINFALL DATA NO-TREND TEST (Spearman's Rho)", fmt_title2)
        row += 1
        row = _write_df_table(ws, dq["trend"]["detail"], row, 0, fmt_hdr, fmt_num_new, fmt_border_center)
        row = _write_summary_kv(ws, dq["trend"]["summary"], row, 0, fmt_label, fmt_num_new, fmt_value_text)
        ws.write(row - 2, 1, dq["trend"]["summary"]["Conclusion"],
                 fmt_accept if "No Trend" in dq["trend"]["summary"]["Conclusion"] else fmt_reject)
        row += 1

        ws.write(row, 0, "3. RAINFALL DATA VARIANCE STATIONARITY TEST (F-Test)", fmt_title2)
        row += 1
        ws.write(row, 0, "Group 1 (First Half)", fmt_border_left)
        row += 1
        row = _write_df_table(ws, dq["var"]["detail1"], row, 0, fmt_hdr, fmt_num_new, fmt_border_center)
        ws.write(row, 0, "Group 2 (Second Half)", fmt_border_left)
        row += 1
        row = _write_df_table(ws, dq["var"]["detail2"], row, 0, fmt_hdr, fmt_num_new, fmt_border_center)
        row = _write_summary_kv(ws, dq["var"]["summary"], row, 0, fmt_label, fmt_num_new, fmt_value_text)
        ws.write(row - 2, 1, dq["var"]["summary"]["Conclusion"],
                 fmt_accept if "Homogeneous" in dq["var"]["summary"]["Conclusion"] and "Not" not in dq["var"]["summary"]["Conclusion"] else fmt_reject)
        row += 1

        ws.write(row, 0, "4. RAINFALL DATA MEAN STATIONARITY TEST (T-Test)", fmt_title2)
        row += 1
        ws.write(row, 0, "Group 1 (First Half)", fmt_border_left)
        row += 1
        row = _write_df_table(ws, dq["mean"]["detail1"], row, 0, fmt_hdr, fmt_num_new, fmt_border_center)
        ws.write(row, 0, "Group 2 (Second Half)", fmt_border_left)
        row += 1
        row = _write_df_table(ws, dq["mean"]["detail2"], row, 0, fmt_hdr, fmt_num_new, fmt_border_center)
        row = _write_summary_kv(ws, dq["mean"]["summary"], row, 0, fmt_label, fmt_num_new, fmt_value_text)
        ws.write(row - 2, 1, dq["mean"]["summary"]["Conclusion"],
                 fmt_accept if "Homogeneous" in dq["mean"]["summary"]["Conclusion"] and "Not" not in dq["mean"]["summary"]["Conclusion"] else fmt_reject)
        row += 1

        ws.write(row, 0, "5. RAINFALL DATA INDEPENDENCE TEST (Lag-1 Serial Correlation)", fmt_title2)
        row += 1
        row = _write_df_table(ws, dq["indep"]["detail"], row, 0, fmt_hdr, fmt_num_new, fmt_border_center)
        row = _write_summary_kv(ws, dq["indep"]["summary"], row, 0, fmt_label, fmt_num_new, fmt_value_text)
        ws.write(row - 2, 1, dq["indep"]["summary"]["Conclusion"],
                 fmt_accept if "Independent" in dq["indep"]["summary"]["Conclusion"] and "Not" not in dq["indep"]["summary"]["Conclusion"] else fmt_reject)
        row += 1

        row += 1
        ws.write(row, 0, "RAINFALL DATA QUALITY TEST - QUICK RECAP", fmt_title)
        row += 2

        ws.write(row, 0, "1. RAINFALL DATA OUTLIER TEST (Grubbs / WRC)", fmt_title2)
        row += 1
        recap_outlier = [
            ("Upper Threshold (Xh)", dq["outlier"]["summary"]["Upper Limit Xh"]),
            ("Lower Threshold (Xl)", dq["outlier"]["summary"]["Lower Limit Xl"]),
            ("Kn Value", dq["outlier"]["summary"]["Kn (Table)"]),
        ]
        for label, val in recap_outlier:
            ws.write(row, 0, label, fmt_border_left)
            ws.write_number(row, 1, float(val), fmt_num_new)
            row += 1
        ws.write(row, 0, "Outlier Status", fmt_label)
        status = "Safe (No Outliers)" if "No Outliers" in dq["outlier"]["summary"]["Conclusion"] else "Outliers Detected"
        ws.write(row, 1, status, fmt_accept if "Safe" in status else fmt_reject)
        row += 2

        ws.write(row, 0, "2. RAINFALL DATA NO-TREND TEST (Spearman's Rho)", fmt_title2)
        row += 1
        recap_trend = [
            ("Correlation Value (Rho)", dq["trend"]["summary"]["Rho (Spearman)"]),
            ("Calculated T", dq["trend"]["summary"]["t Calculated"]),
            ("Critical T", dq["trend"]["summary"]["t Critical"]),
        ]
        for label, val in recap_trend:
            ws.write(row, 0, label, fmt_border_left)
            ws.write_number(row, 1, float(val), fmt_num_new)
            row += 1
        ws.write(row, 0, "Trend Status", fmt_label)
        ws.write(row, 1, dq["trend"]["summary"]["Conclusion"],
                 fmt_accept if "No Trend" in dq["trend"]["summary"]["Conclusion"] else fmt_reject)
        row += 2

        ws.write(row, 0, "3. RAINFALL DATA VARIANCE / F-TEST", fmt_title2)
        row += 1
        recap_var = [
            ("Calculated F", dq["var"]["summary"]["F Calculated"]),
            ("Critical F", dq["var"]["summary"]["F Critical"]),
        ]
        for label, val in recap_var:
            ws.write(row, 0, label, fmt_border_left)
            ws.write_number(row, 1, float(val), fmt_num_new)
            row += 1
        ws.write(row, 0, "Variance Status", fmt_label)
        var_concl = dq["var"]["summary"]["Conclusion"]
        var_short = "Homogeneous Variance (Pass)" if ("Homogeneous" in var_concl and "Not" not in var_concl) else "Non-Homogeneous Variance (Fail)"
        ws.write(row, 1, var_short, fmt_accept if "Pass" in var_short else fmt_reject)
        row += 2

        ws.write(row, 0, "4. RAINFALL DATA MEAN / T-TEST", fmt_title2)
        row += 1
        recap_mean = [
            ("Calculated T", dq["mean"]["summary"]["T Calculated"]),
            ("Critical T", dq["mean"]["summary"]["T Critical"]),
        ]
        for label, val in recap_mean:
            ws.write(row, 0, label, fmt_border_left)
            ws.write_number(row, 1, float(val), fmt_num_new)
            row += 1
        ws.write(row, 0, "Mean Status", fmt_label)
        mean_concl = dq["mean"]["summary"]["Conclusion"]
        mean_short = "Homogeneous Mean (Pass)" if ("Homogeneous" in mean_concl and "Not" not in mean_concl) else "Non-Homogeneous Mean (Fail)"
        ws.write(row, 1, mean_short, fmt_accept if "Pass" in mean_short else fmt_reject)
        row += 2

        ws.write(row, 0, "5. RAINFALL DATA INDEPENDENCE TEST (Lag-1 Serial Correlation)", fmt_title2)
        row += 1
        recap_indep = [
            ("Correlation (r1)", dq["indep"]["summary"]["r1 (Correlation)"]),
            ("Lower Limit", dq["indep"]["summary"]["Lower Limit (Ll)"]),
            ("Upper Limit", dq["indep"]["summary"]["Upper Limit (Ul)"]),
        ]
        for label, val in recap_indep:
            ws.write(row, 0, label, fmt_border_left)
            ws.write_number(row, 1, float(val), fmt_num_new)
            row += 1
        ws.write(row, 0, "Independence Status", fmt_label)
        ws.write(row, 1, dq["indep"]["summary"]["Conclusion"],
                 fmt_accept if ("Independent" in dq["indep"]["summary"]["Conclusion"] and "Not" not in dq["indep"]["summary"]["Conclusion"]) else fmt_reject)
        row += 1

        ws.set_column(0, 0, 26)
        ws.set_column(1, 9, 16)

        # --- Basic Statistics SHEET ---
        sd = result["stats_dasar"]
        lm = result["lmoments"]
        ps = result["pandas_stats"]
        ws = wb.add_worksheet("Basic Statistics")
        writer.sheets["Basic Statistics"] = ws

        ws.write(1, 0, "Data Statistics", fmt_title)
        rows_new = [
            ("Data Name", "Rainfall (mm/day)"), ("Data Count", ps["n"]),
            ("Minimum", ps["min"]), ("Maximum", ps["max"]), ("Mean", ps["mean"]),
            ("Standard Deviation", ps["std"]), ("Kurtosis", ps["ck"]), ("Skewness", ps["cs"])
        ]
        for i, (label, val) in enumerate(rows_new, start=2):
            ws.write(i, 0, label, fmt_border_left)
            if isinstance(val, str):
                ws.write(i, 1, val, fmt_border_left)
            else:
                ws.write_number(i, 1, float(val), fmt_num_new)

        ws.write(12, 0, "BASIC STATISTICS (Method of Moments)", fmt_title)
        rows_old = [("Data Count (n)", sd["n"]), ("Mean", sd["mean"]),
                    ("Standard Deviation (Std Dev)", sd["std"]), ("Coef. of Variation (Cv)", sd["cv"]),
                    ("Coef. of Skewness (Cs)", sd["cs"]), ("Coef. of Kurtosis (Ck)", sd["ck"])]
        for i, (label, val) in enumerate(rows_old, start=14):
            ws.write(i, 0, label, fmt_border)
            ws.write_number(i, 1, float(val), fmt_num)

        ws.write(22, 0, "L-MOMENTS", fmt_title)
        rows_lm = [("L1", lm["L1"]), ("L2", lm["L2"]), ("L3", lm["L3"]), ("L4", lm["L4"]),
                   ("t3 (L-Skewness)", lm["t3"]), ("t4 (L-Kurtosis)", lm["t4"])]
        for i, (label, val) in enumerate(rows_lm, start=24):
            ws.write(i, 0, label, fmt_border)
            ws.write_number(i, 1, float(val), fmt_num)

        ws.set_column(0, 0, 26)
        ws.set_column(1, 1, 20)

        req_dist = ["Normal", "LogNormal", "3-Parameter Log Normal", "Gumbel", "Pearson III", "Log Pearson III", "GEV"]
        disp_dist = ["Normal", "Log Normal", "3-Parameter Log Normal", "Gumbel", "Pearson III", "Log Pearson III", "GEV"]

        # --- Chi-Square Test SHEET ---
        chi_df = result["chi_df"].drop(columns=["Chi-Square Statistic", "dof"]).rename(columns={"Attained alpha": "Attained a"})
        chi_full = result["chi_df"]
        chi_df.to_excel(writer, sheet_name="Chi-Square Test", index=False, startrow=3)
        ws = writer.sheets["Chi-Square Test"]

        ws.write(1, 0, f"ORIGINAL CHI-SQUARE TEST (number of classes k = {result['k_classes']})", fmt_title)
        for c, name in enumerate(chi_df.columns):
            ws.write(3, c, name, fmt_hdr)
        for r in range(len(chi_df)):
            for c in [1, 2, 3]:
                val = chi_df.iat[r, c]
                fmt = fmt_accept if val == "ACCEPT" else (fmt_reject if val == "REJECT" else fmt_border)
                ws.write(r + 4, c, val, fmt)
            ws.write_number(r + 4, 4, 0 if pd.isna(chi_df.iat[r, 4]) else float(chi_df.iat[r, 4]), fmt_num) if not pd.isna(chi_df.iat[r, 4]) else ws.write_blank(r + 4, 4, None, fmt_border)
            pearson_val = chi_full["Chi-Square Statistic"].iat[r]
            if pd.isna(pearson_val):
                ws.write_blank(r + 4, 5, None, fmt_border)
            else:
                ws.write_number(r + 4, 5, float(pearson_val), fmt_num)
            ws.write(r + 4, 0, chi_df.iat[r, 0], fmt_border)
        ws.write(3, 5, "Chi-Square Statistic", fmt_hdr)

        row_offset = len(chi_df) + 7
        for alpha in ALPHAS:
            pct = int(round(alpha * 100))
            ws.write(row_offset, 0, f"Goodness of fit test (Chi-Square) at \u03b1 = {pct}%", fmt_title)
            ws.write(row_offset + 1, 0, "Distribution", fmt_hdr)
            for c, dname in enumerate(disp_dist, start=1):
                ws.write(row_offset + 1, c, dname, fmt_hdr)
            ws.write(row_offset + 2, 0, "Maximum Delta", fmt_border_left)
            ws.write(row_offset + 3, 0, "Critical Delta", fmt_border_left)
            ws.write(row_offset + 4, 0, "Test Result", fmt_border_left)

            for c, r_name in enumerate(req_dist, start=1):
                row2 = chi_full[chi_full["Distribution"] == r_name].iloc[0]
                stat = row2["Chi-Square Statistic"]
                dof = row2["dof"]
                if pd.isna(stat) or pd.isna(dof) or dof <= 0:
                    ws.write(row_offset + 2, c, "-", fmt_border_center)
                    ws.write(row_offset + 3, c, "-", fmt_border_center)
                    ws.write(row_offset + 4, c, "Fail", fmt_reject)
                else:
                    crit = stats.chi2.ppf(1 - alpha, dof)
                    hasil = "Pass" if stat < crit else "Fail"
                    ws.write_number(row_offset + 2, c, stat, fmt_num_new)
                    ws.write_number(row_offset + 3, c, crit, fmt_num_new)
                    ws.write(row_offset + 4, c, hasil, fmt_accept if hasil == "Pass" else fmt_reject)

            row_offset += 7

        ws.set_column(0, 0, 30)
        ws.set_column(1, 7, 16)

        # --- Kolmogorov-Smirnov Test SHEET ---
        ks_df = result["ks_df"].rename(columns={"Attained alpha": "Attained a"})
        ks_df.to_excel(writer, sheet_name="Kolmogorov-Smirnov Test", index=False, startrow=3)
        ws = writer.sheets["Kolmogorov-Smirnov Test"]

        ws.write(1, 0, "ORIGINAL KOLMOGOROV-SMIRNOV TEST", fmt_title)
        for c, name in enumerate(ks_df.columns):
            ws.write(3, c, name, fmt_hdr)
        for r in range(len(ks_df)):
            ws.write(r + 4, 0, ks_df.iat[r, 0], fmt_border)
            for c in [1, 2, 3]:
                val = ks_df.iat[r, c]
                fmt = fmt_accept if val == "ACCEPT" else (fmt_reject if val == "REJECT" else fmt_border)
                ws.write(r + 4, c, val, fmt)
            for c in [4, 5]:
                val = ks_df.iat[r, c]
                if pd.isna(val):
                    ws.write_blank(r + 4, c, None, fmt_border)
                else:
                    ws.write_number(r + 4, c, float(val), fmt_num)
            rank_val = ks_df.iat[r, 6]
            ws.write(r + 4, 6, "" if pd.isna(rank_val) else int(rank_val), fmt_border)

        row_offset = len(ks_df) + 7
        n_data = ps["n"]
        for alpha in ALPHAS:
            pct = int(round(alpha * 100))
            ws.write(row_offset, 0, f"Goodness of fit test (Smirnov-Kolmogorov) at \u03b1 = {pct}%", fmt_title)
            ws.write(row_offset + 1, 0, "Distribution", fmt_hdr)
            for c, dname in enumerate(disp_dist, start=1):
                ws.write(row_offset + 1, c, dname, fmt_hdr)
            ws.write(row_offset + 2, 0, "Maximum Delta", fmt_border_left)
            ws.write(row_offset + 3, 0, "Critical Delta", fmt_border_left)
            ws.write(row_offset + 4, 0, "Test Result", fmt_border_left)

            try:
                ks_crit = stats.kstwo.ppf(1 - alpha, n_data)
            except AttributeError:
                ks_crit = stats.ksone.ppf(1 - alpha / 2, n_data)

            for c, r_name in enumerate(req_dist, start=1):
                row2 = ks_df[ks_df["Distribution"] == r_name].iloc[0]
                dmax = row2["Dmax"]
                if pd.isna(dmax):
                    ws.write(row_offset + 2, c, "-", fmt_border_center)
                    ws.write(row_offset + 3, c, ks_crit, fmt_num_new)
                    ws.write(row_offset + 4, c, "Fail", fmt_reject)
                else:
                    hasil = "Pass" if dmax < ks_crit else "Fail"
                    ws.write_number(row_offset + 2, c, dmax, fmt_num_new)
                    ws.write_number(row_offset + 3, c, ks_crit, fmt_num_new)
                    ws.write(row_offset + 4, c, hasil, fmt_accept if hasil == "Pass" else fmt_reject)

            row_offset += 7

        ws.set_column(0, 0, 30)
        ws.set_column(1, 7, 16)

        # --- GEV - PWM & KS SHEET ---
        ws = wb.add_worksheet("GEV - PWM & KS")
        writer.sheets["GEV - PWM & KS"] = ws
        df_pwm = result["df_pwm"]
        df_ks_gev = result["df_ks_gev"]

        ws.write(1, 0, "Probability in % (GEV Calculation using PWM Method)", fmt_title)
        for c, col_name in enumerate(df_pwm.columns):
            ws.write(2, c, col_name, fmt_hdr)
        for r in range(len(df_pwm)):
            ws.write(r + 3, 0, df_pwm.iat[r, 0], fmt_border_center)
            for c in range(1, len(df_pwm.columns)):
                ws.write_number(r + 3, c, df_pwm.iat[r, c], fmt_num)
        ws.set_column(0, 6, 15)

        row_ks_offset = len(df_pwm) + 7
        ws.write(row_ks_offset, 0, "Kolmogorov-Smirnov GEV", fmt_title)
        for c, col_name in enumerate(df_ks_gev.columns):
            ws.write(row_ks_offset + 2, c, col_name, fmt_hdr)
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
        for c, name in enumerate(quant_df.columns):
            ws.write(3, c, name, fmt_hdr)

        for r in range(len(quant_df)):
            for c, col_name in enumerate(quant_df.columns):
                val = quant_df.iat[r, c]
                if c == 0:
                    ws.write(r + 4, c, val, fmt_border)
                elif col_name == "R100/R2 OK?":
                    ws.write(r + 4, c, val, fmt_accept if val == "Yes" else fmt_reject)
                elif pd.isna(val):
                    ws.write_blank(r + 4, c, None, fmt_border)
                else:
                    ws.write_number(r + 4, c, float(val), fmt_num)

        ws.set_column(0, 0, 25)
        ws.set_column(1, len(quant_df.columns) - 1, 12)

        helper_row = len(quant_df) + 6
        ws.write(helper_row, 0, "T (chart helper)", fmt_border)
        for c, T in enumerate(RETURN_PERIODS, start=1):
            ws.write_number(helper_row, c, float(T))

        chart = wb.add_chart({"type": "scatter", "subtype": "straight_with_markers"})
        for r in range(len(quant_df)):
            chart.add_series({
                "name": [ws.get_name(), r + 3, 0],
                "categories": [ws.get_name(), helper_row, 1, helper_row, n_periods],
                "values": [ws.get_name(), r + 3, 1, r + 3, n_periods],
                "marker": {"type": "circle", "size": 4},
                "line": {"width": 1},
            })
        chart.set_title({"name": "Return Period Curve - All Methods"})
        chart.set_x_axis({"name": "Return Period T (years)"})
        chart.set_y_axis({"name": "Value"})
        chart.set_size({"width": 900, "height": 500})
        chart.set_legend({"font": {"size": 8}})
        ws.insert_chart(3, len(quant_df.columns) + 1, chart)

        # --- Graph Data (hidden helper sheet) ---
        ws_data = wb.add_worksheet("Graph Data")
        writer.sheets["Graph Data"] = ws_data

        T_emp, x_emp, Tline, curves = build_return_period_curve_data(
            result["x"], result["fits"], RETURN_PERIODS, n_curve_points=100
        )
        n_emp = len(T_emp)
        n_curve = len(Tline)

        ws_data.write(0, 0, "T (Empirical)", fmt_hdr)
        ws_data.write(0, 1, "X (Empirical)", fmt_hdr)
        for i in range(n_emp):
            ws_data.write_number(i + 1, 0, float(T_emp[i]))
            ws_data.write_number(i + 1, 1, float(x_emp[i]))

        curve_col = {}
        col = 3
        names_ordered = list(result["fits"].keys())
        for name in names_ordered:
            ws_data.write(0, col, f"T_{name}")
            ws_data.write(0, col + 1, f"X_{name}")
            y = curves[name]
            for i in range(n_curve):
                ws_data.write_number(i + 1, col, float(Tline[i]))
                val = y[i]
                if np.isfinite(val):
                    ws_data.write_number(i + 1, col + 1, float(val))
                else:
                    ws_data.write_blank(i + 1, col + 1, None)
            curve_col[name] = col
            col += 3

        ws_data.hide()

        # --- Graph SHEET (native Excel charts) ---
        ws = wb.add_worksheet("Graph")
        writer.sheets["Graph"] = ws
        ws.write(1, 0, "Frequency Curve - All Methods (native Excel charts, linear T axis)", fmt_title)

        ncols_grid = 3
        chart_width = 480
        chart_height = 300
        row_start = 3
        row_step = 16
        col_step = 9

        data_sheet_name = ws_data.get_name()

        for idx, name in enumerate(names_ordered):
            r_idx, c_idx = divmod(idx, ncols_grid)
            col = curve_col[name]

            chart = wb.add_chart({"type": "scatter", "subtype": "straight_with_markers"})
            chart.add_series({
                "name": "Empirical Data",
                "categories": [data_sheet_name, 1, 0, n_emp, 0],
                "values": [data_sheet_name, 1, 1, n_emp, 1],
                "marker": {"type": "circle", "size": 5, "fill": {"color": "black"}, "border": {"color": "black"}},
                "line": {"none": True},
            })
            chart.add_series({
                "name": f"{name} Fit",
                "categories": [data_sheet_name, 1, col, n_curve, col],
                "values": [data_sheet_name, 1, col + 1, n_curve, col + 1],
                "marker": {"type": "none"},
                "line": {"color": "red", "width": 1.5},
            })
            chart.set_title({"name": name})
            chart.set_x_axis({"name": "Return Period T (years)"})
            chart.set_y_axis({"name": "Value"})
            chart.set_size({"width": chart_width, "height": chart_height})
            chart.set_legend({"font": {"size": 7}})

            insert_row = row_start + r_idx * row_step
            insert_col = c_idx * col_step
            ws.insert_chart(insert_row, insert_col, chart)

        _stamp_all_sheets(wb)


def pilih_distribusi_terbaik_ks(result):
    ks_df_all = result["ks_df"].dropna(subset=["Rank"]).sort_values("Rank")
    quant_df = result["quant_df"]

    if len(ks_df_all) == 0:
        return None, None, "No", None, "Fail"

    def jarak(ratio):
        if ratio < RATIO_LOW:
            return RATIO_LOW - ratio
        if ratio > RATIO_HIGH:
            return ratio - RATIO_HIGH
        return 0.0

    def kumpulkan(ks_subset):
        candidates = []
        for _, row in ks_subset.iterrows():
            name = row["Distribution"]
            rank_val = row["Rank"]
            match = quant_df[quant_df["Distribution"] == name]
            if len(match) == 0:
                continue
            qrow = match.iloc[0]
            ratio = qrow.get("R100/R2", np.nan)
            if pd.isna(ratio):
                continue
            candidates.append((rank_val, name, qrow, ratio))
        return candidates

    def lulus_ks5(name):
        row = ks_df_all[ks_df_all["Distribution"] == name]
        if len(row) == 0:
            return "Fail"
        return row.iloc[0]["Pass KS 5%"]

    ks_lulus = ks_df_all[ks_df_all["Pass KS 5%"] == "Pass"]
    candidates = kumpulkan(ks_lulus)

    if not candidates:
        candidates = kumpulkan(ks_df_all)

    if not candidates:
        best_name = ks_df_all.iloc[0]["Distribution"]
        best_rank = ks_df_all.iloc[0]["Rank"]
        match = quant_df[quant_df["Distribution"] == best_name]
        qrow = match.iloc[0] if len(match) else None
        return best_name, qrow, "No", best_rank, lulus_ks5(best_name)

    for rank_val, name, qrow, ratio in candidates:
        if RATIO_LOW <= ratio <= RATIO_HIGH:
            return name, qrow, "Yes", rank_val, lulus_ks5(name)

    candidates_sorted = sorted(candidates, key=lambda c: (jarak(c[3]), c[0]))
    rank_val, name, qrow, ratio = candidates_sorted[0]
    return name, qrow, "Approaching", rank_val, lulus_ks5(name)


def export_rekap_excel(rekap_rows, out_path, dq_rekap_rows=None):
    rekap_df = pd.DataFrame(rekap_rows)
    with pd.ExcelWriter(out_path, engine="xlsxwriter") as writer:
        rekap_df.to_excel(writer, sheet_name="Best Return Period Summary", index=False, startrow=3)
        wb = writer.book
        ws = writer.sheets["Best Return Period Summary"]

        fmt_title = wb.add_format({"bold": True, "font_size": 13})
        fmt_hdr = wb.add_format({"bold": True, "bg_color": "#D9E1F2", "border": 1, "align": "center", "valign": "vcenter"})
        fmt_border = wb.add_format({"border": 1})
        fmt_num = wb.add_format({"num_format": "0.0000", "border": 1, "align": "center"})
        fmt_num3 = wb.add_format({"num_format": "0.000", "border": 1, "align": "center"})
        fmt_accept = wb.add_format({"bg_color": "#C6EFCE", "font_color": "#006100", "align": "center", "border": 1})
        fmt_reject = wb.add_format({"bg_color": "#FFC7CE", "font_color": "#9C0006", "align": "center", "border": 1})
        fmt_warning = wb.add_format({"bg_color": "#FFEB9C", "font_color": "#9C6500", "align": "center", "border": 1})
        fmt_hdr_out = wb.add_format({"bold": True, "bg_color": "#FCE4D6", "border": 1, "align": "center", "valign": "vcenter", "text_wrap": True})
        fmt_hdr_trn = wb.add_format({"bold": True, "bg_color": "#E2EFDA", "border": 1, "align": "center", "valign": "vcenter", "text_wrap": True})
        fmt_hdr_var = wb.add_format({"bold": True, "bg_color": "#DDEBF7", "border": 1, "align": "center", "valign": "vcenter", "text_wrap": True})
        fmt_hdr_mea = wb.add_format({"bold": True, "bg_color": "#EDD9F5", "border": 1, "align": "center", "valign": "vcenter", "text_wrap": True})
        fmt_hdr_ind = wb.add_format({"bold": True, "bg_color": "#FFF2CC", "border": 1, "align": "center", "valign": "vcenter", "text_wrap": True})

        ws.write(1, 0, "BEST DISTRIBUTION SUMMARY (KS TEST RANKING = 1) & RETURN PERIOD VALUES", fmt_title)
        for c, name in enumerate(rekap_df.columns):
            ws.write(3, c, name, fmt_hdr)

        for r in range(len(rekap_df)):
            for c, col in enumerate(rekap_df.columns):
                val = rekap_df.iat[r, c]
                if col in ("Series", "Selected Distribution", "Selected KS Rank"):
                    ws.write(r + 4, c, val, fmt_border)
                elif col == "Pass KS 5%?":
                    ws.write(r + 4, c, val, fmt_accept if val == "Pass" else fmt_reject)
                elif col == "R100/R2 In Range?":
                    if val == "Yes":
                        fmt_s = fmt_accept
                    elif val == "Approaching":
                        fmt_s = fmt_warning
                    else:
                        fmt_s = fmt_reject
                    ws.write(r + 4, c, val, fmt_s)
                elif pd.isna(val):
                    ws.write_blank(r + 4, c, None, fmt_border)
                else:
                    ws.write_number(r + 4, c, float(val), fmt_num)

        ws.set_column(0, 1, 28)
        ws.set_column(2, len(rekap_df.columns) - 1, 12)

        if dq_rekap_rows:
            ws2 = wb.add_worksheet("Data Quality Test Summary")
            writer.sheets["Data Quality Test Summary"] = ws2

            COL_DEFS = [
                ("No", "No", fmt_hdr, False, None),
                ("Series / Station", "Series", fmt_hdr, False, None),
                ("Upper Threshold\nXh", "out_Xh", fmt_hdr_out, False, None),
                ("Lower Threshold\nXl", "out_Xl", fmt_hdr_out, False, None),
                ("Kn Value", "out_Kn", fmt_hdr_out, False, None),
                ("Outlier Status", "out_Status", fmt_hdr_out, True, "No Outliers"),
                ("Rho\n(Spearman)", "trn_Rho", fmt_hdr_trn, False, None),
                ("T Calculated\n(Trend)", "trn_T_calc", fmt_hdr_trn, False, None),
                ("T Critical\n(Trend)", "trn_T_crit", fmt_hdr_trn, False, None),
                ("Trend Status", "trn_Status", fmt_hdr_trn, True, "No Trend"),
                ("F Calculated", "var_F_calc", fmt_hdr_var, False, None),
                ("F Critical", "var_F_crit", fmt_hdr_var, False, None),
                ("Variance Status", "var_Status", fmt_hdr_var, True, "Homogeneous"),
                ("T Calculated\n(Mean)", "mea_T_calc", fmt_hdr_mea, False, None),
                ("T Critical\n(Mean)", "mea_T_crit", fmt_hdr_mea, False, None),
                ("Mean Status", "mea_Status", fmt_hdr_mea, True, "Homogeneous"),
                ("r1\n(Lag-1)", "ind_r1", fmt_hdr_ind, False, None),
                ("Lower Limit\n(Ll)", "ind_Ll", fmt_hdr_ind, False, None),
                ("Upper Limit\n(Ul)", "ind_Ul", fmt_hdr_ind, False, None),
                ("Independence Status", "ind_Status", fmt_hdr_ind, True, "Independent"),
            ]

            ws2.write(1, 0, "DATA QUALITY TEST SUMMARY - ALL SERIES", fmt_title)
            ws2.write(2, 0, "One row per series/station. Green = Pass, Red = Fail.", fmt_border)

            groups = [
                (0, 1, ""),
                (2, 5, "OUTLIER TEST (Grubbs/WRC)"),
                (6, 9, "NO-TREND TEST (Spearman's Rho)"),
                (10, 12, "VARIANCE TEST (F-Test)"),
                (13, 15, "MEAN TEST (T-Test)"),
                (16, 19, "INDEPENDENCE TEST (Lag-1)"),
            ]
            group_fmt = {
                "OUTLIER TEST (Grubbs/WRC)": wb.add_format({"bold": True, "bg_color": "#FCE4D6", "border": 1, "align": "center", "valign": "vcenter"}),
                "NO-TREND TEST (Spearman's Rho)": wb.add_format({"bold": True, "bg_color": "#E2EFDA", "border": 1, "align": "center", "valign": "vcenter"}),
                "VARIANCE TEST (F-Test)": wb.add_format({"bold": True, "bg_color": "#DDEBF7", "border": 1, "align": "center", "valign": "vcenter"}),
                "MEAN TEST (T-Test)": wb.add_format({"bold": True, "bg_color": "#EDD9F5", "border": 1, "align": "center", "valign": "vcenter"}),
                "INDEPENDENCE TEST (Lag-1)": wb.add_format({"bold": True, "bg_color": "#FFF2CC", "border": 1, "align": "center", "valign": "vcenter"}),
            }
            blank_fmt = wb.add_format({"border": 1})
            for (c_start, c_end, label) in groups:
                if label == "":
                    for c in range(c_start, c_end + 1):
                        ws2.write(4, c, "", blank_fmt)
                elif c_start == c_end:
                    ws2.write(4, c_start, label, group_fmt.get(label, fmt_hdr))
                else:
                    ws2.merge_range(4, c_start, 4, c_end, label, group_fmt.get(label, fmt_hdr))

            for c, (hdr, key, hfmt, is_status, pass_kw) in enumerate(COL_DEFS):
                ws2.write(5, c, hdr, hfmt)

            for r, drow in enumerate(dq_rekap_rows):
                for c, (hdr, key, hfmt, is_status, pass_kw) in enumerate(COL_DEFS):
                    val = drow.get(key, "")
                    if key == "No":
                        ws2.write_number(r + 6, c, r + 1, fmt_num3)
                    elif is_status:
                        if isinstance(val, str) and pass_kw and pass_kw in val and "Not" not in val and "Non" not in val:
                            ws2.write(r + 6, c, val, fmt_accept)
                        else:
                            ws2.write(r + 6, c, val, fmt_reject)
                    elif isinstance(val, (int, float, np.integer, np.floating)) and not isinstance(val, bool):
                        ws2.write_number(r + 6, c, float(val), fmt_num3)
                    else:
                        ws2.write(r + 6, c, str(val) if val is not None else "", fmt_border)

            ws2.set_column(0, 0, 5)
            ws2.set_column(1, 1, 28)
            ws2.set_column(2, 3, 14)
            ws2.set_column(4, 4, 10)
            ws2.set_column(5, 5, 22)
            ws2.set_column(6, 8, 13)
            ws2.set_column(9, 9, 22)
            ws2.set_column(10, 11, 13)
            ws2.set_column(12, 12, 26)
            ws2.set_column(13, 14, 13)
            ws2.set_column(15, 15, 26)
            ws2.set_column(16, 18, 13)
            ws2.set_column(19, 19, 24)
            ws2.set_row(4, 20)
            ws2.set_row(5, 35)

        _stamp_all_sheets(wb)
    print(f"Best return period summary saved in: {out_path}")


def _safe_sheet_filename(name):
    name = str(name).strip()
    name = re.sub(r'[\\/*?:"<>|]', "_", name)
    return name if name else "Series"


# =======================================================================
#                       BATCH PROCESSING (UNIFIED)
# =======================================================================

def run_batch_from_series(series_dict, out_dir, log, progress):
    """Run the full analysis pipeline for every series in `series_dict`
    (a dict {series_name: DataFrame(Year, Rainfall)}), regardless of
    whether the data was imported from Excel or typed manually in the
    app. This is the single processing engine used by the GUI."""
    os.makedirs(out_dir, exist_ok=True)
    log(f"Output folder ready: {out_dir}")

    names = list(series_dict.keys())
    total = len(names)
    log(f"Found {total} series to process:")
    for s in names:
        log(f"   - {s}")

    success, failed = [], []
    rekap_rows = []
    dq_rekap_rows = []

    for idx, name in enumerate(names, start=1):
        log(f"\n=== Processing series: '{name}' ({idx}/{total}) ===")
        try:
            df = series_dict[name]
            ok, msg = validate_series_df(df)
            if not ok:
                raise ValueError(msg)
            log(f"  -> {len(df)} data rows read.")

            log("  Computing statistics & fitting distributions (incl. GEV via PWM/Hosking)...")
            result = run_analysis(df)

            out_name = f"{_safe_sheet_filename(name)}.xlsx"
            out_path = os.path.join(out_dir, out_name)

            log("  Building native Excel charts & writing report file...")
            export_excel(df, result, out_path)
            log(f"  Done! Result saved to: {out_path}")

            best_name, quant_row, in_range, rank_val, lulus_ks5_status = pilih_distribusi_terbaik_ks(result)
            entry = {"Series": name}
            if best_name is None or quant_row is None:
                entry["Selected Distribution"] = "-"
                entry["Selected KS Rank"] = "-"
                entry["Pass KS 5%?"] = "-"
                entry["R100/R2"] = np.nan
                entry["R100/R2 In Range?"] = "-"
                for T in RETURN_PERIODS:
                    entry[f"T={T}"] = np.nan
            else:
                entry["Selected Distribution"] = best_name
                entry["Selected KS Rank"] = int(rank_val) if pd.notna(rank_val) else "-"
                entry["Pass KS 5%?"] = lulus_ks5_status
                entry["R100/R2"] = quant_row.get("R100/R2", np.nan)
                entry["R100/R2 In Range?"] = in_range
                for T in RETURN_PERIODS:
                    entry[f"T={T}"] = quant_row[f"T={T}"]
            rekap_rows.append(entry)

            dq = result["dq_tests"]
            dq_entry = {
                "Series": name,
                "out_Xh": dq["outlier"]["summary"]["Upper Limit Xh"],
                "out_Xl": dq["outlier"]["summary"]["Lower Limit Xl"],
                "out_Kn": dq["outlier"]["summary"]["Kn (Table)"],
                "out_Status": "Safe (No Outliers)" if "No Outliers" in dq["outlier"]["summary"]["Conclusion"] else "Outliers Detected",
                "trn_Rho": dq["trend"]["summary"]["Rho (Spearman)"],
                "trn_T_calc": dq["trend"]["summary"]["t Calculated"],
                "trn_T_crit": dq["trend"]["summary"]["t Critical"],
                "trn_Status": dq["trend"]["summary"]["Conclusion"],
                "var_F_calc": dq["var"]["summary"]["F Calculated"],
                "var_F_crit": dq["var"]["summary"]["F Critical"],
                "var_Status": "Homogeneous Variance (Pass)" if ("Homogeneous" in dq["var"]["summary"]["Conclusion"] and "Not" not in dq["var"]["summary"]["Conclusion"]) else "Non-Homogeneous Variance (Fail)",
                "mea_T_calc": dq["mean"]["summary"]["T Calculated"],
                "mea_T_crit": dq["mean"]["summary"]["T Critical"],
                "mea_Status": "Homogeneous Mean (Pass)" if ("Homogeneous" in dq["mean"]["summary"]["Conclusion"] and "Not" not in dq["mean"]["summary"]["Conclusion"]) else "Non-Homogeneous Mean (Fail)",
                "ind_r1": dq["indep"]["summary"]["r1 (Correlation)"],
                "ind_Ll": dq["indep"]["summary"]["Lower Limit (Ll)"],
                "ind_Ul": dq["indep"]["summary"]["Upper Limit (Ul)"],
                "ind_Status": dq["indep"]["summary"]["Conclusion"],
            }
            dq_rekap_rows.append(dq_entry)

            success.append(name)

        except Exception as e:
            log(f"  FAILED to process series '{name}': {e}")
            failed.append((name, str(e)))
        finally:
            progress(idx, total)

    if rekap_rows:
        rekap_path = os.path.join(out_dir, "Best_Return_Period_Summary.xlsx")
        log(f"\nWriting final summary: {rekap_path}")
        export_rekap_excel(rekap_rows, rekap_path, dq_rekap_rows=dq_rekap_rows)
        log(f"Summary saved to: {rekap_path}")

    log("\n=================== SUMMARY ===================")
    log(f"Succeeded : {len(success)} series -> {success}")
    if failed:
        log(f"Failed    : {len(failed)} series")
        for s, err in failed:
            log(f"   - {s}: {err}")
    log("=================================================")

    return success, failed


def _open_in_explorer(path):
    try:
        if sys.platform.startswith("win"):
            os.startfile(path)
        elif sys.platform == "darwin":
            subprocess.Popen(["open", path])
        else:
            subprocess.Popen(["xdg-open", path])
    except Exception:
        pass


# =======================================================================
#                             GUI  -  THEME
# =======================================================================

COLOR_BG = "#EEF2F9"
COLOR_SIDEBAR = "#122A4D"
COLOR_SIDEBAR_ACTIVE = "#1C3F73"
COLOR_PRIMARY = "#1B4F91"
COLOR_ACCENT = "#0EA5A5"
COLOR_ACCENT_DARK = "#0B8484"
COLOR_CARD = "#FFFFFF"
COLOR_TEXT = "#1E293B"
COLOR_MUTED = "#64748B"
COLOR_SUCCESS = "#15803D"
COLOR_SUCCESS_BG = "#DCFCE7"
COLOR_DANGER = "#B91C1C"
COLOR_DANGER_BG = "#FEE2E2"
COLOR_WARN = "#B45309"
COLOR_WARN_BG = "#FEF3C7"
COLOR_BORDER = "#D9E1EC"

FONT_FAMILY = "Segoe UI" if sys.platform.startswith("win") else "Helvetica"


class Card(tk.Frame):
    """A simple rounded-look card container (flat bg + border) used to
    group related controls in the dashboard."""

    def __init__(self, parent, title=None, subtitle=None, **kwargs):
        super().__init__(parent, bg=COLOR_CARD, highlightbackground=COLOR_BORDER,
                          highlightthickness=1, bd=0, **kwargs)
        self._pad = tk.Frame(self, bg=COLOR_CARD)
        self._pad.pack(fill="both", expand=True, padx=16, pady=14)
        if title:
            tk.Label(self._pad, text=title, bg=COLOR_CARD, fg=COLOR_TEXT,
                      font=(FONT_FAMILY, 12, "bold")).pack(anchor="w")
        if subtitle:
            tk.Label(self._pad, text=subtitle, bg=COLOR_CARD, fg=COLOR_MUTED,
                      font=(FONT_FAMILY, 9)).pack(anchor="w", pady=(0, 6))

    @property
    def body(self):
        return self._pad


def make_pill_button(parent, text, command, bg=COLOR_PRIMARY, fg="white",
                      hover=None, font_size=10, padx=14, pady=7, state="normal"):
    hover = hover or bg
    btn = tk.Label(parent, text=text, bg=bg, fg=fg, font=(FONT_FAMILY, font_size, "bold"),
                    padx=padx, pady=pady, cursor="hand2")
    btn._enabled = (state == "normal")
    btn._bg = bg
    btn._hover = hover

    def on_click(e):
        if btn._enabled:
            command()

    def on_enter(e):
        if btn._enabled:
            btn.configure(bg=btn._hover)

    def on_leave(e):
        if btn._enabled:
            btn.configure(bg=btn._bg)

    btn.bind("<Button-1>", on_click)
    btn.bind("<Enter>", on_enter)
    btn.bind("<Leave>", on_leave)

    def set_state(state):
        btn._enabled = (state == "normal")
        btn.configure(bg=btn._bg if btn._enabled else "#B7C3D6",
                      fg="white" if btn._enabled else "#EEF2F9",
                      cursor="hand2" if btn._enabled else "arrow")

    btn.set_state = set_state
    if state != "normal":
        set_state(state)
    return btn


def status_badge(parent, text, kind="neutral"):
    palette = {
        "success": (COLOR_SUCCESS_BG, COLOR_SUCCESS),
        "danger": (COLOR_DANGER_BG, COLOR_DANGER),
        "warn": (COLOR_WARN_BG, COLOR_WARN),
        "neutral": ("#E2E8F0", COLOR_MUTED),
    }
    bg, fg = palette.get(kind, palette["neutral"])
    return tk.Label(parent, text=text, bg=bg, fg=fg, font=(FONT_FAMILY, 9, "bold"),
                     padx=10, pady=3)


# =======================================================================
#                    MANUAL SERIES EDITOR (spreadsheet-like)
# =======================================================================

class ManualSeriesEditor(tk.Toplevel):
    """A small spreadsheet-like editor for typing rainfall data directly
    into the app: one series = one editable table of (Year, Rainfall)
    rows, exactly like working on one sheet of an Excel workbook."""

    def __init__(self, master, series_name, df, on_save):
        super().__init__(master)
        self.title(f"Manual Data Entry - {series_name}")
        self.geometry("620x620")
        self.minsize(520, 480)
        self.configure(bg=COLOR_BG)
        self.transient(master)
        self.grab_set()

        self.series_name = series_name
        self.on_save = on_save
        self.df = df.copy().reset_index(drop=True)

        self._build_ui()
        self._refresh_table()
        self._refresh_stats()

    def _build_ui(self):
        header = tk.Frame(self, bg=COLOR_PRIMARY)
        header.pack(fill="x")
        tk.Label(header, text=f"\u270E  {self.series_name}", bg=COLOR_PRIMARY, fg="white",
                  font=(FONT_FAMILY, 13, "bold")).pack(anchor="w", padx=16, pady=10)

        toolbar = tk.Frame(self, bg=COLOR_BG)
        toolbar.pack(fill="x", padx=14, pady=(10, 4))
        make_pill_button(toolbar, "+ Add Row", self._add_row, bg=COLOR_ACCENT,
                          hover=COLOR_ACCENT_DARK, font_size=9).pack(side="left", padx=(0, 6))
        make_pill_button(toolbar, "Edit Selected", self._edit_row, bg="#475569",
                          hover="#334155", font_size=9).pack(side="left", padx=6)
        make_pill_button(toolbar, "Delete Selected", self._delete_row, bg=COLOR_DANGER,
                          hover="#991B1B", font_size=9).pack(side="left", padx=6)
        make_pill_button(toolbar, "Sort by Year", self._sort_rows, bg="#475569",
                          hover="#334155", font_size=9).pack(side="left", padx=6)
        make_pill_button(toolbar, "\u2398 Paste / Bulk Add", self._open_bulk_paste, bg="#7C3AED",
                          hover="#6D28D9", font_size=9).pack(side="left", padx=6)

        table_card = Card(self, title="Data Table (Year / Rainfall)",
                           subtitle="Double-click a row to edit it.")
        table_card.pack(fill="both", expand=True, padx=14, pady=6)

        cols = ("no", "year", "value")
        self.tree = ttk.Treeview(table_card.body, columns=cols, show="headings", height=14)
        self.tree.heading("no", text="No")
        self.tree.heading("year", text="Year")
        self.tree.heading("value", text=COL_DATA)
        self.tree.column("no", width=50, anchor="center")
        self.tree.column("year", width=110, anchor="center")
        self.tree.column("value", width=140, anchor="center")
        vsb = ttk.Scrollbar(table_card.body, orient="vertical", command=self.tree.yview)
        self.tree.configure(yscrollcommand=vsb.set)
        self.tree.pack(side="left", fill="both", expand=True)
        vsb.pack(side="right", fill="y")
        self.tree.bind("<Double-1>", lambda e: self._edit_row())

        stats_card = Card(self, title="Live Preview Statistics")
        stats_card.pack(fill="x", padx=14, pady=(0, 6))
        self.lbl_stats = tk.Label(stats_card.body, text="-", bg=COLOR_CARD, fg=COLOR_MUTED,
                                    font=(FONT_FAMILY, 9), justify="left")
        self.lbl_stats.pack(anchor="w")

        footer = tk.Frame(self, bg=COLOR_BG)
        footer.pack(fill="x", padx=14, pady=(0, 14))
        make_pill_button(footer, "Save & Close", self._save_and_close, bg=COLOR_PRIMARY,
                          hover="#123A6B").pack(side="right")
        make_pill_button(footer, "Cancel", self.destroy, bg="#94A3B8",
                          hover="#64748B").pack(side="right", padx=(0, 8))

    def _refresh_table(self):
        for i in self.tree.get_children():
            self.tree.delete(i)
        for i, row in self.df.iterrows():
            self.tree.insert("", "end", iid=str(i), values=(i + 1, row[COL_TAHUN], row[COL_DATA]))

    def _refresh_stats(self):
        n = len(self.df)
        if n == 0:
            self.lbl_stats.configure(text="No data yet - use '+ Add Row' or 'Paste / Bulk Add' to begin.")
            return
        ok, msg = validate_series_df(self.df)
        txt = f"N = {n}"
        if n >= 2:
            x = self.df[COL_DATA].to_numpy(dtype=float)
            txt += f"   |   Mean = {x.mean():.3f}   |   Std Dev = {x.std(ddof=1):.3f}"
        if n >= 3:
            m = basic_moments(x)
            txt += f"   |   Cv = {m['cv']:.3f}   |   Cs = {m['cs']:.3f}"
        txt += f"\nStatus: {msg}"
        self.lbl_stats.configure(text=txt)

    def _row_dialog(self, year=None, value=None):
        dlg = tk.Toplevel(self)
        dlg.title("Row")
        dlg.configure(bg=COLOR_BG)
        dlg.geometry("300x180")
        dlg.transient(self)
        dlg.grab_set()

        tk.Label(dlg, text="Year", bg=COLOR_BG, fg=COLOR_TEXT, font=(FONT_FAMILY, 10)).pack(pady=(16, 2))
        e_year = tk.Entry(dlg, font=(FONT_FAMILY, 11), justify="center")
        e_year.pack(padx=30, fill="x")
        if year is not None:
            e_year.insert(0, str(year))

        tk.Label(dlg, text=COL_DATA, bg=COLOR_BG, fg=COLOR_TEXT, font=(FONT_FAMILY, 10)).pack(pady=(12, 2))
        e_val = tk.Entry(dlg, font=(FONT_FAMILY, 11), justify="center")
        e_val.pack(padx=30, fill="x")
        if value is not None:
            e_val.insert(0, str(value))

        result = {}

        def on_ok():
            try:
                y = float(e_year.get().strip().replace(",", "."))
                v = float(e_val.get().strip().replace(",", "."))
                y = int(y) if float(y).is_integer() else y
                result["year"] = y
                result["value"] = v
                dlg.destroy()
            except Exception:
                messagebox.showerror("Invalid Input", "Please enter valid numbers for Year and Value.", parent=dlg)

        btns = tk.Frame(dlg, bg=COLOR_BG)
        btns.pack(pady=16)
        make_pill_button(btns, "OK", on_ok, bg=COLOR_PRIMARY).pack(side="left", padx=6)
        make_pill_button(btns, "Cancel", dlg.destroy, bg="#94A3B8").pack(side="left", padx=6)

        e_year.focus_set()
        dlg.wait_window()
        return result.get("year"), result.get("value")

    def _add_row(self):
        year, value = self._row_dialog()
        if year is None:
            return
        self.df.loc[len(self.df)] = [year, value]
        self._refresh_table()
        self._refresh_stats()

    def _selected_index(self):
        sel = self.tree.selection()
        if not sel:
            messagebox.showinfo("No Selection", "Please select a row first.", parent=self)
            return None
        return int(sel[0])

    def _edit_row(self):
        idx = self._selected_index()
        if idx is None:
            return
        cur_year = self.df.at[idx, COL_TAHUN]
        cur_val = self.df.at[idx, COL_DATA]
        year, value = self._row_dialog(cur_year, cur_val)
        if year is None:
            return
        self.df.at[idx, COL_TAHUN] = year
        self.df.at[idx, COL_DATA] = value
        self._refresh_table()
        self._refresh_stats()

    def _delete_row(self):
        idx = self._selected_index()
        if idx is None:
            return
        self.df = self.df.drop(index=idx).reset_index(drop=True)
        self._refresh_table()
        self._refresh_stats()

    def _sort_rows(self):
        try:
            self.df = self.df.sort_values(COL_TAHUN).reset_index(drop=True)
            self._refresh_table()
        except Exception:
            pass

    def _open_bulk_paste(self):
        dlg = tk.Toplevel(self)
        dlg.title("Paste / Bulk Add Data")
        dlg.configure(bg=COLOR_BG)
        dlg.geometry("460x420")
        dlg.transient(self)
        dlg.grab_set()

        tk.Label(dlg, text="Paste rows below (one per line).", bg=COLOR_BG, fg=COLOR_TEXT,
                  font=(FONT_FAMILY, 10, "bold")).pack(anchor="w", padx=16, pady=(14, 0))
        tk.Label(dlg, text="Accepted formats:  2010,120.5   /   2010\\t120.5   /   2010 120.5",
                  bg=COLOR_BG, fg=COLOR_MUTED, font=(FONT_FAMILY, 9)).pack(anchor="w", padx=16, pady=(0, 8))

        txt = scrolledtext.ScrolledText(dlg, height=14, font=("Consolas", 10))
        txt.pack(fill="both", expand=True, padx=16)

        def on_add():
            raw = txt.get("1.0", "end")
            rows, bad = parse_bulk_text(raw)
            if not rows:
                messagebox.showwarning("Nothing to Add", "No valid Year / Value pairs were found.", parent=dlg)
                return
            for y, v in rows:
                self.df.loc[len(self.df)] = [y, v]
            self._refresh_table()
            self._refresh_stats()
            msg = f"Added {len(rows)} row(s)."
            if bad:
                msg += f"\n{len(bad)} line(s) could not be parsed and were skipped."
            messagebox.showinfo("Bulk Add Complete", msg, parent=dlg)
            dlg.destroy()

        btns = tk.Frame(dlg, bg=COLOR_BG)
        btns.pack(pady=12)
        make_pill_button(btns, "Add These Rows", on_add, bg=COLOR_ACCENT,
                          hover=COLOR_ACCENT_DARK).pack(side="left", padx=6)
        make_pill_button(btns, "Cancel", dlg.destroy, bg="#94A3B8").pack(side="left", padx=6)

    def _save_and_close(self):
        self.on_save(self.series_name, self.df)
        self.destroy()


# =======================================================================
#                          MAIN APPLICATION
# =======================================================================

APP_TITLE = "Rainfall Frequency Analysis"
APP_VERSION = "3.0"


class RainfallApp(tk.Tk):
    """Main application window with a sidebar-driven, dashboard-style
    workspace:
        1) Data Source   - import an Excel workbook and/or create /
                            edit manually typed series
        2) Series Manager - unified list of every series currently
                            loaded (regardless of source), with quick
                            stats and status
        3) Process & Log  - choose an output folder, run the analysis
                            for every loaded series, and watch progress
    """

    def __init__(self):
        super().__init__()
        self.title(f"{APP_TITLE}  v{APP_VERSION}")
        self.geometry("1180x760")
        self.minsize(980, 640)
        self.configure(bg=COLOR_BG)

        # unified data store: name -> DataFrame(Year, Rainfall)
        self.series_store = {}
        # name -> "Excel: <file>" or "Manual"
        self.series_source = {}

        self.out_dir = None
        self.worker_thread = None
        self.msg_queue = queue.Queue()

        self._build_style()
        self._build_layout()
        self.protocol("WM_DELETE_WINDOW", self._on_close)
        self.after(100, self._poll_queue)
        self._show_page("data")

    # ------------------------------------------------------------------
    # STYLE
    # ------------------------------------------------------------------
    def _build_style(self):
        style = ttk.Style(self)
        try:
            style.theme_use("clam")
        except Exception:
            pass
        style.configure("Treeview", rowheight=26, font=(FONT_FAMILY, 9),
                          background="white", fieldbackground="white", bordercolor=COLOR_BORDER)
        style.configure("Treeview.Heading", font=(FONT_FAMILY, 9, "bold"),
                          background="#DCE6F5", foreground=COLOR_TEXT, relief="flat")
        style.map("Treeview", background=[("selected", COLOR_ACCENT)], foreground=[("selected", "white")])
        style.configure("TProgressbar", background=COLOR_ACCENT, troughcolor="#DCE3EE",
                          bordercolor=COLOR_BG, lightcolor=COLOR_ACCENT, darkcolor=COLOR_ACCENT)

    # ------------------------------------------------------------------
    # LAYOUT
    # ------------------------------------------------------------------
    def _build_layout(self):
        # ---- Sidebar ----
        self.sidebar = tk.Frame(self, bg=COLOR_SIDEBAR, width=230)
        self.sidebar.pack(side="left", fill="y")
        self.sidebar.pack_propagate(False)

        brand = tk.Frame(self.sidebar, bg=COLOR_SIDEBAR)
        brand.pack(fill="x", pady=(22, 26), padx=18)
        tk.Label(brand, text="\U0001F327", bg=COLOR_SIDEBAR, fg="white",
                  font=(FONT_FAMILY, 26)).pack(anchor="w")
        tk.Label(brand, text="Rainfall Frequency", bg=COLOR_SIDEBAR, fg="white",
                  font=(FONT_FAMILY, 14, "bold")).pack(anchor="w", pady=(4, 0))
        tk.Label(brand, text="Analysis Studio", bg=COLOR_SIDEBAR, fg="#9FB4D8",
                  font=(FONT_FAMILY, 10)).pack(anchor="w")

        self.nav_buttons = {}
        nav_items = [
            ("data", "\U0001F4C2  1. Data Source"),
            ("manager", "\U0001F5C2  2. Series Manager"),
            ("process", "\u25B6  3. Process & Results"),
        ]
        for key, label in nav_items:
            b = tk.Label(self.sidebar, text=label, bg=COLOR_SIDEBAR, fg="#D6E2F5",
                          font=(FONT_FAMILY, 11), anchor="w", padx=18, pady=12, cursor="hand2")
            b.pack(fill="x")
            b.bind("<Button-1>", lambda e, k=key: self._show_page(k))
            b.bind("<Enter>", lambda e, b=b, k=key: b.configure(bg=COLOR_SIDEBAR_ACTIVE) if self.current_page != k else None)
            b.bind("<Leave>", lambda e, b=b, k=key: b.configure(bg=COLOR_SIDEBAR) if self.current_page != k else None)
            self.nav_buttons[key] = b

        tk.Frame(self.sidebar, bg=COLOR_SIDEBAR).pack(fill="both", expand=True)

        about_btn = tk.Label(self.sidebar, text="\u2139  About", bg=COLOR_SIDEBAR, fg="#9FB4D8",
                               font=(FONT_FAMILY, 10), anchor="w", padx=18, pady=10, cursor="hand2")
        about_btn.pack(fill="x", side="bottom")
        about_btn.bind("<Button-1>", lambda e: self._show_about())

        # ---- Main area ----
        main = tk.Frame(self, bg=COLOR_BG)
        main.pack(side="left", fill="both", expand=True)

        self.header = tk.Frame(main, bg=COLOR_BG)
        self.header.pack(fill="x", padx=24, pady=(20, 6))
        self.lbl_page_title = tk.Label(self.header, text="", bg=COLOR_BG, fg=COLOR_TEXT,
                                          font=(FONT_FAMILY, 18, "bold"))
        self.lbl_page_title.pack(anchor="w")
        self.lbl_page_subtitle = tk.Label(self.header, text="", bg=COLOR_BG, fg=COLOR_MUTED,
                                             font=(FONT_FAMILY, 10))
        self.lbl_page_subtitle.pack(anchor="w")

        self.page_container = tk.Frame(main, bg=COLOR_BG)
        self.page_container.pack(fill="both", expand=True, padx=24, pady=(6, 16))

        self.pages = {}
        self.pages["data"] = self._build_page_data(self.page_container)
        self.pages["manager"] = self._build_page_manager(self.page_container)
        self.pages["process"] = self._build_page_process(self.page_container)

        self.current_page = None

    def _show_page(self, key):
        titles = {
            "data": ("Data Source", "Import a multi-sheet Excel workbook and/or create series manually - just like adding a new Excel sheet, without leaving the app."),
            "manager": ("Series Manager", "Every series currently loaded, from any source, in one unified list."),
            "process": ("Process & Results", "Choose an output folder and run the full 6-distribution analysis for every loaded series."),
        }
        for k, b in self.nav_buttons.items():
            b.configure(bg=COLOR_SIDEBAR_ACTIVE if k == key else COLOR_SIDEBAR,
                        fg="white" if k == key else "#D6E2F5")
        for k, frame in self.pages.items():
            frame.pack_forget()
        self.pages[key].pack(fill="both", expand=True)
        title, subtitle = titles[key]
        self.lbl_page_title.configure(text=title)
        self.lbl_page_subtitle.configure(text=subtitle)
        self.current_page = key
        if key == "manager":
            self._refresh_manager_table()
        if key == "process":
            self._refresh_process_summary()

    # ------------------------------------------------------------------
    # PAGE 1: DATA SOURCE
    # ------------------------------------------------------------------
    def _build_page_data(self, parent):
        page = tk.Frame(parent, bg=COLOR_BG)

        cols_frame = tk.Frame(page, bg=COLOR_BG)
        cols_frame.pack(fill="both", expand=True)
        cols_frame.columnconfigure(0, weight=1)
        cols_frame.columnconfigure(1, weight=1)
        cols_frame.rowconfigure(0, weight=1)

        # --- Left card: Import from Excel ---
        left = Card(cols_frame, title="\U0001F4E5  Import from Excel",
                    subtitle="Each sheet in the workbook becomes one data series (columns 'Year' and 'Rainfall (mm)').")
        left.grid(row=0, column=0, sticky="nsew", padx=(0, 10))

        make_pill_button(left.body, "Browse Excel File...", self._on_import_excel,
                          bg=COLOR_PRIMARY, hover="#123A6B").pack(anchor="w", pady=(4, 10))

        self.lbl_import_file = tk.Label(left.body, text="No file imported yet.", bg=COLOR_CARD,
                                          fg=COLOR_MUTED, font=(FONT_FAMILY, 9))
        self.lbl_import_file.pack(anchor="w", pady=(0, 8))

        preview_frame = tk.Frame(left.body, bg=COLOR_CARD)
        preview_frame.pack(fill="both", expand=True)
        cols = ("sheet", "n", "mean", "std", "status")
        self.import_tree = ttk.Treeview(preview_frame, columns=cols, show="headings", height=10)
        headers = {"sheet": "Sheet", "n": "N", "mean": "Mean", "std": "Std Dev", "status": "Status"}
        widths = {"sheet": 150, "n": 50, "mean": 80, "std": 80, "status": 150}
        for c in cols:
            self.import_tree.heading(c, text=headers[c])
            self.import_tree.column(c, width=widths[c], anchor="center" if c != "sheet" else "w")
        vsb = ttk.Scrollbar(preview_frame, orient="vertical", command=self.import_tree.yview)
        self.import_tree.configure(yscrollcommand=vsb.set)
        self.import_tree.pack(side="left", fill="both", expand=True)
        vsb.pack(side="right", fill="y")

        # --- Right card: Manual entry ---
        right = Card(cols_frame, title="\u270E  Manual Entry",
                    subtitle="Type or paste your own rainfall series directly - each one behaves like a separate Excel sheet.")
        right.grid(row=0, column=1, sticky="nsew", padx=(10, 0))

        row_btns = tk.Frame(right.body, bg=COLOR_CARD)
        row_btns.pack(fill="x", pady=(4, 10))
        make_pill_button(row_btns, "+ New Manual Series", self._on_new_manual_series,
                          bg=COLOR_ACCENT, hover=COLOR_ACCENT_DARK).pack(side="left")

        tk.Label(right.body, text="Manual series created in this session:", bg=COLOR_CARD,
                  fg=COLOR_MUTED, font=(FONT_FAMILY, 9)).pack(anchor="w", pady=(4, 4))

        manual_frame = tk.Frame(right.body, bg=COLOR_CARD)
        manual_frame.pack(fill="both", expand=True)
        cols2 = ("name", "n", "status")
        self.manual_tree = ttk.Treeview(manual_frame, columns=cols2, show="headings", height=10)
        self.manual_tree.heading("name", text="Series Name")
        self.manual_tree.heading("n", text="N")
        self.manual_tree.heading("status", text="Status")
        self.manual_tree.column("name", width=160, anchor="w")
        self.manual_tree.column("n", width=50, anchor="center")
        self.manual_tree.column("status", width=150, anchor="center")
        vsb2 = ttk.Scrollbar(manual_frame, orient="vertical", command=self.manual_tree.yview)
        self.manual_tree.configure(yscrollcommand=vsb2.set)
        self.manual_tree.pack(side="left", fill="both", expand=True)
        vsb2.pack(side="right", fill="y")
        self.manual_tree.bind("<Double-1>", lambda e: self._on_edit_manual_series())

        manual_btns2 = tk.Frame(right.body, bg=COLOR_CARD)
        manual_btns2.pack(fill="x", pady=(8, 0))
        make_pill_button(manual_btns2, "Edit Selected", self._on_edit_manual_series, bg="#475569",
                          hover="#334155", font_size=9).pack(side="left", padx=(0, 6))
        make_pill_button(manual_btns2, "Delete Selected", self._on_delete_manual_series, bg=COLOR_DANGER,
                          hover="#991B1B", font_size=9).pack(side="left")

        return page

    def _on_import_excel(self):
        path = filedialog.askopenfilename(
            title="Select Input Excel File",
            filetypes=[("Excel files", "*.xlsx *.xls"), ("All files", "*.*")],
        )
        if not path:
            return
        self.lbl_import_file.configure(text=f"Reading: {os.path.basename(path)} ...")
        self.update_idletasks()
        try:
            series, errors = read_excel_all_sheets(path)
        except Exception as e:
            messagebox.showerror("Import Failed", f"Could not read this Excel file:\n\n{e}")
            self.lbl_import_file.configure(text="Import failed.")
            return

        for i in self.import_tree.get_children():
            self.import_tree.delete(i)

        added = 0
        for name, df in series.items():
            unique_name = self._unique_series_name(name)
            self.series_store[unique_name] = df
            self.series_source[unique_name] = f"Excel: {os.path.basename(path)}"
            ok, msg = validate_series_df(df)
            m = basic_moments(df[COL_DATA].to_numpy(dtype=float)) if len(df) >= 2 else None
            mean_txt = f"{m['mean']:.2f}" if m else "-"
            std_txt = f"{m['std']:.2f}" if m else "-"
            self.import_tree.insert("", "end", values=(unique_name, len(df), mean_txt, std_txt, msg))
            added += 1

        for name, err in errors.items():
            self.import_tree.insert("", "end", values=(name, "-", "-", "-", f"Skipped: {err}"))

        self.lbl_import_file.configure(
            text=f"Imported {added} series from '{os.path.basename(path)}' "
                 f"({len(errors)} sheet(s) skipped)." if errors else
            f"Imported {added} series from '{os.path.basename(path)}'.")

        if self.current_page == "manager":
            self._refresh_manager_table()

    def _unique_series_name(self, base_name):
        name = str(base_name).strip() or "Series"
        if name not in self.series_store:
            return name
        i = 2
        while f"{name} ({i})" in self.series_store:
            i += 1
        return f"{name} ({i})"

    def _on_new_manual_series(self):
        dlg = tk.Toplevel(self)
        dlg.title("New Manual Series")
        dlg.configure(bg=COLOR_BG)
        dlg.geometry("340x160")
        dlg.transient(self)
        dlg.grab_set()

        tk.Label(dlg, text="Series / Station Name", bg=COLOR_BG, fg=COLOR_TEXT,
                  font=(FONT_FAMILY, 10, "bold")).pack(pady=(20, 4))
        e_name = tk.Entry(dlg, font=(FONT_FAMILY, 11), justify="center")
        e_name.pack(padx=30, fill="x")
        e_name.focus_set()

        def on_ok():
            name = e_name.get().strip()
            if not name:
                messagebox.showerror("Invalid Name", "Please enter a series name.", parent=dlg)
                return
            unique_name = self._unique_series_name(name)
            dlg.destroy()
            self.series_store[unique_name] = empty_series_df()
            self.series_source[unique_name] = "Manual"
            self._refresh_manual_tree()
            self._open_manual_editor(unique_name)

        btns = tk.Frame(dlg, bg=COLOR_BG)
        btns.pack(pady=18)
        make_pill_button(btns, "Create & Open Editor", on_ok, bg=COLOR_ACCENT,
                          hover=COLOR_ACCENT_DARK).pack(side="left", padx=6)
        make_pill_button(btns, "Cancel", dlg.destroy, bg="#94A3B8").pack(side="left", padx=6)
        dlg.bind("<Return>", lambda e: on_ok())

    def _open_manual_editor(self, name):
        df = self.series_store.get(name, empty_series_df())

        def on_save(series_name, new_df):
            self.series_store[series_name] = new_df.reset_index(drop=True)
            self._refresh_manual_tree()
            if self.current_page == "manager":
                self._refresh_manager_table()
            if self.current_page == "process":
                self._refresh_process_summary()

        ManualSeriesEditor(self, name, df, on_save)

    def _refresh_manual_tree(self):
        for i in self.manual_tree.get_children():
            self.manual_tree.delete(i)
        for name, src in self.series_source.items():
            if src != "Manual":
                continue
            df = self.series_store[name]
            ok, msg = validate_series_df(df)
            self.manual_tree.insert("", "end", values=(name, len(df), msg))

    def _selected_manual_name(self):
        sel = self.manual_tree.selection()
        if not sel:
            messagebox.showinfo("No Selection", "Please select a manual series first.")
            return None
        return self.manual_tree.item(sel[0])["values"][0]

    def _on_edit_manual_series(self):
        name = self._selected_manual_name()
        if name is None:
            return
        self._open_manual_editor(str(name))

    def _on_delete_manual_series(self):
        name = self._selected_manual_name()
        if name is None:
            return
        if messagebox.askyesno("Confirm Delete", f"Remove manual series '{name}'?"):
            self.series_store.pop(name, None)
            self.series_source.pop(name, None)
            self._refresh_manual_tree()
            if self.current_page == "manager":
                self._refresh_manager_table()

    # ------------------------------------------------------------------
    # PAGE 2: SERIES MANAGER
    # ------------------------------------------------------------------
    def _build_page_manager(self, parent):
        page = tk.Frame(parent, bg=COLOR_BG)

        summary = tk.Frame(page, bg=COLOR_BG)
        summary.pack(fill="x", pady=(0, 10))
        self.manager_summary_labels = []
        for i in range(3):
            card = tk.Frame(summary, bg=COLOR_CARD, highlightbackground=COLOR_BORDER, highlightthickness=1)
            card.pack(side="left", fill="x", expand=True, padx=(0 if i == 0 else 8, 0))
            lbl_val = tk.Label(card, text="0", bg=COLOR_CARD, fg=COLOR_PRIMARY, font=(FONT_FAMILY, 20, "bold"))
            lbl_val.pack(anchor="w", padx=14, pady=(10, 0))
            lbl_name = tk.Label(card, text="", bg=COLOR_CARD, fg=COLOR_MUTED, font=(FONT_FAMILY, 9))
            lbl_name.pack(anchor="w", padx=14, pady=(0, 10))
            self.manager_summary_labels.append((lbl_val, lbl_name))
        self.manager_summary_labels[0][1].configure(text="Total series loaded")
        self.manager_summary_labels[1][1].configure(text="Ready to process")
        self.manager_summary_labels[2][1].configure(text="Not ready (needs data)")

        card = Card(page, title="All Loaded Series",
                    subtitle="Combined list from Excel import and manual entry. Select a row to rename, edit, or remove it.")
        card.pack(fill="both", expand=True)

        toolbar = tk.Frame(card.body, bg=COLOR_CARD)
        toolbar.pack(fill="x", pady=(0, 8))
        make_pill_button(toolbar, "\u21BB Refresh", self._refresh_manager_table, bg="#475569",
                          hover="#334155", font_size=9).pack(side="left", padx=(0, 6))
        make_pill_button(toolbar, "Edit / View", self._manager_edit_selected, bg=COLOR_PRIMARY,
                          hover="#123A6B", font_size=9).pack(side="left", padx=6)
        make_pill_button(toolbar, "Rename", self._manager_rename_selected, bg="#475569",
                          hover="#334155", font_size=9).pack(side="left", padx=6)
        make_pill_button(toolbar, "Remove", self._manager_remove_selected, bg=COLOR_DANGER,
                          hover="#991B1B", font_size=9).pack(side="left", padx=6)
        make_pill_button(toolbar, "Remove All", self._manager_remove_all, bg=COLOR_DANGER,
                          hover="#991B1B", font_size=9).pack(side="left", padx=6)

        table_frame = tk.Frame(card.body, bg=COLOR_CARD)
        table_frame.pack(fill="both", expand=True)
        cols = ("name", "source", "n", "mean", "std", "cv", "cs", "status")
        headers = {"name": "Series / Station", "source": "Source", "n": "N", "mean": "Mean",
                   "std": "Std Dev", "cv": "Cv", "cs": "Cs", "status": "Status"}
        widths = {"name": 190, "source": 170, "n": 50, "mean": 80, "std": 80, "cv": 60, "cs": 60, "status": 170}
        self.manager_tree = ttk.Treeview(table_frame, columns=cols, show="headings", height=14)
        for c in cols:
            self.manager_tree.heading(c, text=headers[c])
            self.manager_tree.column(c, width=widths[c], anchor="w" if c in ("name", "source") else "center")
        vsb = ttk.Scrollbar(table_frame, orient="vertical", command=self.manager_tree.yview)
        self.manager_tree.configure(yscrollcommand=vsb.set)
        self.manager_tree.pack(side="left", fill="both", expand=True)
        vsb.pack(side="right", fill="y")
        self.manager_tree.bind("<Double-1>", lambda e: self._manager_edit_selected())

        return page

    def _refresh_manager_table(self):
        for i in self.manager_tree.get_children():
            self.manager_tree.delete(i)
        ready_count = 0
        for name, df in self.series_store.items():
            src = self.series_source.get(name, "-")
            ok, msg = validate_series_df(df)
            if ok:
                ready_count += 1
            if len(df) >= 2:
                x = df[COL_DATA].to_numpy(dtype=float)
                mean_txt = f"{x.mean():.2f}"
                std_txt = f"{x.std(ddof=1):.2f}"
            else:
                mean_txt = std_txt = "-"
            if len(df) >= 3:
                m = basic_moments(df[COL_DATA].to_numpy(dtype=float))
                cv_txt = f"{m['cv']:.2f}"
                cs_txt = f"{m['cs']:.2f}"
            else:
                cv_txt = cs_txt = "-"
            self.manager_tree.insert("", "end", values=(name, src, len(df), mean_txt, std_txt, cv_txt, cs_txt, msg))

        total = len(self.series_store)
        not_ready = total - ready_count
        self.manager_summary_labels[0][0].configure(text=str(total))
        self.manager_summary_labels[1][0].configure(text=str(ready_count))
        self.manager_summary_labels[2][0].configure(text=str(not_ready))

    def _selected_manager_name(self):
        sel = self.manager_tree.selection()
        if not sel:
            messagebox.showinfo("No Selection", "Please select a series first.")
            return None
        return self.manager_tree.item(sel[0])["values"][0]

    def _manager_edit_selected(self):
        name = self._selected_manager_name()
        if name is None:
            return
        name = str(name)
        self._open_manual_editor(name)

    def _manager_rename_selected(self):
        name = self._selected_manager_name()
        if name is None:
            return
        name = str(name)
        new_name = tk.simpledialog.askstring("Rename Series", "New name:", initialvalue=name, parent=self) \
            if hasattr(tk, "simpledialog") else None
        if not new_name:
            return
        new_name = new_name.strip()
        if not new_name or new_name == name:
            return
        if new_name in self.series_store:
            messagebox.showerror("Name In Use", "A series with that name already exists.")
            return
        self.series_store[new_name] = self.series_store.pop(name)
        self.series_source[new_name] = self.series_source.pop(name, "-")
        self._refresh_manager_table()
        self._refresh_manual_tree()

    def _manager_remove_selected(self):
        name = self._selected_manager_name()
        if name is None:
            return
        name = str(name)
        if messagebox.askyesno("Confirm Remove", f"Remove series '{name}' from the workspace?"):
            self.series_store.pop(name, None)
            self.series_source.pop(name, None)
            self._refresh_manager_table()
            self._refresh_manual_tree()

    def _manager_remove_all(self):
        if not self.series_store:
            return
        if messagebox.askyesno("Confirm Remove All", "Remove ALL loaded series from the workspace?"):
            self.series_store.clear()
            self.series_source.clear()
            self._refresh_manager_table()
            self._refresh_manual_tree()

    # ------------------------------------------------------------------
    # PAGE 3: PROCESS & RESULTS
    # ------------------------------------------------------------------
    def _build_page_process(self, parent):
        page = tk.Frame(parent, bg=COLOR_BG)

        top = tk.Frame(page, bg=COLOR_BG)
        top.pack(fill="x")

        card_out = Card(top, title="Output Folder", subtitle="Where the Excel reports will be saved.")
        card_out.pack(side="left", fill="both", expand=True, padx=(0, 10))
        row = tk.Frame(card_out.body, bg=COLOR_CARD)
        row.pack(fill="x")
        make_pill_button(row, "Choose Folder...", self._on_choose_outdir, bg=COLOR_PRIMARY,
                          hover="#123A6B", font_size=9).pack(side="left")
        self.lbl_outdir = tk.Label(card_out.body, text="(not selected)", bg=COLOR_CARD, fg=COLOR_MUTED,
                                     font=(FONT_FAMILY, 9), wraplength=380, justify="left")
        self.lbl_outdir.pack(anchor="w", pady=(8, 0))

        card_run = Card(top, title="Run Analysis", subtitle="Processes every series currently loaded.")
        card_run.pack(side="left", fill="both", expand=True, padx=(10, 0))
        self.btn_process = make_pill_button(card_run.body, "\u25B6  Process All Series", self._on_process,
                                              bg=COLOR_ACCENT, hover=COLOR_ACCENT_DARK, state="disabled")
        self.btn_process.pack(anchor="w")
        self.lbl_process_summary = tk.Label(card_run.body, text="0 series ready.", bg=COLOR_CARD,
                                              fg=COLOR_MUTED, font=(FONT_FAMILY, 9))
        self.lbl_process_summary.pack(anchor="w", pady=(8, 0))

        progress_card = Card(page, title="Progress")
        progress_card.pack(fill="x", pady=(10, 10))
        self.progress = ttk.Progressbar(progress_card.body, mode="determinate")
        self.progress.pack(fill="x", pady=(0, 6))
        self.lbl_status = tk.Label(progress_card.body, text="Waiting to start.", bg=COLOR_CARD,
                                     fg=COLOR_MUTED, font=(FONT_FAMILY, 9), anchor="w")
        self.lbl_status.pack(anchor="w")

        log_card = Card(page, title="Processing Log")
        log_card.pack(fill="both", expand=True)
        self.txt_log = scrolledtext.ScrolledText(
            log_card.body, height=14, state="disabled", font=("Consolas", 9),
            wrap="word", bg="#0f172a", fg="#e2e8f0", insertbackground="#e2e8f0", relief="flat",
        )
        self.txt_log.pack(fill="both", expand=True)

        bottom = tk.Frame(page, bg=COLOR_BG)
        bottom.pack(fill="x", pady=(10, 0))
        self.btn_open_result = make_pill_button(bottom, "\U0001F5C2  Open Result Folder",
                                                   self._on_open_result, bg="#475569", hover="#334155",
                                                   state="disabled")
        self.btn_open_result.pack(side="left")

        return page

    def _refresh_process_summary(self):
        total = len(self.series_store)
        ready = sum(1 for df in self.series_store.values() if validate_series_df(df)[0])
        self.lbl_process_summary.configure(
            text=f"{ready} of {total} loaded series are ready to process." if total else
            "No series loaded yet - go to 'Data Source' first.")
        self._update_process_button()

    def _on_choose_outdir(self):
        d = filedialog.askdirectory(title="Select Output Folder")
        if not d:
            return
        self.out_dir = d
        self.lbl_outdir.configure(text=d)
        self._update_process_button()

    def _update_process_button(self):
        ready = any(validate_series_df(df)[0] for df in self.series_store.values())
        if self.out_dir and ready:
            self.btn_process.set_state("normal")
        else:
            self.btn_process.set_state("disabled")

    def _append_log_direct(self, msg):
        self.txt_log.configure(state="normal")
        self.txt_log.insert("end", str(msg) + "\n")
        self.txt_log.see("end")
        self.txt_log.configure(state="disabled")

    def _clear_log(self):
        self.txt_log.configure(state="normal")
        self.txt_log.delete("1.0", "end")
        self.txt_log.configure(state="disabled")

    def _log(self, msg):
        self.msg_queue.put(("log", msg))

    def _set_progress(self, done, total):
        self.msg_queue.put(("progress", (done, total)))

    def _poll_queue(self):
        try:
            while True:
                kind, payload = self.msg_queue.get_nowait()
                if kind == "log":
                    self._append_log_direct(payload)
                elif kind == "progress":
                    done, total = payload
                    self.progress["maximum"] = max(total, 1)
                    self.progress["value"] = done
                    self.lbl_status.configure(text=f"Processing... ({done}/{total} series)")
                elif kind == "done":
                    success, failed, out_dir = payload
                    self._on_process_finished(success, failed, out_dir)
        except queue.Empty:
            pass
        self.after(100, self._poll_queue)

    def _on_process(self):
        if not self.out_dir:
            return
        ready_series = {name: df for name, df in self.series_store.items() if validate_series_df(df)[0]}
        if not ready_series:
            messagebox.showwarning(APP_TITLE, "No series are ready to process (each series needs at least 5 rows).")
            return
        if self.worker_thread and self.worker_thread.is_alive():
            return

        self.btn_process.set_state("disabled")
        self.btn_open_result.set_state("disabled")
        self.progress["value"] = 0
        self.lbl_status.configure(text="Processing...")
        self._clear_log()
        self._append_log_direct(">>> Starting analysis for all ready series...\n")

        def worker():
            try:
                success, failed = run_batch_from_series(ready_series, self.out_dir, self._log, self._set_progress)
                self.msg_queue.put(("done", (success, failed, self.out_dir)))
            except Exception:
                err_text = traceback.format_exc()
                self._log(err_text)
                self.msg_queue.put(("done", ([], [("FATAL", "See the Processing Log panel for details")], self.out_dir)))

        self.worker_thread = threading.Thread(target=worker, daemon=True)
        self.worker_thread.start()

    def _on_process_finished(self, success, failed, out_dir):
        self.btn_process.set_state("normal")
        self.btn_open_result.set_state("normal")

        if failed and not success:
            self.lbl_status.configure(text="Finished with errors.")
            messagebox.showerror(
                APP_TITLE,
                f"Processing finished with errors.\n\nSucceeded: {len(success)}\nFailed: {len(failed)}\n\n"
                "See the 'Processing Log' panel for details."
            )
        elif failed:
            self.lbl_status.configure(text="Finished (some series failed).")
            messagebox.showwarning(
                APP_TITLE,
                f"Processing finished.\n\nSucceeded: {len(success)}\nFailed: {len(failed)}\n\n"
                "See the 'Processing Log' panel for details."
            )
        else:
            self.lbl_status.configure(text="Finished. All series processed successfully.")
            messagebox.showinfo(
                APP_TITLE,
                f"Processing complete!\n\nSucceeded: {len(success)} series\n\nResults saved to:\n{out_dir}"
            )

    def _on_open_result(self):
        if self.out_dir:
            _open_in_explorer(self.out_dir)

    # ------------------------------------------------------------------
    def _show_about(self):
        messagebox.showinfo(
            APP_TITLE,
            f"{APP_TITLE} v{APP_VERSION}\n\n"
            "Rainfall / annual discharge frequency analysis using 6 distributions:\n"
            "Normal, Log Normal, 3-Parameter Log Normal, Gumbel,\n"
            "Pearson III, Log Pearson III, and GEV.\n\n"
            "Data can be imported from a multi-sheet Excel workbook, typed\n"
            "manually inside the app, or both combined.\n\n"
            f"{CREDENTIAL_TEXT}"
        )

    def _on_close(self):
        if self.worker_thread and self.worker_thread.is_alive():
            if not messagebox.askyesno(
                APP_TITLE,
                "A process is still running. Are you sure you want to close the application?\n"
                "(The unfinished process will be force-stopped.)"
            ):
                return
        self.destroy()


def main():
    import tkinter.simpledialog  # noqa: F401 - registers tk.simpledialog
    app = RainfallApp()
    app.mainloop()


if __name__ == "__main__":
    main()
