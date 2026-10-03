"""Fit mass versus elapsed time; summarize independent repeat runs."""

from __future__ import annotations

import csv
import hashlib
import io
import math
from pathlib import Path

import numpy as np


def _regularized_beta(x: float, a: float, b: float) -> float:
    """Regularized incomplete beta using a continued fraction (Numerical Recipes)."""
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0

    def fraction(aa: float, bb: float, xx: float) -> float:
        qab, qap, qam = aa + bb, aa + 1.0, aa - 1.0
        c = 1.0
        d = 1.0 - qab * xx / qap
        d = 1e-300 if abs(d) < 1e-300 else d
        d = 1.0 / d
        result = d
        for m in range(1, 401):
            m2 = 2 * m
            term = m * (bb - m) * xx / ((qam + m2) * (aa + m2))
            d = 1.0 + term * d
            d = 1e-300 if abs(d) < 1e-300 else d
            c = 1.0 + term / c
            c = 1e-300 if abs(c) < 1e-300 else c
            d = 1.0 / d
            result *= d * c
            term = -(aa + m) * (qab + m) * xx / ((aa + m2) * (qap + m2))
            d = 1.0 + term * d
            d = 1e-300 if abs(d) < 1e-300 else d
            c = 1.0 + term / c
            c = 1e-300 if abs(c) < 1e-300 else c
            d = 1.0 / d
            delta = d * c
            result *= delta
            if abs(delta - 1.0) <= 3e-14:
                return result
        raise ArithmeticError("Student-t critical-value continued fraction did not converge")

    bt = math.exp(math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
                  + a * math.log(x) + b * math.log1p(-x))
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * fraction(a, b, x) / a
    return 1.0 - bt * fraction(b, a, 1.0 - x) / b


def _student_t_quantile(probability: float, degrees_of_freedom: int) -> float:
    if not 0.5 < probability < 1 or degrees_of_freedom < 1:
        raise ValueError("Student-t quantile requires probability in (0.5, 1) and positive degrees of freedom")
    low, high = 0.0, 1.0

    def cdf(value: float) -> float:
        x = degrees_of_freedom / (degrees_of_freedom + value * value)
        return 1.0 - 0.5 * _regularized_beta(x, degrees_of_freedom / 2.0, 0.5)

    while cdf(high) < probability:
        high *= 2.0
        if not math.isfinite(high):
            raise ArithmeticError("could not bracket Student-t quantile")
    for _ in range(100):
        middle = (low + high) / 2.0
        if cdf(middle) < probability:
            low = middle
        else:
            high = middle
    return (low + high) / 2.0


def _grubbs_critical_value(n: int, alpha: float) -> float:
    if n < 3:
        raise ValueError("Grubbs' test requires at least three observations")
    if isinstance(alpha, bool) or not isinstance(alpha, (int, float)) or not np.isfinite(alpha) or not 0 < alpha < 1:
        raise ValueError("Grubbs alpha must be strictly between 0 and 1")
    degrees_of_freedom = n - 2
    t_critical = _student_t_quantile(1.0 - alpha / (2.0 * n), degrees_of_freedom)
    return ((n - 1.0) / np.sqrt(n)) * np.sqrt(t_critical**2 / (degrees_of_freedom + t_critical**2))


def detect_outliers_iqr(residuals: np.ndarray, factor: float = 1.5) -> list[int]:
    """Identify indices of residual outliers using Tukey's fences (IQR method)."""
    residuals = np.asarray(residuals, dtype=float)
    if residuals.ndim != 1 or not np.isfinite(residuals).all():
        raise ValueError("residuals must be a finite one-dimensional array")
    if isinstance(factor, bool) or not isinstance(factor, (int, float)) or not np.isfinite(factor) or factor < 0:
        raise ValueError("IQR factor must be a finite nonnegative number")
    if len(residuals) < 4:
        return []
    q25, q75 = np.percentile(residuals, [25, 75])
    iqr = q75 - q25
    if iqr == 0:
        std = float(np.std(residuals))
        if std == 0:
            return []
        z = np.abs(residuals - np.median(residuals)) / std
        return [int(i) for i in np.where(z > 3.0)[0]]
    lower = q25 - factor * iqr
    upper = q75 + factor * iqr
    outlier_mask = (residuals < lower) | (residuals > upper)
    return [int(i) for i in np.where(outlier_mask)[0]]


def detect_outliers_grubbs(residuals: np.ndarray, alpha: float = 0.05) -> list[int]:
    """Identify indices of residual outliers using Grubbs' maximum normalized residual test."""
    residuals = np.asarray(residuals, dtype=float)
    if residuals.ndim != 1 or not np.isfinite(residuals).all():
        raise ValueError("residuals must be a finite one-dimensional array")
    if isinstance(alpha, bool) or not isinstance(alpha, (int, float)) or not np.isfinite(alpha) or not 0 < alpha < 1:
        raise ValueError("Grubbs alpha must be strictly between 0 and 1")
    n = len(residuals)
    if n < 4:
        return []
    std = float(np.std(residuals, ddof=1))
    if std == 0:
        return []
    mean = float(np.mean(residuals))
    deviations = np.abs(residuals - mean)
    max_idx = int(np.argmax(deviations))
    g_stat = deviations[max_idx] / std
    g_crit = _grubbs_critical_value(n, alpha)
    if g_stat > g_crit:
        return [max_idx]
    return []


def detect_outliers(residuals: np.ndarray, method: str = "iqr", threshold: float | None = None) -> list[int]:
    if not isinstance(method, str):
        raise ValueError("outlier method must be 'iqr' or 'grubbs'")
    method = method.lower()
    if method == "iqr":
        factor = 1.5 if threshold is None else threshold
        return detect_outliers_iqr(residuals, factor=factor)
    elif method == "grubbs":
        alpha = 0.05 if threshold is None else threshold
        return detect_outliers_grubbs(residuals, alpha=alpha)
    else:
        raise ValueError(f"unknown outlier method: '{method}'; choose 'iqr' or 'grubbs'")


def fit_line(time_s, mass_mg, outlier_method: str = "iqr", outlier_threshold: float | None = None):
    time_s = np.asarray(time_s, dtype=float)
    mass_mg = np.asarray(mass_mg, dtype=float)
    if time_s.ndim != 1 or mass_mg.shape != time_s.shape or len(time_s) < 3:
        raise ValueError("each fitted run needs at least three time/mass readings")
    if not np.isfinite(time_s).all() or not np.isfinite(mass_mg).all():
        raise ValueError("time and mass must be finite")
    if np.any(time_s < 0) or np.any(np.diff(time_s) <= 0):
        raise ValueError("times must be nonnegative and strictly increasing within a run")
    centered = time_s - time_s.mean()
    slope = float(centered @ (mass_mg - mass_mg.mean()) / (centered @ centered))
    intercept = float(mass_mg.mean() - slope * time_s.mean())
    residual = mass_mg - (intercept + slope * time_s)
    ss_total = float(np.sum((mass_mg - mass_mg.mean()) ** 2))
    outlier_indices = detect_outliers(residual, method=outlier_method, threshold=outlier_threshold)
    return {
        "mass_rate_mg_s": slope,
        "intercept_mg": intercept,
        "r_squared": None if ss_total == 0 else float(1 - (residual @ residual) / ss_total),
        "residual_rmse_mg": float(np.sqrt(np.mean(residual ** 2))),
        "n_readings": len(time_s),
        "duration_s": float(time_s[-1] - time_s[0]),
        "outlier_indices": outlier_indices,
        "outlier_count": len(outlier_indices),
    }



def flow_if_constant_evaporation(apparent_flow_ul_min, evaporation_mg_s, density_mg_ul):
    """Shift an apparent flow by a constant mass-loss rate.

    This is a hypothetical sensitivity. It is not a measured evaporation rate
    and it is not part of the bootstrap interval.
    """
    values = (apparent_flow_ul_min, evaporation_mg_s, density_mg_ul)
    if any(isinstance(value, bool) or not isinstance(value, (int, float)) or not np.isfinite(value) for value in values):
        raise ValueError("flow, evaporation, and density must be finite numbers")
    if density_mg_ul <= 0:
        raise ValueError("density_mg_ul must be positive")
    correction_ul_min = float(evaporation_mg_s) * 60.0 / float(density_mg_ul)
    return {
        "apparent_flow_ul_min": float(apparent_flow_ul_min),
        "evaporation_mg_s": float(evaporation_mg_s),
        "density_mg_ul": float(density_mg_ul),
        "correction_ul_min": correction_ul_min,
        "implied_delivered_flow_ul_min": float(apparent_flow_ul_min) + correction_ul_min,
        "note": "Constant evaporation is a hypothetical sensitivity, not a measured loss or a calibration.",
    }


def analyze(path, *, density_mg_ul, discard_seconds=0.0, seed=0, bootstrap_draws=2000,
            outlier_method="iqr", outlier_threshold=None):
    if isinstance(density_mg_ul, bool) or not isinstance(density_mg_ul, (int, float)) or not np.isfinite(density_mg_ul) or density_mg_ul <= 0:
        raise ValueError("density_mg_ul must be finite and positive")
    if isinstance(discard_seconds, bool) or not isinstance(discard_seconds, (int, float)) or not np.isfinite(discard_seconds) or discard_seconds < 0:
        raise ValueError("discard_seconds must be finite and nonnegative")
    if outlier_threshold is not None:
        if isinstance(outlier_threshold, bool) or not isinstance(outlier_threshold, (int, float)) or not np.isfinite(outlier_threshold):
            raise ValueError("outlier_threshold must be a finite number")
        valid_threshold = (outlier_threshold >= 0 if outlier_method == "iqr"
                           else 0 < outlier_threshold < 1 if outlier_method == "grubbs" else False)
        if not valid_threshold:
            raise ValueError("IQR threshold must be nonnegative; Grubbs alpha must lie strictly between 0 and 1")
    if not isinstance(bootstrap_draws, int) or bootstrap_draws < 100:
        raise ValueError("use at least 100 bootstrap draws")
    path = Path(path)
    # Parse and hash one snapshot so the reported digest identifies the fitted data.
    input_bytes = path.read_bytes()
    runs = {}
    with io.StringIO(input_bytes.decode("utf-8-sig"), newline="") as handle:
        reader = csv.DictReader(handle)
        expected = {"run_id", "time_s", "mass_mg", "target_flow_ul_min"}
        header = reader.fieldnames or []
        if len(header) != len(set(header)) or not expected.issubset(header):
            raise ValueError("CSV needs unique run_id,time_s,mass_mg,target_flow_ul_min columns")
        for line, row in enumerate(reader, 2):
            if None in row or any(v is None for v in row.values()) or not row["run_id"].strip():
                raise ValueError(f"row {line}: missing fields or run_id")
            try:
                values = [float(row[c]) for c in ["time_s", "mass_mg", "target_flow_ul_min"]]
            except ValueError as exc:
                raise ValueError(f"row {line}: invalid number") from exc
            if not np.isfinite(values).all() or values[0] < 0 or values[2] <= 0:
                raise ValueError(f"row {line}: require finite numbers, nonnegative time, positive target")
            runs.setdefault(row["run_id"].strip(), []).append(values)
    if not runs:
        raise ValueError("no runs in CSV")
    result = []
    for run_id, values in sorted(runs.items()):
        array = np.array(values)
        if np.any(np.diff(array[:, 0]) <= 0):
            raise ValueError(f"{run_id}: duplicate or nonincreasing times")
        if len(set(array[:, 2])) != 1:
            raise ValueError(f"{run_id}: target flow changes within a run; separate runs first")
        selected = array[array[:, 0] - array[0, 0] >= discard_seconds]
        fit = fit_line(selected[:, 0], selected[:, 1], outlier_method=outlier_method, outlier_threshold=outlier_threshold)
        target = float(array[0, 2])
        flow = fit["mass_rate_mg_s"] * 60.0 / density_mg_ul
        drift = None
        if len(selected) >= 6 and abs(fit["mass_rate_mg_s"]) > 1e-12:
            half = len(selected) // 2
            early = fit_line(selected[:half, 0], selected[:half, 1])["mass_rate_mg_s"]
            late = fit_line(selected[half:, 0], selected[half:, 1])["mass_rate_mg_s"]
            drift = 100 * (late - early) / abs(fit["mass_rate_mg_s"])
        flags = []
        if flow <= 0:
            flags.append("nonpositive_measured_flow")
        if fit["r_squared"] is None or fit["r_squared"] < 0.95:
            flags.append("inspect_nonlinearity_or_noise")
        if drift is not None and abs(drift) > 20:
            flags.append("inspect_within_run_drift")
        if fit["outlier_count"] > 0:
            flags.append(f"outlier_readings_detected({fit['outlier_count']})")
        outlier_times = [float(selected[i, 0]) for i in fit["outlier_indices"]]
        outlier_residuals = [round(float(selected[i, 1] - (fit["intercept_mg"] + fit["mass_rate_mg_s"] * selected[i, 0])), 4)
                             for i in fit["outlier_indices"]]
        result.append({"run_id": run_id, **fit, "target_flow_ul_min": target,
                       "measured_flow_ul_min": flow, "error_percent": 100 * (flow - target) / target,
                       "late_vs_early_slope_change_percent": drift,
                       "outlier_times_s": outlier_times, "outlier_residuals_mg": outlier_residuals,
                       "flags": flags})
    summaries = []
    rng = np.random.default_rng(seed)
    for target in sorted({r["target_flow_ul_min"] for r in result}):
        subset = [r for r in result if r["target_flow_ul_min"] == target]
        flows = np.array([r["measured_flow_ul_min"] for r in subset])
        interval = None
        if len(flows) >= 3:
            means = np.array([rng.choice(flows, size=len(flows), replace=True).mean()
                              for _ in range(bootstrap_draws)])
            interval = np.quantile(means, [0.025, 0.975]).tolist()
        summaries.append({"target_flow_ul_min": target, "n_runs": len(flows),
                          "mean_flow_ul_min": float(flows.mean()),
                          "sd_flow_ul_min": None if len(flows) < 2 else float(flows.std(ddof=1)),
                          "mean_error_percent": 100 * float(flows.mean() - target) / target,
                          "mean_flow_ci95": interval,
                          "flagged_runs": [r["run_id"] for r in subset if r["flags"]],
                          "outlier_runs": [r["run_id"] for r in subset if r["outlier_count"] > 0]})
    return {"schema_version": 1, "input_sha256": hashlib.sha256(input_bytes).hexdigest(),
            "configuration": {"density_mg_ul": density_mg_ul, "discard_seconds": discard_seconds,
                              "seed": seed, "bootstrap_draws": bootstrap_draws,
                              "outlier_method": outlier_method, "outlier_threshold": outlier_threshold},
            "runs": result, "targets": summaries,
            "notes": ["Offline research analysis. No hardware commands are generated.",
                      "This report alone is not a physical pump calibration or acceptance decision.",
                      "Intervals resample independent run slopes, not serially correlated readings.",
                      "Fewer than three runs: no interval. Small repeat counts give unstable intervals.",
                      "Density, balance calibration, evaporation and collection losses are not included in uncertainty.",
                      "R-squared <0.95 and >20% half-run drift are diagnostic flags, not acceptance standards.",
                      "IQR and Grubbs outlier flags do not alter regression weights; Grubbs' independent-normal assumptions may not hold for correlated linear-fit residuals.",
                      "Flagged runs remain in summaries; inspect them before interpreting mean flow."]}


def markdown(report):
    lines = ["# Perfusion calibration analysis", "", f"Input SHA-256: `{report['input_sha256']}`", "",
             f"Density used: {report['configuration']['density_mg_ul']:g} mg/µL. "
             f"Startup interval discarded: {report['configuration']['discard_seconds']:g} s. "
             f"Outlier method: {report['configuration'].get('outlier_method', 'iqr')}.", "",
             "| Target (µL/min) | Runs | Mean measured (µL/min) | Error (%) | Repeat SD (µL/min) | 95% run-bootstrap interval (µL/min) |",
             "|---:|---:|---:|---:|---:|---:|"]
    for row in report["targets"]:
        sd = "unavailable" if row["sd_flow_ul_min"] is None else f"{row['sd_flow_ul_min']:.3f}"
        interval = row["mean_flow_ci95"]
        ci = "unavailable (<3 runs)" if interval is None else f"[{interval[0]:.3f}, {interval[1]:.3f}]"
        lines.append(f"| {row['target_flow_ul_min']:.3f} | {row['n_runs']} | {row['mean_flow_ul_min']:.3f} | "
                     f"{row['mean_error_percent']:.2f} | {sd} | {ci} |")
    lines.extend(["", "## Run diagnostics", ""])
    for run in report["runs"]:
        flag_str = ', '.join(run['flags']) or 'no diagnostic flags'
        lines.append(f"- {run['run_id']}: {flag_str}")
    lines.extend(["", "## Interpretation", "", *[f"- {n}" for n in report["notes"]], ""])
    return "\n".join(lines)
