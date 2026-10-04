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


_UNCERTAINTY_COMPONENTS = {
    "density": "mg/uL",
    "balance_gain": "relative_fraction",
    "balance_slope_drift": "mg/s",
    "evaporation_rate": "mg/s",
    "timing_scale": "relative_fraction",
}


def validate_uncertainty_budget(budget):
    """Validate a source-linked, independent-standard-uncertainty budget."""
    if not isinstance(budget, dict) or set(budget) != {
        "schema_version", "combination_assumption", "components"
    }:
        raise ValueError("uncertainty budget needs schema_version, combination_assumption, and components")
    if type(budget["schema_version"]) is not int or budget["schema_version"] != 1:
        raise ValueError("unsupported uncertainty budget schema_version")
    if budget["combination_assumption"] != "independent_standard_uncertainties":
        raise ValueError("only independent_standard_uncertainties can be combined")
    components = budget["components"]
    if not isinstance(components, dict) or set(components) != set(_UNCERTAINTY_COMPONENTS):
        raise ValueError("uncertainty budget must declare all five named components")

    normalized = {}
    for name, unit in _UNCERTAINTY_COMPONENTS.items():
        component = components[name]
        if not isinstance(component, dict):
            raise ValueError(f"uncertainty component {name} must be an object")
        status = component.get("status")
        if status == "not_available":
            if set(component) != {"status", "reason"}:
                raise ValueError(f"unavailable component {name} needs only status and reason")
            reason = component.get("reason")
            if not isinstance(reason, str) or not reason.strip():
                raise ValueError(f"unavailable component {name} needs a nonblank reason")
            normalized[name] = {"status": status, "unit": unit, "reason": reason.strip()}
            continue
        expected_keys = {"status", "standard_uncertainty", "source_record"}
        if name != "density":
            expected_keys.add("correction")
        if status != "measured" or set(component) != expected_keys:
            raise ValueError(f"component {name} must be measured with uncertainty and source, or not_available")
        value = component["standard_uncertainty"]
        source = component["source_record"]
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not np.isfinite(value) or value <= 0:
            raise ValueError(f"component {name} standard_uncertainty must be finite and positive")
        if not isinstance(source, str) or not source.strip():
            raise ValueError(f"measured component {name} needs a nonblank source_record")
        correction = component.get("correction")
        if name != "density":
            if (isinstance(correction, bool) or not isinstance(correction, (int, float))
                    or not np.isfinite(correction)):
                raise ValueError(f"component {name} correction must be finite")
            if name in {"balance_gain", "timing_scale"} and correction <= -1:
                raise ValueError(f"component {name} correction must keep its scale factor positive")
            if name == "evaporation_rate" and correction < 0:
                raise ValueError("evaporation_rate correction must be nonnegative")
        normalized[name] = {
            "status": status,
            "unit": unit,
            "standard_uncertainty": float(value),
            "source_record": source.strip(),
        }
        if name != "density":
            normalized[name]["correction"] = float(correction)
    return {
        "schema_version": 1,
        "combination_assumption": "independent_standard_uncertainties",
        "components": normalized,
    }


def _flow_after_corrections(mass_rate_mg_s, density_mg_ul, budget_components):
    balance_bias = budget_components["balance_slope_drift"].get("correction", 0.0)
    evaporation = budget_components["evaporation_rate"].get("correction", 0.0)
    balance_gain = 1.0 + budget_components["balance_gain"].get("correction", 0.0)
    timing_scale = 1.0 + budget_components["timing_scale"].get("correction", 0.0)
    apparent_flow = (mass_rate_mg_s - balance_bias) * 60.0 / density_mg_ul
    evaporation_result = flow_if_constant_evaporation(apparent_flow, evaporation, density_mg_ul)
    # balance_gain is indicated mass / true mass, so remove its bias by dividing;
    # timing_scale is indicated elapsed time / true elapsed time, so multiply.
    corrected_flow = evaporation_result["implied_delivered_flow_ul_min"] * timing_scale / balance_gain
    applied = [
        name for name in ("balance_slope_drift", "evaporation_rate", "balance_gain", "timing_scale")
        if budget_components[name]["status"] == "measured"
    ]
    return float(corrected_flow), applied


def _flow_uncertainty_budget(corrected_flow_ul_min, density_mg_ul, budget_components):
    """Propagate listed components without mixing them into the run bootstrap."""
    contributions = {}
    complete = True
    balance_gain = 1.0 + budget_components["balance_gain"].get("correction", 0.0)
    timing_scale = 1.0 + budget_components["timing_scale"].get("correction", 0.0)
    flow_without_gain = corrected_flow_ul_min / balance_gain
    flow_without_timing = corrected_flow_ul_min / timing_scale
    rate_sensitivity = 60.0 * balance_gain * timing_scale / density_mg_ul
    for name, component in budget_components.items():
        if component["status"] != "measured":
            complete = False
            contributions[name] = None
            continue
        standard_uncertainty = component["standard_uncertainty"]
        if name == "density":
            contribution = abs(corrected_flow_ul_min) * standard_uncertainty / density_mg_ul
        elif name == "balance_gain":
            contribution = abs(flow_without_gain) * standard_uncertainty
        elif name == "timing_scale":
            contribution = abs(flow_without_timing) * standard_uncertainty
        else:
            contribution = rate_sensitivity * standard_uncertainty
        contributions[name] = float(contribution)
    combined = None
    if complete:
        combined = float(math.sqrt(sum(value * value for value in contributions.values())))
    if not any(component["status"] == "measured" for component in budget_components.values()):
        status = "unavailable"
    elif complete:
        status = "quantified_for_declared_components"
    else:
        status = "partial"
    return {
        "status": status,
        "combined_standard_uncertainty_ul_min": combined,
        "component_contributions_ul_min": contributions,
        "unquantified_components": [
            name for name, component in budget_components.items() if component["status"] != "measured"
        ],
        "combination_assumption": "independent_standard_uncertainties",
        "note": (
            "First-order contributions around the flow after available source-linked corrections. "
            "This is separate from run resampling and excludes unlisted loss mechanisms."
        ),
        "flow_after_corrections_ul_min": float(corrected_flow_ul_min),
        "corrections_applied": [
            name for name in ("balance_slope_drift", "evaporation_rate", "balance_gain", "timing_scale")
            if budget_components[name]["status"] == "measured"
        ],
    }


def analyze(path, *, density_mg_ul, discard_seconds=0.0, seed=0, bootstrap_draws=2000,
            outlier_method="iqr", outlier_threshold=None, uncertainty_budget=None,
            uncertainty_budget_sha256=None):
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
    normalized_budget = (
        None if uncertainty_budget is None else validate_uncertainty_budget(uncertainty_budget)
    )
    if uncertainty_budget_sha256 is not None and (
        not isinstance(uncertainty_budget_sha256, str)
        or len(uncertainty_budget_sha256) != 64
        or any(char not in "0123456789abcdef" for char in uncertainty_budget_sha256.lower())
    ):
        raise ValueError("uncertainty_budget_sha256 must be a SHA-256 hex digest")
    budget_components = (
        normalized_budget["components"] if normalized_budget is not None else {
            name: {"status": "not_available", "unit": unit,
                   "reason": "No measurement uncertainty budget was supplied."}
            for name, unit in _UNCERTAINTY_COMPONENTS.items()
        }
    )
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
        corrected_flow, corrections_applied = _flow_after_corrections(
            fit["mass_rate_mg_s"], density_mg_ul, budget_components
        )
        drift = None
        if len(selected) >= 6 and abs(fit["mass_rate_mg_s"]) > 1e-12:
            half = len(selected) // 2
            early = fit_line(selected[:half, 0], selected[:half, 1])["mass_rate_mg_s"]
            late = fit_line(selected[half:, 0], selected[half:, 1])["mass_rate_mg_s"]
            drift = 100 * (late - early) / abs(fit["mass_rate_mg_s"])
        flags = []
        if flow <= 0:
            flags.append("nonpositive_measured_flow")
        if corrected_flow <= 0 < flow:
            flags.append("nonpositive_flow_after_available_corrections")
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
                       "flow_after_available_corrections_ul_min": corrected_flow,
                       "corrections_applied": corrections_applied,
                       "late_vs_early_slope_change_percent": drift,
                       "outlier_times_s": outlier_times, "outlier_residuals_mg": outlier_residuals,
                       "flags": flags})
    summaries = []
    rng = np.random.default_rng(seed)
    for target in sorted({r["target_flow_ul_min"] for r in result}):
        subset = [r for r in result if r["target_flow_ul_min"] == target]
        flows = np.array([r["measured_flow_ul_min"] for r in subset])
        corrected_flows = np.array([r["flow_after_available_corrections_ul_min"] for r in subset])
        interval = None
        corrected_interval = None
        if len(flows) >= 3:
            means = np.array([rng.choice(flows, size=len(flows), replace=True).mean()
                              for _ in range(bootstrap_draws)])
            interval = np.quantile(means, [0.025, 0.975]).tolist()
            corrected_means = np.array([rng.choice(corrected_flows, size=len(corrected_flows), replace=True).mean()
                                        for _ in range(bootstrap_draws)])
            corrected_interval = np.quantile(corrected_means, [0.025, 0.975]).tolist()
        mean_corrected_flow = float(corrected_flows.mean())
        correction_count = sum(
            budget_components[name]["status"] == "measured"
            for name in ("balance_gain", "balance_slope_drift", "evaporation_rate", "timing_scale")
        )
        correction_status = (
            "none_available" if correction_count == 0
            else "all_declared_corrections_available" if correction_count == 4
            else "partial"
        )
        summaries.append({"target_flow_ul_min": target, "n_runs": len(flows),
                          "mean_flow_ul_min": float(flows.mean()),
                          "sd_flow_ul_min": None if len(flows) < 2 else float(flows.std(ddof=1)),
                          "mean_error_percent": 100 * float(flows.mean() - target) / target,
                          "mean_flow_ci95": interval,
                          "mean_flow_after_available_corrections_ul_min": mean_corrected_flow,
                          "mean_flow_after_available_corrections_ci95": corrected_interval,
                          "correction_status": correction_status,
                          "measurement_system_uncertainty": _flow_uncertainty_budget(
                              mean_corrected_flow, density_mg_ul, budget_components
                          ),
                          "flagged_runs": [r["run_id"] for r in subset if r["flags"]],
                          "outlier_runs": [r["run_id"] for r in subset if r["outlier_count"] > 0]})
    return {"schema_version": 1, "input_sha256": hashlib.sha256(input_bytes).hexdigest(),
            "configuration": {"density_mg_ul": density_mg_ul, "discard_seconds": discard_seconds,
                              "seed": seed, "bootstrap_draws": bootstrap_draws,
                              "outlier_method": outlier_method, "outlier_threshold": outlier_threshold},
            "measurement_uncertainty_budget": {
                "status": "not_provided" if normalized_budget is None else "provided",
                "input_sha256": uncertainty_budget_sha256,
                "combination_assumption": (
                    None if normalized_budget is None
                    else normalized_budget["combination_assumption"]
                ),
                "components": budget_components,
            },
            "runs": result, "targets": summaries,
            "notes": ["Offline research analysis. No hardware commands are generated.",
                      "This report alone is not a physical pump calibration or acceptance decision.",
                      "Intervals resample independent run slopes, not serially correlated readings.",
                      "Fewer than three runs: no interval. Small repeat counts give unstable intervals.",
                      "The run-bootstrap interval excludes density, balance, evaporation and timing uncertainty.",
                      "The optional uncertainty budget applies only source-linked corrections that are explicitly marked measured.",
                      "Correlated uncertainty components require covariance propagation; this tool only combines components declared independent.",
                      "Retained droplets, collection losses and other unlisted effects are not quantified.",
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
    budget = report["measurement_uncertainty_budget"]
    lines.extend(["", "## Measurement-system uncertainty", ""])
    lines.append(f"Budget status: {budget['status']}. Input SHA-256: `{budget['input_sha256'] or 'not supplied'}`.")
    if budget["combination_assumption"]:
        lines.append(f"Combination assumption: `{budget['combination_assumption']}`.")
    lines.append("Component contributions are separate from run-bootstrap intervals. The combined value is a standard uncertainty, not a 95% interval.")
    lines.append("")
    component_names = list(budget["components"])
    lines.append("| Target (µL/min) | " + " | ".join(component_names) + " | Combined standard uncertainty (µL/min) | Status | Unquantified |")
    lines.append("|---:|" + "---:|" * len(component_names) + "---:|---|---|")
    for row in report["targets"]:
        uncertainty = row["measurement_system_uncertainty"]
        combined = uncertainty["combined_standard_uncertainty_ul_min"]
        combined_text = "unavailable" if combined is None else f"{combined:.4f}"
        missing = ", ".join(uncertainty["unquantified_components"]) or "none"
        component_values = [
            "unavailable" if uncertainty["component_contributions_ul_min"][name] is None
            else f"{uncertainty['component_contributions_ul_min'][name]:.4f}"
            for name in component_names
        ]
        lines.append(
            f"| {row['target_flow_ul_min']:.3f} | " + " | ".join(component_values)
            + f" | {combined_text} | {uncertainty['status']} | {missing} |"
        )
    lines.extend(["", "## Flow after available corrections", "",
                  "Uncorrected measured flow remains in the main summary above. This table applies only source-linked corrections marked `measured`.",
                  "",
                  "| Target (µL/min) | Mean after available corrections (µL/min) | Correction status | 95% run-bootstrap interval after corrections (µL/min) |",
                  "|---:|---:|---|---:|"])
    for row in report["targets"]:
        interval = row["mean_flow_after_available_corrections_ci95"]
        ci = "unavailable (<3 runs)" if interval is None else f"[{interval[0]:.3f}, {interval[1]:.3f}]"
        lines.append(
            f"| {row['target_flow_ul_min']:.3f} | "
            f"{row['mean_flow_after_available_corrections_ul_min']:.3f} | "
            f"{row['correction_status']} | {ci} |"
        )
    lines.extend(["", "## Run diagnostics", ""])
    for run in report["runs"]:
        flag_str = ', '.join(run['flags']) or 'no diagnostic flags'
        lines.append(f"- {run['run_id']}: {flag_str}")
    lines.extend(["", "## Interpretation", "", *[f"- {n}" for n in report["notes"]], ""])
    return "\n".join(lines)
