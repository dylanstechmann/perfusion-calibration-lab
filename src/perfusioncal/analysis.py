"""Fit mass versus elapsed time; summarize independent repeat runs."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path

import numpy as np


def fit_line(time_s, mass_mg):
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
    return {"mass_rate_mg_s": slope, "intercept_mg": intercept,
            "r_squared": None if ss_total == 0 else float(1 - (residual @ residual) / ss_total),
            "residual_rmse_mg": float(np.sqrt(np.mean(residual ** 2))),
            "n_readings": len(time_s), "duration_s": float(time_s[-1] - time_s[0])}


def analyze(path, *, density_mg_ul, discard_seconds=0.0, seed=0, bootstrap_draws=2000):
    if not np.isfinite(density_mg_ul) or density_mg_ul <= 0:
        raise ValueError("density_mg_ul must be finite and positive")
    if not np.isfinite(discard_seconds) or discard_seconds < 0:
        raise ValueError("discard_seconds must be finite and nonnegative")
    if not isinstance(bootstrap_draws, int) or bootstrap_draws < 100:
        raise ValueError("use at least 100 bootstrap draws")
    path = Path(path)
    runs = {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
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
        fit = fit_line(selected[:, 0], selected[:, 1])
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
        result.append({"run_id": run_id, **fit, "target_flow_ul_min": target,
                       "measured_flow_ul_min": flow, "error_percent": 100 * (flow - target) / target,
                       "late_vs_early_slope_change_percent": drift, "flags": flags})
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
                          "flagged_runs": [r["run_id"] for r in subset if r["flags"]]})
    return {"schema_version": 1, "input_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "configuration": {"density_mg_ul": density_mg_ul, "discard_seconds": discard_seconds,
                              "seed": seed, "bootstrap_draws": bootstrap_draws},
            "runs": result, "targets": summaries,
            "notes": ["Offline research analysis. No hardware commands are generated.",
                      "Intervals resample independent run slopes, not serially correlated readings.",
                      "Fewer than three runs: no interval. Small repeat counts give unstable intervals.",
                      "Density, balance calibration, evaporation and collection losses are not included in uncertainty.",
                      "R-squared <0.95 and >20% half-run drift are diagnostic flags, not acceptance standards.",
                      "Flagged runs remain in summaries; inspect them before interpreting mean flow."]}


def markdown(report):
    lines = ["# Perfusion calibration analysis", "", f"Input SHA-256: `{report['input_sha256']}`", "",
             "| Target (µL/min) | Runs | Mean measured (µL/min) | Error (%) | Repeat SD |",
             "|---:|---:|---:|---:|---:|"]
    for row in report["targets"]:
        sd = "unavailable" if row["sd_flow_ul_min"] is None else f"{row['sd_flow_ul_min']:.3f}"
        lines.append(f"| {row['target_flow_ul_min']:.3f} | {row['n_runs']} | {row['mean_flow_ul_min']:.3f} | "
                     f"{row['mean_error_percent']:.2f} | {sd} |")
    lines.extend(["", "## Run diagnostics", ""])
    for run in report["runs"]:
        lines.append(f"- {run['run_id']}: {', '.join(run['flags']) or 'no diagnostic flags'}")
    lines.extend(["", "## Interpretation", "", *[f"- {n}" for n in report["notes"]], ""])
    return "\n".join(lines)
