import copy
import csv
import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

from perfusioncal.analysis import (
    analyze, flow_if_constant_evaporation, fit_line, markdown,
    validate_uncertainty_budget,
)
from perfusioncal.cli import main, write_demo


class CalibrationTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.path = Path(self.tmp.name) / "measurements.csv"

    def write_rows(self, rows):
        with self.path.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["run_id", "time_s", "mass_mg", "target_flow_ul_min"])
            writer.writerows(rows)

    @staticmethod
    def complete_budget():
        return {
            "schema_version": 1,
            "combination_assumption": "independent_standard_uncertainties",
            "components": {
                "density": {"status": "measured", "standard_uncertainty": 0.01,
                            "source_record": "density-certificate"},
                "balance_gain": {"status": "measured", "standard_uncertainty": 0.001,
                                 "source_record": "balance-certificate", "correction": 0.002},
                "balance_slope_drift": {"status": "measured", "standard_uncertainty": 0.001,
                                        "source_record": "blank-study", "correction": 0.1},
                "evaporation_rate": {"status": "measured", "standard_uncertainty": 0.002,
                                     "source_record": "matched-blank", "correction": 0.2},
                "timing_scale": {"status": "measured", "standard_uncertainty": 0.001,
                                 "source_record": "timebase-certificate", "correction": 0.003},
            },
        }

    def test_uncertainty_budget_requires_source_linked_finite_components(self):
        root = Path(__file__).resolve().parents[1]
        example = json.loads((root / "examples/measurement_uncertainty_budget.example.json").read_text())
        normalized = validate_uncertainty_budget(example)
        self.assertEqual(normalized["schema_version"], 1)
        self.assertEqual(normalized["components"]["density"]["status"], "not_available")

        bad = self.complete_budget()
        bad["components"]["balance_gain"]["source_record"] = "  "
        with self.assertRaisesRegex(ValueError, "source_record"):
            validate_uncertainty_budget(bad)
        bad = self.complete_budget()
        bad["components"]["timing_scale"]["correction"] = -1
        with self.assertRaisesRegex(ValueError, "positive"):
            validate_uncertainty_budget(bad)
        bad = self.complete_budget()
        bad["components"]["density"]["standard_uncertainty"] = True
        with self.assertRaisesRegex(ValueError, "finite and positive"):
            validate_uncertainty_budget(bad)
        bad = self.complete_budget()
        bad["combination_assumption"] = "correlated"
        with self.assertRaisesRegex(ValueError, "only independent"):
            validate_uncertainty_budget(bad)

    def test_uncertainty_budget_applies_corrections_and_combines_complete_components(self):
        self.write_rows([["one", t, 2 * t, 60] for t in [0, 1, 2, 3]])
        budget = self.complete_budget()
        report = analyze(self.path, density_mg_ul=2, uncertainty_budget=budget)
        target = report["targets"][0]
        # (2.0 - 0.1 + 0.2) * 60/2 * 1.003 / 1.002
        expected = (2.0 - 0.1 + 0.2) * 30 * 1.003 / 1.002
        self.assertAlmostEqual(target["mean_flow_after_available_corrections_ul_min"], expected)
        uncertainty = target["measurement_system_uncertainty"]
        self.assertEqual(uncertainty["status"], "quantified_for_declared_components")
        self.assertAlmostEqual(
            uncertainty["combined_standard_uncertainty_ul_min"],
            np.sqrt(sum(value ** 2 for value in uncertainty["component_contributions_ul_min"].values())),
        )
        self.assertEqual(len(uncertainty["corrections_applied"]), 4)
        self.assertIsNone(target["mean_flow_ci95"])

    def test_uncertainty_contributions_match_independent_finite_differences(self):
        # A substantial non-unit gain makes multiply/divide errors visible.
        self.write_rows([["one", t, 14 + 2 * t, 60] for t in [0, 1, 2, 3]])
        density = 2.0
        epsilon = 1e-6
        for gain_correction in [0.25, -0.4]:
            budget = self.complete_budget()
            budget["components"]["balance_gain"]["correction"] = gain_correction
            report = analyze(self.path, density_mg_ul=density, uncertainty_budget=budget)
            contributions = report["targets"][0]["measurement_system_uncertainty"]["component_contributions_ul_min"]
            for name, component in budget["components"].items():
                flows = []
                for direction in [-1, 1]:
                    perturbed = copy.deepcopy(budget)
                    perturbed_density = density
                    if name == "density":
                        perturbed_density += direction * epsilon
                    else:
                        perturbed["components"][name]["correction"] += direction * epsilon
                    result = analyze(self.path, density_mg_ul=perturbed_density, uncertainty_budget=perturbed)
                    flows.append(result["targets"][0]["mean_flow_after_available_corrections_ul_min"])
                sensitivity = abs(flows[1] - flows[0]) / (2 * epsilon)
                expected = sensitivity * component["standard_uncertainty"]
                with self.subTest(gain=gain_correction, component=name):
                    self.assertAlmostEqual(contributions[name], expected, delta=expected * 1e-7)

    def test_corrected_bootstrap_uses_the_same_run_resamples(self):
        rows = [[f"run-{i}", t, 5 + slope * t, 60]
                for i, slope in enumerate([0.8, 1.0, 1.3, 1.7, 2.1])
                for t in [0, 1, 2, 3]]
        self.write_rows(rows)
        budget = self.complete_budget()
        budget["components"]["balance_gain"]["correction"] = 0.25
        report = analyze(self.path, density_mg_ul=2, uncertainty_budget=budget,
                         bootstrap_draws=137, seed=31)
        target = report["targets"][0]
        expected = [(bound + (0.2 - 0.1) * 60 / 2) * 1.003 / 1.25
                    for bound in target["mean_flow_ci95"]]
        np.testing.assert_allclose(target["mean_flow_after_available_corrections_ci95"], expected,
                                   rtol=1e-13, atol=1e-13)
        uncorrected = analyze(self.path, density_mg_ul=2, bootstrap_draws=137, seed=31)["targets"][0]
        self.assertEqual(uncorrected["mean_flow_ci95"], uncorrected["mean_flow_after_available_corrections_ci95"])

    def test_partial_and_unavailable_budgets_do_not_claim_combined_uncertainty(self):
        self.write_rows([["one", t, 2 * t, 60] for t in [0, 1, 2, 3]])
        partial = self.complete_budget()
        partial["components"]["density"] = {
            "status": "not_available", "reason": "no traceable density measurement"
        }
        result = analyze(self.path, density_mg_ul=2, uncertainty_budget=partial)
        uncertainty = result["targets"][0]["measurement_system_uncertainty"]
        self.assertEqual(uncertainty["status"], "partial")
        self.assertIsNone(uncertainty["combined_standard_uncertainty_ul_min"])
        self.assertIn("density", uncertainty["unquantified_components"])
        no_budget = analyze(self.path, density_mg_ul=2)
        self.assertEqual(no_budget["targets"][0]["measurement_system_uncertainty"]["status"], "unavailable")

    def test_uncertainty_budget_cli_binds_bytes_and_renders_report(self):
        self.write_rows([["one", t, 2 * t, 60] for t in [0, 1, 2, 3]])
        budget_path = Path(self.tmp.name) / "budget.json"
        budget_path.write_text(json.dumps(self.complete_budget()), encoding="utf-8")
        output = Path(self.tmp.name) / "cli-output"
        self.assertEqual(main(["analyze", str(self.path), "--density-mg-ul", "2",
                               "--uncertainty-budget", str(budget_path), "--out", str(output)]), 0)
        report = json.loads((output / "report.json").read_text(encoding="utf-8"))
        self.assertEqual(report["measurement_uncertainty_budget"]["input_sha256"],
                         hashlib.sha256(budget_path.read_bytes()).hexdigest())
        text = (output / "REPORT.md").read_text(encoding="utf-8")
        self.assertIn("standard uncertainty, not a 95% interval", text)
        self.assertIn("quantified_for_declared_components", text)
        self.assertEqual(main(["analyze", str(self.path), "--density-mg-ul", "2",
                               "--uncertainty-budget", str(Path(__file__).resolve().parents[1] / "examples/measurement_uncertainty_budget.example.json"),
                               "--out", str(Path(self.tmp.name) / "example-output")]), 0)

    def test_known_slope_with_tare_offset_and_density(self):
        self.write_rows([["one", t, 14 + 2 * t, 60] for t in [0, 1, 2, 3]])
        result = analyze(self.path, density_mg_ul=2)
        run = result["runs"][0]
        self.assertAlmostEqual(run["measured_flow_ul_min"], 60)
        self.assertAlmostEqual(run["intercept_mg"], 14)
        self.assertAlmostEqual(run["error_percent"], 0)
        self.assertIsNone(result["targets"][0]["mean_flow_ci95"])
        self.assertIn("unavailable (<3 runs)", markdown(result))

    def test_checked_in_fixture_report_identifies_the_fixture_bytes(self):
        root = Path(__file__).resolve().parents[1]
        fixture = root / "examples/synthetic_measurements.csv"
        example = json.loads((root / "examples/results/report.json").read_text(encoding="utf-8"))
        self.assertEqual(example["input_sha256"], hashlib.sha256(fixture.read_bytes()).hexdigest())

    def test_constant_evaporation_is_a_signed_shift_only(self):
        shifted = flow_if_constant_evaporation(50, 0.01, 1)
        self.assertAlmostEqual(shifted["correction_ul_min"], 0.6)
        self.assertAlmostEqual(shifted["implied_delivered_flow_ul_min"], 50.6)
        self.assertIn("hypothetical", shifted["note"])
        with self.assertRaises(ValueError):
            flow_if_constant_evaporation(50, 0.01, 0)

    def test_demo_exposes_injected_under_delivery(self):
        write_demo(self.path)
        report = analyze(self.path, density_mg_ul=1, bootstrap_draws=100)
        self.assertEqual(len(report["runs"]), 15)
        for target in report["targets"]:
            self.assertLess(target["mean_error_percent"], -5)
            self.assertGreater(target["mean_error_percent"], -10)
            self.assertEqual(target["n_runs"], 5)
            self.assertIsNotNone(target["mean_flow_ci95"])
        self.assertEqual(report, analyze(self.path, density_mg_ul=1, bootstrap_draws=100))
        readable = markdown(report)
        self.assertIn("Density used: 1 mg/µL", readable)
        self.assertIn("95% run-bootstrap interval", readable)
        first_interval = report["targets"][0]["mean_flow_ci95"]
        self.assertIn(f"[{first_interval[0]:.3f}, {first_interval[1]:.3f}]", readable)
        self.assertIn("not a physical pump calibration", readable)

    def test_reject_duplicate_time_and_nonfinite_values(self):
        for rows in [[["one", 0, 0, 50], ["one", 0, 2, 50], ["one", 1, 3, 50]],
                     [["one", 0, "nan", 50]], [["one", -1, 0, 50]],
                     [["one", 0, 1, 0]]]:
            self.write_rows(rows)
            with self.assertRaises(ValueError):
                analyze(self.path, density_mg_ul=1)

    def test_target_change_requires_separate_run(self):
        self.write_rows([["one", 0, 0, 50], ["one", 1, 1, 50], ["one", 2, 2, 60]])
        with self.assertRaisesRegex(ValueError, "changes"):
            analyze(self.path, density_mg_ul=1)

    def test_warmup_filter_is_relative_to_run_start(self):
        self.write_rows([["one", t, mass, 60] for t, mass in [(100, -20), (101, 0), (102, 1), (103, 2)]])
        report = analyze(self.path, density_mg_ul=1, discard_seconds=1)
        self.assertAlmostEqual(report["runs"][0]["measured_flow_ul_min"], 60)
        with self.assertRaises(ValueError):
            analyze(self.path, density_mg_ul=1, discard_seconds=3)

    def test_drift_flag_and_constant_trace(self):
        self.write_rows([["drift", t, t if t < 5 else 5 + 2 * (t - 5), 60] for t in range(10)])
        self.assertIn("inspect_within_run_drift", analyze(self.path, density_mg_ul=1)["runs"][0]["flags"])
        self.write_rows([["stalled", t, 5, 60] for t in range(4)])
        report = analyze(self.path, density_mg_ul=1)
        self.assertIn("nonpositive_measured_flow", report["runs"][0]["flags"])
        self.assertEqual(report["targets"][0]["flagged_runs"], ["stalled"])

    def test_fit_does_not_require_zero_start_mass(self):
        fit = fit_line(np.arange(5), np.arange(5) * 0.5 - 100)
        self.assertAlmostEqual(fit["mass_rate_mg_s"], 0.5)
        self.assertAlmostEqual(fit["r_squared"], 1)

    def test_hash_identifies_the_parsed_snapshot_if_input_changes(self):
        self.write_rows([["one", t, 14 + 2 * t, 60] for t in [0, 1, 2, 3]])
        original = self.path.read_bytes()
        read_bytes = Path.read_bytes
        reads = []

        def replace_after_read(path):
            data = read_bytes(path)
            reads.append(path)
            path.write_text("changed after reading\n", encoding="utf-8")
            return data

        with patch.object(Path, "read_bytes", replace_after_read):
            report = analyze(self.path, density_mg_ul=2)
        self.assertEqual(reads, [self.path])
        self.assertEqual(report["input_sha256"], hashlib.sha256(original).hexdigest())
        self.assertAlmostEqual(report["runs"][0]["measured_flow_ul_min"], 60)
        self.assertNotEqual(self.path.read_bytes(), original)

    def test_failed_report_publication_can_retry_without_partial_output(self):
        self.write_rows([["one", t, 14 + 2 * t, 60] for t in [0, 1, 2, 3]])
        output = Path(self.tmp.name) / "nested" / "report"
        args = ["analyze", str(self.path), "--density-mg-ul", "2", "--out", str(output)]
        rename = Path.rename

        def fail_second_move(path, target):
            if path.name == "REPORT.md":
                self.assertTrue((output / "report.json").exists())
                raise OSError("simulated publication failure")
            return rename(path, target)

        with patch.object(Path, "rename", fail_second_move):
            with self.assertRaises(SystemExit):
                main(args)
        self.assertFalse(output.exists())
        self.assertEqual(main(args), 0)
        report = json.loads((output / "report.json").read_text(encoding="utf-8"))
        self.assertEqual(report["input_sha256"], hashlib.sha256(self.path.read_bytes()).hexdigest())
        self.assertIn("# Perfusion calibration analysis", (output / "REPORT.md").read_text(encoding="utf-8"))

    def test_existing_report_directory_is_preserved(self):
        self.write_rows([["one", t, t, 60] for t in [0, 1, 2]])
        output = Path(self.tmp.name) / "report"
        output.mkdir()
        marker = output / "keep.txt"
        marker.write_text("original", encoding="utf-8")
        with self.assertRaises(SystemExit):
            main(["analyze", str(self.path), "--density-mg-ul", "1", "--out", str(output)])
        self.assertEqual(marker.read_text(encoding="utf-8"), "original")
        self.assertEqual(list(output.iterdir()), [marker])

    def test_outlier_detection_iqr_and_grubbs(self):
        from perfusioncal.analysis import detect_outliers, detect_outliers_grubbs, detect_outliers_iqr, _grubbs_critical_value
        # Baseline residuals with one prominent spike
        clean = np.zeros(20)
        self.assertEqual(detect_outliers_iqr(clean), [])
        spiked = np.zeros(20)
        spiked[10] = 50.0  # Spike at index 10
        self.assertEqual(detect_outliers_iqr(spiked, factor=1.5), [10])
        self.assertEqual(detect_outliers_grubbs(spiked), [10])
        self.assertEqual(detect_outliers(spiked, method="iqr"), [10])
        self.assertEqual(detect_outliers(spiked, method="grubbs"), [10])
        with self.assertRaises(ValueError):
            detect_outliers(spiked, method="unknown")
        self.assertAlmostEqual(_grubbs_critical_value(20, 0.05), 2.708246, places=5)
        self.assertAlmostEqual(_grubbs_critical_value(100, 0.05), 3.384083, places=5)

    def test_grubbs_critical_value_detects_residual_above_exact_cutoff(self):
        times = np.arange(20, dtype=float)
        error = np.tile([-1.0, 1.0], 10)
        error[10] = 4.0
        fit = fit_line(times, times + error, outlier_method="grubbs", outlier_threshold=0.05)
        self.assertEqual(fit["outlier_count"], 1)
        self.assertEqual(fit["outlier_indices"], [10])

    def test_outlier_thresholds_reject_invalid_values(self):
        for method, threshold in [("grubbs", 0), ("grubbs", 1), ("iqr", -1), ("iqr", float("nan"))]:
            with self.subTest(method=method, threshold=threshold), self.assertRaises(ValueError):
                analyze(self.path, density_mg_ul=1, outlier_method=method, outlier_threshold=threshold)

    def test_outlier_flagging_in_analyze_and_markdown(self):
        # Create a run where one mass reading is an outlier
        times = list(range(0, 100, 10))  # 10 readings: 0, 10, ..., 90
        masses = [5.0 + 1.0 * t for t in times]
        masses[5] += 25.0  # Spike at time 50s
        rows = [["spike_run", t, m, 60] for t, m in zip(times, masses)]
        self.write_rows(rows)

        report = analyze(self.path, density_mg_ul=1.0)
        run = report["runs"][0]
        self.assertEqual(run["outlier_count"], 1)
        self.assertEqual(run["outlier_indices"], [5])
        self.assertEqual(run["outlier_times_s"], [50.0])
        self.assertIn("outlier_readings_detected(1)", run["flags"])
        self.assertIn("spike_run", report["targets"][0]["outlier_runs"])

        md = markdown(report)
        self.assertIn("outlier_readings_detected(1)", md)
        self.assertIn("Outlier method: iqr", md)

    def test_cli_outlier_options(self):
        times = list(range(0, 100, 10))
        masses = [5.0 + 1.0 * t for t in times]
        masses[4] += 30.0
        self.write_rows([["run_a", t, m, 50] for t, m in zip(times, masses)])
        output = Path(self.tmp.name) / "cli_outlier_out"
        rc = main(["analyze", str(self.path), "--density-mg-ul", "1.0",
                   "--outlier-method", "grubbs", "--out", str(output)])
        self.assertEqual(rc, 0)
        rep = json.loads((output / "report.json").read_text(encoding="utf-8"))
        self.assertEqual(rep["configuration"]["outlier_method"], "grubbs")
        self.assertEqual(rep["runs"][0]["outlier_count"], 1)


if __name__ == "__main__":
    unittest.main()
