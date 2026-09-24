import csv
import tempfile
import unittest
from pathlib import Path

import numpy as np

from perfusioncal.analysis import analyze, fit_line
from perfusioncal.cli import write_demo


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

    def test_known_slope_with_tare_offset_and_density(self):
        self.write_rows([["one", t, 14 + 2 * t, 60] for t in [0, 1, 2, 3]])
        result = analyze(self.path, density_mg_ul=2)
        run = result["runs"][0]
        self.assertAlmostEqual(run["measured_flow_ul_min"], 60)
        self.assertAlmostEqual(run["intercept_mg"], 14)
        self.assertAlmostEqual(run["error_percent"], 0)
        self.assertIsNone(result["targets"][0]["mean_flow_ci95"])

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


if __name__ == "__main__":
    unittest.main()
