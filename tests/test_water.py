from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from perfusioncal.analysis import analyze, markdown
from perfusioncal.cli import main
from perfusioncal.water import (
    MAXIMUM_TEMPERATURE_C,
    MINIMUM_TEMPERATURE_C,
    water_density_kg_m3,
    water_density_mg_ul,
    water_density_record,
)

EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "synthetic_measurements.csv"


class WaterDensityTests(unittest.TestCase):
    def test_matches_published_table_values(self):
        # Tanaka et al. 2001 recommended densities, kg/m^3. The closed form reproduces
        # the table to a few parts per million, so a 0.01 kg/m^3 tolerance is generous.
        for temperature, expected in ((0, 999.8428), (5, 999.9668), (10, 999.7027),
                                      (20, 998.2071), (25, 997.0479), (30, 995.6502),
                                      (37, 993.3316), (40, 992.2187)):
            with self.subTest(temperature=temperature):
                self.assertAlmostEqual(water_density_kg_m3(temperature), expected, delta=0.01)

    def test_density_peaks_near_four_degrees_and_falls_above_it(self):
        peak = max((water_density_kg_m3(value / 100) for value in range(0, 1000)),)
        self.assertAlmostEqual(water_density_kg_m3(3.98), peak, delta=1e-4)
        self.assertGreater(water_density_kg_m3(3.98), water_density_kg_m3(0.0))
        temperatures = [5, 10, 15, 20, 25, 30, 35, 40]
        densities = [water_density_kg_m3(value) for value in temperatures]
        self.assertEqual(densities, sorted(densities, reverse=True))

    def test_milligrams_per_microlitre_conversion(self):
        self.assertAlmostEqual(water_density_mg_ul(20), water_density_kg_m3(20) * 1e-3, places=12)
        self.assertAlmostEqual(water_density_mg_ul(20), 0.998207, places=5)
        # The familiar 1 mg/uL is wrong by about 0.18% at 20 C.
        self.assertAlmostEqual(
            water_density_record(20)["relative_error_of_unit_density_approximation"],
            0.0018, delta=0.0002)

    def test_temperatures_outside_the_published_range_are_refused(self):
        for bad in (MINIMUM_TEMPERATURE_C - 0.1, MAXIMUM_TEMPERATURE_C + 0.1, -10, 60):
            with self.subTest(bad=bad):
                with self.assertRaisesRegex(ValueError, "outside the published"):
                    water_density_kg_m3(bad)
        for bad in (float("nan"), float("inf"), float("-inf")):
            with self.subTest(bad=bad):
                with self.assertRaises(ValueError):
                    water_density_kg_m3(bad)
        for bad in ("20", None, True):
            with self.subTest(bad=bad):
                with self.assertRaises(ValueError):
                    water_density_kg_m3(bad)

    def test_record_carries_citation_scope_and_limits(self):
        record = water_density_record(37)
        self.assertEqual(record["source"], "tanaka_2001_pure_water")
        self.assertIn("Metrologia 38", record["citation"])
        self.assertIn("doi.org/10.1088/0026-1394/38/4/3", record["citation"])
        self.assertEqual(record["fluid"], "pure air-free water at standard atmospheric pressure")
        self.assertEqual(record["valid_temperature_range_c"], [0.0, 40.0])
        self.assertIn("Not valid for culture medium", record["scope"])
        joined = " ".join(record["limitations"])
        self.assertIn("air-buoyancy", joined)
        self.assertIn("does not make a flow estimate a physical pump calibration", joined)
        json.dumps(record, allow_nan=False)


class DensityProvenanceTests(unittest.TestCase):
    def test_supplied_density_is_recorded_as_unexplained(self):
        report = analyze(EXAMPLE, density_mg_ul=1.0)
        self.assertEqual(report["fluid_density"]["source"], "supplied_by_caller")
        self.assertEqual(report["fluid_density"]["density_mg_ul"], 1.0)
        self.assertIn("does not record how it was determined",
                      " ".join(report["fluid_density"]["limitations"]))

    def test_water_temperature_changes_flow_by_the_density_ratio(self):
        unit = analyze(EXAMPLE, density_mg_ul=1.0)
        water = water_density_record(25)
        corrected = analyze(EXAMPLE, density_mg_ul=water["density_mg_ul"],
                            density_provenance=water)
        self.assertEqual(corrected["fluid_density"]["source"], "tanaka_2001_pure_water")
        self.assertEqual(corrected["configuration"]["density_mg_ul"], water["density_mg_ul"])
        for before, after in zip(unit["targets"], corrected["targets"]):
            self.assertAlmostEqual(after["mean_flow_ul_min"],
                                   before["mean_flow_ul_min"] / water["density_mg_ul"],
                                   places=6)
        self.assertIn("Fluid density", markdown(corrected))
        self.assertIn("tanaka_2001_pure_water", markdown(corrected))

    def test_cli_requires_exactly_one_density_source(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "out"
            with self.assertRaises(SystemExit):
                main(["analyze", str(EXAMPLE), "--out", str(output)])
            with self.assertRaises(SystemExit):
                main(["analyze", str(EXAMPLE), "--density-mg-ul", "1.0",
                      "--water-temperature-c", "20", "--out", str(output)])
            self.assertFalse(output.exists())

    def test_cli_water_temperature_publishes_the_citation(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "water"
            self.assertEqual(
                main(["analyze", str(EXAMPLE), "--water-temperature-c", "22.5",
                      "--out", str(output)]), 0)
            report = json.loads((output / "report.json").read_text(encoding="utf-8"))
            self.assertEqual(report["fluid_density"]["temperature_c"], 22.5)
            self.assertAlmostEqual(report["configuration"]["density_mg_ul"],
                                   water_density_mg_ul(22.5), places=12)
            text = (output / "REPORT.md").read_text(encoding="utf-8")
            self.assertIn("Metrologia 38", text)
            self.assertIn("would bias flow by", text)

    def test_cli_rejects_an_out_of_range_temperature(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "hot"
            with self.assertRaises(SystemExit):
                main(["analyze", str(EXAMPLE), "--water-temperature-c", "95", "--out", str(output)])
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
