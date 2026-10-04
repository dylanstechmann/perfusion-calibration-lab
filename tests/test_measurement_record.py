"""Constructed schema fixtures test intake behavior; these are not measurements."""

import copy
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from perfusioncal.analysis import analyze, markdown
from perfusioncal.cli import main
from perfusioncal.measurement_record import qualify_measurements


class MeasurementRecordTests(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.csv = self.root / "test.csv"
        self.data = b"run_id,time_s,mass_mg,target_flow_ul_min\na,0,7,45\na,10,22,45\na,20,37,45\nb,0,8,45\nb,10,23,45\nb,20,38,45\n"
        self.csv.write_bytes(self.data)
        (self.root / "native.txt").write_bytes(b"constructed-test-export")
        self.record = self.root / "record.json"
        common = {key: "constructed-test-only" for key in ("independence_record", "pump_identity", "firmware_identity", "syringe_or_tubing", "fluid", "balance_identity", "timebase_record", "collection_setup", "density_record")}
        common.update(complete_trace=True, independent_acquisition=True, source_file="native", fluid_temperature_degC=20.0, density_mg_ul=2.0)
        self.document = {"schema_version": 1, "origin": "direct_measurement", "scope": "public_third_party",
                         "source_url": "https://example.invalid/test-only", "license": "test fixture", "attribution": "constructed schema test", "acquisition_record": "test-only",
                         "units": {"time": "s", "mass": "mg", "target_flow": "ul/min", "density": "mg/ul", "temperature": "degC"},
                         "csv_sha256": hashlib.sha256(self.data).hexdigest(),
                         "source_files": {"native": {"path": "native.txt", "sha256": hashlib.sha256(b"constructed-test-export").hexdigest()}},
                         "runs": {key: dict(common, source_trace_id=key) for key in ("a", "b")}}

    def qualify(self, document=None):
        self.record.write_text(json.dumps(document or self.document), encoding="utf-8")
        return qualify_measurements(self.data, self.record, density_mg_ul=2.0)

    def test_source_is_preserved_and_density_known_slope_survives_intake(self):
        self.qualify()
        report = analyze(self.csv, density_mg_ul=2.0, measurement_record=self.record)
        self.assertEqual(report["measurement_source_qualification"]["n_complete_source_traces"], 2)
        self.assertAlmostEqual(report["targets"][0]["mean_flow_ul_min"], 45.0, places=12)
        self.assertIn("Third-party measurements do not calibrate", markdown(report))
        self.assertIn("not a physical pump calibration", " ".join(report["notes"]))

    def test_ordinary_numeric_csv_has_no_measurement_provenance(self):
        report = analyze(self.csv, density_mg_ul=2.0)
        self.assertEqual(report["measurement_source_qualification"]["status"], "not_provided")

    def test_rejects_unsupported_data_origins(self):
        for origin in ("synthetic", "figure_digitized", "averaged", None):
            with self.subTest(origin=origin):
                document = copy.deepcopy(self.document); document["origin"] = origin
                with self.assertRaisesRegex(ValueError, "direct_measurement"):
                    self.qualify(document)

    def test_rejects_missing_runs_and_incomplete_or_nonindependent_traces(self):
        document = copy.deepcopy(self.document); del document["runs"]["b"]
        with self.assertRaisesRegex(ValueError, "exactly all"):
            self.qualify(document)
        for field in ("complete_trace", "independent_acquisition"):
            document = copy.deepcopy(self.document); document["runs"]["b"][field] = False
            with self.assertRaisesRegex(ValueError, "complete independently"):
                self.qualify(document)

    def test_rejects_one_source_trace_under_multiple_run_and_file_ids(self):
        document = copy.deepcopy(self.document)
        document["source_files"]["alias"] = copy.deepcopy(document["source_files"]["native"])
        document["runs"]["b"].update(source_file="alias", source_trace_id="a")
        with self.assertRaisesRegex(ValueError, "split into replicate"):
            self.qualify(document)

    def test_rejects_modified_csv_and_source_exports(self):
        document = copy.deepcopy(self.document); document["csv_sha256"] = "0" * 64
        with self.assertRaisesRegex(ValueError, "CSV hash mismatch"):
            self.qualify(document)
        (self.root / "native.txt").write_bytes(b"changed")
        with self.assertRaisesRegex(ValueError, "source_file hash mismatch"):
            self.qualify()

    def test_rejects_wrong_units_density_and_nonfinite_temperature(self):
        for field, value in (("density_mg_ul", 1.0), ("density_mg_ul", True), ("fluid_temperature_degC", float("nan"))):
            document = copy.deepcopy(self.document); document["runs"]["a"][field] = value
            with self.assertRaises(ValueError):
                self.qualify(document)
        document = copy.deepcopy(self.document); document["units"]["mass"] = "g"
        with self.assertRaisesRegex(ValueError, "canonical units"):
            self.qualify(document)

    def test_rejects_missing_acquisition_context_duplicate_keys_and_traversal(self):
        document = copy.deepcopy(self.document); document["runs"]["a"]["density_record"] = ""
        with self.assertRaisesRegex(ValueError, "density_record"):
            self.qualify(document)
        self.record.write_text('{"schema_version": 1, "schema_version": 1}')
        with self.assertRaisesRegex(ValueError, "duplicate"):
            qualify_measurements(self.data, self.record, density_mg_ul=2)
        document = copy.deepcopy(self.document); document["source_files"]["native"]["path"] = "../escape.txt"
        with self.assertRaisesRegex(ValueError, "inside"):
            self.qualify(document)

    def test_cli_refuses_synthetic_measurement_claim_before_output(self):
        document = copy.deepcopy(self.document); document["origin"] = "synthetic"
        self.record.write_text(json.dumps(document))
        output = self.root / "report"
        with self.assertRaises(SystemExit) as error:
            main(["analyze", str(self.csv), "--density-mg-ul", "2", "--measurement-record", str(self.record), "--out", str(output)])
        self.assertEqual(error.exception.code, 2)
        self.assertFalse(output.exists())
