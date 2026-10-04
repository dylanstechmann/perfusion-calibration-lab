"""Qualify declared direct measurements without certifying acquisition or hardware."""

import csv
import hashlib
import io
import json
import math
from pathlib import Path


def _object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate measurement-record key: {key}")
        result[key] = value
    return result


def _text(record, key):
    value = record.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"measurement record needs nonempty {key}")
    return value


def _number(record, key):
    value = record.get(key)
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"measurement record needs finite {key}")
    return value


def qualify_measurements(input_bytes, record_path, *, density_mg_ul):
    """Check source bytes, units, complete-run declarations and source trace identity.

    The provenance declarations remain the provider's responsibility. This gate
    cannot independently prove that a physical run happened or that IDs are true.
    Synthetic fixtures remain usable through analyze without a measurement record.
    """
    path = Path(record_path)
    record_bytes = path.read_bytes()
    document = json.loads(record_bytes.decode("utf-8-sig"), object_pairs_hook=_object)
    if not isinstance(document, dict) or type(document.get("schema_version")) is not int or document["schema_version"] != 1:
        raise ValueError("measurement record schema_version must be 1")
    if document.get("origin") != "direct_measurement":
        raise ValueError("measurement record requires direct_measurement; synthetic, digitized and averaged traces are ineligible")
    if document.get("scope") not in {"public_third_party", "local_instrument"}:
        raise ValueError("scope must be public_third_party or local_instrument")
    for key in ("source_url", "license", "attribution", "acquisition_record"):
        _text(document, key)
    if document.get("units") != {"time": "s", "mass": "mg", "target_flow": "ul/min", "density": "mg/ul", "temperature": "degC"}:
        raise ValueError("measurement record must explicitly use the canonical units")
    if document.get("csv_sha256") != hashlib.sha256(input_bytes).hexdigest():
        raise ValueError("measurement CSV hash mismatch")
    files = document.get("source_files")
    if not isinstance(files, dict) or not files:
        raise ValueError("measurement record needs source_files")
    verified = {}
    source_root = path.parent.resolve()
    for file_id, entry in files.items():
        if not isinstance(file_id, str) or not file_id.strip() or not isinstance(entry, dict):
            raise ValueError("invalid source_file entry")
        relative = Path(_text(entry, "path"))
        if relative.is_absolute():
            raise ValueError("source_file path must be relative to measurement record")
        source_path = (source_root / relative).resolve()
        if not source_path.is_relative_to(source_root):
            raise ValueError("source_file path must stay inside measurement-record directory")
        content = source_path.read_bytes()
        digest = hashlib.sha256(content).hexdigest()
        if entry.get("sha256") != digest:
            raise ValueError(f"source_file hash mismatch: {file_id}")
        verified[file_id] = digest
    reader = csv.DictReader(io.StringIO(input_bytes.decode("utf-8-sig")))
    required = {"run_id", "time_s", "mass_mg", "target_flow_ul_min"}
    header = reader.fieldnames or []
    if len(header) != len(set(header)) or not required.issubset(header):
        raise ValueError("measurement CSV needs unique canonical columns")
    run_ids = set()
    for line, row in enumerate(reader, 2):
        if None in row or any(value is None for value in row.values()) or not row["run_id"].strip():
            raise ValueError(f"row {line}: incomplete measurement CSV")
        run_ids.add(row["run_id"].strip())
    if not run_ids:
        raise ValueError("no runs in measurement CSV")
    runs = document.get("runs")
    if not isinstance(runs, dict) or set(runs) != run_ids:
        raise ValueError("measurement record must cover exactly all CSV run_ids")
    identities = set()
    for run_id, entry in runs.items():
        if not isinstance(entry, dict):
            raise ValueError(f"invalid measurement run: {run_id}")
        if entry.get("complete_trace") is not True or entry.get("independent_acquisition") is not True:
            raise ValueError("each run must be a complete independently acquired trace")
        for key in ("source_trace_id", "independence_record", "pump_identity", "firmware_identity", "syringe_or_tubing", "fluid", "balance_identity", "timebase_record", "collection_setup", "density_record"):
            _text(entry, key)
        file_id = _text(entry, "source_file")
        if file_id not in verified:
            raise ValueError(f"unknown source_file for {run_id}")
        # Use bytes identity too: aliasing one original file under two IDs cannot
        # turn slices of the same source trace into independent acquisitions.
        identity = (verified[file_id], entry["source_trace_id"])
        if identity in identities:
            raise ValueError("one source trace cannot be split into replicate run_ids")
        identities.add(identity)
        _number(entry, "fluid_temperature_degC")
        measured_density = _number(entry, "density_mg_ul")
        if measured_density <= 0 or measured_density != density_mg_ul:
            raise ValueError("run density must be positive and match --density-mg-ul")
    return {"status": "qualified_declared_direct_measurement", "scope": document["scope"],
            "record_sha256": hashlib.sha256(record_bytes).hexdigest(),
            "source_url": document["source_url"], "license": document["license"],
            "attribution": document["attribution"], "source_files_sha256": verified,
            "runs": runs, "n_complete_source_traces": len(identities),
            "limits": ["Source declarations are checked for consistency, not independently authenticated.",
                       "Third-party measurements do not calibrate the user's instrument.",
                       "This intake gate does not establish hardware accuracy or measurement-system uncertainty."]}
