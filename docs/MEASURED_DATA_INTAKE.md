# Source qualification and measured-trace intake

The public-source audit on 2026-10-04 acquired **zero eligible raw gravimetric
trace datasets**. [The structured catalog](public-source-qualification.json)
records a bounded review, not a claim that suitable public data does not exist.

[Jung, Shara and Bruns (2025)](https://doi.org/10.1371/journal.pone.0326938)
measured dispensed water every 30 seconds for 5 minutes, with three repeats per
syringe/speed condition. The article analyzes averaged readings and assumes
1 g/mL density. Its motor settings are not volumetric flow targets, and we did
not acquire individual machine-readable traces. Those averages must not become
invented independent runs or density measurements.

[Ha et al., Zenodo version v1.0.3](https://zenodo.org/records/21704252) advertises
raw gravimetric validation data, but its files are explicitly restricted and
its associated unauthenticated GitHub URL returned 404. No raw measurements
were downloaded. The record describes CC BY-SA 4.0 data licensing. Access and
an applicable license are both needed before importing files.

[Page access results and hashes](public-source-access-20261004.json) identify
local snapshots in `artifacts/public-source-audit-20261004/`. These are source
pages, not measurement data. Recheck access with a fresh output directory:

```bash
.venv/bin/python scripts/probe_public_sources.py --out artifacts/new-source-audit
```

A repeat HTML hash can change because pages are dynamic. HTTP success alone
can include a robot-check page and never establishes raw-data eligibility.

## Use the existing recording and uncertainty formats

Continue using [the physical CSV recording guide](../examples/PHYSICAL_MEASUREMENT_RECORDING.md)
and [the existing uncertainty budget](../examples/measurement_uncertainty_budget.example.json).
No new acquisition process, pump driver or uncertainty-budget model is added.

For direct measurements, add a JSON sidecar and invoke:

```bash
perfusioncal analyze /path/to/measurements.csv --density-mg-ul ACTUAL_DENSITY --measurement-record /path/to/record.json --out artifacts/new-measured-analysis
```

`ACTUAL_DENSITY` must be the recorded value for the fluid and temperature.
A third-party trace checks the offline workflow for that source; it cannot
calibrate the user's pump. A numeric CSV without a sidecar remains analyzable,
but its report marks measured-data provenance as `not_provided`.

## JSON sidecar contract (schema version 1)

The top-level object needs these fields:

| Field | Required value or meaning |
|---|---|
| `schema_version` | Integer `1` |
| `origin` | `direct_measurement`; synthetic, averaged and figure-digitized inputs are rejected |
| `scope` | `public_third_party` or `local_instrument` |
| `source_url`, `license`, `attribution`, `acquisition_record` | Nonempty provenance and applicable use permission; use a local record identifier for private measurements |
| `units` | Exactly `{"time":"s","mass":"mg","target_flow":"ul/min","density":"mg/ul","temperature":"degC"}` |
| `csv_sha256` | SHA-256 of the exact supplied normalized CSV bytes |
| `source_files` | Object keyed by file ID, each with relative `path` and exact native-file `sha256` |
| `runs` | Object keyed by exactly the CSV's complete `run_id` set |

Each run object needs `complete_trace: true`, `independent_acquisition: true`,
`source_file` (a key in `source_files`), `source_trace_id` (the original trace's
identifier, not an invented replicate label), and a nonempty
`independence_record` identifying the acquisition evidence. Also record
`pump_identity`, `firmware_identity`, `syringe_or_tubing`, `fluid`,
`balance_identity`, `timebase_record`, `collection_setup`, `density_record`,
finite `fluid_temperature_degC` and positive `density_mg_ul`.
The density must match the analysis argument for every run. Analyze separate
fluid/temperature/density conditions separately.

Native source paths must remain inside the sidecar directory. A source trace
is identified by **native bytes hash plus original trace ID**. Reusing that
identity under different run IDs or file-ID aliases is rejected. An export
containing several truly independent runs may use one file with distinct,
source-supported trace IDs. Keep the complete native export; do not crop one
continuous trace into purported replicates. Unit conversions must be documented
in the acquisition record and must preserve the native bytes and trace IDs.

This gate checks hashes, declared identities and required context. It cannot
independently prove that a physical run occurred, detect deliberately invented
source identities, or quantify density/balance/evaporation/timing uncertainty.
The separate existing uncertainty budget is still required for those effects.
The analyzer retains flagged runs. Passing intake or observing high R-squared
is not instrument acceptance or a physical calibration certificate. Constructed
unit-test fixtures establish software behavior only.
