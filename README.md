# Perfusion Calibration Lab

This is a personal hobby and learning project, developed with substantial assistance from AI coding tools.

**Measure delivered flow from balance readings and compare independent runs.**
An offline companion to [open-perfusion-rig](https://github.com/dylanstechmann/open-perfusion-rig):
the pump repo predicts displacement from geometry; this repo analyzes what a
balance recorded. No controller or hardware connection is required.

## Quickstart

Python 3.10+ and NumPy; CPU only.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e .
python -m unittest discover -s tests -v
perfusioncal analyze examples/synthetic_measurements.csv \
  --density-mg-ul 1 --out artifacts/demo
```

The included fixture assumes density **exactly 1 mg/µL** and deliberately
injects about 7% underdelivery with noise. It is not a physical calibration or
a claim about water density at your temperature. See the generated
[example report](examples/results/REPORT.md).

`report.json` contains all fitted runs, flags, settings, input hash and software
versions. `REPORT.md` shows the density, repeat standard deviation and the
run-bootstrap interval where at least three independent runs exist. That
interval excludes density and balance bias. Use a new output directory for
each run so previous analyses are preserved.

The input hash is calculated from the same CSV bytes used for fitting, so it
identifies the analyzed snapshot even if the source file changes during a run.
The checked-in synthetic CSV is pinned to LF line endings so its example report
has the same byte hash on Windows and Linux. User CSVs retain their own byte hash.
Reports are prepared beside the destination and then published to a new output
directory. Existing outputs are preserved, and ordinary write failures remove
the incomplete new directory so the command can be retried. This does not
guarantee an atomic report after a process crash or power loss.

## Your measurements

```csv
run_id,time_s,mass_mg,target_flow_ul_min
run01,0,7.0,50
run01,30,31.9,50
run01,60,56.8,50
```

- One row per balance reading; units are explicit in the header.
- Times must increase within each run; duplicate times are rejected.
- A run has one positive target flow and at least three usable readings.
- Use a new `run_id` for each independent repeat and each target change.
- Enter the fluid density for the measurement conditions with `--density-mg-ul`.
- Optional `--discard-seconds 30` excludes a predefined startup interval,
  measured relative to the first timestamp of each run. Choose it before analysis.

Do not label slices of a single trace as independent repeated runs.
[`examples/physical_measurements.template.csv`](examples/physical_measurements.template.csv)
is a blank header-only capture template; see
[`examples/PHYSICAL_MEASUREMENT_RECORDING.md`](examples/PHYSICAL_MEASUREMENT_RECORDING.md)
for source records to retain alongside an actual trace.

## Calculation

An ordinary least-squares line with a free intercept fits mass versus time:

```text
mass_mg = intercept_mg + mass_rate_mg_s * time_s
flow_ul_min = 60 * mass_rate_mg_s / density_mg_ul
error_percent = 100 * (flow_ul_min - target_flow_ul_min) / target_flow_ul_min
```

### Where the density comes from

Flow is a mass rate divided by a density, so the density is part of the result.
Pass one of two things:

- `--density-mg-ul` — the density of the fluid you actually pumped, at its
  temperature. The report records that it was supplied and that nothing is known
  about how it was determined.
- `--water-temperature-c` — the published density of pure air-free water at that
  temperature, from Tanaka et al., *Metrologia* **38** (2001) 301–309
  ([doi:10.1088/0026-1394/38/4/3](https://doi.org/10.1088/0026-1394/38/4/3)),
  valid 0–40 °C. The report carries the citation, the temperature, the density in
  both mg/µL and kg/m³, and how much the familiar 1 mg/µL would have biased the
  result — about +0.18% at 20 °C and +0.30% at 37 °C.

The formula is for **pure water only**. Culture medium, buffered saline and
gas-saturated water all differ, no air-buoyancy correction is applied (weighing
in air biases apparent mass by roughly 0.1% for water), and a bath thermometer
reading is not necessarily the temperature of fluid leaving the needle. Using the
right density does not turn this analysis into a physical pump calibration.

The free intercept accommodates a constant tare offset. Repeat runs at the same
target are weighted equally. The report includes their sample standard deviation
and a percentile bootstrap interval for mean flow, resampling entire runs.
There is no interval with fewer than three runs. A small repeat count still gives
an unstable interval; repeated balance readings are not independent replicates.

## Diagnostics and limits

Low R² (<0.95), nonpositive flow, large changes between early/late half-run
slopes (>20%), and detected mass reading outliers (via IQR Tukey fences or
Grubbs' test on linear residuals) are inspection flags. They are transparent
software heuristics, not calibration acceptance standards. At least six
readings are needed for the half-run check. Flagged runs and outlier counts
remain in the summaries and are listed there for operator review without
silently discarding readings.

Options:
- `--outlier-method {iqr, grubbs}` (default: `iqr`)
- `--outlier-threshold` (default: 1.5 for IQR, 0.05 for Grubbs)
- `--uncertainty-budget FILE` (optional source-linked standard-uncertainty budget)

Run-bootstrap intervals do not include systematic uncertainty from fluid density,
the balance, evaporation, retained droplets, collection losses or the timing instrument.
`flow_if_constant_evaporation` only shifts an apparent flow by a constant
mass-loss rate you supply. It does not estimate that rate.
Define a measurement procedure and an uncertainty budget before making a
hardware accuracy claim. This software has been tested on constructed traces;
no physical pump has been calibrated in this repository.

The optional JSON budget uses the format in
[`examples/measurement_uncertainty_budget.example.json`](examples/measurement_uncertainty_budget.example.json).
It records standard uncertainty for density, balance gain, balance slope drift,
evaporation rate, and timing scale. Components without supporting measurements
must say `not_available` with a reason. For measured components, `source_record`
identifies the certificate, blank run, or other measurement record. The density
point value still comes from `--density-mg-ul`; the budget records its standard
uncertainty.

Measured components other than density also require a point `correction`:

- `balance_gain`: fractional excess/deficit in indicated mass divided by true
  mass. The flow calculation divides by `1 + correction`.
- `balance_slope_drift`: additive blank-trace bias in mg/s. The calculation
  subtracts it from the fitted mass rate.
- `evaporation_rate`: nonnegative matched-blank mass loss in mg/s. The
  calculation adds it using `flow_if_constant_evaporation`.
- `timing_scale`: fractional excess/deficit in indicated elapsed time divided
  by true elapsed time. The flow calculation multiplies by `1 + correction`.

The correction model is
`Q = 60 * (mass_rate - drift + evaporation) * timing_scale / (density * balance_gain)`,
where each scale is `1 + correction`. Drift and evaporation rates must use the
same indicated mass/time units as the trace. Their first-order sensitivity is
`60 * timing_scale / (density * balance_gain)`; the gain divides both the point
flow and those uncertainty contributions. Regression tests compare every named
sensitivity against independent finite differences at non-unit gain values.
Apparent and corrected bootstrap intervals use the same sampled runs, preserving
the exact affine relation between them without adding Monte Carlo differences.

For the two scale terms, `standard_uncertainty` and `correction` are fractional
values (for example, `0.001` means 0.1%). Density uncertainty is in mg/µL;
slope-drift and evaporation uncertainty are in mg/s. The report shows
flow-equivalent contributions and combines them only when all five components
are measured and declared independent. This first-order combined value is a
standard uncertainty, not a 95% interval. Correlated components need covariance
propagation, which this tool does not implement. Corrections are applied only
when the corresponding component is marked `measured`; the original apparent
flow and its run-bootstrap interval remain visible. Neither interval covers
unlisted collection losses.

For example, a `balance_gain` entry with `correction: 0.002` represents a
measured indication gain of 1.002, which the calculation divides out. A
`standard_uncertainty` of `0.0001` is the uncertainty in that dimensionless
gain. Replace `source_record` with the identifier for the actual calibration
record. Numeric values in this example description are illustrative, not
calibration data.

This is analysis for research instruments, not medical infusion control. It
produces no motor commands and does not automatically change a pump calibration.

## Reproduce the fixture

```bash
perfusioncal demo --seed 0 --out artifacts/new_synthetic_measurements.csv
```

The fixture and code are MIT licensed; no personal or biological data are used.

## Measured data provenance

See the [public-source audit and executable intake contract](docs/MEASURED_DATA_INTAKE.md). The 2026-10-04 audit acquired zero eligible raw datasets. Source-linked intake checks preserve hashes and reject unsupported measurement origins or split-source replicates; they do not establish hardware or biological validation.
