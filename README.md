# Perfusion Calibration Lab

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

## Calculation

An ordinary least-squares line with a free intercept fits mass versus time:

```text
mass_mg = intercept_mg + mass_rate_mg_s * time_s
flow_ul_min = 60 * mass_rate_mg_s / density_mg_ul
error_percent = 100 * (flow_ul_min - target_flow_ul_min) / target_flow_ul_min
```

The free intercept accommodates a constant tare offset. Repeat runs at the same
target are weighted equally. The report includes their sample standard deviation
and a percentile bootstrap interval for mean flow, resampling entire runs.
There is no interval with fewer than three runs. A small repeat count still gives
an unstable interval; repeated balance readings are not independent replicates.

## Diagnostics and limits

Low R² (<0.95), nonpositive flow, and large changes between early/late half-run
slopes (>20%) are inspection flags. They are transparent software heuristics,
not calibration acceptance standards. At least six readings are needed for
the half-run check. Flagged runs remain in the summaries and are listed there.

Intervals do not include systematic uncertainty from fluid density, the balance,
evaporation, retained droplets, collection losses or the timing instrument.
`flow_if_constant_evaporation` only shifts an apparent flow by a constant
mass-loss rate you supply. It does not estimate that rate.
Define a measurement procedure and an uncertainty budget before making a
hardware accuracy claim. This software has been tested on constructed traces;
no physical pump has been calibrated in this repository.

This is analysis for research instruments, not medical infusion control. It
produces no motor commands and does not automatically change a pump calibration.

## Reproduce the fixture

```bash
perfusioncal demo --seed 0 --out artifacts/new_synthetic_measurements.csv
```

The fixture and code are MIT licensed; no personal or biological data are used.
