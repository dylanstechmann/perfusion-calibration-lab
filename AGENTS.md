# Agent instructions — perfusion-calibration-lab

Work only in this repository. Offline gravimetric analysis: mass vs time, a
free intercept, flow from an explicit density. The synthetic fixture is about
7% underdelivery on purpose. No physical pump has been calibrated here.

## Do not

- Split one continuous trace into several `run_id`s and call them replicates.
- Hide flagged runs (low R², nonpositive flow, early/late slope change).
- Emit motor commands or write a calibration back to `open-perfusion-rig`.
- Treat the bootstrap interval as including density, evaporation, or balance bias. It only resamples runs.
- Claim water is exactly 1 mg/µL at the user's temperature.

## First commands

```bash
python -m pip install -e .
python -m unittest discover -s tests -v
perfusioncal analyze examples/synthetic_measurements.csv --density-mg-ul 1 --out artifacts/agent-demo
```

Use a new `--out` directory.

## Improve, in this order

1. If the estimator changes, add a fixture with a known slope and density and assert flow within a tight tolerance.
2. `flow_if_constant_evaporation` is the evaporation sensitivity. Do not add another copy of it, and do not fold it into the bootstrap.
3. Do not add a hardware driver.

## Done when

Tests pass and the report still says this is not a physical calibration.
