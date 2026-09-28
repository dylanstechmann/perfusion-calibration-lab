# Perfusion calibration analysis

Input SHA-256: `0878d191bc76b93435b2e5a337fd62f9cb3c641d0ce25ca6fef0034d6dc7a030`

Density used: 1 mg/µL. Startup interval discarded: 0 s.

| Target (µL/min) | Runs | Mean measured (µL/min) | Error (%) | Repeat SD (µL/min) | 95% run-bootstrap interval (µL/min) |
|---:|---:|---:|---:|---:|---:|
| 25.000 | 5 | 23.238 | -7.05 | 0.100 | [23.151, 23.309] |
| 50.000 | 5 | 46.476 | -7.05 | 0.399 | [46.184, 46.814] |
| 100.000 | 5 | 92.855 | -7.15 | 0.572 | [92.427, 93.282] |

## Run diagnostics

- synthetic-100-0: no diagnostic flags
- synthetic-100-1: no diagnostic flags
- synthetic-100-2: no diagnostic flags
- synthetic-100-3: no diagnostic flags
- synthetic-100-4: no diagnostic flags
- synthetic-25-0: no diagnostic flags
- synthetic-25-1: no diagnostic flags
- synthetic-25-2: no diagnostic flags
- synthetic-25-3: no diagnostic flags
- synthetic-25-4: no diagnostic flags
- synthetic-50-0: no diagnostic flags
- synthetic-50-1: no diagnostic flags
- synthetic-50-2: no diagnostic flags
- synthetic-50-3: no diagnostic flags
- synthetic-50-4: no diagnostic flags

## Interpretation

- Offline research analysis. No hardware commands are generated.
- This report alone is not a physical pump calibration or acceptance decision.
- Intervals resample independent run slopes, not serially correlated readings.
- Fewer than three runs: no interval. Small repeat counts give unstable intervals.
- Density, balance calibration, evaporation and collection losses are not included in uncertainty.
- R-squared <0.95 and >20% half-run drift are diagnostic flags, not acceptance standards.
- Flagged runs remain in summaries; inspect them before interpreting mean flow.
