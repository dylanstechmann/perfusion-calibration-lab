# Physical measurement recording template

`physical_measurements.template.csv` contains only the input header. It is a
blank format template, not a measurement record. Add one row for each balance
reading:

- `run_id`: one identifier for a complete, independently repeated pump run.
  Give each new run a new ID; do not label slices of one trace as independent.
- `time_s`: elapsed seconds from that run's start, increasing strictly within
  each `run_id`.
- `mass_mg`: balance reading in milligrams, including the original tare or
  collection-vessel offset; the fit estimates a free intercept.
- `target_flow_ul_min`: commanded setpoint recorded for that run, in µL/min.

Keep the instrument export and contemporaneous run notes as source records.
Record the pump and firmware identity, syringe and fluid, fluid temperature,
balance identity and certificate, timebase source, collection setup, run date,
and operator in the lab's record system. Preserve the original instrument
files and note their hashes when available. Use one fresh output directory per
analysis.

Measure density under the actual fluid and temperature conditions. Keep
balance-gain, balance-drift, matched evaporation-blank, and timebase records
separate and enter them in the uncertainty budget only when supported by those
records. This analyzer does not calculate evaporation from a blank trace or
operate a pump. It analyzes the supplied CSV offline. No physical values belong
in this template until measured.
