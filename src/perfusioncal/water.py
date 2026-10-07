"""Temperature-dependent density of pure air-free water.

A gravimetric flow estimate divides a mass rate by a density. Using 1 mg/µL is a
convenient approximation that is wrong by about 0.2% at room temperature, and it
hides the fact that density was assumed rather than determined. This module
offers the published formula instead, for pure water only.

Reference: M. Tanaka, G. Girard, R. Davis, A. Peuto and N. Bignell, "Recommended
table for the density of water between 0 °C and 40 °C based on recent
experimental reports", Metrologia 38 (2001) 301-309,
https://doi.org/10.1088/0026-1394/38/4/3. The CIPM-recommended form used here is

    rho(t) = a5 * [1 - (t + a1)^2 * (t + a2) / (a3 * (t + a4))]

for standard atmospheric pressure, air-free water, with t in degrees Celsius and
rho in kg/m^3.

What this does **not** cover: culture medium, buffered saline, water carrying
dissolved gases, or any solution with solutes; pressure other than standard
atmospheric; and the air-buoyancy correction that applies when weighing on a
balance in air. For any of those, measure the density of the fluid you actually
pumped and pass it explicitly.
"""

from __future__ import annotations

# Tanaka et al. 2001, Table 1 coefficients.
A1_C = -3.983035
A2_C = 301.797
A3_C2 = 522528.9
A4_C = 69.34881
A5_KG_M3 = 999.974950

MINIMUM_TEMPERATURE_C = 0.0
MAXIMUM_TEMPERATURE_C = 40.0
CITATION = (
    "Tanaka et al., Metrologia 38 (2001) 301-309, "
    "https://doi.org/10.1088/0026-1394/38/4/3"
)
SCOPE = (
    "Pure air-free water at standard atmospheric pressure, 0-40 degrees Celsius. "
    "Not valid for culture medium, buffered saline, gas-saturated water or any solution "
    "with solutes, and it does not include the air-buoyancy correction for weighing in air."
)


def water_density_kg_m3(temperature_c: float) -> float:
    """Density of pure air-free water at ``temperature_c``, in kg/m^3."""
    if isinstance(temperature_c, bool) or not isinstance(temperature_c, (int, float)):
        raise ValueError("temperature must be a number in degrees Celsius")
    temperature = float(temperature_c)
    if temperature != temperature or temperature in (float("inf"), float("-inf")):
        raise ValueError("temperature must be finite")
    if not MINIMUM_TEMPERATURE_C <= temperature <= MAXIMUM_TEMPERATURE_C:
        raise ValueError(
            f"temperature {temperature} C is outside the published {MINIMUM_TEMPERATURE_C}-"
            f"{MAXIMUM_TEMPERATURE_C} C range of the Tanaka formula; measure the density instead"
        )
    numerator = (temperature + A1_C) ** 2 * (temperature + A2_C)
    denominator = A3_C2 * (temperature + A4_C)
    return A5_KG_M3 * (1.0 - numerator / denominator)


def water_density_mg_ul(temperature_c: float) -> float:
    """Density of pure air-free water at ``temperature_c``, in mg/µL.

    1 kg/m^3 is 1e-3 mg/µL, so the familiar 1 mg/µL corresponds to 1000 kg/m^3.
    """
    return water_density_kg_m3(temperature_c) * 1e-3


def water_density_record(temperature_c: float) -> dict:
    """The density plus everything a reviewer needs to judge whether it applies."""
    density_kg_m3 = water_density_kg_m3(temperature_c)
    density_mg_ul = density_kg_m3 * 1e-3
    approximation_error = (1.0 - density_mg_ul) / density_mg_ul
    return {
        "source": "tanaka_2001_pure_water",
        "citation": CITATION,
        "fluid": "pure air-free water at standard atmospheric pressure",
        "temperature_c": float(temperature_c),
        "valid_temperature_range_c": [MINIMUM_TEMPERATURE_C, MAXIMUM_TEMPERATURE_C],
        "density_kg_m3": density_kg_m3,
        "density_mg_ul": density_mg_ul,
        "relative_error_of_unit_density_approximation": approximation_error,
        "scope": SCOPE,
        "limitations": [
            "This is a published formula for pure water, not a measurement of the fluid that was "
            "actually pumped. A medium, buffer or gas-saturated solution has a different density.",
            "No air-buoyancy correction is applied. Weighing in air biases an apparent mass by "
            "roughly 0.1% for water unless the balance was calibrated to compensate.",
            "A recorded temperature is only as good as the thermometer and its placement; a bath "
            "reading is not necessarily the temperature of fluid leaving the needle.",
            "Using the correct density does not make a flow estimate a physical pump calibration.",
        ],
    }
