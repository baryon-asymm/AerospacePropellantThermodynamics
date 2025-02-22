"""Module providing fundamental physical and thermodynamic constants.

This module defines key constants used in thermodynamic and combustion calculations.
These constants are provided in SI units for consistency and ease of use in scientific computations.

Constants:
    GAS_CONSTANT (float): The universal gas constant in joules per mole-kelvin [J/(mol·K)].
        Value: 8.31446261815324 J/(mol·K).
    STANDARD_PRESSURE (float): Standard atmospheric pressure in pascals [Pa].
        Value: 101325.0 Pa.
    CALORIE_TO_JOULES (float): Conversion factor from calories to joules [J/cal].
        Value: 4.184 J/cal.

Usage Example:
    >>> from constants import GAS_CONSTANT, STANDARD_PRESSURE, CALORIE_TO_JOULES
    >>> print(GAS_CONSTANT)  # Universal gas constant
    8.31446261815324
    >>> print(STANDARD_PRESSURE)  # Standard atmospheric pressure
    101325.0
    >>> print(CALORIE_TO_JOULES)  # Conversion factor
    4.184

Note:
    - All constants are provided in SI units to ensure compatibility with thermodynamic equations.
    - These constants are widely used in calculations involving ideal gas laws, energy conversions,
      and thermodynamic properties.
"""

GAS_CONSTANT = 8.31446261815324

STANDARD_PRESSURE = 101325.0

CALORIE_TO_JOULES = 4.184
