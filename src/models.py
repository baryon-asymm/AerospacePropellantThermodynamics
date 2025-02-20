"""Module for thermodynamic equilibrium calculations in combustion processes.

This module provides classes to represent chemical formulas, reaction products,
and temperature ranges for use in thermodynamic equilibrium computations based on
thermodynamic and thermophysical properties of combustion products as described in:
"Thermodynamic and Thermophysical Properties of Combustion Products. Volume 1"
by V.P. Glushko.

"""

from dataclasses import dataclass

@dataclass(frozen=True)
class TemperatureRange:
    """Represents a valid temperature range for thermodynamic properties.

    Attributes:
        min (float): Minimum temperature in the range (K).
        max (float): Maximum temperature in the range (K).
    """
    min: float
    max: float


@dataclass(frozen=True)
class PropellantComposition:
    """Represents the chemical composition of a rocket propellant.

    Attributes:
        enthaly (float): Enthalpy of the propellant (J/kg).
        composition (dict[str, float]): Mapping of element symbols to their
            mass fractions or stoichiometric quantities. Example: {'C': 0.5, 'H': 0.2}.
    """
    enthaly: float
    composition: dict[str, float]


@dataclass(frozen=True)
class ReactionProduct:
    """Represents a product species in a combustion reaction, using thermodynamic
    and thermophysical properties from "Thermodynamic and Thermophysical Properties
    of Combustion Products. Volume 1" by V.P. Glushko.

    Attributes:
        formula (str): Chemical formula of the product (e.g., 'CO2', 'H2O').
        coefficients (list[float]): Thermodynamic coefficients for property calculations.
            These coefficients are derived from tabulated data in Glushko's work and
            are used to compute properties such as enthalpy, entropy, and heat capacity.
        phase (str): Physical phase of the product. Common values: 'gas', 'condensed'.
        temperature_range (TemperatureRange): Temperature validity range for coefficients.
        is_condensed (bool): Flag indicating if the product is in a condensed phase.
    """
    formula: str
    coefficients: list[float]
    phase: str
    temperature_range: TemperatureRange
    is_condensed: bool
