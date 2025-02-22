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

    Args:
        min (float): Minimum temperature in Kelvin.
        max (float): Maximum temperature in Kelvin.

    Note:
        This class is immutable due to the `frozen=True` decorator.
    """
    min: float
    max: float


@dataclass(frozen=True)
class PropellantComposition:
    """Represents the chemical composition of a rocket propellant, expressed as a conditional
    chemical formula for 1 kg of propellant. The composition specifies the molar quantities
    of each element in the propellant, normalized such that the total mass (computed as the
    sum of molar quantities multiplied by their respective molar masses) equals 1 kg.

    Attributes:
        enthalpy (float): Enthalpy of the propellant in joules per kilogram (J/kg).
        composition (dict[str, float]): Mapping of element symbols to their molar quantities
            (in moles) in 1 kg of propellant. The molar quantities are normalized such that
            the sum of the products of these quantities and their respective molar masses
            equals 1 kg.

    Args:
        enthalpy (float): Enthalpy of the propellant in J/kg.
        composition (dict[str, float]): Mapping of element symbols to their molar quantities
            (in moles) in 1 kg of propellant. For example:
            {
                "C": 11.92,
                "H": 29.58,
                "O": 22.85,
                "N": 13.25,
                "Cl": 2.12,
                "Al": 7.41
            }
            The values represent the number of moles of each element in 1 kg of propellant.

    Note:
        - The molar quantities in `composition` must satisfy the condition:
          sum(moles * molar_mass for element, moles in composition.items()) == 1 kg.
        - This class is immutable due to the `frozen=True` decorator.
    """
    enthalpy: float
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
            are typically used in polynomial fits for enthalpy, entropy, and heat capacity.
        phase (str): Physical phase of the product. Possible values: 'gas', 'condensed'.
        temperature_range (TemperatureRange): Temperature validity range for coefficients.
        is_condensed (bool): Flag indicating if the product is in a condensed phase.

    Args:
        formula (str): Chemical formula of the product.
        coefficients (list[float]): Thermodynamic coefficients for property calculations.
        phase (str): Physical phase of the product ('gas' or 'condensed').
        temperature_range (TemperatureRange): Valid temperature range for the coefficients.
        is_condensed (bool): Indicates if the product is in a condensed phase.

    Note:
        This class is immutable due to the `frozen=True` decorator.
    """
    formula: str
    coefficients: list[float]
    phase: str
    temperature_range: TemperatureRange
    is_condensed: bool
