"""Module for thermodynamic property calculations in combustion systems.

This module provides classes and methods for calculating thermodynamic properties such as enthalpy,
entropy, heat capacity, and Gibbs free energy for individual substances and entire systems. The
calculations are based on thermodynamic coefficients derived from:

"Thermodynamic and Thermophysical Properties of Combustion Products. Volume 1"
by V.P. Glushko.

These coefficients are used to compute thermodynamic properties for both gas-phase and condensed-phase
species.

Key functionalities include:
- Calculating thermodynamic properties (enthalpy, entropy, heat capacity, Gibbs free energy) for
  individual substances using coefficients from Glushko's work.
- Aggregating thermodynamic properties for a system of multiple substances, considering partial pressures
  for gas-phase species and normalization for condensed-phase species.

Classes:
    ThermodynamicIndividualCalculator: Provides static methods for thermodynamic property calculations
        for individual substances.
    ThermodynamicSystemContext: Represents the context for thermodynamic system calculations, including
        temperature, pressure, substance amounts, coefficients, and phase information.
    ThermodynamicSystemCalculator: Provides static methods for aggregating thermodynamic properties
        across a system of multiple substances.
"""

import numpy as np

from dataclasses import dataclass

from constants import CALORIE_TO_JOULES, GAS_CONSTANT, STANDARD_PRESSURE

class ThermodynamicIndividualCalculator:
    """Provides static methods for thermodynamic property calculations for individual substances.

    These methods use thermodynamic coefficients derived from "Thermodynamic and Thermophysical
    Properties of Combustion Products. Volume 1" by V.P. Glushko to compute enthalpy, entropy,
    heat capacity, and Gibbs free energy for a single substance at a given temperature and pressure.
    """

    @staticmethod
    def calculate_enthalpy(coefficients: np.ndarray, temperature: float) -> float:
        """Calculate molar enthalpy using thermodynamic coefficients from Glushko's work.

        Args:
            coefficients (np.ndarray): Array of 9 thermodynamic coefficients derived from Glushko's work.
            temperature (float): Current temperature in Kelvin.

        Returns:
            float: Enthalpy in J/mol.

        Notes:
            - The coefficients are assumed to be in calorie-based units and are converted to joules
              using the `CALORIE_TO_JOULES` constant.
            - The coefficients correspond to polynomial fits described in Glushko's work.
        """
        t = temperature * 1e-3
        return CALORIE_TO_JOULES * (
            coefficients[1]
            + coefficients[2] * t
            + coefficients[3] * t**2
            + coefficients[4] * t**3
            + coefficients[5] * t**4
            + coefficients[6] * t**5
            + coefficients[7] * t**6
            + coefficients[8] * t**7
        )

    @staticmethod
    def calculate_entropy(coefficients: np.ndarray, temperature: float, partial_pressure: float = 0.0) -> float:
        """Calculate molar entropy with optional pressure correction for gas-phase species.

        Args:
            coefficients (np.ndarray): Array of 9 thermodynamic coefficients derived from Glushko's work.
            temperature (float): Current temperature in Kelvin.
            partial_pressure (float): Partial pressure in Pa (for gases). Default is 0.0.

        Returns:
            float: Entropy in J/(mol·K).

        Notes:
            - If `partial_pressure` is provided, the entropy of mixing is subtracted using the ideal gas law.
            - The coefficients correspond to polynomial fits described in Glushko's work.
        """
        t = temperature * 1e-3
        std_entropy = CALORIE_TO_JOULES * (
            coefficients[0]
            + 1e-3 * coefficients[2] * np.log(t)
            + 1e-3 * (
                2 * coefficients[3] * t
                + 1.5 * coefficients[4] * t**2
                + (4/3) * coefficients[5] * t**3
                + 1.25 * coefficients[6] * t**4
                + 1.2 * coefficients[7] * t**5
                + (7/6) * coefficients[8] * t**6
            )
        )

        if partial_pressure > 0:
            std_entropy -= GAS_CONSTANT * np.log(partial_pressure / STANDARD_PRESSURE)

        return std_entropy

    @staticmethod
    def calculate_heat_capacity(coefficients: np.ndarray, temperature: float) -> float:
        """Calculate heat capacity using thermodynamic coefficients from Glushko's work.

        Args:
            coefficients (np.ndarray): Array of 9 thermodynamic coefficients derived from Glushko's work.
            temperature (float): Current temperature in Kelvin.

        Returns:
            float: Heat capacity in J/(mol·K).

        Notes:
            - The coefficients correspond to polynomial fits described in Glushko's work.
        """
        t = temperature * 1e-3
        return CALORIE_TO_JOULES * 1e-3 * (
            coefficients[2]
            + 2 * coefficients[3] * t
            + 3 * coefficients[4] * t**2
            + 4 * coefficients[5] * t**3
            + 5 * coefficients[6] * t**4
            + 6 * coefficients[7] * t**5
            + 7 * coefficients[8] * t**6
        )

    @staticmethod
    def calculate_gibbs_energy(enthalpy: float, entropy: float, temperature: float) -> float:
        """Calculate Gibbs free energy using enthalpy and entropy.

        Args:
            enthalpy (float): Enthalpy in J/mol.
            entropy (float): Entropy in J/(mol·K).
            temperature (float): Current temperature in Kelvin.

        Returns:
            float: Gibbs free energy in J/mol.
        """
        return enthalpy - temperature * entropy

@dataclass
class ThermodynamicSystemContext:
    """Represents the context for thermodynamic system calculations.

    Attributes:
        temperature (float): System temperature in Kelvin.
        pressure (float): System pressure in Pascals (Pa).
        substance_amounts (np.ndarray): Array of substance amounts (in moles) for each species.
        coefficients (np.ndarray): Array of thermodynamic coefficients derived from Glushko's work
            for each species.
        is_condensed (np.ndarray): Boolean array indicating whether each species is in the condensed phase.
    """
    temperature: float
    pressure: float
    substance_amounts: np.ndarray
    coefficients: np.ndarray
    is_condensed: np.ndarray

class ThermodynamicSystemCalculator:
    """Provides static methods for aggregating thermodynamic properties across a system of multiple substances.

    These methods calculate total enthalpy, entropy, heat capacity, and Gibbs free energy for a system,
    accounting for both gas-phase and condensed-phase species. The thermodynamic coefficients are derived
    from "Thermodynamic and Thermophysical Properties of Combustion Products. Volume 1" by V.P. Glushko.
    """

    @staticmethod
    def calculate_enthalpy(context: ThermodynamicSystemContext) -> float:
        """Calculate total enthalpy of the system.

        Args:
            context (ThermodynamicSystemContext): Context containing system parameters.

        Returns:
            float: Total enthalpy in J.
        """
        substance_amounts = context.substance_amounts
        coefficients = context.coefficients
        temperature = context.temperature
        total_enthalpy = 0.0

        for i in range(len(substance_amounts)):
            q = substance_amounts[i]
            coeffs = coefficients[i]
            h = ThermodynamicIndividualCalculator.calculate_enthalpy(coeffs, temperature)
            total_enthalpy += q * h

        return total_enthalpy

    @staticmethod
    def calculate_entropy(context: ThermodynamicSystemContext) -> float:
        """Calculate total entropy of the system.

        Args:
            context (ThermodynamicSystemContext): Context containing system parameters.

        Returns:
            float: Total entropy in J/K.
        """
        substance_amounts = context.substance_amounts
        coefficients = context.coefficients
        is_condensed = context.is_condensed
        temperature = context.temperature
        pressure = context.pressure
        total_entropy = 0.0
        total_gas_moles = np.sum(substance_amounts * ~is_condensed)

        for i in range(len(substance_amounts)):
            q = substance_amounts[i]
            coeffs = coefficients[i]
            condensed_flag = is_condensed[i]
            partial_pressure = ThermodynamicSystemCalculator._calculate_partial_pressure(
                pressure, q, total_gas_moles, condensed_flag
            )
            s = ThermodynamicIndividualCalculator.calculate_entropy(coeffs, temperature, partial_pressure)
            total_entropy += q * s

        return total_entropy

    @staticmethod
    def calculate_heat_capacity(context: ThermodynamicSystemContext) -> float:
        """Calculate total heat capacity of the system.

        Args:
            context (ThermodynamicSystemContext): Context containing system parameters.

        Returns:
            float: Total heat capacity in J/K.
        """
        substance_amounts = context.substance_amounts
        coefficients = context.coefficients
        is_condensed = context.is_condensed
        temperature = context.temperature
        total_heat_capacity = 0.0
        total_gas_moles = np.sum(substance_amounts * ~is_condensed)

        for i in range(len(substance_amounts)):
            q = substance_amounts[i]
            coeffs = coefficients[i]
            condensed_flag = is_condensed[i]
            n = ThermodynamicSystemCalculator._calculate_substance_amount_over_total_gas_moles(
                q, total_gas_moles, condensed_flag
            )
            cp = ThermodynamicIndividualCalculator.calculate_heat_capacity(coeffs, temperature)
            total_heat_capacity += n * cp

        return total_heat_capacity

    @staticmethod
    def calculate_gibbs_energy(context: ThermodynamicSystemContext) -> float:
        """Calculate total Gibbs free energy of the system.

        Args:
            context (ThermodynamicSystemContext): Context containing system parameters.

        Returns:
            float: Total Gibbs free energy in J.
        """
        substance_amounts = context.substance_amounts
        coefficients = context.coefficients
        is_condensed = context.is_condensed
        temperature = context.temperature
        pressure = context.pressure
        total_gibbs_energy = 0.0
        total_gas_moles = np.sum(substance_amounts * ~is_condensed)

        for i in range(len(substance_amounts)):
            q = substance_amounts[i]
            coeffs = coefficients[i]
            condensed_flag = is_condensed[i]
            partial_pressure = ThermodynamicSystemCalculator._calculate_partial_pressure(
                pressure, q, total_gas_moles, condensed_flag
            )
            h = ThermodynamicIndividualCalculator.calculate_enthalpy(coeffs, temperature)
            s = ThermodynamicIndividualCalculator.calculate_entropy(coeffs, temperature, partial_pressure)
            gibbs_energy = ThermodynamicIndividualCalculator.calculate_gibbs_energy(h, s, temperature)
            total_gibbs_energy += q * gibbs_energy

        return total_gibbs_energy

    @staticmethod
    def _calculate_partial_pressure(total_pressure: float, substance_amount: float, total_gas_moles: float, is_condensed: bool) -> float:
        """Calculate the partial pressure of a gas-phase species.

        Args:
            total_pressure (float): Total system pressure in Pascals (Pa).
            substance_amount (float): Amount of the species in moles.
            total_gas_moles (float): Total moles of gas-phase species.
            is_condensed (bool): Flag indicating whether the species is in the condensed phase.

        Returns:
            float: Partial pressure in Pascals (Pa).
        """
        return total_pressure * (substance_amount / total_gas_moles) if not is_condensed else 0.0

    @staticmethod
    def _calculate_substance_amount_over_total_gas_moles(substance_amount: float, total_gas_moles: float, is_condensed: bool) -> float:
        """Calculate the normalized amount of a species relative to total gas-phase moles.

        Args:
            substance_amount (float): Amount of the species in moles.
            total_gas_moles (float): Total moles of gas-phase species.
            is_condensed (bool): Flag indicating whether the species is in the condensed phase.

        Returns:
            float: Normalized amount (dimensionless).
        """
        return substance_amount / total_gas_moles if not is_condensed else 0.0
