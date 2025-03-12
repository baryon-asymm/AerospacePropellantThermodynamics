"""Module for calculating thermodynamic properties of an optimized combustion system.

This module provides the `ThermodynamicPropertiesCalculator` class, which computes and displays
various thermodynamic properties of a combustion system after optimization. These properties include
total moles, condensed and gas-phase contributions, Gibbs free energy, enthalpy, entropy, heat capacity,
and derived quantities such as specific gas constant, heat ratio, and volume heat capacity.

The calculations are based on the context provided by the `ThermodynamicSystemContext` class.
"""

import numpy as np
from dataclasses import dataclass

from calculators import ThermodynamicSystemCalculator, ThermodynamicSystemContext
from utils import parse_chemical_formula, compute_molar_mass
from constants import GAS_CONSTANT

@dataclass(frozen=True)
class ThermodynamicPropertiesContext:
    """Data class for storing thermodynamic properties of the optimized system.

    This data class stores various thermodynamic properties of a combustion system after optimization.
    The properties include the thermodynamic system context, heat capacity, and average molar mass of gas products.

    Attributes:
        thermodynamic_system_context (ThermodynamicSystemContext): Context containing system parameters such as
            temperature, pressure, substance amounts, coefficients, and phase information.
        specific_heat_capacity_volumetric (float): Volumetric heat capacity of the optimized system.
        gas_average_molar_mass (float): Average molar mass of gas products.
    """
    thermodynamic_system_context: ThermodynamicSystemContext
    specific_heat_capacity_volumetric: float
    gas_average_molar_mass: float

class ThermodynamicPropertiesCalculator:
    """Class for calculating and displaying thermodynamic properties of the optimized system.

    This class computes various thermodynamic properties of a combustion system, including total moles,
    condensed and gas-phase contributions, Gibbs free energy, enthalpy, entropy, heat capacity, and derived
    quantities such as specific gas constant, heat ratio, and volume heat capacity.

    Attributes:
        context (ThermodynamicSystemContext): Context containing system parameters such as temperature,
            pressure, substance amounts, coefficients, and phase information.
        filtered_products (list): List of filtered reaction products used in the optimization.

    Methods:
        calculate_and_display_properties(): Calculates and displays all thermodynamic properties.
        _calculate_weight_fraction_of_condensed_products() -> float: Calculates the total weight fraction
            of condensed products.
    """

    def __init__(self,
                 context,
                 filtered_products: list):
        """
        Args:
            context (ThermodynamicSystemContext): Context containing system parameters such as temperature,
                pressure, substance amounts, coefficients, and phase information.
            filtered_products (list): List of filtered reaction products used in the optimization.
        """
        self.context = context
        self.filtered_products = filtered_products

    def calculate_and_display_properties(self) -> ThermodynamicPropertiesContext:
        """Calculate and display all thermodynamic properties of the optimized system.

        This method computes and prints the following thermodynamic properties:
        - Total moles of substances.
        - Total condensed moles and gas-phase moles.
        - Gibbs free energy, enthalpy, entropy, and heat capacity.
        - Heat capacity derived from enthalpy.
        - Weight fraction of condensed products.
        - Average molar mass of gas products.
        - Specific gas constant, heat capacity, heat ratio, and volume heat capacity.

        Notes:
            - The results are printed directly to the console.
            - The context's `substance_amounts` attribute contains the optimized amounts.
        """

        # Basic properties
        total_moles = self.context.substance_amounts.sum()
        total_condensed_moles = np.sum(self.context.substance_amounts * self.context.is_condensed)
        total_gas_moles = total_moles - total_condensed_moles

        print("Total moles:", total_moles)
        print("Total condensed moles:", total_condensed_moles)

        # Gibbs energy, enthalpy, entropy, and heat capacity
        final_gibbs_energy = ThermodynamicSystemCalculator.calculate_gibbs_energy(self.context)
        final_enthalpy = ThermodynamicSystemCalculator.calculate_enthalpy(self.context)
        final_entropy = ThermodynamicSystemCalculator.calculate_entropy(self.context)
        final_heat_capacity = ThermodynamicSystemCalculator.calculate_heat_capacity(self.context)

        print("Final Gibbs Energy:", final_gibbs_energy)
        print("Final Enthalpy:", final_enthalpy)
        print("Final Entropy:", final_entropy)
        print("Final Heat Capacity (molar):", final_heat_capacity)

        # Heat capacity from enthalpy
        delta_temperature = 1e-3
        self.context.temperature += delta_temperature
        c_p_from_enthalpy = (ThermodynamicSystemCalculator.calculate_enthalpy(self.context) - final_enthalpy) / delta_temperature
        print("Final Heat Capacity (from enthalpy):", c_p_from_enthalpy)

        # Weight fraction of condensed products
        total_weight_fraction = self._calculate_weight_fraction_of_condensed_products()
        print("Total weight fraction of condensed products:", total_weight_fraction)

        # Average molar mass of gas products
        average_molar_mass = (1 - total_weight_fraction) / total_gas_moles
        print("Average molar mass of gas products:", average_molar_mass)

        # Specific gas constant and heat capacity
        specific_gas_constant = GAS_CONSTANT / average_molar_mass
        specific_heat_capacity = final_heat_capacity / average_molar_mass
        specific_heat_ratio = specific_heat_capacity / (specific_heat_capacity - specific_gas_constant)
        volume_heat_capacity = specific_heat_capacity - specific_gas_constant

        print("Specific gas constant:", specific_gas_constant)
        print("Heat capacity:", specific_heat_capacity)
        print("Specific heat ratio:", specific_heat_ratio)
        print("Volume heat capacity:", volume_heat_capacity)

        return ThermodynamicPropertiesContext(
            thermodynamic_system_context=self.context,
            specific_heat_capacity_volumetric=volume_heat_capacity,
            gas_average_molar_mass=average_molar_mass
        )

    def _calculate_weight_fraction_of_condensed_products(self) -> float:
        """Calculate the total weight fraction of condensed products.

        This method computes the total weight fraction of condensed products by summing the contributions
        of each condensed species. The contribution of each species is calculated as the product of its
        molar mass and its optimized amount.

        Returns:
            float: Total weight fraction of condensed products.

        Notes:
            - The molar mass of each species is computed using its chemical formula and the `compute_molar_mass`
              utility function.
            - Only condensed-phase species (as indicated by `is_condensed`) are included in the calculation.
        """
        total_weight_fraction = 0.0
        for product, amount, is_condensed_flag in zip(self.filtered_products, self.context.substance_amounts, self.context.is_condensed):
            if is_condensed_flag:
                elements = parse_chemical_formula(product.formula)
                molar_mass = compute_molar_mass(elements)
                total_weight_fraction += molar_mass * amount
        return total_weight_fraction
