"""Module for writing thermodynamic optimization results to a JSON file.

This module provides utility functions to structure and write the results of
thermodynamic optimization into a JSON file. The output includes input substance parameters,
combustion products (formula, phase, and moles), and the resulting optimal temperature.
"""

import json

from typing import List

from models import PropellantComposition, ReactionProduct
from calculators import ThermodynamicSystemContext
from thermodynamic_properties import ThermodynamicPropertiesContext


def prepare_combustion_products(
    context: ThermodynamicSystemContext,
    filtered_products: List[ReactionProduct]
) -> List[dict]:
    """
    Prepares a list of combustion products with their formula, phase, and moles.

    Args:
        context (ThermodynamicSystemContext): Optimized thermodynamic system context.
        filtered_products (List[ReactionProduct]): List of filtered reaction products.

    Returns:
        List[dict]: List of dictionaries representing combustion products.
    """
    combustion_products = []
    for product, amount in zip(filtered_products, context.substance_amounts):
        combustion_products.append({
            "formula": product.formula,
            "phase": product.phase,
            "moles": float(amount)  # Ensure amount is converted to float for JSON serialization
        })
    return combustion_products


def write_to_json(
    file_path: str,
    propellant: PropellantComposition,
    total_mass: float,
    context: ThermodynamicPropertiesContext,
    filtered_products: List[ReactionProduct]
):
    """
    Writes the optimization result to a JSON file.

    Args:
        file_path (str): Path to the output JSON file.
        propellant (PropellantComposition): Input propellant composition.
        total_mass (float): Total mass of the propellant in kilograms (kg).
        context (ThermodynamicSystemContext): Optimized thermodynamic system context.
        filtered_products (List[ReactionProduct]): List of filtered reaction products.
    """
    # Prepare combustion products
    combustion_products_data = prepare_combustion_products(context.thermodynamic_system_context, filtered_products)

    # Prepare the overall result
    result = {
        "pressure": context.thermodynamic_system_context.pressure,
        "temperature": float(context.thermodynamic_system_context.temperature),  # Ensure temperature is converted to float
        "specific_heat_capacity_volumetric": float(context.specific_heat_capacity_volumetric),  # Ensure heat capacity is converted to float
        "gas_average_molar_mass": float(context.gas_average_molar_mass),  # Ensure molar mass is converted to float
        "propellant": {
            "enthalpy": propellant.enthalpy,
            "composition": propellant.composition,
            "total_mass_kg": float(total_mass)  # Add total mass to the output
        },
        "combustion_products": combustion_products_data
    }

    # Write to JSON file
    with open(file_path, "w", encoding="utf-8") as file:
        json.dump(result, file, indent=4)
