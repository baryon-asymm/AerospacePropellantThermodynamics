"""Module for reading and parsing JSON input files.

This module provides functions to load combustion products and the input substance data
from JSON files and convert them into domain objects. The thermodynamic properties are based on
"Thermodynamic and Thermophysical Properties of Combustion Products. Volume 1" by V.P. Glushko.
"""

import json
from typing import List
from models import ReactionProduct, TemperatureRange, PropellantComposition

def load_combustion_products(filepath: str) -> List[ReactionProduct]:
    """Loads combustion products from a JSON file and returns a list of ReactionProduct objects.
    
    The JSON file must contain an array of objects, each with the following structure:
        {
          "formula": "O",
          "coefficients": [
            45.168916,
            58008.607,
            5353.7423,
            -412.44632,
            246.19247,
            -86.140481,
            17.415382,
            -1.8288189,
            0.077299666
          ],
          "phase": "gas",
          "temperature_range": {
            "min": 1000,
            "max": 5000
          }
        }
    These coefficients correspond to thermodynamic properties derived from tabulated data
    in Glushko's work and are used to compute properties such as enthalpy, entropy, and heat capacity.
    
    Args:
        filepath (str): The path to the combustion products JSON file.
    
    Returns:
        List[ReactionProduct]: A list of ReactionProduct domain objects.
    
    Raises:
        FileNotFoundError: If the specified file does not exist.
        json.JSONDecodeError: If the JSON file is malformed.
        ValueError: If required keys are missing in the JSON data.
    """
    products: List[ReactionProduct] = []
    try:
        with open(filepath, "r", encoding="utf-8") as file:
            data = json.load(file)
    except FileNotFoundError:
        raise FileNotFoundError(f"The file '{filepath}' was not found.")
    except json.JSONDecodeError:
        raise json.JSONDecodeError("The file contains invalid JSON.", "", 0)
    for item in data:
        try:
            temp_range = TemperatureRange(
                min=item["temperature_range"]["min"],
                max=item["temperature_range"]["max"]
            )
            is_condensed = (item["phase"] == "condensed")
            product = ReactionProduct(
                formula=item["formula"],
                coefficients=item["coefficients"],
                phase=item["phase"],
                temperature_range=temp_range,
                is_condensed=is_condensed
            )
            products.append(product)
        except KeyError as error:
            raise ValueError(f"Missing key in combustion product data: {error}")
    return products

def load_input_substance(filepath: str) -> PropellantComposition:
    """Loads input substance data from a JSON file and returns a PropellantComposition object.
    
    The JSON file must contain an object with the following structure:
        {
            "enthalpy": 1.23e-10,
            "composition": {
                "C": 45.67,
                "H": 67.98,
                ...
            }
        }
    The enthalpy value corresponds to the thermodynamic property of the substance, and the composition
    represents its chemical composition.
    
    Args:
        filepath (str): Path to the input substance JSON file.
    
    Returns:
        PropellantComposition: Object containing enthalpy and chemical composition.
    
    Raises:
        FileNotFoundError: If specified file does not exist.
        json.JSONDecodeError: If JSON file is malformed.
        ValueError: If required keys are missing in JSON data.
    """
    try:
        with open(filepath, "r", encoding="utf-8") as file:
            data = json.load(file)
    except FileNotFoundError:
        raise FileNotFoundError(f"The file '{filepath}' was not found.")
    except json.JSONDecodeError:
        raise json.JSONDecodeError(f"Invalid JSON in {filepath}", e.doc, e.pos) from e
    
    try:
        return PropellantComposition(
            enthaly=data["enthalpy"],
            composition=data["composition"]
        )
    except KeyError as e:
        raise ValueError(f"Missing required field: {e}") from e
