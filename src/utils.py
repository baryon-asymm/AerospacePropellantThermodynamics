"""Module for thermodynamic calculations and stoichiometric matrix construction.

This module provides utility functions for processing chemical formulas, filtering reaction products,
and constructing matrices required for thermodynamic equilibrium calculations. It leverages data from
`models`, `molar_masses`, and other constants to ensure accurate and consistent computations.

Key functionalities include:
- Parsing chemical formulas into elemental compositions.
- Filtering reaction products based on temperature and elemental compatibility.
- Constructing stoichiometric matrices and optimization matrices for thermodynamic calculations.
- Calculating the total mass of a propellant based on its composition.

Functions:
    create_elemental_vector(propellant_composition: PropellantComposition) -> np.ndarray:
        Creates a column vector representation of element quantities from propellant composition.
    parse_chemical_formula(formula: str) -> Dict[str, int]:
        Parses a chemical formula into its constituent elements with stoichiometric counts.
    compute_molar_mass(elements: Dict[str, int]) -> float:
        Calculates the molar mass of a compound from its elemental composition.
    filter_products_by_elements(reaction_products: List[ReactionProduct],
                               propellant_composition: PropellantComposition) -> List[ReactionProduct]:
        Filters reaction products to include only those composed of elements present in the propellant.
    construct_stoichiometric_matrix(filtered_products: List[ReactionProduct],
                                    propellant_composition: PropellantComposition) -> np.ndarray:
        Constructs a stoichiometric matrix aligned with the elements in the propellant composition.
    filter_and_construct_matrices(reaction_products: List[ReactionProduct],
                                 propellant: PropellantComposition,
                                 temperature: float) -> Tuple[List[ReactionProduct], np.ndarray]:
        Filters products by temperature and constructs the corresponding stoichiometric matrix.
    prepare_optimization_matrices(reaction_products: List[ReactionProduct],
                                  temperature_kelvin: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        Prepares matrices for thermodynamic optimization, including initial guesses and condensed phase flags.
    compute_total_mass(propellant_composition: PropellantComposition) -> float:
        Calculates the total mass of the propellant based on its elemental composition.
"""

import numpy as np

from typing import List, Dict, Tuple

from models import ReactionProduct, PropellantComposition
from molar_masses import ELEMENT_MOLAR_MASSES

def create_elemental_vector(propellant_composition: PropellantComposition) -> np.ndarray:
    """Creates a column vector representation of element quantities from propellant composition.

    The row order corresponds to the order of element symbols in `propellant_composition.composition.keys()`.

    Args:
        propellant_composition (PropellantComposition): Contains chemical composition data.

    Returns:
        np.ndarray: A 2D column vector (n × 1) with element quantities in moles.
    """
    element_quantities = list(propellant_composition.composition.values())
    return np.vstack(element_quantities)

def parse_chemical_formula(formula: str) -> Dict[str, int]:
    """Parses a chemical formula into its constituent elements with stoichiometric counts.

    Correctly handles:
    - Multi-character element symbols (e.g., Pb, Fe, Na).
    - Implicit counts (e.g., O = O1).
    - Complex formulas (e.g., Fe3O4, C6H5OH).

    Args:
        formula (str): Chemical formula following standard notation.

    Returns:
        Dict[str, int]: Dictionary mapping element symbols to their integer counts.

    Raises:
        ValueError: If the formula contains invalid characters or structure.
    """
    elements = {}
    i = 0
    n = len(formula)

    while i < n:
        # Element symbol must start with an uppercase letter
        if not formula[i].isupper():
            raise ValueError(f"Invalid character '{formula[i]}' at position {i} - elements must start with uppercase.")

        # Capture element symbol (1 uppercase + 0+ lowercase)
        element = [formula[i]]
        i += 1
        while i < n and formula[i].islower():
            element.append(formula[i])
            i += 1
        element = ''.join(element)

        # Capture numeric count
        count = []
        while i < n and formula[i].isdigit():
            count.append(formula[i])
            i += 1

        # Store element with count (default 1 if unspecified)
        elements[element] = int(''.join(count)) if count else 1

    return elements

def compute_molar_mass(elements: Dict[str, int]) -> float:
    """Calculates the molar mass of a compound from its elemental composition.

    Args:
        elements (Dict[str, int]): Dictionary of elements and their counts (from `parse_chemical_formula`).

    Returns:
        float: Molar mass in kg/mol.

    Raises:
        KeyError: If any element is not found in `ELEMENT_MOLAR_MASSES`.
    """
    molar_mass = 0.0

    for element, count in elements.items():
        if element not in ELEMENT_MOLAR_MASSES:
            raise KeyError(f"Element '{element}' not found in molar mass database.")
        molar_mass += ELEMENT_MOLAR_MASSES[element] * count

    return molar_mass

def filter_products_by_elements(reaction_products: List[ReactionProduct],
                               propellant_composition: PropellantComposition) -> List[ReactionProduct]:
    """Filters reaction products to include only those composed of elements present in the propellant.

    Args:
        reaction_products (List[ReactionProduct]): Candidate reaction products.
        propellant_composition (PropellantComposition): Reference composition for filtering.

    Returns:
        List[ReactionProduct]: Products whose ALL constituent elements exist in the propellant composition.

    Examples:
        Propellant contains ['C', 'H', 'O']:
        - Kept: H2O (H, O), CO2 (C, O)
        - Removed: NaCl (Na, Cl), CH3Cl (C, H, Cl)
    """
    propellant_elements = set(propellant_composition.composition.keys())

    return [
        product for product in reaction_products
        if propellant_elements.issuperset(parse_chemical_formula(product.formula))
    ]

def construct_stoichiometric_matrix(filtered_products: List[ReactionProduct],
                                    propellant_composition: PropellantComposition) -> np.ndarray:
    """Constructs a stoichiometric matrix aligned with the elements in the propellant composition.

    Matrix dimensions: (num_propellant_elements × num_filtered_products)
    Row order matches `propellant_composition.composition.keys()` order.

    Args:
        filtered_products (List[ReactionProduct]): Reaction products containing propellant elements.
        propellant_composition (PropellantComposition): Defines row order and relevant elements.

    Returns:
        np.ndarray: Stoichiometric coefficients matrix (n_elements × n_products).
    """
    elements = list(propellant_composition.composition.keys())
    matrix = np.zeros((len(elements), len(filtered_products)), dtype=np.float64)

    for col_idx, product in enumerate(filtered_products):
        counts = parse_chemical_formula(product.formula)
        for row_idx, element in enumerate(elements):
            matrix[row_idx, col_idx] = counts.get(element, 0.0)

    return matrix

def filter_and_construct_matrices(
    reaction_products: List[ReactionProduct],
    propellant: PropellantComposition,
    temperature: float
) -> Tuple[List[ReactionProduct], np.ndarray]:
    """Filters products by temperature and constructs the corresponding stoichiometric matrix.

    Args:
        reaction_products (List[ReactionProduct]): List of all reaction products.
        propellant (PropellantComposition): Propellant composition reference.
        temperature (float): Target temperature in Kelvin.

    Returns:
        Tuple[List[ReactionProduct], np.ndarray]: Filtered products and stoichiometric matrix.
    """
    # First filter by element compatibility
    element_filtered = filter_products_by_elements(reaction_products, propellant)

    # Then filter by temperature validity
    temp_filtered = [
        p for p in element_filtered
        if p.temperature_range.min <= temperature < p.temperature_range.max
    ]

    # Create stoichiometric matrix
    matrix = construct_stoichiometric_matrix(temp_filtered, propellant)

    return temp_filtered, matrix

def prepare_optimization_matrices(
    reaction_products: List[ReactionProduct],
    temperature_kelvin: float
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Prepares matrices for thermodynamic optimization, including initial guesses and condensed phase flags.

    Args:
        reaction_products (List[ReactionProduct]): List of ReactionProduct objects.
        temperature_kelvin (float): Target temperature in Kelvin.

    Returns:
        Tuple[np.ndarray, np.ndarray, np.ndarray]:
            - initial_guess: Initial substance quantities [mol] (1 × N).
            - coeff_matrix: Thermodynamic coefficients (N × 9).
            - is_condensed: Boolean flags for condensed phases (1 × N).

    Raises:
        ValueError: For invalid temperature or missing data.
    """
    if temperature_kelvin <= 0:
        raise ValueError(f"Invalid temperature: {temperature_kelvin} K")

    # Filter valid products
    valid_products = [
        p for p in reaction_products
        if p.temperature_range.min <= temperature_kelvin < p.temperature_range.max
    ]

    if not valid_products:
        raise ValueError(f"No valid products at {temperature_kelvin} K")

    # Initialize matrices
    n_products = len(valid_products)
    initial_guess = np.ones((n_products, 1), dtype=np.float64)
    is_condensed = np.zeros(n_products, dtype=bool)

    # 9 coefficients per NASA polynomial
    coeff_matrix = np.zeros((n_products, 9), dtype=np.float64)
    for idx, product in enumerate(valid_products):
        if len(product.coefficients) != 9:
            raise ValueError(f"Invalid coefficients for {product.formula}")

        coeff_matrix[idx] = product.coefficients
        is_condensed[idx] = product.is_condensed

    return initial_guess, coeff_matrix, is_condensed

def compute_total_mass(propellant_composition: PropellantComposition) -> float:
    """
    Calculates the total mass of the propellant based on its composition.

    The total mass is computed as the sum of the products of the molar quantities
    of each element and their respective molar masses.

    Args:
        propellant_composition (PropellantComposition): The propellant composition
            containing the molar quantities of each element.

    Returns:
        float: Total mass of the propellant in kilograms (kg).

    Raises:
        KeyError: If any element in the composition is not found in the molar mass database.
    """
    total_mass = 0.0

    for element, moles in propellant_composition.composition.items():
        if element not in ELEMENT_MOLAR_MASSES:
            raise KeyError(f"Element '{element}' not found in molar mass database.")
        total_mass += moles * ELEMENT_MOLAR_MASSES[element]

    return total_mass
