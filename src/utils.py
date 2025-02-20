import numpy as np
from typing import List, Dict, Tuple
from models import ReactionProduct, PropellantComposition

def create_matrix_vector(propellant_composition: PropellantComposition) -> np.ndarray:
    """Creates a column vector representation of element quantities from propellant composition.
    
    The row order corresponds to element symbols in propellant_composition.composition.keys() order.
    
    Args:
        propellant_composition: Contains chemical composition data
        
    Returns:
        2D numpy array: Column vector with element quantities
    """
    values = list(propellant_composition.composition.values())
    return np.array(values).reshape(-1, 1)

def extract_elements_from_formula(formula: str) -> Dict[str, int]:
    """Parses chemical formula into constituent elements with stoichiometric counts.
    
    Correctly handles:
    - Multi-character element symbols (e.g., Pb, Fe, Na)
    - Implicit counts (e.g., O = O1)
    - Complex formulas (e.g., Fe3O4, C6H5OH)

    Args:
        formula: Chemical formula following standard notation

    Returns:
        Dictionary mapping element symbols to their integer counts

    Raises:
        ValueError: For invalid characters or formula structure
    """
    elements = {}
    i = 0
    n = len(formula)
    
    while i < n:
        # Element symbol must start with uppercase
        if not formula[i].isupper():
            raise ValueError(f"Invalid character '{formula[i]}' at position {i} - elements must start with uppercase")
        
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

def filter_reaction_products(reaction_products: List[ReactionProduct],
                             propellant_composition: PropellantComposition) -> List[ReactionProduct]:
    """Filters products to those containing ONLY elements present in propellant composition.
    
    Args:
        reaction_products: Candidate reaction products
        propellant_composition: Reference composition for filtering
        
    Returns:
        Products whose ALL constituent elements exist in propellant composition
    
    Examples:
        Propellant contains ['C', 'H', 'O']
        Kept: H2O (H,O), CO2 (C,O)
        Removed: NaCl (Na,Cl), CH3Cl (C,H,Cl)
    """
    propellant_elements = set(propellant_composition.composition.keys())
    
    return [
        product for product in reaction_products
        if propellant_elements.issuperset(extract_elements_from_formula(product.formula))
    ]

def create_stoichiometric_matrix(filtered_products: List[ReactionProduct],
                                 propellant_composition: PropellantComposition) -> np.ndarray:
    """Constructs matrix of stoichiometric coefficients aligned with propellant elements.
    
    Matrix dimensions: (num_propellant_elements × num_filtered_products)
    Row order matches propellant_composition.composition.keys() order
    
    Args:
        filtered_products: Reaction products containing propellant elements
        propellant_composition: Defines row order and relevant elements
        
    Returns:
        2D numpy array: Stoichiometric coefficients matrix
    """
    elements = list(propellant_composition.composition.keys())
    matrix = np.zeros((len(elements), len(filtered_products)), dtype=np.float64)
    
    for col_idx, product in enumerate(filtered_products):
        counts = extract_elements_from_formula(product.formula)
        for row_idx, element in enumerate(elements):
            matrix[row_idx, col_idx] = counts.get(element, 0.0)
    
    return matrix

def create_temperature_filtered_matrices(
    reaction_products: List[ReactionProduct],
    propellant: PropellantComposition,
    temperature: float
) -> Tuple[List[ReactionProduct], np.ndarray]:
    """Filter products by temperature and create stoichiometric matrix.
    
    Args:
        reaction_products: List of all reaction products
        propellant: Propellant composition reference
        temperature: Target temperature in Kelvin
        
    Returns:
        Tuple of (filtered products, stoichiometric matrix)
    """
    # First filter by element compatibility
    element_filtered = filter_reaction_products(reaction_products, propellant)
    
    # Then filter by temperature validity
    temp_filtered = [
        p for p in element_filtered
        if p.temperature_range.min <= temperature <= p.temperature_range.max
    ]
    
    # Create stoichiometric matrix
    matrix = create_stoichiometric_matrix(temp_filtered, propellant)
    
    return temp_filtered, matrix

def create_optimization_matrices(
    reaction_products: List[ReactionProduct],
    temperature_kelvin: float
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Creates optimization matrices with condensed phase flags.
    
    Args:
        reaction_products: List of ReactionProduct objects
        temperature_kelvin: Target temperature in Kelvin
        
    Returns:
        Tuple containing:
        - initial_guess: Initial substance quantities [mol] (1 × N)
        - coeff_matrix: Thermodynamic coefficients (N × 9)
        - is_condensed: Boolean flags for condensed phases (1 × N)
        
    Raises:
        ValueError: For invalid temperature or missing data
    """
    if temperature_kelvin <= 0:
        raise ValueError(f"Invalid temperature: {temperature_kelvin} K")

    # Filter valid products
    valid_products = [
        p for p in reaction_products
        if p.temperature_range.min <= temperature_kelvin <= p.temperature_range.max
    ]
    
    if not valid_products:
        raise ValueError(f"No valid products at {temperature_kelvin} K")

    # Initialize matrices
    n_products = len(valid_products)
    initial_guess = np.ones((n_products, 1), dtype=np.float64)
    is_condensed = np.zeros((1, n_products), dtype=bool)
    
    # 9 coefficients per NASA polynomial
    coeff_matrix = np.zeros((n_products, 9), dtype=np.float64)  

    for idx, product in enumerate(valid_products):
        if len(product.coefficients) != 9:
            raise ValueError(f"Invalid coefficients for {product.formula}")
            
        coeff_matrix[idx] = product.coefficients
        is_condensed[0, idx] = product.is_condensed

    return initial_guess, coeff_matrix, is_condensed
