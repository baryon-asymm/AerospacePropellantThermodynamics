import json
from models import ReactionProduct, TemperatureRange, PropellantComposition
from json_reader import load_combustion_products, load_input_substance
from utils import (
    create_matrix_vector,
    create_temperature_filtered_matrices,
    create_optimization_matrices
)
from scipy.optimize import minimize
import numpy as np
from gibbs import GibbsEnergyCalculator

# Load data
combustion_products = load_combustion_products(r"C:\Projects\RocketPropellantTEM\src\combustion_products.json")
input_substance = load_input_substance(r"C:\Projects\RocketPropellantTEM\src\input_substance.json")

# Create matrix-vector for input substance
matrix_vector = create_matrix_vector(input_substance)
print("Input Substance Vector:\n", matrix_vector)

# Filter products and create matrices
temperature = 3594.0
filtered_products, stoichiometric_matrix = create_temperature_filtered_matrices(
    combustion_products,
    input_substance,
    temperature
)
print(f"Found {len(filtered_products)} temperature-compatible products")

# Create optimization matrices
initial_guess, coeffs, is_condensed = create_optimization_matrices(filtered_products, temperature)

# Initialize calculator with correct parameters
calculator = GibbsEnergyCalculator(
    pressure=100*101325,  # Standard pressure in Pa
    temperature=temperature
)

# Set up optimization constraints
constraints = {
    'type': 'eq',
    'fun': lambda x: (stoichiometric_matrix @ x.reshape(-1, 1)).flatten() - matrix_vector.flatten()
}

# Run optimization
result = minimize(
    fun=lambda x: calculator.calculate(x, coeffs, is_condensed),  # Reshape for matrix ops
    x0=initial_guess.flatten(),  # Flatten initial guess
    method='trust-constr',
    constraints=constraints,
    bounds=[(0, None)] * len(filtered_products),
    options={
        'disp': True,
        'maxiter': 10000,
    }
)

# Process results
if result.success:
    optimized_quantities = result.x.reshape(1, -1)
    print("Optimization successful!")
    print("Final quantities:", optimized_quantities)
else:
    print("Optimization failed:", result.message)

# Calculate final Gibbs energy
final_gibbs = calculator.calculate(optimized_quantities.flatten(), coeffs, is_condensed)
print("Final Gibbs energy:", final_gibbs)

from thermodynamics import ThermodynamicCalculator

# Calculate final enthalpy
total_enthalpy = 0.0
for i in range(len(result.x)):
    q = optimized_quantities[0, i]
    cs = coeffs[i]
    # Calculate thermodynamic properties
    h = ThermodynamicCalculator.calculate_enthalpy(cs, temperature)
    total_enthalpy += q * h
print("Total enthalpy:", np.sum(total_enthalpy))
