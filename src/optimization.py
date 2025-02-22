"""Module for optimizing combustion temperature and product composition.

This module provides classes for optimizing the temperature and composition of combustion products
to achieve thermodynamic equilibrium. The optimization process minimizes Gibbs free energy while
ensuring mass balance constraints are satisfied. The thermodynamic properties are based on:

"Thermodynamic and Thermophysical Properties of Combustion Products. Volume 1"
by V.P. Glushko.

Key functionalities include:
- Optimizing the combustion temperature using root-finding methods.
- Optimizing the composition of combustion products for a given temperature by minimizing Gibbs free energy.
- Ensuring mass balance constraints between the propellant and combustion products.

Classes:
    CombustionCompositionOptimizer: Optimizes the composition of combustion products for a given temperature.
    TemperatureOptimizer: Optimizes the combustion temperature by minimizing the difference between
        the enthalpy of the combustion products and the initial propellant.
"""

import numpy as np

from scipy.optimize import minimize, LinearConstraint, root_scalar
from typing import List

from models import PropellantComposition, ReactionProduct
from calculators import ThermodynamicSystemCalculator, ThermodynamicSystemContext
from utils import (
    create_elemental_vector,
    filter_and_construct_matrices,
    prepare_optimization_matrices
)

class CombustionCompositionOptimizer:
    """Optimizes the composition of combustion products for a given temperature.

    This class handles the optimization of substance amounts (moles) for combustion products at a
    specified temperature. It minimizes Gibbs free energy while ensuring mass balance constraints
    between the propellant and combustion products.

    Attributes:
        pressure (float): System pressure in Pascals (Pa).
        temperature (float): Current temperature in Kelvin.
        propellant_vector (np.ndarray): Column vector representing the elemental composition of the propellant.
        stoichiometric_matrix (np.ndarray): Matrix of stoichiometric coefficients for the combustion products.
        initial_guess (np.ndarray): Initial guess for substance amounts (moles).
        coefficients (np.ndarray): Thermodynamic coefficients for the combustion products.
        is_condensed (np.ndarray): Boolean array indicating whether each product is in the condensed phase.
    """

    def __init__(self,
                 pressure: float,
                 temperature: float,
                 propellant_vector: np.ndarray,
                 stoichiometric_matrix: np.ndarray,
                 initial_guess: np.ndarray,
                 coefficients: np.ndarray,
                 is_condensed: np.ndarray):
        self.pressure = pressure
        self.temperature = temperature
        self.propellant_vector = propellant_vector
        self.stoichiometric_matrix = stoichiometric_matrix
        self.initial_guess = initial_guess
        self.coefficients = coefficients
        self.is_condensed = is_condensed

    def optimize(self) -> np.ndarray:
        """Optimize the composition of combustion products by minimizing Gibbs free energy.

        Returns:
            np.ndarray: Optimized substance amounts (moles) for the combustion products.
        """
        # Set up optimization context
        context = ThermodynamicSystemContext(
            substance_amounts=self.initial_guess,
            coefficients=self.coefficients,
            is_condensed=self.is_condensed,
            temperature=self.temperature,
            pressure=self.pressure
        )

        def calculate_gibbs_energy(x):
            context.substance_amounts = x
            return ThermodynamicSystemCalculator.calculate_gibbs_energy(context)

        # Define linear constraint for mass balance
        linear_constraint = LinearConstraint(
            self.stoichiometric_matrix,
            self.propellant_vector.flatten(),
            self.propellant_vector.flatten()
        )

        # Define bounds for variables (non-negative)
        bounds = [(0, None)] * len(self.initial_guess)

        # Run optimization
        result = minimize(
            fun=calculate_gibbs_energy,
            x0=self.initial_guess.flatten(),
            method='trust-constr',
            constraints=linear_constraint,
            bounds=bounds,
            options={
                'verbose': 0,
                'maxiter': 5000
            }
        )

        if result.success:
            return result.x.reshape(1, -1).flatten()
        else:
            raise RuntimeError("Optimization failed to converge")

class TemperatureOptimizer:
    """Optimizes the combustion temperature to achieve thermodynamic equilibrium.

    This class optimizes the combustion temperature by minimizing the difference between the enthalpy
    of the combustion products and the initial propellant. It uses a root-finding method (`brentq`) to
    find the optimal temperature within a specified range.

    Attributes:
        pressure (float): System pressure in Pascals (Pa).
        min_temperature (float): Minimum allowable temperature in Kelvin.
        max_temperature (float): Maximum allowable temperature in Kelvin.
        propellant (PropellantComposition): Propellant composition data.
        products (List[ReactionProduct]): List of candidate reaction products.
        propellant_vector (np.ndarray): Column vector representing the elemental composition of the propellant.
        cached_initial_guess (Optional[np.ndarray]): Cached initial guess for substance amounts (moles).
    """

    def __init__(self,
                 pressure: float,
                 min_temperature: float,
                 max_temperature: float,
                 propellant: PropellantComposition,
                 products: List[ReactionProduct]):
        self.pressure = pressure
        self.min_temperature = min_temperature
        self.max_temperature = max_temperature
        self.propellant = propellant
        self.products = products
        self.propellant_vector = create_elemental_vector(propellant)
        self.cached_initial_guess = None

    def optimize(self) -> float:
        """Optimize the combustion temperature using root-finding.

        Returns:
            float: Optimal combustion temperature in Kelvin.
        """
        left_bound = self.min_temperature
        right_bound = self.max_temperature

        # Use root_scalar with the 'brentq' method
        result = root_scalar(
            f=self._calculate_error,
            bracket=[left_bound, right_bound],
            method='brentq'
        )

        if result.converged:
            return result.root
        else:
            raise RuntimeError("Optimization failed to converge")

    def optimize_context_at_temperature(self, temperature: float) -> ThermodynamicSystemContext:
        """Optimize the composition of combustion products at a given temperature and return the context.

        Args:
            temperature (float): The temperature in Kelvin for which to optimize the context.

        Returns:
            ThermodynamicSystemContext: Optimized thermodynamic system context for the given temperature.
        """
        # Filter products by temperature and create stoichiometric matrix
        temp_filtered, stoichiometric_matrix = filter_and_construct_matrices(
            self.products, self.propellant, temperature
        )

        print(f"Optimizing context at temperature {temperature} K, found {len(temp_filtered)} temperature-compatible products.")

        # Create optimization matrices
        initial_guess, coefficients, is_condensed = prepare_optimization_matrices(
            temp_filtered, temperature
        )

        # Check if initial guess is cached
        if self.cached_initial_guess is not None and len(self.cached_initial_guess) == len(initial_guess):
            initial_guess = self.cached_initial_guess

        # Optimize combustion product composition
        composition_optimizer = CombustionCompositionOptimizer(
            pressure=self.pressure,
            temperature=temperature,
            propellant_vector=self.propellant_vector,
            stoichiometric_matrix=stoichiometric_matrix,
            initial_guess=initial_guess,
            coefficients=coefficients,
            is_condensed=is_condensed
        )

        optimized_substance_amounts = composition_optimizer.optimize()
        self.cached_initial_guess = optimized_substance_amounts

        # Create and return the optimized context
        context = ThermodynamicSystemContext(
            substance_amounts=optimized_substance_amounts,
            coefficients=coefficients,
            is_condensed=is_condensed,
            temperature=temperature,
            pressure=self.pressure
        )
        return context

    def _calculate_error(self, temperature: float) -> float:
        """Calculate the error (enthalpy difference) at a given temperature.

        Args:
            temperature (float): Current temperature in Kelvin.

        Returns:
            float: Error (difference between enthalpy of products and propellant) in J.
        """
        # Use the public method to get the optimized context
        context = self.optimize_context_at_temperature(temperature)

        # Calculate enthalpy error
        return (
            ThermodynamicSystemCalculator.calculate_enthalpy(context) - self.propellant.enthalpy
        )
