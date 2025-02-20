from typing import Callable
import numpy as np
from thermodynamics import ThermodynamicCalculator

class GibbsEnergyCalculator:
    def __init__(self, pressure: float, temperature: float):
        self.R = 8.314  # J/(mol·K)
        self.P0 = pressure
        self.T = temperature

    def calculate(self, 
                quantities: np.ndarray, 
                coefficients: np.ndarray,
                is_condensed: np.ndarray) -> float:
        """Calculate total Gibbs free energy of the system."""
        total_gibbs = 0.0
        total_gas_moles = np.sum(quantities * ~is_condensed)
        
        for i in range(len(quantities)):
            q = quantities[i]
            coeffs = coefficients[i]
            
            # Calculate thermodynamic properties
            h = ThermodynamicCalculator.calculate_enthalpy(coeffs, self.T)
            s = self._calculate_entropy(coeffs, q, total_gas_moles, is_condensed[0, i])
            
            total_gibbs += q * (h - self.T * s)
            
        return total_gibbs

    def _calculate_entropy(self, 
                          coeffs: np.ndarray,
                          quantity: float,
                          total_gas: float,
                          is_condensed: bool) -> float:
        """Calculate entropy with pressure correction."""
        pressure = self.P0 * (quantity / total_gas) if not is_condensed else 0.0
        return ThermodynamicCalculator.calculate_entropy(
            coeffs, self.T, pressure
        )
