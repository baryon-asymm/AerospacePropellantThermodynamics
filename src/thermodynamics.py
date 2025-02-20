import numpy as np

class ThermodynamicCalculator:
    @staticmethod
    def calculate_enthalpy(coefficients: np.ndarray, temperature: float) -> float:
        """Calculate molar enthalpy using NASA polynomial coefficients.
        
        Args:
            coefficients: Array of 9 polynomial coefficients
            temperature: Current temperature in Kelvin
            
        Returns:
            Enthalpy in J/mol
        """
        t = temperature * 1e-3  # Reduced temperature
        return 4.184 * (
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
    def calculate_entropy(coefficients: np.ndarray, 
                         temperature: float, 
                         partial_pressure: float = 0.0) -> float:
        """Calculate molar entropy with optional pressure correction.
        
        Args:
            coefficients: Array of 9 polynomial coefficients
            temperature: Current temperature in Kelvin
            partial_pressure: Partial pressure in Pa (for gases)
            
        Returns:
            Entropy in J/(mol·K)
        """
        R = 8.314  # Universal gas constant
        std_pressure = 101325  # Standard pressure in Pa
        t = temperature * 1e-3
        
        std_entropy = 4.184 * (
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
            std_entropy -= R * np.log(partial_pressure / std_pressure)
            
        return std_entropy

    @staticmethod
    def calculate_heat_capacity(coeffs: np.ndarray, T: float) -> float:
        """Calculate heat capacity using NASA polynomial coefficients."""
        Tr = T * 1e-3
        return 4.184 * 1e-3 * (coeffs[2] + 2*coeffs[3]*Tr + 3*coeffs[4]*Tr**2 +
                             4*coeffs[5]*Tr**3)
