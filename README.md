# Thermodynamic Combustion Optimization

This project performs thermodynamic optimization for combustion systems to determine the equilibrium composition and properties of combustion products. It uses thermodynamic coefficients from *Thermodynamic and Thermophysical Properties of Combustion Products. Volume 1* by V.P. Glushko to calculate properties such as enthalpy, entropy, heat capacity, and Gibbs free energy.

## Features

- **Thermodynamic Property Calculations:**  
  Computes enthalpy, entropy, heat capacity, and Gibbs free energy for individual substances and entire systems.

- **Combustion Optimization:**  
  Optimizes the combustion temperature and product composition by minimizing Gibbs free energy while enforcing mass balance constraints.

- **JSON Input/Output:**  
  Reads propellant and combustion product data from JSON files and writes the optimization results to a JSON file.

- **Mass Balance:**  
  Ensures mass conservation between the propellant and combustion products.

- **Phase Support:**  
  Handles both gas-phase and condensed-phase species.

- **Derived Properties:**  
  Calculates additional properties such as the specific gas constant, heat ratio, and volumetric heat capacity.

## Installation

1. **Clone the Repository:**

   ```bash
   git clone <repository_url>
   cd thermodynamic-combustion-optimization
   ```

2. **(Optional) Create a Virtual Environment:**

   ```bash
   python3 -m venv venv
   source venv/bin/activate  # On Linux/macOS
   # venv\Scripts\activate    # On Windows
   ```

3. **Install Dependencies:**

   ```bash
   pip install numpy scipy
   ```

## Usage

### 1. Prepare Input JSON Files

Sample JSON files are provided in the `data` directory:

- **Propellant (`propellant.json`):**  
  Contains data about the propellant, including its enthalpy and elemental composition.  
  **Example:**

  ```json
  {
      "enthalpy": -8.23e+6,
      "composition": {
          "C": 45.67,
          "H": 67.98,
          "O": 22.85,
          "N": 13.25,
          "Cl": 2.12,
          "Al": 7.41
      }
  }
  ```

- **Combustion Products (`combustion_products.json`):**  
  Contains potential combustion products along with their thermodynamic coefficients, phase information, and temperature ranges.  
  **Example:**

  ```json
  [
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
      },
      {
          "formula": "H2O",
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
  ]
  ```

### 2. Run the Main Script

All Python files are located in the `src` directory. Run the optimization by executing the following command:

```bash
python src/main.py --propellant data/propellant.json --combustion-products data/combustion_products.json --pressure <chamber_pressure_pa> --output-json <path_to_output.json>
```

- **`data/propellant.json`:** Path to the propellant JSON file.
- **`data/combustion_products.json`:** Path to the combustion products JSON file.
- **`<chamber_pressure_pa>`:** Chamber pressure in Pascals (Pa). Must be a positive value.
- **`<path_to_output.json>` (optional):** Path to the output JSON file. If omitted, results are printed to the console only.

**Example:**

```bash
python src/main.py --propellant data/propellant.json --combustion-products data/combustion_products.json --pressure 1013250 --output-json output/results.json
```

## Output

The script outputs the optimal combustion temperature and various thermodynamic properties to the console. If an output JSON file path is provided, the results are saved in JSON format, including:

- Propellant data
- Combustion products (formula, phase, and moles)
- Optimal combustion temperature

**Example JSON Output Structure:**

```json
{
    "pressure": 1013250.0,
    "temperature": 3000.0,
    "propellant": {
        "enthalpy": -8230000.0,
        "composition": {
            "C": 45.67,
            "H": 67.98,
            "O": 22.85,
            "N": 13.25,
            "Cl": 2.12,
            "Al": 7.41
        },
        "total_mass_kg": 1.0
    },
    "combustion_products": [
        {
            "formula": "CO2",
            "phase": "gas",
            "moles": 1.234
        },
        {
            "formula": "H2O",
            "phase": "gas",
            "moles": 2.345
        }
    ]
}
```

## Project Structure

```
AerospacePropellantThermodynamics/
├── data/                   # Sample JSON files: propellant.json and combustion_products.json
├── output/                 # Directory for output JSON files (created automatically if it doesn't exist)
└── src/                    # All Python source files
    ├── calculators.py      # Classes for calculating thermodynamic properties.
    ├── constants.py        # Defines physical and thermodynamic constants.
    ├── json_reader.py      # Loads data from JSON files.
    ├── json_writer.py      # Writes results to a JSON file.
    ├── main.py             # Main script for running the optimization.
    ├── models.py           # Data classes representing the propellant and combustion products.
    ├── molar_masses.py     # Provides molar masses of chemical elements.
    ├── optimization.py     # Optimization of combustion temperature and product composition.
    ├── thermodynamic_properties.py  # Calculates and displays system thermodynamic properties.
    └── utils.py            # Utility functions for processing chemical formulas and constructing matrices.
```

## Error Handling

The script handles several types of errors:

- **Invalid Pressure:**  
  Checks that the specified chamber pressure is a positive number.

- **File Issues:**  
  Detects missing or malformed JSON files and required field omissions.

- **Optimization Errors:**  
  Catches optimization convergence failures and invalid temperature values.

In the event of an error, an appropriate error message and help information are displayed.

## Contributing

Contributions are welcome! Please submit a pull request with your improvements or bug fixes.

## License

This project is licensed under the [MIT License](LICENSE).

## Acknowledgements

This project uses thermodynamic coefficients from:

*V.P. Glushko, "Thermodynamic and Thermophysical Properties of Combustion Products. Volume 1"*
