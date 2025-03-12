import argparse
import os
import sys
import numpy as np

from json_reader import load_combustion_products, load_input_substance
from json_writer import write_to_json
from optimization import TemperatureOptimizer
from thermodynamic_properties import ThermodynamicPropertiesCalculator
from utils import filter_and_construct_matrices, compute_total_mass

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Perform thermodynamic optimization for combustion products.")
    
    # Required arguments
    parser.add_argument(
        "--propellant",
        required=True,
        help="Path to the input propellant JSON file (e.g., propellant.json)."
    )
    parser.add_argument(
        "--combustion-products",
        required=True,
        help="Path to the combustion products JSON file (e.g., combustion_products.json)."
    )
    parser.add_argument(
        "--pressure",
        type=float,
        required=True,
        help="Chamber pressure in Pascals (Pa). Must be a positive value."
    )
    
    # Optional argument
    parser.add_argument(
        "--output-json",
        default=None,
        help="Path to the output JSON file. If not provided, no output JSON will be generated."
    )
    
    # Parse arguments
    args = parser.parse_args()

    try:
        # Validate pressure
        if args.pressure <= 0:
            raise ValueError("Pressure must be a positive value.")

        # Resolve absolute paths for input files
        input_propellant_path = os.path.abspath(args.propellant)
        combustion_products_path = os.path.abspath(args.combustion_products)

        # Load data
        print("Loading data...")
        combustion_products = load_combustion_products(combustion_products_path)
        input_substance = load_input_substance(input_propellant_path)
        print(f"Loaded {len(combustion_products)} combustion products")

        # Define temperature bounds
        min_temperature = min(product.temperature_range.min for product in combustion_products)
        max_temperature = max(product.temperature_range.max for product in combustion_products)
        delta_temperature = 1e-3

        # Optimize temperature
        print("Optimizing temperature...")
        optimizer = TemperatureOptimizer(
            pressure=args.pressure,
            min_temperature=min_temperature + delta_temperature,
            max_temperature=max_temperature - delta_temperature,
            propellant=input_substance,
            products=combustion_products
        )
        optimal_temperature = optimizer.optimize()
        print("Optimal Temperature:", optimal_temperature)

        # Get optimized context and filtered products
        context = optimizer.optimize_context_at_temperature(optimal_temperature)
        filtered_products, _ = filter_and_construct_matrices(
            optimizer.products, optimizer.propellant, optimal_temperature
        )

        # Calculate additional thermodynamic properties
        print("Calculating thermodynamic properties...")
        thermo_calculator = ThermodynamicPropertiesCalculator(
            context=context,
            filtered_products=filtered_products
        )
        thermo_properties_context = thermo_calculator.calculate_and_display_properties()

        # Write results to JSON if output path is provided
        if args.output_json:
            print("Writing results to JSON...")
            total_mass = compute_total_mass(input_substance)
            output_json_path = os.path.abspath(args.output_json)  # Resolve absolute path for output
            
            # Ensure the output directory exists
            output_directory = os.path.dirname(output_json_path)
            os.makedirs(output_directory, exist_ok=True)  # Create directory if it doesn't exist

            write_to_json(output_json_path, input_substance, total_mass, thermo_properties_context, filtered_products)
            print(f"Results successfully written to {output_json_path}")
        else:
            print("Output JSON file path not provided. Skipping JSON generation.")

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
