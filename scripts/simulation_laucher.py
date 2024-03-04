import argparse
import os
from simulation import simulation

if __name__ == "__main__":
    current_dir = os.getcwd() # Pour moi à changer après.
    parser= argparse.ArgumentParser()
    parser.add_argument('--rate', '-r', help='Rates of extinction, speciation, colonisation', type=float, nargs='+', required=True, action='store')
    parser.add_argument('--proportion', '-p', help='Proportion of demes mass extincted', type=float, nargs='+', required=True, action='store')
    parser.add_argument('--name', '-n', help='Name of simulation', type=str, required=True, action='store')
    parser.add_argument('--output_dir', '-o', help='Name of simulaiton', type=str, default=f'{current_dir}/results')
    args = parser.parse_args()

    simulation(args.rate, args.proportion, args.name, args.output_dir)
