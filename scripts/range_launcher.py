import argparse
from setup_range import setup_range

if __name__ == "__main__":
    parser= argparse.ArgumentParser()
    parser.add_argument('--simulation_dir', '-s', help='Directory with all simulations', type=str, required=True, action='store')
    parser.add_argument('--number_simulation', '-n', help='Number of simulations', type=int, required=True, action='store')
    parser.add_argument('--scenario', '-sc', help='The number of the scenario', type=int, required=True, action='store')
    parser.add_argument('--output_file', '-o', help='path and name of the output file', type=str, required=True, action='store')
    args = parser.parse_args()

    setup_range(args.simulation_dir, args.number_simulation, args.scenario, args.output_file)