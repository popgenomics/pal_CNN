import argparse
from data_to_plots import data_to_dataframes

if __name__ == "__main__":
    parser= argparse.ArgumentParser()
    parser.add_argument('--file', '-f', help='Input file', type=str, required=True, action='store')
    parser.add_argument('--times', '-t', help='Max times in the simulation', type=int, required=True, action='store')
    parser.add_argument('--simulations', '-s', help='simulation id', type=int, required=True, action='store')
    parser.add_argument('--output_dir', '-o', help='path and name of the output directory', type=str, required=True, action='store')
    parser.add_argument('--demes', '-d', help='Number of demes in the simulation', type=int, required=True, action='store')
    args = parser.parse_args()

    data_to_dataframes(args.file, args.times, args.simulations, args.output_dir, args.demes)