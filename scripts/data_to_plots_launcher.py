import argparse
from data_to_plots import data_to_dataframes

if __name__ == "__main__":
    parser= argparse.ArgumentParser()
    parser.add_argument('--file', '-f', help='Input file', type=str, required=True, action='store')
    parser.add_argument('--times', '-t', help='Max times in the simulation', type=int, required=True, action='store')
    parser.add_argument('--output_file', '-o', help='path and name of the output file', type=str, required=True, action='store')
    args = parser.parse_args()

    data_to_dataframes(args.file, args.times, args.output_file)