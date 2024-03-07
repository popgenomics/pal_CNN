import argparse
from data_to_plots import data_to_dataframes

if __name__ == "__main__":
    parser= argparse.ArgumentParser()
    parser.add_argument('--file', '-f', help='Input file', type=str, required=True, action='store')
    parser.add_argument('--range', '-r', help='File with the max value of each layers', type=str, required=True, action='store')
    parser.add_argument('--output_file', '-o', help='path and name of the output file', type=str, required=True, action='store')
    args = parser.parse_args()

    data_to_dataframes(args.file, args.range, args.output_file)