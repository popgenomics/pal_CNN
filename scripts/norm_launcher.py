import argparse
from normalize import normalize_array
from normalize import to_image

if __name__ == "__main__":
    parser= argparse.ArgumentParser()
    parser.add_argument('--csv', '-c', help='Path to input csv file', type=str, required=True, action='store')
    parser.add_argument('--csv_norm', '-cn', help='Path to output normalize csv file', type=str, required=True, action='store')
    parser.add_argument('--range', '-r', help='Path to the file with max values of each stats', type=str, required=True, action='store')
    #parser.add_argument('--demes', '-d', help='Number of localities in simulation', type=int, required=True, action='store')
    parser.add_argument('--greyout', '-g', help='Path to pixels representation in greyscale.\nIf not specified the images not will be produce.', type=str, required=False, action='store')
    parser.add_argument('--rgbout', '-rgb', help='Path to pixels representation in rgba.\nIf not specified the images not will be produce.', type=str, required=False, action='store')
    args = parser.parse_args()

    normalize_array(csv_file=args.csv, range_file=args.range, outfile=args.csv_norm)
    to_image(array_norm=args.csv_norm, greyout=args.greyout, rgbout=args.rgbout)

