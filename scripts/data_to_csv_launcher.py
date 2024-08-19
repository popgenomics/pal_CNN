import argparse
from data_to_csv import add_holes, data_to_dataframes, pbdb_to_dataframes

if __name__ == "__main__":
    parser= argparse.ArgumentParser()
    parser.add_argument('--file', '-f', help='Input file', type=str, required=True, action='store')
    parser.add_argument('--times', '-t', help='Max times in the simulation', type=int, required=True, action='store')
    parser.add_argument('--simulations', '-s', help='simulation id', type=str, required=True, action='store')
    parser.add_argument('--output_dir', '-o', help='path and name of the output directory', type=str, required=True, action='store')
    parser.add_argument('--demes', '-d', help='Number of demes in the simulation', type=int, required=True, action='store')
    parser.add_argument('--proportion_holes', '-p', help='Proportions of holes in the simulation', type=float, action='store')
    parser.add_argument('--global_holes', '-g', help='The holes are add in global not in localities', action='store_true')
    parser.add_argument('--add_holes', '-ah', help='add holes', action='store_true')
    parser.add_argument('--pbdb', '-pbdb', help='Use pbdb', action='store_true')
    parser.add_argument('--sep', '-sp', action='store', default=',')
    parser.add_argument('--size_holes', '-si', action='store', choices=[1,10], type=int)
    args = parser.parse_args()

    if args.add_holes:
        add_holes(data=args.file, nb_time=args.times, simul=args.simulations, fileout=args.output_dir, nb_deme=args.demes, prop_holes=args.proportion_holes, global_holes=args.global_holes, size_holes=args.size_holes)
    elif args.pbdb:
        if args.sep == 't':
            current_sep = '\t'
        else:
            current_sep = args.sep
        pbdb_to_dataframes(data=args.file, nb_time=args.times, simul=args.simulations, fileout=args.output_dir, nb_deme=args.demes, sep=current_sep)
    else:
        data_to_dataframes(data=args.file, nb_time=args.times, simul=args.simulations, fileout=args.output_dir, nb_deme=args.demes)