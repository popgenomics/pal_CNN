import argparse
from simulation import simulation

if __name__ == "__main__":
    parser= argparse.ArgumentParser()
    parser.add_argument('--proba_E1', '-e1', help='Rates of extinction E1', type=float, required=True, action='store')
    parser.add_argument('--proba_E2', '-e2', help='Rates of extinction E2', type=float, required=True, action='store')
    parser.add_argument('--proba_E3', '-e3', help='Rates of extinction E3', type=float, required=True, action='store')
    parser.add_argument('--proba_D1', '-d1', help='Rates of speciation D1', type=float, required=True, action='store')
    parser.add_argument('--proba_D2', '-d2', help='Rates of speciation D2', type=float, required=True, action='store')
    parser.add_argument('--proba_C1', '-c1', help='Rates of colonisation C1', type=float, required=True, action='store')
    parser.add_argument('--proportion_E2', '-pe2', help='Proportion of genus mass extincted', type=float,required=True, action='store')
    parser.add_argument('--proportion_E3', '-pe3', help='Proportion of demes mass extincted', type=float,required=True, action='store')
    parser.add_argument('--name', '-n', help='Name of simulation', type=str, required=True, action='store')
    parser.add_argument('--output_dir', '-o', help='Output directory', type=str, required=True, action='store')
    parser.add_argument('--nb_localities', '-nl', help='Number of localities', type=int, default=10, action='store')
    parser.add_argument('--nb_generations', '-ng', help='Number of generations', type=int, default=600, action='store')
    parser.add_argument('--nb_max_genus', '-nmg', help='Max number of genus', type=int, default=2000, action='store')
    parser.add_argument('--nb_max_lineages', '-ngl', help='Max number of lineages in genus', type=int, default=50, action='store')
    parser.add_argument('--verbose', '-v', help='Display the generations', action='store_true')
    args = parser.parse_args()

    rates = [args.proba_E1, args.proba_E2, args.proba_E3, args.proba_D1, args.proba_D2, args.proba_C1]
    proportions = [args.proportion_E2, args.proportion_E3]
    simulation(rate=rates, proportion=proportions, simulation_name=args.name, result_repertory=args.output_dir, verbose=args.verbose, nb_localities=args.nb_localities, nb_generations=args.nb_generations, nb_max_genus=args.nb_max_genus, nb_max_lineages_in_genus=args.nb_max_lineages)