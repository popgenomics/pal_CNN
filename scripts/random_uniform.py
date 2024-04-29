from numpy.random import uniform
import argparse

def random_parameter(type, simulation_number):
    if type == 'E1':
        parameter = uniform(low=0.01, high=0.1, size=1)
    elif type == 'E2':
        parameter = uniform(low=0.01, high=0.2, size=1)
    elif type == 'E3':
        parameter = uniform(low=0.01, high=0.2, size=1)
    elif type == 'pE2':
        parameter = uniform(low=0.2, high=1.0, size=1)
    elif type == 'pE3':
        parameter = uniform(low=0.2, high=1.0, size=1)
    else:
        raise ValueError
    return parameter[0]

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--type', '-t', help='The parameter at define', type=str, required=True, action='store')
    parser.add_argument('--simul_nb', '-s', help='The number of the simulation. Just use to be sure to have new random value', type=str, required=True, action='store')
    args = parser.parse_args()
    print(random_parameter(type=args.type, simulation_number=args.simul_nb))