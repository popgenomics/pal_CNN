# pal_CNN
# execution (for the moment)

## Simulation part
### Start a simulation.

Use simulation_launcher.py.

Command line with required arguments : 
```
$ python simulation_launcher.py -e1 -e2 -3 -d1 -d2 -d3 -c1 -pe2 -pe3 -n -o -s -k -t
```
Required arguments formats :

-e1 : float, rate of extinction E1.

-e2 : float, rate of extinction E2.

-e3 : float, rate of extinction E3.

-d1 : float, rates of speciation D1.

-d2 : float, rates of speciation D2.

-d3 : float, rates of speciation D3.

-c1 : float, rates of colonisation C1.

-pe2 : float, proportion of genus mass extincted with E2 events.

-pe3 : float, proportion of genus mass extincted with E3 events.

-n : str, the name of the simulation, use to create the output file.

-o : str, path of the directory where the files will be saved.

-s : list of int, value of sigma (the strength of the selection) associate at each trait separate by space. 
Respect the same order than parameter -k and -sdeme.

-k : list of int, value of k (the curvature of the curve) associate at each trait separate by space.
Respect the same order than parameter -s and -sdeme.

-t : int, number of traits.

Optionals arguments are :  

-sdeme : list of float, value of s_deme for each locality and each trait. With 1 trait and 2 localities the command line
will be ```-sdeme 0.1 0.3```. 0.1 is the s_deme for deme 1 and 0.3 s_deme for deme 2.

With 2 traits and 2 localities the command line will be ```-sdeme 0.1 0.5 0.3 0.2```. 0.1 is the s_deme for trait 1 in deme 1,
0.5 is the s_deme for trait 2 in deme 1, 0.3 and s_deme for trait 1 in deme 2 and 0.2 is s_deme for trait 2 in deme 2.
If not specify the s_deme will be random, 0 or 1.

-nl : int, number of localities. 10 by default.

-ng : int, number of generations. 600 by default.

-nmg : int, max number of genus. 2,000 by default.

-ngl : int, max number if lineages in genus. 50 by default.

-v : Display the number of generations and the number of extincted lineages.

### Plots
1. Prepare the data to be plots

Use data_to_plots_launcher.py.
```
$ python3 data_to_plots_launcher.py -f {input file} -t {max number of generations} -o {output file}
```
Arguments :

-f : str, input file (produce by the simulator).

-n : int, number of generation in the simulation.

-o : str, output file that will be create.

2. Setup the scale of y axis of each plots.

Use range_launcher.py.
This scripts work only if the output file of the simulation have the form S{scenario number}_{simulation number}... 
```
$ python3 range_launcher.py -s {path of directory with all the simulation} -n {number of simulation} -sc {the number of this scenario} -o {path of output file}
```
Arguments :

-s : str, path of the directory with all the simulation.

-n : int, number of simulations.

-sc : int, scenario number. Use to find file.

-o : str, output file

3. Plots

Use plot_fossils.R
```
Rscript plots_fossils.R {simulation name} {images directory} {output file} {size of image}
```
Arguments :

- The name of the simulation
- The path of the images directory
- The name of output file
- The image size in x and y axis (in cm).

## Convolutional neural network (CNN) part
In dev