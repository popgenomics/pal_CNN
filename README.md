# pal_CNN
# execution

## Simulation part
For the proof of concept use python3 but for the density dependent part you need to use pypy.

The command lines below are made to launch density-dependent simulations.
### Start a simulation.

Use simulation_density_pypy.py.

Command line with required arguments : 
```
$ pypy simulation_density_pypy.py -e1 {rate E1} -e2 {rate E2} -e3 {rate E3} -d1 {rate D1} -d2 {rate D2} -c1 {rate C1} -pe2 {rate pE2} -pe3 {rate pE3} -n {simulation name} -o {output directory} -k {response to density} -c {all carrying capacity}
```
#### Required arguments formats :

-e1 : float, rate of extinction E1.

-e2 : float, rate of extinction E2.

-e3 : float, rate of extinction E3.

-d1 : float, rates of speciation D1.

-d2 : float, rates of speciation D2.

-c1 : float, rates of colonisation C1.

-pe2 : float, proportion of genus mass extincted with E2 events.

-pe3 : float, proportion of genus mass extincted with E3 events.

-n : str, the name of the simulation, use to create the output file.

-o : str, path of the directory where the files will be saved.

-k : int, values of k, the response to the density.

-c : list of int, values of carrying capacity. Must be the same length as the number of localities. 

#### Optionals arguments are :  

-cs : int between 0 (default), 1 and 2, number to activate the density dependent speciation:
- 0: don't use density dependence.
- 1: apply density dependence on D2.
- 2: apply density dependence on D1 and D2.

-ce : int between 0 (default), 1 and 2, number to activate the density dependent extinction:
- 0: don\'t use density dependence.
- 1: apply density dependence on E1.
- 2: apply density dependence on E2.

-nl : int, number of localities. 10 by default.

-ng : int, number of generations. 600 by default.

-nmg : int, max number of living genus. 2,000 by default.

-ngl : int, max number if lineages in genus. 50 by default.

-v : Display the number of generations and the number of extincted lineages.

### Product images
1. Processes the data to produce a matrix in csv format.

Use data_to_csv_launcher.py.
```
$ python3 data_to_csv_launcher.py -f {input file} -t {max number of generations} -s -o {output file} -d {number of demes}
-sp {separator}
```
Arguments :

-f : str, input file (produce by the simulator).

-t : int, number of generation in simulation.

-s : str, simulation name.

-o : str, output file that will be create.

-d : int, number of localities in simulation

-g : Use to delete a proportion (with the argument -p) of data in simulation. 

-p : float, between 0 and 1, proportion of data at deleted.

-pbdb : Used to create a matrix from data in the Paleobiology Database.

-sp : str, separator use to delimit column. By default use ','.

-si : int 1 or 10, size of holes added in the data, 1 Ma or 10 Ma.

2. Retrieve the maximum values for each statistic in the matrix.

Use range_launcher.py.
```
$ python3 range_launcher.py -s {path of directory with all the simulation} -o {path of output file}
```
Arguments :

-s : str, path of the directory with all the simulation.

-o : str, output file

3. Produce the images.

Use norm_launcher.py.

```
$ python3 norm_launcher.py -c {input csv path} -cn {output normalize csv path} -r {range file path} -g {greyscale picture path} -rgb {rgb picture path}
```
Arguments :

-c : str, path to input csv file with the matrix of current simulation.

-cn : str, path to output normalize csv file with the normalize matrix.

-r : str, path to the file with max values of each statistics

-g : str, path to pixels representation in greyscale. If not specified the images not will be produce.

-rgb : str path to pixels representation in rgba. If not specified the images not will be produce.

4. create dataset

Use dataset_pixel.py

```
$ python3 dataset_pixel.py -o {output dataset} -i {directory input} -c {number of channel(s)}
```
Arguments : 

-o : str, path and name of the output, specify a h5py format extension.

-i : str, path of the directory with all the pictures.

-c : int, 1, 3 or 4, number of channels use in the pictures :
- 1 : greyscale
- 3 : RGB
- 4 : RGBA

### Old plots
Use plot_fossils.R

```
$ Rscript plots_fossils.R {simulation name} {images directory} {output file} {size of image}
```
Arguments :

- The name of the simulation
- The path of the images directory
- The name of output file
- The image size in x and y axis (in cm).
