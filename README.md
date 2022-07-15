# Ant Colony Optimization 2DHP

Implementation of the Ant Colony Optimization meta-heuristic for the Protein Structure Prediction problem using 2D HP model. This approach includes the pull move heuristic for local search.

## Requirements
* gcc
* g++

## Compilation

To compile, execute each one of the following commands from the root of the project:
```
gcc -Wall -O2  -c ./file.c -o ./file.o
```
```
gcc -Wall -O2  -c ./main.c -o ./main.o
```
```
gcc -Wall -O2  -c ./aco.c -o ./aco.o
```
```
g++  -o ./aco-2dhp ./file.o ./main.o ./aco.o  -s
```

## How to run

After compiling it, you can run the aco-2dhp executable by executing the following command from the root of the project.
```
./aco-2dhp input sequence8
```

You can open the `input` file, located at the root of the project, in a text editor to change the parameters and see other sequences available.


The output will be printed following the template bellow
```
iterations|protein|directions|energy|final population average|final population standard deviation|convergence|found on iteration|time|iteration,energy/.../iteration,energy/
```

## Running the script.sh
The `script.sh` script can be used to run the aco-2dhp executable for all the proteins from the `input`. The script will generate latex tables for all the results obtained. The script will also plot figures representing the best resultant conformation of each protein, as long the `2dhp-plot` executable is avaiable in the root of the project (The `2dhp-plot` executable avaiable in the root of the project was compiled in a linux x64. To compile a new `2dhp-plot` executable, check out this other repository: https://github.com/douglas444/2dhp-plot )

To run the `script.sh` script, execute the following commands from the root of the project
```
chmod +x ./script.sh
```
```
./script.sh
```
The results are going to be generated as follows:
```
./results/outputlog
./results/summary
./results/figures/sequence[1..8]
./results/tables/sequence[1..8]
./results/energy_evolution/sequence[1..8]/run[1..20]
```
To generate pdf/png from the latex files, execute the following command, changing `filename` to the path to the latex file you want to convert to pdf/png
```
pdflatex --shell-escape filename
```
