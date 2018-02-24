# Ant Colony Optimization 2DHP

### Executable aco-2dhp
###### How to compile
```
$ gcc -Wall -O2  -c ./file.c -o ./obj/Release/file.o
$ gcc -Wall -O2  -c ./main.c -o ./obj/Release/main.o
$ gcc -Wall -O2  -c ./aco.c -o ./obj/Release/aco.o
$ g++  -o ./bin/Release/aco-2dhp ./obj/Release/file.o ./obj/Release/main.o ./obj/Release/aco.o  -s
```
###### How to run
```
$ [sudo] ./bin/Release/aco-2dhp input_file_path sequence_key_on_input_file
```
###### Output
```
iterations|protein|directions|energy|final population average|final population standard deviation|convergence|found on iteration|time|iteration,energy/.../iteration,energy/
```

### Bash script.sh
###### How to run
```
$ [sudo] chmod +x ./script.sh
$ [sudo] ./script.sh
```
###### Output
```
./results/outputlog
./results/summary
./results/figures/sequence[1..8]
./results/tables/sequence[1..8]
./results/energy_evolution/sequence[1..8]/run[1..20]
```
###### Generate pdf/png from latex files:
```
$ pdflatex --shell-escape filename
```
