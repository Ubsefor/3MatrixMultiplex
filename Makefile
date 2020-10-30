

all:
	gcc -fopenmp -Wall -Wpedantic -Ofast 3MatrixMultiplex/*.c -o 3MatrixMultiplex

no-opt:
	gcc -fopenmp -Wall -Wpedantic -O0 3MatrixMultiplex/*.c -o 3MatrixMultiplex
	
