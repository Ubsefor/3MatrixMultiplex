

all:
	gcc -fopenmp -Wall -Wpedantic -Ofast 3MatrixMultiplex/*.c -o 3MatrixMultiplex-exe
