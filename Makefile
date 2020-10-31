

all:
	gcc -fopenmp -std=gnu11 -Wall -Wpedantic -Ofast 3MatrixMultiplex/*.c -o 3MatrixMultiplex-exe

no-opt:
	gcc -fopenmp -std=gnu11 -Wall -Wpedantic -O0 3MatrixMultiplex/*.c -o 3MatrixMultiplex-exe
	
clean:
	rm *-exe

