

all: mini

mini:
	gcc -fopenmp -DMINI_DATASET -std=gnu11 -Wall -Wpedantic -Ofast 3MatrixMultiplex/*.c -o 3MatrixMultiplex-exe

small:
	gcc -fopenmp -DSMALL_DATASET -std=gnu11 -Wall -Wpedantic -Ofast 3MatrixMultiplex/*.c -o 3MatrixMultiplex-exe

medium:
	gcc -fopenmp -DMEDIUM_DATASET -std=gnu11 -Wall -Wpedantic -Ofast 3MatrixMultiplex/*.c -o 3MatrixMultiplex-exe
	
large:
	gcc -fopenmp -DLARGE_DATASET -std=gnu11 -Wall -Wpedantic -Ofast 3MatrixMultiplex/*.c -o 3MatrixMultiplex-exe

extra:
	gcc -fopenmp -DEXTRALARGE_DATASET -std=gnu11 -Wall -Wpedantic -Ofast 3MatrixMultiplex/*.c -o 3MatrixMultiplex-exe

no-opt: noopt-mini

noopt-mini:
	gcc -fopenmp -DMINI_DATASET -std=gnu11 -Wall -Wpedantic -O0 3MatrixMultiplex/*.c -o 3MatrixMultiplex-exe
	
noopt-small:
	gcc -fopenmp -DSMALL_DATASET -std=gnu11 -Wall -Wpedantic -O0 3MatrixMultiplex/*.c -o 3MatrixMultiplex-exe
	
noopt-medium:
	gcc -fopenmp -DMEDIUM_DATASET -std=gnu11 -Wall -Wpedantic -O0 3MatrixMultiplex/*.c -o 3MatrixMultiplex-exe
	
clean:
	rm *-exe

