
help:
	echo "Available makes: [noopt-]user, [noopt-]<dataset size>, clean, distclean"
	echo "But bevare, dataset builds are less friendly"
	echo "Clean option removes executable"
	echo "Distclean also removes benchmarks"
	
user:
	gcc -fopenmp -std=gnu11 -Wall -Wpedantic -Ofast 3MatrixMultiplex/*.c -o 3MatrixMultiplex-exe

mini:
	gcc -fopenmp -DMINI_DATASET -DBENCH -std=gnu11 -Wall -Wpedantic -Ofast 3MatrixMultiplex/*.c -o 3MatrixMultiplex-exe

small:
	gcc -fopenmp -DSMALL_DATASET -DBENCH -std=gnu11 -Wall -Wpedantic -Ofast 3MatrixMultiplex/*.c -o 3MatrixMultiplex-exe

medium:
	gcc -fopenmp -DMEDIUM_DATASET -DBENCH -std=gnu11 -Wall -Wpedantic -Ofast 3MatrixMultiplex/*.c -o 3MatrixMultiplex-exe
	
large:
	gcc -fopenmp -DLARGE_DATASET -DBENCH -std=gnu11 -Wall -Wpedantic -Ofast 3MatrixMultiplex/*.c -o 3MatrixMultiplex-exe

extra:
	gcc -fopenmp -DEXTRALARGE_DATASET -DBENCH -std=gnu11 -Wall -Wpedantic -Ofast 3MatrixMultiplex/*.c -o 3MatrixMultiplex-exe

noopt-user:
	gcc -fopenmp -std=gnu11 -Wall -Wpedantic -O0 3MatrixMultiplex/*.c -o 3MatrixMultiplex-exe

noopt-mini:
	gcc -fopenmp -DMINI_DATASET -DBENCH -std=gnu11 -Wall -Wpedantic -O0 3MatrixMultiplex/*.c -o 3MatrixMultiplex-exe
	
noopt-small:
	gcc -fopenmp -DSMALL_DATASET -DBENCH -std=gnu11 -Wall -Wpedantic -O0 3MatrixMultiplex/*.c -o 3MatrixMultiplex-exe
	
noopt-medium:
	gcc -fopenmp -DMEDIUM_DATASET -DBENCH -std=gnu11 -Wall -Wpedantic -O0 3MatrixMultiplex/*.c -o 3MatrixMultiplex-exe
	
noopt-large:
	gcc -fopenmp -DLARGE_DATASET -DBENCH -std=gnu11 -Wall -Wpedantic -O0 3MatrixMultiplex/*.c -o 3MatrixMultiplex-exe
	
noopt-extra:
	gcc -fopenmp -DEXTRALARGE_DATASET -DBENCH -std=gnu11 -Wall -Wpedantic -O0 3MatrixMultiplex/*.c -o 3MatrixMultiplex-exe
	
clean:
	rm *-exe
	
distclean:
	rm *-exe
	rm -rf benchmarks

