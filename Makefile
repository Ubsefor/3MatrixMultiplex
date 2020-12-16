
help:
	echo "Available makes: [noopt-]user, [noopt-]<dataset size>, clean, distclean"
	echo "But bevare, dataset builds are less friendly"
	echo "Clean option removes executable"
	echo "Distclean also removes benchmarks"
	
user:
	mpicc -std=c99 -Wall -O3 3MatrixMultiplex/*.c -o 3MatrixMultiplex-exe

mini:
	mpicc -DMINI_DATASET -DBENCH -std=c99 -Wall -O3 3MatrixMultiplex/*.c -o 3MatrixMultiplex-exe

small:
	mpicc -DSMALL_DATASET -DBENCH -std=c99 -Wall  -O3 3MatrixMultiplex/*.c -o 3MatrixMultiplex-exe

medium:
	mpicc -DMEDIUM_DATASET -DBENCH -std=c99 -Wall  -O3 3MatrixMultiplex/*.c -o 3MatrixMultiplex-exe
	
large:
	mpicc -DLARGE_DATASET -DBENCH -std=c99 -Wall  -O3 3MatrixMultiplex/*.c -o 3MatrixMultiplex-exe

extra:
	mpicc -DEXTRALARGE_DATASET -DBENCH -std=c99 -Wall  -O3 3MatrixMultiplex/*.c -o 3MatrixMultiplex-exe

noopt-user:
	mpicc -std=c99 -Wall  -O0 3MatrixMultiplex/*.c -o 3MatrixMultiplex-exe

noopt-mini:
	mpicc -DMINI_DATASET -DBENCH -std=c99 -Wall  -O0 3MatrixMultiplex/*.c -o 3MatrixMultiplex-exe
	
noopt-small:
	mpicc -DSMALL_DATASET -DBENCH -std=c99 -Wall  -O0 3MatrixMultiplex/*.c -o 3MatrixMultiplex-exe
	
noopt-medium:
	mpicc -DMEDIUM_DATASET -DBENCH -std=c99 -Wall  -O0 3MatrixMultiplex/*.c -o 3MatrixMultiplex-exe
	
noopt-large:
	mpicc -DLARGE_DATASET -DBENCH -std=c99 -Wall  -O0 3MatrixMultiplex/*.c -o 3MatrixMultiplex-exe
	
noopt-extra:
	mpicc -DEXTRALARGE_DATASET -DBENCH -std=c99 -Wall  -O0 3MatrixMultiplex/*.c -o 3MatrixMultiplex-exe
	
clean:
	rm *-exe
	
distclean:
	rm *-exe
	rm -rf benchmarks

