
help:
	echo "Available makes: [noopt-]user, [noopt-]<dataset size>, clean, distclean"
	echo "But bevare, dataset builds are less friendly"
	echo "Clean option removes executable"
	echo "Distclean also removes benchmarks"
	
user:
	mpicc -std=c99 -Wall -O3 3MatrixMultiplex/*.c -o 3MatrixMultiplex.unix

mini:
	mpicc -DMINI_DATASET -DBENCH -std=c99 -Wall -O3 3MatrixMultiplex/*.c -o 3MatrixMultiplex.unix

small:
	mpicc -DSMALL_DATASET -DBENCH -std=c99 -Wall  -O3 3MatrixMultiplex/*.c -o 3MatrixMultiplex.unix

medium:
	mpicc -DMEDIUM_DATASET -DBENCH -std=c99 -Wall  -O3 3MatrixMultiplex/*.c -o 3MatrixMultiplex.unix
	
large:
	mpicc -DLARGE_DATASET -DBENCH -std=c99 -Wall  -O3 3MatrixMultiplex/*.c -o 3MatrixMultiplex.unix

extra:
	mpicc -DEXTRALARGE_DATASET -DBENCH -std=c99 -Wall  -O3 3MatrixMultiplex/*.c -o 3MatrixMultiplex.unix

noopt-user:
	mpicc -std=c99 -Wall  -O0 3MatrixMultiplex/*.c -o 3MatrixMultiplex.unix

noopt-mini:
	mpicc -DMINI_DATASET -DBENCH -std=c99 -Wall  -O0 3MatrixMultiplex/*.c -o 3MatrixMultiplex.unix
	
noopt-small:
	mpicc -DSMALL_DATASET -DBENCH -std=c99 -Wall  -O0 3MatrixMultiplex/*.c -o 3MatrixMultiplex.unix
	
noopt-medium:
	mpicc -DMEDIUM_DATASET -DBENCH -std=c99 -Wall  -O0 3MatrixMultiplex/*.c -o 3MatrixMultiplex.unix
	
noopt-large:
	mpicc -DLARGE_DATASET -DBENCH -std=c99 -Wall  -O0 3MatrixMultiplex/*.c -o 3MatrixMultiplex.unix
	
noopt-extra:
	mpicc -DEXTRALARGE_DATASET -DBENCH -std=c99 -Wall  -O0 3MatrixMultiplex/*.c -o 3MatrixMultiplex.unix
	
clean:
	rm *.unix
	
distclean:
	rm *.unix
	rm -rf benchmarks

