#!/bin/sh

#  bench-runner.sh
#  3MatrixMultiplex
#
#  Created by Alexander Makhov on 24/11/20.
#  

DATASET=("mini" "small" "medium" "large" "extra")

echo "Preparing to run benchmarks..."

if make noopt-user ; then
    echo "Seems to compile noopt successfully"
else
    echo "Compilation of noopt returned error!"
    exit -1;
fi

BENCHPATH=./benchmarks/noopt/

mkdir -p $BENCHPATH

for bench in ${DATASET[*]}; do
    make noopt-$bench ;
    for i in 0..10; do
        ./3MatrixMultiplex-exe $((2**i)) >> $BENCHPATH/$bench ;
        echo "" >> $BENCHPATH/$bench
    done
done

echo "Done noopt benchmarks!"

if make user ; then
    echo "Seems to compile noopt successfully"
else
    echo "Compilation of noopt returned error!"
    exit -1;
fi

BENCHPATH=benchmarks/ofast/
mkdir -p $BENCHPATH

for bench in ${DATASET[*]}; do
    make $bench ;
    for i in 0..10; do
        ./3MatrixMultiplex-exe $((2**i)) >> $BENCHPATH/$bench ;
        echo "" >> $BENCHPATH/$bench
    done
done

echo "Done benchmarks!"






