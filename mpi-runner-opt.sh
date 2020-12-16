#!/bin/sh

#  bench-runner.sh
#  3MatrixMultiplex
#
#  Created by Alexander Makhov on 24/11/20.
#  

DATASET=("mini" "small" "medium" "large" "extra")

echo "Preparing to submit tasks..."

if make user ; then
    echo "Seems to compile optimized successfully"
else
    echo "Compilation of opt returned error!"
    exit -1;
fi

BENCHPATH=benchmarks/o3/
mkdir -p $BENCHPATH

for bench in ${DATASET[*]}; do
    make $bench ;
    for i in {0..10}; do
        echo "Submitting opt $bench for $((2**$i)) threads"
        mpisubmit.bg -np $((2**$i)) -w 00:03:00 3MatrixMultiplex-exe -stdout $BENCHPATH/$bench-$((2**$i)).out ;
        touch $BENCHPATH/$bench-$((2**$i)).out
        sleep 20;
        clear
    done
done

echo "Done benchmarks!"






