#!/bin/bash

# numero de rodagens
# nao pode ter espaço n_runs = X
n_runs=5
l_runs=50
n_threads=4

# compilar
g++ -O2 -std=c++20 -o r r.cpp

for run in $(seq 1 $n_runs); do

  for i in $(seq 1 $l_runs); do

    ./r $run $l_runs $i &
    (( $(jobs -r | wc -l) >= n_threads )) && wait -n

  done

done

wait

# apaga executavel
rm r
