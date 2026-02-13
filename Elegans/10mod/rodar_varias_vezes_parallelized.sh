#!/bin/bash

# numero de rodagens
# nao pode ter espaço n_runs = X
l_runs=50
n_threads=4

# compilar
g++ -O2 -std=c++20 -o r r.cpp

for i in $(seq 1 $l_runs); do

  ./r $l_runs $i &
  (( $(jobs -r | wc -l) >= n_threads )) && wait -n

done

wait

# apaga executavel
rm r
