#! /bin/bash
#
gfortran -c -Wall r_lambda.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran r_lambda.o $HOME/Documents/R/lib/rkf45.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm r_lambda.o
#
mv a.out r_lambda
./r_lambda
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm r_lambda
#
echo "Normal end of execution."
