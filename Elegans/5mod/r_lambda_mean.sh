#! /bin/bash
#
gfortran -c -Wall r_lambda_mean.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran r_lambda_mean.o $HOME/Documents/Mestrado/lib/rkf45.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm r_lambda_mean.o
#
mv a.out r_lambda_mean
./r_lambda_mean
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm r_lambda_mean
#
echo "Normal end of execution."
