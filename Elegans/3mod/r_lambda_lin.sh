#! /bin/bash
#
gfortran -c -Wall r_lambda_lin.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran r_lambda_lin.o $HOME/Documents/Mestrado/lib/rkf45.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm r_lambda_lin.o
#
mv a.out r_lambda_lin
./r_lambda_lin
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm r_lambda_lin
#
echo "Normal end of execution."
