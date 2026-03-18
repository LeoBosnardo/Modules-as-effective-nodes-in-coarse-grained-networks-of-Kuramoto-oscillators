#! /bin/bash
#
gfortran -c -Wall epsilon_lambda.f90
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
gfortran epsilon_lambda.o $HOME/Documents/Mestrado/lib/rkf45.o
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm epsilon_lambda.o
#
mv a.out epsilon_lambda
./epsilon_lambda
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm epsilon_lambda
#
echo "Normal end of execution."
