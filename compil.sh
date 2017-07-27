#!/bin/bash

mkdir outputs
mkdir ./outputs/density
mkdir ./outputs/grain_d
mkdir ./outputs/tempera
mkdir ./outputs/theta_w

cd src
make
rm -f *.o *.mod
mv SMPP ../SMPP
