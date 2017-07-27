#!/bin/bash

cd outputs/theta_w
rm -f *.vts *.pvd
cd ../..

cd outputs/tempera 
rm -f *.vts *.pvd
cd ../..

cd outputs/density 
rm -f *.vts *.pvd
cd ../.. 

cd outputs/grain_d
rm -f *.vts *.pvd
cd ../..

cd outputs/ 
rm -f *.dat
