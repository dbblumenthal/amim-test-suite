#!/bin/bash

cd src
g++ -o ../clustex2 -I Eigen/ -fopenmp classes.cpp main.cpp liblapack.a libblas.a libf2c.a libtmg.a
cd ..

#end
