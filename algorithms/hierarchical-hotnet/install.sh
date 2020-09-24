#!/bin/bash

cd src
f2py -c fortran_module.f95 -m fortran_module > /dev/null

#end
