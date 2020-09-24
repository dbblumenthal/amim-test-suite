#!/bin/bash

cd src/CLAPACK-3.2.1/
cp make.inc.example make.inc
make f2clib
cp F2CLIBS/libf2c.a ../
make tmglib
cp tmglib_LINUX.a ../libtmg.a
make blaslib
cp blas_LINUX.a ../libblas.a
make lapacklib
cp lapack_LINUX.a ../liblapack.a
cd ../..

#end
