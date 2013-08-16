#!/bin/sh
cp -p -R `ls --ignore=test_matrices` /cygdrive/c/MinGW/mingw/bin
cd /cygdrive/c/MinGW/mingw/bin
make clean
make
cd -
