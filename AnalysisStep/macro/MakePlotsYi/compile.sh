#!/bin/sh

g++ `root-config --cflags` LeptSfProvider.cpp -c -o LeptSfProvider.o
g++ `root-config --cflags` XSecProvider.cpp -c -o XSecProvider.o
g++ `root-config --glibs` `root-config --cflags` TestRun.cpp -o Run LeptSfProvider.o XSecProvider.o
g++ `root-config --glibs` `root-config --cflags` -I. Plotting/Plotting.cpp -o Plotting/Run

