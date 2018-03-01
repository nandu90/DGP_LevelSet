#!/bin/bash

PROCS=$1
rm -rf RUN
echo "Running Job on $PROCS processors"
echo "See RUN folder"
mkdir RUN
cp control.txt RUN/
cd RUN
mpirun -np $PROCS ../bin/DGPLS


