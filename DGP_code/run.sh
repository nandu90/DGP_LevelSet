#!/bin/bash

PROCS=$1
rm -rf RUN
echo "Running Job on $PROCS processors"
echo "See RUN folder"
mkdir RUN
cp control.txt RUN/
cd RUN
mpirun -np $PROCS ../bin/DGPLS

cd ..
rm -r RUN0
mkdir RUN0
cp -r RUN/* RUN0/
