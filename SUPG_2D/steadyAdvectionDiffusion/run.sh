#!/bin/bash

PROCS=$1

echo "Running Job on $PROCS processors"
echo "See RUN folder"
mkdir RUN
cd RUN
cp ../control.txt .
mpirun -np $PROCS ../bin/DGPLS

