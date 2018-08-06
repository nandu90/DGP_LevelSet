#!/bin/bash

PROCS=$1
DELETE=$2

if [ $DELETE -gt 0 ]; then
    echo "New case. Deleting old directories"
    rm -rf RUN
    mkdir RUN;
fi
echo "Running Job on $PROCS processors"
echo "See RUN folder"
cp control.txt RUN/
cd RUN
mpirun -np $PROCS ../bin/DGPLS

