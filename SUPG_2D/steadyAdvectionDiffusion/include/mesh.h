/*Routine to partition mesh to various processors*/



#ifndef MESH_H
#define MESH_H
#include "common.h"

void partition();
void gridgen(double **, double **, struct elemdata **);
void gridread(double **, double **, double ****, double **, double **, double **);
void meshOperations(double **, double **, double ****, double **, double **, double **);

#endif 
