/*Routine to partition mesh to various processors*/



#ifndef MESH_H
#define MESH_H

void partition();
void gridgen(double **, double **, double ****, double **, double **, double **);
void gridread(double **, double **, double ****, double **, double **, double **);
void meshOperations(double **, double **, double ****, double **, double **, double **);
void meshShift(double **, double **, double **, double **);
#endif 
