#ifndef DGPFUNC_H
#define DGPFUNC_H

double mappingJacobianDeterminant(int, int, double, double);

void basis1D(double, double *);

void basis2D(double, double, double *);

void massmatrix();
#endif
