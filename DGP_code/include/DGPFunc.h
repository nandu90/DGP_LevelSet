#ifndef DGPFUNC_H
#define DGPFUNC_H

double mappingJacobianDeterminant(int, int, double, double, double **, double **, double **);

void basis1D(double, double *);

void basis2D(double, double, double *);

void massmatrix(double **, double **, double ****);

void mappingFunc(double *, double []);

void basisDiff2D( double, double, double *, int);

void basisdiff1D(double, double *);

double lineJacobian(int, int, double, double **, int);
#endif
