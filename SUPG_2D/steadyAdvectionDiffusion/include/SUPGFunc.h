#ifndef SUPGFUNC_H
#define SUPGFUNC_H

double mappingJacobianDeterminant(int, int, double, double, double **, double **, double *, double*);

void basis1D(double, double *);

void basis2D(double, double, double *);

void mappingFunc(double *, double [], int, double **, double **);

void basisDiff2D( double, double, double *, int);

void basisdiff1D(double, double *);

void GaussPoints1D(double *, double *, int, int);

void GaussPoints2D(double **, double **, int, int);


#endif
