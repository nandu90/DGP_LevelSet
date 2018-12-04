/***************************************************************************

Author: nsaini
Created: 2018-07-24

***************************************************************************/

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "common.h"

void basisdiff1D(double, double *);

void basis1D(double, double *);

double mappingJacobianDeterminant(int, double, double *, double *, double *);

void massmatrix(double **, double *);

void solveSystem(double **, double *, double *, int, int);

void naturalToCartesian(double *, double *, int);

void dgesv_ (int *, int *, double *, int *, int *, double *, int *, int *);

void output(double *, double *, int);

double max(double, double);

double min(double, double);

char* getexepath();

char* concat(char [], char []);

double minmod(double, double, double);

void stiffness(double **, double *);

void convection(double **, double *);

void forceVector(double *, double *);

void weight1D(double, double *, double *, int);

void weightdiff1D(double, double *, double *, int);
#endif
