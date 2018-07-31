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

void massmatrix(double ***, double *);

void solveSystem(double **, double *, double *);

void level_setBC(double **);

void initialize(struct elemsclr, double *);

void initializeVel(struct elemsclr, double *);

void initializeLS(struct elemsclr, double *);

void naturalToCartesian(double *, double *, int);

void dgesv_ (int *, int *, double *, int *, int *, double *, int *, int *);

void output(struct elemsclr, double *, int);

double max(double, double);

double min(double, double);

char* getexepath();

char* concat(char [], char []);

void domainIntegral(double *, struct elemsclr, double **);

void getRHS(struct elemsclr, double *, double **);

void Runge_Kutta(struct elemsclr, double *, double, double **);

void eulerIncrement(double **, double ***, double **, double);

void fluxes(double *, double *, struct elemsclr);

double upwind(double, double, double, double);

void boundaryIntegral(double **, double *, double *);

void errorNormL1(double **, double **, double *, double *, double *);

void errorNormL2(double **, double **, double *, double *, double *);

#endif
