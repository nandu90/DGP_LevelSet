/***************************************************************************

Author: nsaini
Created: 2018-03-11

***************************************************************************/

#ifndef SOLVERS_H
#define SOLVERS_H

#include "common.h"

void itrdrv(struct elemsclr ,double **, double **, double **, double **, double ** ,double ****);

void solveSystem( double **, double *, double *, int, int);

void dgesv_ (int *, int *, double *, int *, int *, double *, int *, int *);

void dgels_ (char *, int *, int *, int *, double *, int *, double *, int *, double *, int *, int *);

void euler(double ***, double ****, double ***, double);

void eulerIncrement(double ***, double ****, double ***, double);

void Runge_Kutta(struct elemsclr, double **, double **, double, double ***, double ****);

#endif
