/***************************************************************************

Author: nsaini
Created: 2018-11-11

***************************************************************************/


#ifndef SUPGSOLVER_H
#define SUPGSOLVER_H

#include "common.h"

void stiffness(double ***, double ***, double **, double **);

void convection(double ***, double ***, struct elemsclr, double **, double **);

void massMatrix(double ****, double **, double **);

void forceVector(double ***, double **, double **);
#endif
