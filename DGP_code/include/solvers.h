/***************************************************************************

Author: nsaini
Created: 2018-03-11

***************************************************************************/

#ifndef SOLVERS_H
#define SOLVERS_H

void solveSystem( double **, double *, double *);

void dgesv_ (int *, int *, double *, int *, int *, double *, int *, int *);

void dgels_ (char *, int *, int *, int *, double *, int *, double *, int *, double *, int *, int *);

void euler(double ***, double ****, double ***, double);

#endif
