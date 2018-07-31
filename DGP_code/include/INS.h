/***************************************************************************

Author: nsaini
Created: 2018-04-25

***************************************************************************/


#ifndef INS_H
#define INS_H

#include "common.h"

void find_density_visc(double **, double **, double **);

void heavy_func(double **, double **, double);

void initializeINS(struct elemsclr);

void INSlevel_setBC(double **, int **);

void INSinitialize(struct elemsclr);

void rhscalc(struct elemsclr, double **, double **, double ****, double **, int **, double **, double **);

void flux(double , double , double , double , double , double , double , double , int , double *, double *);

void vel_BC(double **, double **, int **);

void grad_level_setBC(double **, int **);

void surface(struct elemsclr, double **, double **, double ****);

void body(struct elemsclr, double **, double **, double **);

void variable_pressure(double **, double **, double **, double , double **, double **, double **, double ****, double **, int **);

void pointJacobi(double ***, double **, double **, int **, double **);

void pressureBC(double **, int **);

void grad_func(double **, double **, double ****, int **);

void delta_func(double **, double **, double);

void rhs_redist2(double **, double **, double **, double ****, int **);

void hyperbolic(struct elemsclr, double ****);

void heavy_funcDG(double ***, double ***, double);

#endif
