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
#endif
