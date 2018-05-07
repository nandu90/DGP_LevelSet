/***************************************************************************

Author: nsaini
Created: 2018-03-04

***************************************************************************/

#ifndef ICBC_H
#define ICBC_H

void initialize(struct elemsclr, double **, double **);

void initializeVel(struct elemsclr, double **, double **);

void initializeLS(struct elemsclr, double **, double **);


//------------------------------------------------------------------------//
//Boundary Conditions
void level_setBC(double ***, int **);
//------------------------------------------------------------------------//


//------------------------------------------------------------------------//
//Special Cases
double zalesak(double, double);
void GaussianStep(double, double, double*);
//------------------------------------------------------------------------//

#endif
