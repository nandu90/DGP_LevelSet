/***************************************************************************

Author: nsaini
Created: 2018-03-27

***************************************************************************/

#ifndef RHS_H
#define RHS_H

void getRHS(struct elemsclr, double **, double **, double ***, double ****);

void domainIntegral(double **, double **, struct elemsclr, double ***);

void fluxes(double ***, double ***, double **, double **, struct elemsclr);

double upwind(double , double , double , double);

void boundaryIntegral(double ***, double ***, double ***, double **, double **, double ****);

#endif
