/***************************************************************************

Author: nsaini
Created: 2018-03-27

***************************************************************************/

#ifndef RHS_H
#define RHS_H

void getRHS(struct elemsclr, double **, double **, double ***);

void domainIntegral(double **, double **, struct elemsclr, double ***);

void getFluxes(double ***, double ***, struct elemsclr);
#endif
