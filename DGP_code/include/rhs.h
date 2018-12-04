/***************************************************************************

Author: nsaini
Created: 2018-03-27

***************************************************************************/

#ifndef RHS_H
#define RHS_H

void getRHS(struct elemsclr, double **, double **, double ***, double ****);

void domainIntegral(double **, double **, struct elemsclr, double ***);

void fluxes(double ***, double ***, double **, double **, struct elemsclr, double *, double *, int, double * , double *, int);

double upwind(double , double , double , double);

void boundaryIntegral(double ***, double ***, double ***, double **, double **, double ****, double *, double *, int, double * , double *, int);


void sourceIntegral(double **, double **, struct elemsclr, double ***);

double sourceTerm(double, double);

double getterm1(double, double, double);
double getterm2(double, double, double);
double getterm3(double, double, double);
double getphi(double, double, double);
double acos2(double);
double sqrt2(double);

#endif
