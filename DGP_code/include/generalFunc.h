
#ifndef GENERALFUNC_H
#define GENERALFUNC_H

char* getexepath();

char* concat(char *, char *);

double max(double, double);

double min(double, double);

void naturalToCartesian(double **, double **, double **, int, int);

void errorNormL2(double ***, double ***, double *, double *);

void errorNormL1(double ***, double ***, double *, double *);

void calc_vf(double ***, double, double *);
#endif
