/***************************************************************************

Author: nsaini
Created: 2018-07-24

***************************************************************************/

#ifndef COMMON_H
#define COMMON_H

///Standard Libraries to include
#define _GNU_SOURCE
#define PI 3.1415926535897

#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <fenv.h>
#include <stdio.h>
#include <limits.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <time.h>
#include <string.h>
#include <stdbool.h>
#include <dirent.h>
#include <fenv.h>

//------------------------------------------------------------------------//
//Input Variable
int polyorder;
int xelem;
int xnode;
double deltat;
double xlen;

int tgauss;
int ncoeff;
int case_tog;

double xb_in;

/*struct elemsclr
{
    double **u;
    double **phi;
    double ***mass;
    };*/

double *zeta;
double *weights;

//------------------------------------------------------------------------//

#endif
