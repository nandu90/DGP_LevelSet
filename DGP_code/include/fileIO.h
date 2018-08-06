

#ifndef FILEIO_H
#define FILEIO_H
#include "common.h"

void control();

void output_xml(struct elemsclr, int, double **, double **);

void prevfileread(struct elemsclr, double *, double **, double **);

void filewrite(struct elemsclr, double, int);

#endif

