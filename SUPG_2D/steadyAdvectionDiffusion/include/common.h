
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   common.h
 * Author: nsaini3
 *
 * Created on October 27, 2016, 1:44 PM
 */

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
#include <omp.h>
#include <mpi.h>
#include <string.h>
#include <stdbool.h>
#include <dirent.h>
#include <fenv.h>

//------------------------------------------------------------------------//

//All input variables to be read are collected in a separate file
//The input file is included in common.h so that you don't have to call both files - just common.h"
#include "input.h"
//------------------------------------------------------------------------//

//------------------------------------------------------------------------//
//MPI Relevant variables
int master;
int nprocs;
int myrank;

//------------------------------------------------------------------------//

//------------------------------------------------------------------------//
//Mesh related variables
int xelem; //Total elem in x
int yelem; //Total elem in y
int xnode; //Total nodes in x
int ynode; //Total nodes in y

int procm, procn; //Processor matrix m X n
int elemm, elemn; //Mode of number of elements on each processor

int gxnode;
int gynode;

//------------------------------------------------------------------------//


//------------------------------------------------------------------------//
//MPI communication related variables
int bhailog[4];
int per[4];         //Complements bhailog and tells which side is periodic
double **sendptr;
double **recvptr;
int **io_info;

double **INSsendptr;
double **INSrecvptr;

struct bhaiarray
{
    double *sendrbuf;
    double *recvrbuf;
    
    double *sendlbuf;
    double *recvlbuf;
    
    double *sendubuf;
    double *recvubuf;
    
    double *senddbuf;
    double *recvdbuf;
}bhai, INSbhai;


//------------------------------------------------------------------------//

//------------------------------------------------------------------------//
//Element struct. Only the struct is defined here. Everything is allocated inside
//the main program so as to avoid allocating large arrays in common block
struct elemsclr
{
    double ***u;
    double ***v;
    double ***phi;
    double ****mass;
    int **iBC;
};

struct nodesclr
{
    double ***phi;
};


//------------------------------------------------------------------------//

//------------------------------------------------------------------------//
//DGP related Variables
int ncoeff;       //Number of coefficients in the solution expansion
                  //The above is equal to the number of basis functions
//------------------------------------------------------------------------//

//------------------------------------------------------------------------//
//SUPG related Variables
int tdof;
int Tdof;

struct dofdata
{
    int gindex;
    int BC;
    int controlproc;
};


struct elemdata
{
    int edofn;
    int *edofs;
};
//------------------------------------------------------------------------//

#endif /* COMMON_H */

