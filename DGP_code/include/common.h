
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
}bhai;


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
    /*double ***p;
    double ***rho;
    double ***mu;*/
};


//------------------------------------------------------------------------//

//------------------------------------------------------------------------//
//DGP related Variables
double **zeta;
double **weights;
int tgauss;       //Total Gauss Quadrature points in the element
int xgpts;        //Number of Gauss points in x
int ygpts;        //Number of Gauss Points in y
int ncoeff;       //Number of coefficients in the solution expansion
                  //The above is equal to the number of basis functions
//------------------------------------------------------------------------//



/*int debug;

///Global Variable declaration (so that we do not have to pass around information between functions)
double nu;
double cfl;
double tol;
int itermax;


    
int advect_steps;
double advect_deltat;
int solnread;
int bub_conv_scheme;
double rhof;
double rhog;
double muf;
double mug;
double epsilon;
double sf_coeff;
double relax;
double ptol;
double re_time;
int re_loops;
int print_gap;
int startstep;
double gx;
double gy;


//Some Simulation Control variables/
int sf_toggle;
int flow_solve;
int p_solver;
int advect_solve;
int sol_type;
int vf_control;
int time_control;
double max_cfl;
int redist_method;
int case_tog;*/












#endif /* COMMON_H */

