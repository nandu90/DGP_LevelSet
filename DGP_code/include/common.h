
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

//MPI relevant variables
//MPI_STATUS status;

int debug;
int myrank;
int master;
int nprocs;
int procm, procn; //Processor matrix m X n
int elemm, elemn; //Mode of number of elements on each processor
int **iBC;        //Determines if a particular cell is a boundary or interior cell
int bhailog[4];
double **sendptr;
double **recvptr;
int **io_info;
int polyorder;    //Order of the basis functions

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

///Global Variable declaration (so that we do not have to pass around information between functions)
double nu;
double cfl;
double tol;
int itermax;

//Local nodes and elements//
int xelem; //Total elem in x
int yelem; //Total elem in y
int zelem; //Total elem in z
int xnode; //Total nodes in x
int ynode; //Total nodes in y
int znode; //Total nodes in z
/////////////////////////////
//Global nodes and elements//
int gxelem; //Total elem in x
int gyelem; //Total elem in y
int gzelem; //Total elem in z
int gxnode;
int gynode;
int gznode;
/////////////////////////////

double xlen;
double ylen;
double zlen;

double **x;
double **y;
double **xc;
double **yc;
double **vol;
double ****area;
double ****mass;

    
///Variables for bubble
double rb_in;
double xb_in;
double yb_in;
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
int x_bound;
int y_bound;
int advect_solve;
int sol_type;
int vf_control;
int time_control;
double max_cfl;
int redist_method;
int case_tog;






struct elemsclr
{
  double ***p;
  double ***u;
  double ***v;
  double ***phi;
  double ***rho;
  double ***mu;
};





#endif /* COMMON_H */

