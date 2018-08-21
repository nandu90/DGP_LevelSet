/***************************************************************************

Author: nsaini
Created: 2018-03-04

***************************************************************************/

#ifndef INPUT_H
#define INPUT_H

//All input variables to be read from the control file are specified here

//------------------------------------------------------------------------//
//Mesh related
int gxelem;
int gyelem;

double xlen;
double ylen;
//------------------------------------------------------------------------//


//------------------------------------------------------------------------//
// DGP related
int polyorder;
int quadtype;
int basistype;
int limit;
//------------------------------------------------------------------------//


//------------------------------------------------------------------------//
//Level-Set
double xb_in;
double yb_in;
double rb_in;

//------------------------------------------------------------------------//


//------------------------------------------------------------------------//
//Boundary Conditions
int x_bound;
int y_bound;

//------------------------------------------------------------------------//

//------------------------------------------------------------------------//
//Time Control
double advect_deltat;
double max_cfl;
int time_control;
int RKstages;
double totaltime;
//------------------------------------------------------------------------//

//------------------------------------------------------------------------//
//Solution Control
int startstep;
int itermax;
int print_gap;
int print_restart_gap;
//------------------------------------------------------------------------//


//------------------------------------------------------------------------//
//Case Control
int case_tog;
//------------------------------------------------------------------------//


//------------------------------------------------------------------------//
//Below parameters are added for INS
double epsilon;
double rhof;
double rhog;
double muf;
double mug;
double sf_coeff;
//------------------------------------------------------------------------//

//------------------------------------------------------------------------//
//INS control parameters
int flow_solve;
int sf_toggle;
//------------------------------------------------------------------------//

//------------------------------------------------------------------------//
//Body force vector
double gx;
double gy;
//------------------------------------------------------------------------//

//------------------------------------------------------------------------//
//Presure control
double ptol;
//------------------------------------------------------------------------//

//------------------------------------------------------------------------//
//Re-distancing control
double re_time;
int re_loops;
int redist_method;
//------------------------------------------------------------------------//

//------------------------------------------------------------------------//
//mesh control
int meshread;
//------------------------------------------------------------------------//


#endif
