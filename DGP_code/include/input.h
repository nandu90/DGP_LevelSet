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
//------------------------------------------------------------------------//

//------------------------------------------------------------------------//
//Solution Control
int startstep;
int itermax;
int print_gap;
//------------------------------------------------------------------------//


//------------------------------------------------------------------------//
//Case Control
int case_tog;
//------------------------------------------------------------------------//


#endif
