/***************************************************************************

Author: nsaini
Created: 2018-11-29

***************************************************************************/


#include "common.h"
#include "memory.h"
#include "polylib.h"
#include "DGPFunc.h"
#include "rhs.h"
#include "solvers.h"
#include "generalFunc.h"


void sourceIntegral(double **x, double **y, struct elemsclr elem, double ***rhs)
{
    //------------------------------------------------------------------------//
    //Loop Indexes
    int ielem, jelem;
    int igauss;
    int icoeff, icoeff1;
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Temporary variables
    double *basis;
    allocator1(&basis, ncoeff);

    double *inv,  *jacobian;
    allocator1(&inv, 4);
    allocator1(&jacobian, 4);

    double detJ;

    double theta;
    double arm;
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Define quad points and weights here 
    int extra;
    if(quadtype == 1)
    {
	extra = 1;
    }
    else
    {
	extra = 0;
    }
    double **zeta, **weights;
    int tgauss = pow(polyorder + 1 + extra, 2);

    allocator2(&zeta, tgauss,2);
    allocator2(&weights, tgauss,2);
    
    GaussPoints2D(zeta, weights, quadtype, tgauss); 
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Allocate source vector - source w.r.t Cartesian Coordinates
    double *cartSource;
    allocator1(&cartSource, tgauss);

    //Allocate the Vandermonde matrix
    double **vand;
    allocator2(&vand, tgauss, ncoeff);

    //Loop over the quadrature points to fill the Vandermonde Matrix
    for(igauss=0; igauss<tgauss; igauss++)
    {
	//Get the basis vector
	basis2D(zeta[igauss][0], zeta[igauss][1], basis);
	//Fill up row of the Vandermonde matrix
	for(icoeff=0; icoeff<ncoeff; icoeff++)
	{
	    vand[igauss][icoeff] = basis[icoeff];
	}
    }

    //Allocate coordinate matrix corresponding to zs - solution points
    double **xs;
    allocator2(&xs, tgauss, 2);

    double *naturalSource;
    allocator1(&naturalSource, ncoeff);

    double *sourceContribution;
    allocator1(&sourceContribution, ncoeff);

    double recsource;
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Loop over the elements
    for(ielem = 1; ielem<xelem-1; ielem++)
    {
	for(jelem = 1; jelem<yelem-1; jelem++)
	{
	    //Initialize contribution vector to 0
	    for(icoeff=0; icoeff<ncoeff; icoeff++)
	    {
		sourceContribution[icoeff] = 0.0;
	    }
	    
	    //Convert natural coordinates at quadrature points to Cartesian
	    naturalToCartesian(xs, x, y, ielem, jelem, zeta, tgauss);

	    //Get the source values at Cartesian Quadrature points
	    for(igauss=0; igauss<tgauss; igauss++)
	    {
		theta = atan2(xs[igauss][1], xs[igauss][0]);
		arm = sqrt(pow(xs[igauss][0],2.0) + pow(xs[igauss][1],2.0));
		cartSource[igauss] = sourceTerm(theta, arm);
	    }

	    //Solve the system to get the source coefficents in natural coordinates
	    solveSystem(vand, cartSource, naturalSource, tgauss, ncoeff);
	    
	    //Loop over the test functions
	    for(icoeff=0; icoeff<ncoeff; icoeff++)
	    {
		//Loop over the Gauss Quadrature points
		for(igauss =0; igauss<tgauss; igauss++)
		{
		    //Get the Test function value at this Quadrature point
		    basis2D(zeta[igauss][0], zeta[igauss][1], basis);

		    //Get the Determinant
		    detJ = mappingJacobianDeterminant(ielem, jelem, zeta[igauss][0], zeta[igauss][1], x, y, inv, jacobian);

		    //Reconstruct the source term at the Quadrature point
		    recsource = 0.0;
		    for(icoeff1 =0; icoeff1<ncoeff; icoeff1++)
		    {
			recsource += basis[icoeff1]*naturalSource[icoeff1];
		    }
		    
		    //Add to the source contribution
		    sourceContribution[icoeff] += weights[igauss][0]*weights[igauss][1]*basis[icoeff]*detJ*recsource;
		}
	    }

	    //Finally add to the Residual Term
	    for(icoeff=0; icoeff<ncoeff; icoeff++)
	    {
		rhs[ielem][jelem][icoeff] += sourceContribution[icoeff];
	    }
	}
    }
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Deallocators
    deallocator1(&basis, ncoeff);
    deallocator1(&inv, 4);
    deallocator1(&jacobian, 4);
    deallocator2(&zeta, tgauss,2);
    deallocator2(&weights, tgauss,2);
    deallocator1(&cartSource, tgauss);
    deallocator2(&vand, tgauss, ncoeff);
    deallocator2(&xs, tgauss, 2);
    deallocator1(&naturalSource, ncoeff);
    deallocator1(&sourceContribution, ncoeff);
    //------------------------------------------------------------------------//


}


double sourceTerm(double theta, double r)
{
    double t = simtime;
    double source;

    theta = fabs(theta);
    
    if(theta > PI/2) //2nd or 3rd quadrant
    {
	if(r == 0.0 && theta != PI)
	{
	    
	}
	else if(r != 0.0 && theta == PI)
	{

	}
	else if(r == 0.0 && theta == PI)
	{

	}
	else
	{
	    source = ((sin(theta)*(t*(theta*(-5.0/2.0)+3.141592653589793*(5.0/2.0)+pow(theta-3.141592653589793,2.0)*(5.0/2.0)+pow(theta-3.141592653589793,3.0)*4.058888755698753-sqrt(-theta+3.141592653589793)+pow(theta-3.141592653589793,4.0)*1.295846037797472)+1.0)*2.0-t*(atan(theta-3.141592653589793)*1.364494638965387E-2+1.014485301813865/(exp(theta*1.0E1-3.141592653589793*(1.5E1/2.0))+1.0)-3.93673619140503E-4)*3.273239544735159)*(atan(theta-3.141592653589793)*2.233158905420314E-2+1.683771314041424E31/(exp(theta*1.0E1-3.141592653589793*(1.5E1/2.0))*1.014120480182584E31+1.014120480182584E31)-sin(theta)*(theta*(-5.0/2.0)+3.141592653589793*(5.0/2.0)+pow(theta-3.141592653589793,2.0)*(5.0/2.0)+pow(theta-3.141592653589793,3.0)*4.058888755698753-sqrt(-theta+3.141592653589793)+pow(theta-3.141592653589793,4.0)*1.295846037797472)-6.442940289448513E-4)-cos(t)*pow(cos(theta),2.0)*(sin(t)+1.0)*2.0)*1.0/sqrt(pow(cos(theta),2.0)*pow(sin(t)+1.0,2.0)+pow(sin(theta)*(t*(theta*(-5.0/2.0)+3.141592653589793*(5.0/2.0)+pow(theta-3.141592653589793,2.0)*(5.0/2.0)+pow(theta-3.141592653589793,3.0)*4.058888755698753-sqrt(-theta+3.141592653589793)+pow(theta-3.141592653589793,4.0)*1.295846037797472)+1.0)-t*(atan(theta-3.141592653589793)*1.364494638965387E-2+1.014485301813865/(exp(theta*1.0E1-3.141592653589793*(1.5E1/2.0))+1.0)-3.93673619140503E-4)*1.63661977236758,2.0))*(1.0/2.0)+cos(theta)*(-cos(r*cos(theta))*cos(r*sin(theta))*cos(theta)*(sqrt(pow(cos(theta),2.0)*pow(sin(t)+1.0,2.0)+pow(sin(theta)*(t*(theta*(-5.0/2.0)+3.141592653589793*(5.0/2.0)+pow(theta-3.141592653589793,2.0)*(5.0/2.0)+pow(theta-3.141592653589793,3.0)*4.058888755698753-sqrt(-theta+3.141592653589793)+pow(theta-3.141592653589793,4.0)*1.295846037797472)+1.0)-t*(atan(theta-3.141592653589793)*1.364494638965387E-2+1.014485301813865/(exp(theta*1.0E1-3.141592653589793*(1.5E1/2.0))+1.0)-3.93673619140503E-4)*1.63661977236758,2.0))-sqrt((r*r)*pow(cos(theta),2.0)+(r*r)*pow(sin(theta),2.0)))+sin(r*cos(theta))*sin(r*sin(theta))*sin(theta)*(sqrt(pow(cos(theta),2.0)*pow(sin(t)+1.0,2.0)+pow(sin(theta)*(t*(theta*(-5.0/2.0)+3.141592653589793*(5.0/2.0)+pow(theta-3.141592653589793,2.0)*(5.0/2.0)+pow(theta-3.141592653589793,3.0)*4.058888755698753-sqrt(-theta+3.141592653589793)+pow(theta-3.141592653589793,4.0)*1.295846037797472)+1.0)-t*(atan(theta-3.141592653589793)*1.364494638965387E-2+1.014485301813865/(exp(theta*1.0E1-3.141592653589793*(1.5E1/2.0))+1.0)-3.93673619140503E-4)*1.63661977236758,2.0))-sqrt((r*r)*pow(cos(theta),2.0)+(r*r)*pow(sin(theta),2.0)))+cos(r*sin(theta))*sin(r*cos(theta))*1.0/sqrt((r*r)*pow(cos(theta),2.0)+(r*r)*pow(sin(theta),2.0))*(r*pow(cos(theta),2.0)*2.0+r*pow(sin(theta),2.0)*2.0)*(1.0/2.0))-sin(theta)*(-cos(r*cos(theta))*cos(r*sin(theta))*sin(theta)*(sqrt(pow(cos(theta),2.0)*pow(sin(t)+1.0,2.0)+pow(sin(theta)*(t*(theta*(-5.0/2.0)+3.141592653589793*(5.0/2.0)+pow(theta-3.141592653589793,2.0)*(5.0/2.0)+pow(theta-3.141592653589793,3.0)*4.058888755698753-sqrt(-theta+3.141592653589793)+pow(theta-3.141592653589793,4.0)*1.295846037797472)+1.0)-t*(atan(theta-3.141592653589793)*1.364494638965387E-2+1.014485301813865/(exp(theta*1.0E1-3.141592653589793*(1.5E1/2.0))+1.0)-3.93673619140503E-4)*1.63661977236758,2.0))-sqrt((r*r)*pow(cos(theta),2.0)+(r*r)*pow(sin(theta),2.0)))+sin(r*cos(theta))*sin(r*sin(theta))*cos(theta)*(sqrt(pow(cos(theta),2.0)*pow(sin(t)+1.0,2.0)+pow(sin(theta)*(t*(theta*(-5.0/2.0)+3.141592653589793*(5.0/2.0)+pow(theta-3.141592653589793,2.0)*(5.0/2.0)+pow(theta-3.141592653589793,3.0)*4.058888755698753-sqrt(-theta+3.141592653589793)+pow(theta-3.141592653589793,4.0)*1.295846037797472)+1.0)-t*(atan(theta-3.141592653589793)*1.364494638965387E-2+1.014485301813865/(exp(theta*1.0E1-3.141592653589793*(1.5E1/2.0))+1.0)-3.93673619140503E-4)*1.63661977236758,2.0))-sqrt((r*r)*pow(cos(theta),2.0)+(r*r)*pow(sin(theta),2.0)))+cos(r*cos(theta))*sin(r*sin(theta))*1.0/sqrt((r*r)*pow(cos(theta),2.0)+(r*r)*pow(sin(theta),2.0))*(r*pow(cos(theta),2.0)*2.0+r*pow(sin(theta),2.0)*2.0)*(1.0/2.0))+(cos(theta)*(cos(r*cos(theta))*sin(r*sin(theta))*((sin(theta)*(t*(theta*(-5.0/2.0)+3.141592653589793*(5.0/2.0)+pow(theta-3.141592653589793,2.0)*(5.0/2.0)+pow(theta-3.141592653589793,3.0)*4.058888755698753-sqrt(-theta+3.141592653589793)+pow(theta-3.141592653589793,4.0)*1.295846037797472)+1.0)*2.0-t*(atan(theta-3.141592653589793)*1.364494638965387E-2+1.014485301813865/(exp(theta*1.0E1-3.141592653589793*(1.5E1/2.0))+1.0)-3.93673619140503E-4)*3.273239544735159)*(cos(theta)*(t*(theta*(-5.0/2.0)+3.141592653589793*(5.0/2.0)+pow(theta-3.141592653589793,2.0)*(5.0/2.0)+pow(theta-3.141592653589793,3.0)*4.058888755698753-sqrt(-theta+3.141592653589793)+pow(theta-3.141592653589793,4.0)*1.295846037797472)+1.0)-t*(1.364494638965387E-2/(pow(theta-3.141592653589793,2.0)+1.0)-exp(theta*1.0E1-3.141592653589793*(1.5E1/2.0))*1.0/pow(exp(theta*1.0E1-3.141592653589793*(1.5E1/2.0))+1.0,2.0)*1.014485301813865E1)*1.63661977236758+t*sin(theta)*(theta*5.0-3.141592653589793*5.0+pow(theta-3.141592653589793,2.0)*1.217666626709626E1+pow(theta-3.141592653589793,3.0)*5.183384151189888+1.0/sqrt(-theta+3.141592653589793)*(1.0/2.0)-5.0/2.0))-cos(theta)*sin(theta)*pow(sin(t)+1.0,2.0)*2.0)*1.0/sqrt(pow(cos(theta),2.0)*pow(sin(t)+1.0,2.0)+pow(sin(theta)*(t*(theta*(-5.0/2.0)+3.141592653589793*(5.0/2.0)+pow(theta-3.141592653589793,2.0)*(5.0/2.0)+pow(theta-3.141592653589793,3.0)*4.058888755698753-sqrt(-theta+3.141592653589793)+pow(theta-3.141592653589793,4.0)*1.295846037797472)+1.0)-t*(atan(theta-3.141592653589793)*1.364494638965387E-2+1.014485301813865/(exp(theta*1.0E1-3.141592653589793*(1.5E1/2.0))+1.0)-3.93673619140503E-4)*1.63661977236758,2.0))*(1.0/2.0)+r*cos(r*cos(theta))*cos(r*sin(theta))*cos(theta)*(sqrt(pow(cos(theta),2.0)*pow(sin(t)+1.0,2.0)+pow(sin(theta)*(t*(theta*(-5.0/2.0)+3.141592653589793*(5.0/2.0)+pow(theta-3.141592653589793,2.0)*(5.0/2.0)+pow(theta-3.141592653589793,3.0)*4.058888755698753-sqrt(-theta+3.141592653589793)+pow(theta-3.141592653589793,4.0)*1.295846037797472)+1.0)-t*(atan(theta-3.141592653589793)*1.364494638965387E-2+1.014485301813865/(exp(theta*1.0E1-3.141592653589793*(1.5E1/2.0))+1.0)-3.93673619140503E-4)*1.63661977236758,2.0))-sqrt((r*r)*pow(cos(theta),2.0)+(r*r)*pow(sin(theta),2.0)))+r*sin(r*cos(theta))*sin(r*sin(theta))*sin(theta)*(sqrt(pow(cos(theta),2.0)*pow(sin(t)+1.0,2.0)+pow(sin(theta)*(t*(theta*(-5.0/2.0)+3.141592653589793*(5.0/2.0)+pow(theta-3.141592653589793,2.0)*(5.0/2.0)+pow(theta-3.141592653589793,3.0)*4.058888755698753-sqrt(-theta+3.141592653589793)+pow(theta-3.141592653589793,4.0)*1.295846037797472)+1.0)-t*(atan(theta-3.141592653589793)*1.364494638965387E-2+1.014485301813865/(exp(theta*1.0E1-3.141592653589793*(1.5E1/2.0))+1.0)-3.93673619140503E-4)*1.63661977236758,2.0))-sqrt((r*r)*pow(cos(theta),2.0)+(r*r)*pow(sin(theta),2.0)))))/r-(sin(theta)*(cos(r*sin(theta))*sin(r*cos(theta))*((sin(theta)*(t*(theta*(-5.0/2.0)+3.141592653589793*(5.0/2.0)+pow(theta-3.141592653589793,2.0)*(5.0/2.0)+pow(theta-3.141592653589793,3.0)*4.058888755698753-sqrt(-theta+3.141592653589793)+pow(theta-3.141592653589793,4.0)*1.295846037797472)+1.0)*2.0-t*(atan(theta-3.141592653589793)*1.364494638965387E-2+1.014485301813865/(exp(theta*1.0E1-3.141592653589793*(1.5E1/2.0))+1.0)-3.93673619140503E-4)*3.273239544735159)*(cos(theta)*(t*(theta*(-5.0/2.0)+3.141592653589793*(5.0/2.0)+pow(theta-3.141592653589793,2.0)*(5.0/2.0)+pow(theta-3.141592653589793,3.0)*4.058888755698753-sqrt(-theta+3.141592653589793)+pow(theta-3.141592653589793,4.0)*1.295846037797472)+1.0)-t*(1.364494638965387E-2/(pow(theta-3.141592653589793,2.0)+1.0)-exp(theta*1.0E1-3.141592653589793*(1.5E1/2.0))*1.0/pow(exp(theta*1.0E1-3.141592653589793*(1.5E1/2.0))+1.0,2.0)*1.014485301813865E1)*1.63661977236758+t*sin(theta)*(theta*5.0-3.141592653589793*5.0+pow(theta-3.141592653589793,2.0)*1.217666626709626E1+pow(theta-3.141592653589793,3.0)*5.183384151189888+1.0/sqrt(-theta+3.141592653589793)*(1.0/2.0)-5.0/2.0))-cos(theta)*sin(theta)*pow(sin(t)+1.0,2.0)*2.0)*1.0/sqrt(pow(cos(theta),2.0)*pow(sin(t)+1.0,2.0)+pow(sin(theta)*(t*(theta*(-5.0/2.0)+3.141592653589793*(5.0/2.0)+pow(theta-3.141592653589793,2.0)*(5.0/2.0)+pow(theta-3.141592653589793,3.0)*4.058888755698753-sqrt(-theta+3.141592653589793)+pow(theta-3.141592653589793,4.0)*1.295846037797472)+1.0)-t*(atan(theta-3.141592653589793)*1.364494638965387E-2+1.014485301813865/(exp(theta*1.0E1-3.141592653589793*(1.5E1/2.0))+1.0)-3.93673619140503E-4)*1.63661977236758,2.0))*(-1.0/2.0)+r*cos(r*cos(theta))*cos(r*sin(theta))*sin(theta)*(sqrt(pow(cos(theta),2.0)*pow(sin(t)+1.0,2.0)+pow(sin(theta)*(t*(theta*(-5.0/2.0)+3.141592653589793*(5.0/2.0)+pow(theta-3.141592653589793,2.0)*(5.0/2.0)+pow(theta-3.141592653589793,3.0)*4.058888755698753-sqrt(-theta+3.141592653589793)+pow(theta-3.141592653589793,4.0)*1.295846037797472)+1.0)-t*(atan(theta-3.141592653589793)*1.364494638965387E-2+1.014485301813865/(exp(theta*1.0E1-3.141592653589793*(1.5E1/2.0))+1.0)-3.93673619140503E-4)*1.63661977236758,2.0))-sqrt((r*r)*pow(cos(theta),2.0)+(r*r)*pow(sin(theta),2.0)))+r*sin(r*cos(theta))*sin(r*sin(theta))*cos(theta)*(sqrt(pow(cos(theta),2.0)*pow(sin(t)+1.0,2.0)+pow(sin(theta)*(t*(theta*(-5.0/2.0)+3.141592653589793*(5.0/2.0)+pow(theta-3.141592653589793,2.0)*(5.0/2.0)+pow(theta-3.141592653589793,3.0)*4.058888755698753-sqrt(-theta+3.141592653589793)+pow(theta-3.141592653589793,4.0)*1.295846037797472)+1.0)-t*(atan(theta-3.141592653589793)*1.364494638965387E-2+1.014485301813865/(exp(theta*1.0E1-3.141592653589793*(1.5E1/2.0))+1.0)-3.93673619140503E-4)*1.63661977236758,2.0))-sqrt((r*r)*pow(cos(theta),2.0)+(r*r)*pow(sin(theta),2.0)))))/r;
	}
    }

    else
    {
	if(r == 0.0 && theta != 0.0)
	{
	    
	}
	else if(r != 0.0 && theta == 0.0)
	{

	}
	else if(r == 0.0 && theta == 0.0)
	{

	}
	else
	{
	    source  = ((sin(theta)*(t*(theta*(5.0/2.0)+(theta*theta)*(5.0/2.0)-sqrt(theta)-(theta*theta*theta)*4.058888755698753+(theta*theta*theta*theta)*1.295846037797472)+1.0)*2.0+t*(atan(theta)*1.364494638965387E-2-1.014485301813865/(exp(theta*-1.0E1+3.141592653589793*(5.0/2.0))+1.0)+3.93673619140503E-4)*3.273239544735159)*(atan(theta)*2.233158905420314E-2+sin(theta)*(theta*(5.0/2.0)+(theta*theta)*(5.0/2.0)-sqrt(theta)-(theta*theta*theta)*4.058888755698753+(theta*theta*theta*theta)*1.295846037797472)-1.683771314041424E31/(exp(theta*-1.0E1+3.141592653589793*(5.0/2.0))*1.014120480182584E31+1.014120480182584E31)+6.442940289448513E-4)+cos(t)*pow(cos(theta),2.0)*(sin(t)+1.0)*2.0)*1.0/sqrt(pow(cos(theta),2.0)*pow(sin(t)+1.0,2.0)+pow(sin(theta)*(t*(theta*(5.0/2.0)+(theta*theta)*(5.0/2.0)-sqrt(theta)-(theta*theta*theta)*4.058888755698753+(theta*theta*theta*theta)*1.295846037797472)+1.0)+t*(atan(theta)*1.364494638965387E-2-1.014485301813865/(exp(theta*-1.0E1+3.141592653589793*(5.0/2.0))+1.0)+3.93673619140503E-4)*1.63661977236758,2.0))*(-1.0/2.0)+cos(theta)*(-cos(r*cos(theta))*cos(r*sin(theta))*cos(theta)*(sqrt(pow(cos(theta),2.0)*pow(sin(t)+1.0,2.0)+pow(sin(theta)*(t*(theta*(5.0/2.0)+(theta*theta)*(5.0/2.0)-sqrt(theta)-(theta*theta*theta)*4.058888755698753+(theta*theta*theta*theta)*1.295846037797472)+1.0)+t*(atan(theta)*1.364494638965387E-2-1.014485301813865/(exp(theta*-1.0E1+3.141592653589793*(5.0/2.0))+1.0)+3.93673619140503E-4)*1.63661977236758,2.0))-sqrt((r*r)*pow(cos(theta),2.0)+(r*r)*pow(sin(theta),2.0)))+sin(r*cos(theta))*sin(r*sin(theta))*sin(theta)*(sqrt(pow(cos(theta),2.0)*pow(sin(t)+1.0,2.0)+pow(sin(theta)*(t*(theta*(5.0/2.0)+(theta*theta)*(5.0/2.0)-sqrt(theta)-(theta*theta*theta)*4.058888755698753+(theta*theta*theta*theta)*1.295846037797472)+1.0)+t*(atan(theta)*1.364494638965387E-2-1.014485301813865/(exp(theta*-1.0E1+3.141592653589793*(5.0/2.0))+1.0)+3.93673619140503E-4)*1.63661977236758,2.0))-sqrt((r*r)*pow(cos(theta),2.0)+(r*r)*pow(sin(theta),2.0)))+cos(r*sin(theta))*sin(r*cos(theta))*1.0/sqrt((r*r)*pow(cos(theta),2.0)+(r*r)*pow(sin(theta),2.0))*(r*pow(cos(theta),2.0)*2.0+r*pow(sin(theta),2.0)*2.0)*(1.0/2.0))-sin(theta)*(-cos(r*cos(theta))*cos(r*sin(theta))*sin(theta)*(sqrt(pow(cos(theta),2.0)*pow(sin(t)+1.0,2.0)+pow(sin(theta)*(t*(theta*(5.0/2.0)+(theta*theta)*(5.0/2.0)-sqrt(theta)-(theta*theta*theta)*4.058888755698753+(theta*theta*theta*theta)*1.295846037797472)+1.0)+t*(atan(theta)*1.364494638965387E-2-1.014485301813865/(exp(theta*-1.0E1+3.141592653589793*(5.0/2.0))+1.0)+3.93673619140503E-4)*1.63661977236758,2.0))-sqrt((r*r)*pow(cos(theta),2.0)+(r*r)*pow(sin(theta),2.0)))+sin(r*cos(theta))*sin(r*sin(theta))*cos(theta)*(sqrt(pow(cos(theta),2.0)*pow(sin(t)+1.0,2.0)+pow(sin(theta)*(t*(theta*(5.0/2.0)+(theta*theta)*(5.0/2.0)-sqrt(theta)-(theta*theta*theta)*4.058888755698753+(theta*theta*theta*theta)*1.295846037797472)+1.0)+t*(atan(theta)*1.364494638965387E-2-1.014485301813865/(exp(theta*-1.0E1+3.141592653589793*(5.0/2.0))+1.0)+3.93673619140503E-4)*1.63661977236758,2.0))-sqrt((r*r)*pow(cos(theta),2.0)+(r*r)*pow(sin(theta),2.0)))+cos(r*cos(theta))*sin(r*sin(theta))*1.0/sqrt((r*r)*pow(cos(theta),2.0)+(r*r)*pow(sin(theta),2.0))*(r*pow(cos(theta),2.0)*2.0+r*pow(sin(theta),2.0)*2.0)*(1.0/2.0))+(cos(theta)*(cos(r*cos(theta))*sin(r*sin(theta))*((sin(theta)*(t*(theta*(5.0/2.0)+(theta*theta)*(5.0/2.0)-sqrt(theta)-(theta*theta*theta)*4.058888755698753+(theta*theta*theta*theta)*1.295846037797472)+1.0)*2.0+t*(atan(theta)*1.364494638965387E-2-1.014485301813865/(exp(theta*-1.0E1+3.141592653589793*(5.0/2.0))+1.0)+3.93673619140503E-4)*3.273239544735159)*(cos(theta)*(t*(theta*(5.0/2.0)+(theta*theta)*(5.0/2.0)-sqrt(theta)-(theta*theta*theta)*4.058888755698753+(theta*theta*theta*theta)*1.295846037797472)+1.0)+t*(1.364494638965387E-2/(theta*theta+1.0)-exp(theta*-1.0E1+3.141592653589793*(5.0/2.0))*1.0/pow(exp(theta*-1.0E1+3.141592653589793*(5.0/2.0))+1.0,2.0)*1.014485301813865E1)*1.63661977236758+t*sin(theta)*(theta*5.0-(theta*theta)*1.217666626709626E1-1.0/sqrt(theta)*(1.0/2.0)+(theta*theta*theta)*5.183384151189888+5.0/2.0))-cos(theta)*sin(theta)*pow(sin(t)+1.0,2.0)*2.0)*1.0/sqrt(pow(cos(theta),2.0)*pow(sin(t)+1.0,2.0)+pow(sin(theta)*(t*(theta*(5.0/2.0)+(theta*theta)*(5.0/2.0)-sqrt(theta)-(theta*theta*theta)*4.058888755698753+(theta*theta*theta*theta)*1.295846037797472)+1.0)+t*(atan(theta)*1.364494638965387E-2-1.014485301813865/(exp(theta*-1.0E1+3.141592653589793*(5.0/2.0))+1.0)+3.93673619140503E-4)*1.63661977236758,2.0))*(1.0/2.0)+r*cos(r*cos(theta))*cos(r*sin(theta))*cos(theta)*(sqrt(pow(cos(theta),2.0)*pow(sin(t)+1.0,2.0)+pow(sin(theta)*(t*(theta*(5.0/2.0)+(theta*theta)*(5.0/2.0)-sqrt(theta)-(theta*theta*theta)*4.058888755698753+(theta*theta*theta*theta)*1.295846037797472)+1.0)+t*(atan(theta)*1.364494638965387E-2-1.014485301813865/(exp(theta*-1.0E1+3.141592653589793*(5.0/2.0))+1.0)+3.93673619140503E-4)*1.63661977236758,2.0))-sqrt((r*r)*pow(cos(theta),2.0)+(r*r)*pow(sin(theta),2.0)))+r*sin(r*cos(theta))*sin(r*sin(theta))*sin(theta)*(sqrt(pow(cos(theta),2.0)*pow(sin(t)+1.0,2.0)+pow(sin(theta)*(t*(theta*(5.0/2.0)+(theta*theta)*(5.0/2.0)-sqrt(theta)-(theta*theta*theta)*4.058888755698753+(theta*theta*theta*theta)*1.295846037797472)+1.0)+t*(atan(theta)*1.364494638965387E-2-1.014485301813865/(exp(theta*-1.0E1+3.141592653589793*(5.0/2.0))+1.0)+3.93673619140503E-4)*1.63661977236758,2.0))-sqrt((r*r)*pow(cos(theta),2.0)+(r*r)*pow(sin(theta),2.0)))))/r-(sin(theta)*(cos(r*sin(theta))*sin(r*cos(theta))*((sin(theta)*(t*(theta*(5.0/2.0)+(theta*theta)*(5.0/2.0)-sqrt(theta)-(theta*theta*theta)*4.058888755698753+(theta*theta*theta*theta)*1.295846037797472)+1.0)*2.0+t*(atan(theta)*1.364494638965387E-2-1.014485301813865/(exp(theta*-1.0E1+3.141592653589793*(5.0/2.0))+1.0)+3.93673619140503E-4)*3.273239544735159)*(cos(theta)*(t*(theta*(5.0/2.0)+(theta*theta)*(5.0/2.0)-sqrt(theta)-(theta*theta*theta)*4.058888755698753+(theta*theta*theta*theta)*1.295846037797472)+1.0)+t*(1.364494638965387E-2/(theta*theta+1.0)-exp(theta*-1.0E1+3.141592653589793*(5.0/2.0))*1.0/pow(exp(theta*-1.0E1+3.141592653589793*(5.0/2.0))+1.0,2.0)*1.014485301813865E1)*1.63661977236758+t*sin(theta)*(theta*5.0-(theta*theta)*1.217666626709626E1-1.0/sqrt(theta)*(1.0/2.0)+(theta*theta*theta)*5.183384151189888+5.0/2.0))-cos(theta)*sin(theta)*pow(sin(t)+1.0,2.0)*2.0)*1.0/sqrt(pow(cos(theta),2.0)*pow(sin(t)+1.0,2.0)+pow(sin(theta)*(t*(theta*(5.0/2.0)+(theta*theta)*(5.0/2.0)-sqrt(theta)-(theta*theta*theta)*4.058888755698753+(theta*theta*theta*theta)*1.295846037797472)+1.0)+t*(atan(theta)*1.364494638965387E-2-1.014485301813865/(exp(theta*-1.0E1+3.141592653589793*(5.0/2.0))+1.0)+3.93673619140503E-4)*1.63661977236758,2.0))*(-1.0/2.0)+r*cos(r*cos(theta))*cos(r*sin(theta))*sin(theta)*(sqrt(pow(cos(theta),2.0)*pow(sin(t)+1.0,2.0)+pow(sin(theta)*(t*(theta*(5.0/2.0)+(theta*theta)*(5.0/2.0)-sqrt(theta)-(theta*theta*theta)*4.058888755698753+(theta*theta*theta*theta)*1.295846037797472)+1.0)+t*(atan(theta)*1.364494638965387E-2-1.014485301813865/(exp(theta*-1.0E1+3.141592653589793*(5.0/2.0))+1.0)+3.93673619140503E-4)*1.63661977236758,2.0))-sqrt((r*r)*pow(cos(theta),2.0)+(r*r)*pow(sin(theta),2.0)))+r*sin(r*cos(theta))*sin(r*sin(theta))*cos(theta)*(sqrt(pow(cos(theta),2.0)*pow(sin(t)+1.0,2.0)+pow(sin(theta)*(t*(theta*(5.0/2.0)+(theta*theta)*(5.0/2.0)-sqrt(theta)-(theta*theta*theta)*4.058888755698753+(theta*theta*theta*theta)*1.295846037797472)+1.0)+t*(atan(theta)*1.364494638965387E-2-1.014485301813865/(exp(theta*-1.0E1+3.141592653589793*(5.0/2.0))+1.0)+3.93673619140503E-4)*1.63661977236758,2.0))-sqrt((r*r)*pow(cos(theta),2.0)+(r*r)*pow(sin(theta),2.0)))))/r;
	}
    }

    return source;
}
