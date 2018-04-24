/***************************************************************************

Author: nsaini
Created: 2018-03-04

***************************************************************************/

#include "common.h"
#include "icbc.h"
#include "generalFunc.h"
#include "DGPFunc.h"
#include "memory.h"
#include "solvers.h"
#include "commu.h"

double zalesak(double x, double y)
{
    double phi = 0.0;

    int k;
    
    double a,b,c,d;
    a=(0.5-0.03)*100.0;
    b=(0.5+0.03)*100.0;
    c=(0.75-0.15)*100.0;
    d=(0.75+0.15-0.06)*100.0;
    if(x-a == 0.0 || x-b == 0.0 || y-c == 0.0 || y-d == 0.0)
    {
	phi =0.0;
    }
    else
    {
	
	double dist[4];
	dist[0] = sqrt(pow(x-a,2.0) + pow(y-c,2.0));
	dist[1] = sqrt(pow(x-b,2.0) + pow(y-c,2.0));
	dist[2] = sqrt(pow(x-b,2.0) + pow(y-d,2.0));
	dist[3] = sqrt(pow(x-a,2.0) + pow(y-d,2.0));
	
	int min = 0;
	for (k=1; k<4; k++)
	{
	    if(dist[min] > dist[k])
	    {
		min=k;
	    }
	}
	
	double base, height;
	if(min == 0)
	{
	    base = a - x;
	    height = c - y;
	}
	else if(min == 1)
	{
	    base = b - x;
	    height = c - y;
	}
	else if(min == 2)
	{
	    base = b - x;
	    height = d - y;
	}
	else if(min == 3)
	{
	    base = a - x;
	    height = d - y;
	}
	
	double angle;
	angle = atan2(height,base);
	if(min == 0)
	{
	    if(angle < PI/4.0 && angle > -3.0*PI/4.0)
	    {
		phi = base;//(height);
	    }
	    else
	    {
		phi = height;//(base);
	    }
	}
	
	else if(min == 2)
	{
	    if(angle < PI/4.0 && angle > -3.0*PI/4.0)
	    {
		phi = -height;//(base);
	    }
	    else
	    {
		phi = -base;//(height);
	    }
	}
	
	else if(min == 1)
	{
	    if(angle < 3.0*PI/4.0 && angle > -PI/4.0)
	    {
		phi = height;//(base);
	    }
	    else
	    {
		phi = -base;//(height);
	    }
	}
	
	else if (min == 3)
	{
	    if(angle <= 3.0*PI/4.0 && angle >= -PI/4.0)
	    {
		phi = base;//(height);
	    }
	    else
	    {
		phi = -height;//(base);
	    }
	}
    }

    phi = -phi;
    double cir = sqrt(pow(x-50.0,2.0) + pow(y-75.0,2.0)) - 15.0;
    if(phi < 0.0)
    {
	phi = max(cir,phi);
    }
    return phi;
}

void initializeVel(struct elemsclr elem, double **x, double **y)
{
    int i,j,k,l;

    //Allocate solution vector - known soln at Gauss quadrature points
    double *us, *vs;
    allocator1(&us, tgauss);
    allocator1(&vs, tgauss);

    //Allocate basis vector
    double *basis;
    allocator1(&basis, ncoeff);

    //Allocate the Vandermonde matrix
    double **vand;
    allocator2(&vand, tgauss, ncoeff);

    //Loop over the quadrature points to fill the Vandermonde Matrix
    for(k=0; k<tgauss; k++)
    {
	//Get the basis vector
	basis2D(zeta[k][0], zeta[k][1], basis);
	//Fill up row of the Vandermonde matrix
	for(l=0; l<ncoeff; l++)
	{
	    vand[k][l] = basis[l];
	}
    }

    //Allocate coordinate matrix corresponding to zs - solution points
    double **xs;
    allocator2(&xs, tgauss, 2);
    
    for (i=0; i<xelem; i++)
    {
        for (j=0; j<yelem; j++)
        {
	    //Convert natural coordinates at quadrature points to Cartesian
	    naturalToCartesian(xs, x, y, i, j);
	    
	    //Get the vel values at the Cartesian Quadrature points
	    for(k=0; k<tgauss; k++)
	    {
		if(case_tog == 1 || case_tog == 2)
		{
		    us[k] = 1.0;
		    vs[k] = 0.0;
		}
		else if(case_tog == 3)
		{
		    us[k] = PI*(50.0 - xs[k][1])/314.0;
		    vs[k] = PI*(xs[k][0] - 50.0)/314.0;
		}
	    }
	    //Solve the system to get the coefficients
	    solveSystem(vand, us, elem.u[i][j]);
	    solveSystem(vand, vs, elem.v[i][j]);
	}
    }

    //------------------------------------------------------------------------//
    //Deallocators
    deallocator1(&us, tgauss);
    deallocator1(&vs, tgauss);
    deallocator1(&basis, ncoeff);
    deallocator2(&vand, tgauss, ncoeff);
    deallocator2(&xs, tgauss, 2);
    //------------------------------------------------------------------------//

}

void initializeLS(struct elemsclr elem, double **x, double **y)
{
    int i,j,k,l;

    //Allocate solution vector - known soln at Gauss quadrature points
    double *ls;
    allocator1(&ls, tgauss);
    
    //Allocate basis vector
    double *basis;
    allocator1(&basis, ncoeff);

    //Allocate the Vandermonde matrix
    double **vand;
    allocator2(&vand, tgauss, ncoeff);

    //Loop over the quadrature points to fill the Vandermonde Matrix
    for(k=0; k<tgauss; k++)
    {
	//Get the basis vector
	basis2D(zeta[k][0], zeta[k][1], basis);
	//Fill up row of the Vandermonde matrix
	for(l=0; l<ncoeff; l++)
	{
	    vand[k][l] = basis[l];
	}
    }

    //Allocate coordinate matrix corresponding to zs - solution points
    double **xs;
    allocator2(&xs, tgauss, 2);
    
        
    for (i=0; i<xelem; i++)
    {
        for (j=0; j<yelem; j++)
        {
	    //Convert natural coordinates at quadrature points to Cartesian
	    naturalToCartesian(xs, x, y, i, j);
		
	    //Get the LS value at the Cartesian Quadrature points
	    for(k=0; k<tgauss; k++)
	    {
		//------------------------------------------------------------------------//
		//Unique initial conditions for different cases
		//Gaussian Wave
		if(case_tog == 1)
		{
		    	double sigmax = 25.0;
			double sigmay = 25.0;
			double term1 = 0.5*pow((xs[k][0] - xb_in)/sigmax,2.0);
			double term2 = 0.5*pow((xs[k][1] - yb_in)/sigmay,2.0);
			ls[k] = 1.0*exp(-(term1 + term2));
		}
		else if(case_tog == 2)
		{
		    ls[k] = sqrt(pow(xb_in - xs[k][0],2.0) + pow(yb_in - xs[k][1],2.0)) - rb_in;
		}
		else if(case_tog == 3)
		{
		    ls[k] = zalesak(xs[k][0], xs[k][1]);
		}
	       
		//------------------------------------------------------------------------//

	    }
	    //Solve the system to get the coefficients
	    solveSystem(vand, ls, elem.phi[i][j]);

	    /*//------------------------------------------------------------------------//
	    //Reconstruct the solution - check
	    printf("\nReconstructed soln at gauss pts is\n");
	    for(k=0; k<tgauss; k++)
	    {
		ls[k] = 0.0;
		basis2D(zeta[k][0], zeta[k][1], basis);
		for(l=0; l<pow(polyorder+1,2); l++)
		{
		    ls[k] += basis[l]*elem.phi[i][j][l];
		}

		printf("%.4f\n",ls[k]);
	    }
	    exit(1);
	    //------------------------------------------------------------------------//*/

	}
    }
    
    deallocator1(&basis, ncoeff);
    deallocator2(&vand,tgauss, ncoeff);
    deallocator2(&xs,tgauss,2);
    deallocator1(&ls,tgauss);
}

void initialize(struct elemsclr elem, double **x, double **y)
{
    //------------------------------------------------------------------------//
    //Initialize Velocities
    initializeVel(elem, x, y);
    //Apply BC

    //Communicate
    commu2(elem.u);
    commu2(elem.v);
    //------------------------------------------------------------------------//

    
    //------------------------------------------------------------------------//
    //Initialize LevelSet field
    initializeLS(elem, x, y);
    //Apply BC
    level_setBC(elem.phi, elem.iBC);
    //Communicate
    commu2(elem.phi);
    //------------------------------------------------------------------------//

}

