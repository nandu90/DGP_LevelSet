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
    a=0.475;
    b=0.525;
    c=0.6;
    d=0.85;
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
    double cir = sqrt(pow(x-0.5,2.0) + pow(y-0.75,2.0)) - 0.15;
    if(phi < 0.0)
    {
	phi = max(cir,phi);
    }
    return phi;
}

double zalesak2(double x, double y)
{
    double r = sqrt(pow(x-0.5,2.0) + pow(y-0.75,2.0)) - 0.15;

    double phi = 0.0;

    if(r >= 0.0)
    {
	phi = 0.0;
    }
    else
    {
	phi = 1.0;
    }
    
    //Add the notch
    if(x >= 0.475 && x <= 0.525 && y >= 0.6 && y <= 0.85)
    {
	phi = 0.0;
    }
    
    if(x >= 0.0 && x <= 1.0 && y >= 0.0 && y <= 0.5)
    {
	double term1 = 200.0*pow((x - 0.5),2.0);
	double term2 = 200.0*pow((y - 0.25),2.0);
	phi = exp(-(term1 + term2));
    }
    
    return phi;
}

void GaussianStep(double x, double y, double *ls)
{
    double xmin = 90.0;
    double xmax = 110.0;
    double ymin = 15.0;
    double ymax = 35.0;

    if(x >= xmin && x <= xmax)
    {
	if(y >= ymin && y<=ymax)
	{
	    (*ls) = 1.0;
	}
    }
}

void initializeVel(struct elemsclr elem, double **x, double **y)
{
    int i,j,k,l;

    //int ielem, jelem, icoeff;

    //------------------------------------------------------------------------//
    //Define quad points and weights here 
    int extra = 0;
    
    double **zeta, **weights;
    int tgauss = pow(polyorder + 1 + extra, 2);

    allocator2(&zeta, tgauss,2);
    allocator2(&weights, tgauss,2);
    
    GaussPoints2D(zeta, weights, quadtype, tgauss); 
    //------------------------------------------------------------------------//
    
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

    double *inv,  *jacobian;
    allocator1(&inv, 4);
    allocator1(&jacobian, 4);

    double uxy = 0.0, vxy = 0.0, magxy = 0.0,mag = 0.0;
    
    
    //Allocate coordinate matrix corresponding to zs - solution points
    double **xs;
    allocator2(&xs, tgauss, 2);
    
    for (i=0; i<xelem; i++)
    {
        for (j=0; j<yelem; j++)
        {
	    
	    //Convert natural coordinates at quadrature points to Cartesian
	    naturalToCartesian(xs, x, y, i, j, zeta, tgauss);
	    
	    //Get the vel values at the Cartesian Quadrature points
	    for(k=0; k<tgauss; k++)
	    {
		if(case_tog == 1 || case_tog == 2 || case_tog == 5)
		{
		    uxy = 1.0;
		    vxy = 0.0;
		}
		else if(case_tog == 3 || case_tog == 6)
		{
		    uxy =  PI*(0.5 - xs[k][1])/3.14;
		    vxy =  PI*(xs[k][0] - 0.5)/3.14;		    
		}
		else if(case_tog == 7)
		{
		    uxy = 1.0;
		    vxy = 0.0;
		}
		else if(case_tog == 8)
		{
		    uxy = sin(2.0*PI*xs[k][1])*pow(sin(PI*xs[k][0]),2.0);
		    vxy = -sin(2.0*PI*xs[k][0])*pow(sin(PI*xs[k][1]),2.0);
		    //printf("%d %d %d %.4f %.4f\n",i,j,k,us[k],vs[k]);
		}
	    
		//------------------------------------------------------------------------//
		//Transform the velocity vector onto local coordinate
		//Get the magnitude in xy coordinates
		
		magxy = sqrt(pow(uxy,2.0) + pow(vxy,2.0));
		
		mappingJacobianDeterminant(i, j, zeta[k][0], zeta[k][1], x, y, inv, jacobian);
		
		//Use the inverse Jacobian matrix to transform
		us[k] = uxy * inv[0] + vxy * inv[1];
		vs[k] = uxy * inv[2] + vxy * inv[3];
		
		//Rescale the obtained vector to be equal to original vector
		mag = sqrt(pow(us[k],2.0) + pow(vs[k],2.0));
		if(mag == 0.0)
		{
		    us[k] = 0.0;
		    vs[k] = 0.0;
		}
		else
		{
		    us[k] = us[k]*magxy/mag;
		    vs[k] = vs[k]*magxy/mag;
		}
		//------------------------------------------------------------------------//
	    }
	    
	    //Solve the system to get the coefficients
	    solveSystem(vand, us, elem.u[i][j], tgauss, ncoeff);
	    solveSystem(vand, vs, elem.v[i][j], tgauss, ncoeff);
	}
    }

    
    //------------------------------------------------------------------------//
    //Deallocators
    deallocator2(&zeta, tgauss, 2);
    deallocator2(&weights, tgauss, 2);
    deallocator1(&us, tgauss);
    deallocator1(&vs, tgauss);
    deallocator1(&basis, ncoeff);
    deallocator2(&vand, tgauss, ncoeff);
    deallocator2(&xs, tgauss, 2);
    deallocator1(&inv, 4);
    deallocator1(&jacobian, 4);
    //------------------------------------------------------------------------//

}

void initializeLS(struct elemsclr elem, double **x, double **y)
{
    int i,j,k,l;

    //------------------------------------------------------------------------//
    //Define quad points and weights here 
    int extra = 0;
    double **zeta, **weights;
    int tgauss = pow(polyorder + 1 + extra, 2);

    allocator2(&zeta, tgauss,2);
    allocator2(&weights, tgauss,2);
    
    GaussPoints2D(zeta, weights, quadtype, tgauss); 
    //------------------------------------------------------------------------//
    
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
    printf("vandermonde matrix\n");
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

    //------------------------------------------------------------------------//
    double line = 2.5*2.0*rb_in;
    double temp;
    double term1;
    double term2; 
    //------------------------------------------------------------------------//

    
        
    for (i=0; i<xelem; i++)
    {
        for (j=0; j<yelem; j++)
        {
	    //Convert natural coordinates at quadrature points to Cartesian
	    naturalToCartesian(xs, x, y, i, j, zeta, tgauss);

	    
	    
	    //Get the LS value at the Cartesian Quadrature points
	    for(k=0; k<tgauss; k++)
	    {
		//------------------------------------------------------------------------//
		//Unique initial conditions for different cases
		//Gaussian Wave
		if(case_tog == 1)
		{		    	
			term1 = 200.0*pow((xs[k][0] - xb_in),2.0);
			term2 = 200.0*pow((xs[k][1] - yb_in),2.0);
			ls[k] = exp(-(term1));// + term2));
		}
		else if(case_tog == 2 || case_tog == 4)
		{
		    ls[k] = sqrt(pow(xb_in - xs[k][0],2.0) + pow(yb_in - xs[k][1],2.0)) - rb_in;
		}
		else if(case_tog == 3)
		{
		    //ls[k] = zalesak(xs[k][0], xs[k][1]);
		    ls[k] = zalesak2(xs[k][0], xs[k][1]);
		}
		else if(case_tog == 5)
		{
		    if(xs[k][0] <= 0.6)
		    {
			term1 = 200.0*pow((xs[k][0] - 0.3),2.0);
			term2 = 200.0*pow((xs[k][1] - 0.25),2.0);
			ls[k] = exp(-(term1 + term2));
		    }
		    else if(xs[k][0] >= 0.6 && xs[k][0] <= 0.8 && xs[k][1] >= 0.15 && xs[k][1] <= 0.35)
		    {
			ls[k] = 1.0;
		    }
		    else
		    {
			ls[k] = 0.0;
		    }
		    /*double xmin = 100.0;
		    double xmax = 120.0;
		    double ymin = 15.0;
		    double ymax = 35.0;
		    
		    if(xs[k][0] >= xmin && xs[k][0] <= xmax)
		    {
			if(xs[k][1] >= ymin && xs[k][1]<=ymax)
			{
			    ls[k] = 1.0;
			}
		    }
		    //GaussianStep(xs[k][0], xs[k][1], &ls[k]);*/
		}
		else if(case_tog == 6)
		{
		    //ls[k] = 1.0;
		    ls[k] = sqrt(pow(xs[k][0]-xb_in,2.0) + pow(xs[k][1]-yb_in,2.0)) - rb_in;
		}
		else if(case_tog == 7)
		{
		    ls[k] = sin(xs[k][0]);
		}
		else if(case_tog == 8)
		{
		    ls[k] = sqrt(pow(xs[k][0]-xb_in,2.0) + pow(xs[k][1]-yb_in,2.0)) - rb_in;
		}
		else
		{
		    if(myrank == master)printf("Please enter a valid case number. Exiting...\n");
		}

		//Bubble breakup case
		if(case_tog == 4)
		{
		    temp = line - xs[k][1];
		    if(temp < ls[k])
		    {
			ls[k] = temp;
		    }
		}
		//------------------------------------------------------------------------//

	    }

	    
	    //Solve the system to get the coefficients
	    solveSystem(vand, ls, elem.phi[i][j], tgauss, ncoeff);
	    //exit(1);
	    //------------------------------------------------------------------------//
	    //Reconstruct the solution - check
	    /*if(i == 2 && j == 2)
	    {
		printf("\nReconstructed soln at gauss pts is\n");
		for(k=0; k<tgauss; k++)
		{
		    ls[k] = 0.0;
		    basis2D(zeta[k][0], zeta[k][1], basis);
		    for(l=0; l<pow(polyorder+1,2); l++)
		    {
			ls[k] += basis[l]*elem.phi[i][j][l];
			printf("%.4e\n",elem.phi[i][j][l]);
		    }
		    
		    printf("%.4f %.4f %.4f\n",zeta[k][0], zeta[k][1], ls[k]);
		}
		exit(1);
		}*/
	    
	    
	    //------------------------------------------------------------------------//

	}
    }
    deallocator1(&basis, ncoeff);
    deallocator2(&vand,tgauss, ncoeff);
    deallocator2(&xs,tgauss,2);
    deallocator1(&ls,tgauss);
    deallocator2(&zeta, tgauss, 2);
    deallocator2(&weights, tgauss, 2);
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

