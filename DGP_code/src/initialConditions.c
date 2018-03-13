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

    
    
    for (i=2; i<xelem-2; i++)
    {
        for (j=2; j<yelem-2; j++)
        {
	    //Get the vel values at the Cartesian Quadrature points
	    for(k=0; k<tgauss; k++)
	    {
		us[k] = 1.0;
		vs[k] = 0.0;
	    }
	    //Solve the system to get the coefficients
	    solveSystem(vand, us, elem.u[i][j]);
	    solveSystem(vand, vs, elem.v[i][j]);
	}
    }
    
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
    
        
    for (i=2; i<xelem-2; i++)
    {
        for (j=2; j<yelem-2; j++)
        {
	    //Convert natural coordinates at quadrature points to Cartesian
	    naturalToCartesian(xs, x, y, i, j);
		
	    //Get the LS value at the Cartesian Quadrature points
	    for(k=0; k<tgauss; k++)
	    {
		ls[k] = sqrt(pow(xb_in - xs[k][0],2.0) + pow(yb_in - xs[k][1],2.0)) - rb_in;
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

