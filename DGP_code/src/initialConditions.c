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

void initializevel(struct elemsclr elem, double **x, double **y)
{
    int i,j,k;

    double xvertex[4];
    double yvertex[4];

    double uvel[4];
    double vvel[4];

    for (i=0; i<xelem; i++)
    {
        for (j=0; j<yelem; j++)
        {
	    //Assign coordinate value to vertices
	    //Get the LS value at all the vertices
	    for(k=0; k<4; k++)
	    {
		uvel[k] = 1.0;
		vvel[k] = 0.0;
	    }
	    //Now map the value at the quadrature points
	}
    }
}

void initializeLS(struct elemsclr elem, double **x, double **y)
{
    int i,j,k,l;

    //Allocate solution vector
    double *ls;
    allocator1(&ls, tgauss);
    

    
    double *basis;
    allocator1(&basis, pow(polyorder+1,2));

    //Depending on the polyorder determine how many points you need to construct the soln
    //This will be same as the number of basis function(which is equal to number of coefficients)
    double **zs;
    allocator2(&zs, tgauss, 2);
    //Populate the coordinates
    //getSolnNaturalCoord(zs);

    //Allocate the Vandermonde matrix
    double **vand;
    allocator2(&vand, tgauss, pow(polyorder+1,2));

    //Loop over the solution (not quadrature) points
    for(k=0; k<tgauss; k++)
    {
	//Get the basis vector
	basis2D(zeta[k][0], zeta[k][1], basis);
	//Fill up row of the Vandermonde matrix
	for(l=0; l<pow(polyorder+1,2); l++)
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
	    //Assign coordinate value to solution points
	    getvertices(xs, x, y, i, j);
	    //------------------------------------------------------------------------//
	    //check
	    printf("The coordinates are:\n");
	    for(k=0; k<tgauss; k++)
	    {
		printf("%.4f %.4f    ", xs[k][0], xs[k][1]);
		printf("%.4f %.4f\n", zeta[k][0], zeta[k][1]);
	    }
	    printf("\n\n");
	    //------------------------------------------------------------------------//

		
	    //Get the LS value at the solution points
	    for(k=0; k<tgauss; k++)
	    {
		ls[k] = sqrt(pow(xb_in - xs[k][0],2.0) + pow(yb_in - xs[k][1],2.0)) - rb_in;
	    }
	    //Solve the system to get the coefficients
	    solveSystem(vand, ls, elem.phi[i][j]);

	    //------------------------------------------------------------------------//
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
	    //------------------------------------------------------------------------//

	}
    }
    
    deallocator1(&basis,pow(polyorder+1,2));
    deallocator2(&zs,tgauss,2);
    deallocator2(&vand,tgauss, pow(polyorder+1,2));
    deallocator2(&xs,tgauss,2);
    deallocator1(&ls,tgauss);
}

void initialize(struct elemsclr elem, double **x, double **y)
{
    //------------------------------------------------------------------------//
    //Initialize Velocities
    //initializevel(elem, x, y);
    //Apply BC

    //------------------------------------------------------------------------//

    
    //------------------------------------------------------------------------//
    //Initialize LevelSet field
    initializeLS(elem, x, y);
    //Apply BC
    level_setBC(elem.phi, elem.iBC);
    //------------------------------------------------------------------------//

}

