/***************************************************************************

Author: nsaini
Created: 2018-04-23

***************************************************************************/


#include "common.h"
#include "DGPFunc.h"
#include "memory.h"
#include "solvers.h"

double minmod(double a, double b)
{
    double result;
    if(a*b > 0.0)
    {
        if (fabs(a) < fabs(b))
        {
            result= a;
        }
        else if (fabs (a) >= fabs(b))
        {
            result= b;
        }
    }
    else if (a*b <= 0.0)
    {
        result= 0.0;
    }
    return result;
}

void modifiedLimiter(double ***scalar)
{
    //------------------------------------------------------------------------//
    //Loop indexes
    int ielem, jelem, icoeff;
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Temporary variables
    double *basis;
    allocator1(&basis, ncoeff);

    double zeta1, zeta2;
    double **scalarR, **scalarL;
    double **scalarT, **scalarB;
    allocator2(&scalarR, xelem, yelem);
    allocator2(&scalarL, xelem, yelem);
    allocator2(&scalarT, xelem, yelem);
    allocator2(&scalarB, xelem, yelem);

    double **mean;
    allocator2(&mean, xelem, yelem);

    double a,b,c;

    double *soln;
    allocator1(&soln, 4);

    double **matrix;
    allocator2(&matrix, 4, 4);
    //------------------------------------------------------------------------//
    //Get the mean values at the cell centers
    for(ielem=1; ielem<xelem-1; ielem++)
    {
	for(jelem=1; jelem<yelem-1; jelem++)
	{
	    //------------------------------------------------------------------------//
	    //Slope limiting in x-direction
	    //Reconstruct solution at the edge centers
	    zeta1 = 1.0;
	    zeta2 = 0.0;
	    basis2D(zeta1, zeta2, basis);
	    for(icoeff=0; icoeff<ncoeff; icoeff++)
	    {
		scalarR[ielem][jelem] += scalar[ielem][jelem][icoeff]*basis[icoeff];
	    }

	    zeta1 = -1.0;
	    zeta2 = 0.0;
	    basis2D(zeta1, zeta2, basis);
	    for(icoeff=0; icoeff<ncoeff; icoeff++)
	    {
		scalarL[ielem][jelem] += scalar[ielem][jelem][icoeff]*basis[icoeff];
	    }

	    mean[ielem][jelem] = 0.5*(scalarR[ielem][jelem] + scalarL[ielem][jelem]);

	    zeta1 = 0.0;
	    zeta2 = 1.0;
	    basis2D(zeta1, zeta2, basis);
	    for(icoeff=0; icoeff<ncoeff; icoeff++)
	    {
		scalarT[ielem][jelem] += scalar[ielem][jelem][icoeff]*basis[icoeff];
	    }

	    zeta1 = 0.0;
	    zeta2 = -1.0;
	    basis2D(zeta1, zeta2, basis);
	    for(icoeff=0; icoeff<ncoeff; icoeff++)
	    {
		scalarB[ielem][jelem] += scalar[ielem][jelem][icoeff]*basis[icoeff];
	    }
	}
    }
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Reconstruct the value at the cell edge centers
    matrix[0][0] = 1.0;
    matrix[0][1] = 1.0;

    matrix[1][1] = 1.0;
    matrix[1][2] = 1.0;

    matrix[2][2] = 1.0;
    matrix[2][3] = 1.0;

    matrix[3][0] = 1.0;
    matrix[3][3] = 1.0;
    for(ielem=2; ielem<xelem-2; ielem++)
    {
	for(jelem=2; jelem<yelem-2; jelem++)
	{
	    //Left face
	    a = mean[ielem][jelem]-scalarL[ielem][jelem];
	    b = alpha*(mean[ielem][jelem] - mean[ielem-1][jelem]);
	    c = alpha*(mean[ielem+1][jelem] - mean[ielem][jelem]);

	    scalarL[ielem][jelem] = mean[ielem][jelem] - minmod(a,minmod(b,c));

	    //Right face
	    a = scalarR[ielem][jelem] - mean[ielem][jelem];
	    scalarR[ielem][jelem] = mean[ielem][jelem] + minmod(a,minmod(b,c));

	    //Bottom face
	    a = mean[ielem][jelem]-scalarB[ielem][jelem];
	    b = alpha*(mean[ielem][jelem] - mean[ielem][jelem-1]);
	    c = alpha*(mean[ielem][jelem+1] - mean[ielem][jelem]);

	    scalarB[ielem][jelem] =  mean[ielem][jelem] - minmod(a,minmod(b,c));

	    //Top face
	    a = scalarT[ielem][jelem] - mean[ielem][jelem];
	    scalarT[ielem][jelem] = mean[ielem][jelem] + minmod(a,minmod(b,c));

	    //Get velocity at vertices
	    
	    
	    
	    
	}
    }
    
    //------------------------------------------------------------------------//
    //Deallocators
    deallocator1(&basis, ncoeff);
    deallocator2(&mean, xelem, yelem);
    deallocator2(&scalarR, xelem, yelem);
    deallocator2(&scalarL, xelem, yelem);

    deallocator1(&soln, 4);
    deallocator2(&matrix, 4, 4);
    //------------------------------------------------------------------------//

}
