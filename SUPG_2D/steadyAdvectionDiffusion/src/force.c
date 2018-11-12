/***************************************************************************

Author: nsaini
Created: 2018-11-11

***************************************************************************/


#include "common.h"
#include "memory.h"
#include "SUPGFunc.h"
#include "supgSolver.h"

void forceVector(double ***f, double **x, double **y)
{
    //------------------------------------------------------------------------//
    //Temporary variables
    int ielem, jelem;
    int icoeff;
    int igauss;
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Array to store test function
    double *w;
    allocator1(&w, supgcoeff);
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
    int tgauss = pow(supgorder + 1 + extra, 2);

    allocator2(&zeta, tgauss,2);
    allocator2(&weights, tgauss,2);
    
    GaussPoints2D(zeta, weights, quadtype, tgauss); 
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Other temporary variables
    double *inv, *jacobian;
    allocator1(&inv, 4);
    allocator1(&jacobian, 4);

    double detJ;
    double *sum;
    allocator1(&sum, supgcoeff);
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Loop over the elements
    for(ielem=2; ielem<xelem-2; ielem++)
    {
	for(jelem=2; jelem<yelem-2; jelem++)
	{
	    //Loop over the test functions
	    for(icoeff=0; icoeff<supgcoeff; icoeff++)
	    {
		//Initialize sum array
		sum[icoeff] = 0.0;

		//Loop over the Gauss quadrature points
		for(igauss=0; igauss<tgauss; igauss++)
		{
		    //Get the test functions - HAS TO BE CHANGED LATER
		    basis2D(zeta[igauss][0], zeta[igauss][1], w);
		    
		    //Get the Jacobian inverse and determinant
		    detJ = mappingJacobianDeterminant(ielem, jelem, zeta[igauss][0], zeta[igauss][1], x, y, inv, jacobian);

		    //Assemble terms
		    sum[icoeff] += weights[igauss][0]*weights[igauss][1]*w[icoeff]*detJ*g;
		}
	    }

	    //Add to force vector
	    for(icoeff=0; icoeff<supgcoeff; icoeff++)
	    {
		f[ielem][jelem][icoeff] += sum[icoeff];
	    }
	}
    }
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Deallocators
    deallocator1(&w, supgcoeff);
    deallocator2(&zeta, tgauss,2);
    deallocator2(&weights, tgauss,2);
    deallocator1(&inv, 4);
    deallocator1(&jacobian, 4);
    deallocator1(&sum, supgcoeff);
    //------------------------------------------------------------------------//

}
