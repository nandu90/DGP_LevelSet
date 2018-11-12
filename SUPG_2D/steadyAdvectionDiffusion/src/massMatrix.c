/***************************************************************************

Author: nsaini
Created: 2018-11-11

***************************************************************************/


#include "common.h"
#include "memory.h"
#include "SUPGFunc.h"
#include "supgSolver.h"

void massMatrix(double ****M, double **x, double **y)
{
    //------------------------------------------------------------------------//
    //Temporary variables
    int ielem, jelem;
    int igauss;
    int icoeff;
    int icoeff1;
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Array to store basis function
    double *b;
    allocator1(&b, supgcoeff);
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Array to store test functions
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

    double **mass;
    allocator2(&mass, supgcoeff, supgcoeff);
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Loop over elements
    for(ielem=2; ielem<xelem-2; ielem++)
    {
	for(jelem=2; jelem<yelem-2; jelem++)
	{
	    //Loop over the Gauss Quadrature points
	    for(igauss=0; igauss<tgauss; igauss++)
	    {
		//Get the Basis vector
		basis2D(zeta[igauss][0], zeta[igauss][1], b);

		//Get the test functions - HAS TO BE CHANGED LATER
		basis2D(zeta[igauss][0], zeta[igauss][1], w);

		//Get the Jacobian inverse and determinant
		detJ = mappingJacobianDeterminant(ielem, jelem, zeta[igauss][0], zeta[igauss][1], x, y, inv, jacobian);

		//Loop over the basis matrix
		for(icoeff=0; icoeff<supgcoeff; icoeff++)
		{
		    for(icoeff1=0; icoeff1<supgcoeff; icoeff1++)
		    {
			M[ielem][jelem][icoeff][icoeff1] += weights[igauss][0]*weights[igauss][1]*detJ*w[icoeff]*b[icoeff1];
		    }
		}
	    }
	}
    }
    //------------------------------------------------------------------------//
    
    //------------------------------------------------------------------------//
    //Deallocators
    deallocator1(&b, supgcoeff);
    deallocator1(&w, supgcoeff);
    deallocator2(&zeta, tgauss,2);
    deallocator2(&weights, tgauss,2);
    deallocator1(&inv, 4);
    deallocator1(&jacobian, 4);
    //------------------------------------------------------------------------//

}
