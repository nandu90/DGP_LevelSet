/***************************************************************************

Author: nsaini
Created: 2018-11-11

***************************************************************************/


#include "common.h"
#include "memory.h"
#include "SUPGFunc.h"
#include "supgSolver.h"

void convection(double ***rhs, double ***phi, struct elemsclr elem, double **x, double **y)
{
    //------------------------------------------------------------------------//
    //Temporary Variables
    int ielem, jelem;
    int igauss;
    int icoeff;
    int icoeff1;
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Array to store differential of basis functions
    double *bdiff1, *bdiff2;
    allocator1(&bdiff1, supgcoeff);
    allocator1(&bdiff2, supgcoeff);
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Array for weight functions
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
    double term1, term2, term3, term4;
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

		//Loop over the Gauss Quadrature points
		for(igauss=0; igauss<tgauss; igauss++)
		{
		    //Initialize temporary terms
		    term1 = 0.0;
		    term2 = 0.0;
		    term3 = 0.0;
		    term4 = 0.0;

		    //Get the basis differential vector
		    basisDiff2D(zeta[igauss][0], zeta[igauss][1],bdiff1, 1);
		    basisDiff2D(zeta[igauss][1], zeta[igauss][0],bdiff2, 2);

		    //Get the test functions - HAS TO BE CHANGED LATER
		    basis2D(zeta[igauss][0], zeta[igauss][1], w);
		    
		    //Get the Jacobian inverse and determinant
		    detJ = mappingJacobianDeterminant(ielem, jelem, zeta[igauss][0], zeta[igauss][1], x, y, inv, jacobian);

		    //Assemble terms
		    for(icoeff1 = 0; icoeff1<supgcoeff; icoeff1++)
		    {
			term1 += phi[ielem][jelem][icoeff1]*bdiff1[icoeff1];
			term2 += phi[ielem][jelem][icoeff1]*bdiff2[icoeff1];
		    }

		    term3 = inv[0] + inv[1];
		    term4 = inv[2] + inv[3];

		    sum[icoeff] += weights[igauss][0]*weights[igauss][1]*(elem.u[ielem][jelem][icoeff]*term1*term3 + elem.v[ielem][jelem][icoeff]*term2*term4)*detJ*w[icoeff];
		}
	    }

	    //Add to element residual array
	    for(icoeff=0; icoeff<supgcoeff; icoeff++)
	    {
		rhs[ielem][jelem][icoeff] += sum[icoeff];
	    }
	}
    }
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Deallocators
    deallocator1(&bdiff1, supgcoeff);
    deallocator1(&bdiff2, supgcoeff);
    deallocator1(&w, supgcoeff);
    deallocator2(&zeta, tgauss,2);
    deallocator2(&weights, tgauss,2);
    deallocator1(&inv, 4);
    deallocator1(&jacobian, 4);
    deallocator1(&sum, supgcoeff);
    //------------------------------------------------------------------------//

}
