/***************************************************************************

Author: nsaini
Created: 2018-09-06

***************************************************************************/

#include "common.h"
#include "memory.h"
#include "functions.h"
#include "polylib.h"

void convection(double **C, double *x)
{
    //------------------------------------------------------------------------//
    int ielem, icoeff, igauss;
    int icoeff1;
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Define temp array to store element contibutions
    double ***conv;
    allocator3(&conv, xelem, ncoeff, ncoeff);
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Define array to store differential of basis functions
    double *bdiff;
    allocator1(&bdiff, ncoeff);
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Define array to store weight functions
    double *w;
    allocator1(&w, ncoeff);
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Other temporary variables
    double *inv, *jacobian;
    allocator1(&inv, 1);
    allocator1(&jacobian, 1);

    double detJ;
    //------------------------------------------------------------------------//

    
    //------------------------------------------------------------------------//
    //Loop over the elements
    for(ielem=0; ielem<xelem; ielem++)
    {
	//Loop over the Gauss Quadrature points
	for(igauss=0; igauss<tgauss; igauss++)
	{
	    //Get the basis vector
	    basisdiff1D(zeta[igauss], bdiff);
	    
	    //Get the weight vector
	    basis1D(zeta[igauss], w);

	    //Get the determinant
	    detJ = mappingJacobianDeterminant(ielem, zeta[igauss], x, inv, jacobian);
	    
	    //Loop over the weight vector
	    for(icoeff=0; icoeff<ncoeff; icoeff++)
	    {
		//Loop over the basis vector
		for(icoeff1=0; icoeff1<ncoeff; icoeff1++)
		{
		    conv[ielem][icoeff][icoeff1] += weights[igauss]*w[icoeff]*bdiff[icoeff1]*inv[0]*detJ;
		}
	    }
	}
    }
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Assemble the global Convection matrix
    for(ielem=0; ielem<xelem; ielem++)
    {
	//Define mapping to global matrix
	int row[2] = {ielem, ielem+1};
	int col[2] = {ielem, ielem+1};

	for(icoeff=0; icoeff<ncoeff; icoeff++)
	{
	    for(icoeff1=0; icoeff1<ncoeff; icoeff1++)
	    {
		C[row[icoeff]][col[icoeff1]] += conv[ielem][icoeff][icoeff1];
	    }
	}
    }
    //------------------------------------------------------------------------//


    //------------------------------------------------------------------------//
    //Deallocators
    deallocator3(&conv, xelem, ncoeff, ncoeff);
    deallocator1(&bdiff, ncoeff);
    deallocator1(&w, ncoeff);
    deallocator1(&inv, 1);
    deallocator1(&jacobian, 1);
    //------------------------------------------------------------------------//


}
