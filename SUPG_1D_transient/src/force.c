/***************************************************************************

Author: nsaini
Created: 2018-09-06

***************************************************************************/

#include "common.h"
#include "memory.h"
#include "functions.h"
#include "polylib.h"

void forceVector(double *F, double *x)
{
    //------------------------------------------------------------------------//
    int ielem, icoeff, igauss;
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Define temp array to store element contibutions
    double **force;
    allocator2(&force, xelem, ncoeff);
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Define array to store differential of weight functions
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
	    //Get the weight vector
	    weight1D(zeta[igauss], w, x, ielem);

	    //Get the determinant
	    detJ = mappingJacobianDeterminant(ielem, zeta[igauss], x, inv, jacobian);
	    
	    //Loop over the weight vector
	    for(icoeff=0; icoeff<ncoeff; icoeff++)
	    {
		force[ielem][icoeff] += weights[igauss]*w[icoeff]*detJ;
	    }
	}
    }
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Assemble the global Force vector
    for(ielem=0; ielem<xelem; ielem++)
    {
	//Define mapping to global matrix
	int col[2] = {ielem, ielem+1};

	for(icoeff=0; icoeff<ncoeff; icoeff++)
	{
	    F[col[icoeff]] += force[ielem][icoeff];
	}
    }
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Deallocators
    deallocator2(&force, xelem, ncoeff);
    deallocator1(&w, ncoeff);
    deallocator1(&inv, 1);
    deallocator1(&jacobian, 1);
    //------------------------------------------------------------------------//

}
