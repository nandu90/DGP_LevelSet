/***************************************************************************

Author: nsaini
Created: 2018-07-24

***************************************************************************/

#include "common.h"
#include "memory.h"
#include "functions.h"
#include "polylib.h"


void massmatrix(double **M, double *x)
{
    //------------------------------------------------------------------------//
    //Temporary variables
    int ielem;
    int igauss;
    int icoeff;
    int icoeff1;
    //---------------------------------s---------------------------------------//

    //------------------------------------------------------------------------//
    //Define temp array to store element contibutions
    double ***mass;
    allocator3(&mass, xelem, ncoeff, ncoeff);
    //------------------------------------------------------------------------//
    
    //------------------------------------------------------------------------//
    //Define temp array for basis
    double *b;
    allocator1(&b, ncoeff);
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Define temp array for weight functions
    double *w;
    allocator1(&w, ncoeff);
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //other temporary variables
    double *inv, *jacobian;
    allocator1(&inv, 1);
    allocator1(&jacobian, 1);

    double detJ;
    //------------------------------------------------------------------------//

    for(ielem=0; ielem<xelem; ielem++)
    {
	//Loop over the Gauss quadrature points
	for(igauss=0; igauss<tgauss; igauss++)
	{
	    //Get the basis vector
	    basis1D(zeta[igauss], b);

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
		    mass[ielem][icoeff][icoeff1] += b[icoeff]*w[icoeff1]*weights[igauss]*detJ;
		}
	    }
	}
    }

    //------------------------------------------------------------------------//
    //Assemble the global stiffness matrix
    for(ielem=0; ielem<xelem; ielem++)
    {
	//Define mapping to global matrix
	int row[2] = {ielem, ielem+1};
	int col[2] = {ielem, ielem+1};

	for(icoeff=0; icoeff<ncoeff; icoeff++)
	{
	    for(icoeff1=0; icoeff1<ncoeff; icoeff1++)
	    {
		M[row[icoeff]][col[icoeff1]] += mass[ielem][icoeff][icoeff1];
	    }
	}
    }
    //------------------------------------------------------------------------//
    
    //------------------------------------------------------------------------//
    //Deallocators
    deallocator3(&mass, xelem, ncoeff, ncoeff);
    deallocator1(&b, ncoeff);
    deallocator1(&w, ncoeff);
    deallocator1(&inv, 1);
    deallocator1(&jacobian, 1);
    //------------------------------------------------------------------------//

}
