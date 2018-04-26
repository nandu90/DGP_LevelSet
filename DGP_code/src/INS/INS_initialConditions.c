/***************************************************************************

Author: nsaini
Created: 2018-04-25

***************************************************************************/


#include "common.h"
#include "memory.h"
#include "DGPFunc.h"
#include "INS.h"
#include "generalFunc.h"

void INSinitialize(struct elemsclr elem)
{
    //------------------------------------------------------------------------//
    //Loop Indexes
    int ielem,jelem;
    int icoeff;
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Temporary variables
    double recphi;
    double *basis;
    allocator1(&basis, ncoeff);
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //First off assign level set values to phi2 - at centroid
    for(ielem=1; ielem<xelem-1; ielem++)
    {
	for(jelem=1; jelem<yelem-1; jelem++)
	{
	    recphi = 0.0;
	    basis2D(0.0, 0.0, basis);
	    for(icoeff=0; icoeff<ncoeff; icoeff++)
	    {
		recphi += basis[icoeff]*elem.phi[ielem][jelem][icoeff];
	    }
	    elem.phi2[ielem][jelem] = recphi;
	}
    }
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Get the Heavyside function
    double eps=epsilon*max(xlen/(gxelem), ylen/(gyelem));

    
    double **H;
    allocator2(&H, xelem, yelem);

    heavy_func(H, elem.phi2, eps);
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Get the density and viscosity
    find_density_visc(H, elem.rho, elem.mu);
    INSlevel_setBC(elem.rho, elem.iBC);
    INSlevel_setBC(elem.mu, elem.iBC);
    //------------------------------------------------------------------------//
    
    //Deallocators
    deallocator2(&H, xelem, yelem);
    deallocator1(&basis, ncoeff);
    //------------------------------------------------------------------------//

}
