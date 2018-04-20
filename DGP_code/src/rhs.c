/***************************************************************************

Author: nsaini
Created: 2018-03-29

***************************************************************************/

#include "common.h"
#include "rhs.h"
#include "DGPFunc.h"
#include "memory.h"
#include "commu.h"

void getRHS(struct elemsclr elem, double **x, double **y, double ***rhs)
{
    //------------------------------------------------------------------------//
    //Get the domain integral
    double ***domIntegral;
    allocator3(&domIntegral, xelem, yelem, ncoeff);

    domainIntegral(x, y, elem, domIntegral);
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Get the fluxes on the right and top face of each element
    double ***rflux, ***tflux;
    allocator3(&rflux, xelem, yelem, ygpts);
    allocator3(&tflux, xelem, yelem, xgpts);
    
    fluxes(rflux, tflux, elem);

    //Get the boundary Integral
    double ***boundIntegral;
    allocator3(&boundIntegral, xelem, yelem, ncoeff);

    boundaryIntegral(boundIntegral, rflux, tflux);
    //------------------------------------------------------------------------//


    //------------------------------------------------------------------------//
    //Finally construct the RHS vector
    int ielem, jelem, icoeff;
    for(ielem = 0; ielem <xelem; ielem++)
    {
	for(jelem = 0; jelem < yelem; jelem++)
	{
	    for(icoeff = 0; icoeff < ncoeff; icoeff++)
	    {
		rhs[ielem][jelem][icoeff] = domIntegral[ielem][jelem][icoeff] - boundIntegral[ielem][jelem][icoeff];
	    }
	}
    }

    //Communciate the RHS info
    commu2(rhs);
    //------------------------------------------------------------------------//


    
    //------------------------------------------------------------------------//
    //Deallocate
    deallocator3(&domIntegral, xelem, yelem, ncoeff);
    deallocator3(&rflux, xelem, yelem, ygpts);
    deallocator3(&tflux, xelem, yelem, xgpts);
    deallocator3(&boundIntegral, xelem, yelem, ncoeff);
    //------------------------------------------------------------------------//

}
