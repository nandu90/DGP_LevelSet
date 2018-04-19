/***************************************************************************

Author: nsaini
Created: 2018-03-29

***************************************************************************/

#include "common.h"
#include "rhs.h"
#include "DGPFunc.h"
#include "memory.h"

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
    //Deallocate
    deallocator3(&domIntegral, xelem, yelem, ncoeff);
    deallocator3(&rflux, xelem, yelem, ygpts);
    deallocator3(&tflux, xelem, yelem, xgpts);
    deallocator3(&boundIntegral, xelem, yelem, ncoeff);
    //------------------------------------------------------------------------//

}
