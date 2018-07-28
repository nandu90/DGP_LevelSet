/***************************************************************************

Author: nsaini
Created: 2018-07-25

***************************************************************************/

#include "common.h"
#include "memory.h"
#include "functions.h"

void getRHS(struct elemsclr elem, double *x, double **rhs)
{
    //------------------------------------------------------------------------//
    //get the contribution from domain integral
    
    domainIntegral(x, elem, rhs);
    
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    double *flux;
    allocator1(&flux, xelem);

    fluxes(flux, x, elem);
    
    //get the contribution from boundary integral
    boundaryIntegral(rhs, flux, x);
    
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Deallocators
    deallocator1(&flux, xelem);
    //------------------------------------------------------------------------//

}
