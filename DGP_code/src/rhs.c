/***************************************************************************

Author: nsaini
Created: 2018-03-29

***************************************************************************/

#include "common.h"
#include "rhs.h"
#include "DGPFunc.h"
#include "memory.h"
#include "commu.h"

void getRHS(struct elemsclr elem, double **x, double **y, double ***rhs, double ****area)
{
    //------------------------------------------------------------------------//
    //Get the contribution from domain integral
    domainIntegral(x, y, elem, rhs);
    
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Get the fluxes on the right and top face of each element
    double ***rflux, ***tflux;
    allocator3(&rflux, xelem, yelem, ygpts);
    allocator3(&tflux, xelem, yelem, xgpts);
    
    fluxes(rflux, tflux, x, y, elem);
    
    //Get the contribution from boundary integral
    boundaryIntegral(rhs, rflux, tflux, x, y, area);
    //------------------------------------------------------------------------//


    //------------------------------------------------------------------------//
    //Communciate the RHS info
    commu2(rhs);
    //------------------------------------------------------------------------//

    /*printf("RHS \n");
    for(icoeff=0; icoeff<ncoeff; icoeff++)
    {
	printf("%.4e ",rhs[2][2][icoeff]);
    }
    printf("\n");
    //exit(1);*/
    
    //------------------------------------------------------------------------//
    //Deallocate
    deallocator3(&rflux, xelem, yelem, ygpts);
    deallocator3(&tflux, xelem, yelem, xgpts);
    //------------------------------------------------------------------------//

}
