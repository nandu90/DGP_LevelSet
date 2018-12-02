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
    double *zx, *wx;
    double *zy, *wy;
    int extra;
    if(quadtype == 1)
    {
	extra = 1;
    }
    else
    {
	extra = 0;
    }
    int tx = polyorder + 1 + extra;
    int ty = polyorder + 1 + extra;

    allocator1(&zx, tx);
    allocator1(&zy, ty);
    allocator1(&wx, tx);
    allocator1(&wy, ty);

    GaussPoints1D(zx, wx, quadtype, tx);
    GaussPoints1D(zy, wy, quadtype, ty);

    double ***rflux, ***tflux;
    allocator3(&rflux, xelem, yelem, ty);
    allocator3(&tflux, xelem, yelem, tx);
    
    fluxes(rflux, tflux, x, y, elem, zx, wx, tx, zy, wy, ty);

    //Get the contribution from boundary integral
    boundaryIntegral(rhs, rflux, tflux, x, y, area, zx, wx, tx, zy, wy, ty);
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Get the Contribution from source integral for MMS case
    if(case_tog == 10)
    {
	sourceIntegral(x,y,elem,rhs);
    }
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
    deallocator3(&rflux, xelem, yelem, ty);
    deallocator3(&tflux, xelem, yelem, tx);
    deallocator1(&zx, tx);
    deallocator1(&wx, tx);
    deallocator1(&zy, ty);
    deallocator1(&wy, ty);
    //------------------------------------------------------------------------//

}
