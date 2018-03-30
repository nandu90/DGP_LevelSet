/***************************************************************************

Author: nsaini
Created: 2018-03-28

***************************************************************************/


#include "common.h"
#include "DGPFunc.h"
#include "rhs.h"
#include "polylib.h"
#include "memory.h"

void getFluxes(double ***rflux, double ***tflux, struct elemsclr elem)
{
    //------------------------------------------------------------------------//
    //Loop Indexes
    int ielem, jelem;
    int ixgauss, iygauss;
    int icoeff;
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Temporary variables
    double *zx, *zy;
    double *wx, *wy;
    allocator1(&zx, xgpts);
    allocator1(&zy, ygpts);
    allocator1(&wx, xgpts);
    allocator1(&wy, ygpts);
    //Get the quadrature points and weights
    zwgl(zx,wx,xgpts);
    zwgl(zy,wy,ygpts);

    double *basisx, *basisy;
    allocator1(&basisx, ncoeff);
    allocator1(&basisy, ncoeff);

    double recphi, recu, recv;

    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Loop over the elements
    for(ielem = 0; ielem<xelem; ielem++)
    {
	for(jelem = 0; jelem<yelem; jelem++)
	{
	    //Initialize the integral to 0
	    for(icoeff = 0; icoeff < ncoeff; icoeff++)
	    {
		rflux[ielem][jelem][icoeff] = 0.0;
	    }
	    //Loop over the quadrature points on the right face(y-quadrature)
	    for(iygauss=0; iygauss<ygpts; iygauss++)
	    {
		//Get the basis
		basis2D(1.0, zy[iygauss], basisy);

		//Get the flux vector
		//Reconstruct the solution at the quadrature point
		recphi = 0.0;
		recu = 0.0;
		recv = 0.0;
		for(icoeff = 0; icoeff<ncoeff; icoeff++)
		{
		    recphi += basisy[icoeff]*elem.phi[ielem][jelem][icoeff];
		    recu += basisy[icoeff]*elem.u[ielem][jelem][icoeff];
		    recv += basisy[icoeff]*elem.v[ielem][jelem][icoeff];
		}

		//------------------------------------------------------------------------//
		//Now sum to the integral
		//The following is calculated on the assumption that the face is parallel to y-axis.
		//Will have to have unit face normal vectors if operating on a curved grid.
		//Lets not worry about that now 
		for(icoeff=0; icoeff<ncoeff; icoeff++)
		{
		    rflux[ielem][jelem][icoeff] = wy[iygauss]*basisy[icoeff]*recphi*(recu+recv);
		}
		//------------------------------------------------------------------------//

	    }
		//Initialize the integral to 0
	    for(icoeff = 0; icoeff < ncoeff; icoeff++)
	    {
		tflux[ielem][jelem][icoeff] = 0.0;
	    }
	    //Loop over the quadrature points on the right face(y-quadrature)
	    for(ixgauss=0; ixgauss<xgpts; ixgauss++)
	    {
		//Get the basis
		basis2D(zx[ixgauss], 1.0, basisx);
		
		//Get the flux vector
		//Reconstruct the solution at the quadrature point
		recphi = 0.0;
		recu = 0.0;
		recv = 0.0;
		for(icoeff = 0; icoeff<ncoeff; icoeff++)
		{
		    recphi += basisx[icoeff]*elem.phi[ielem][jelem][icoeff];
		    recu += basisx[icoeff]*elem.u[ielem][jelem][icoeff];
		    recv += basisx[icoeff]*elem.v[ielem][jelem][icoeff];
		}
		
		//------------------------------------------------------------------------//
		//Now sum to the integral
		//The following is calculated on the assumption that the face is parallel to y-axis.
		//Will have to have unit face normal vectors if operating on a curved grid.
		//Lets not worry about that now 
		for(icoeff=0; icoeff<ncoeff; icoeff++)
		{
		    tflux[ielem][jelem][icoeff] = wx[ixgauss]*basisx[icoeff]*recphi*(recu+recv);
		}
		//------------------------------------------------------------------------//
	    }
	}
    }

    //------------------------------------------------------------------------//
    //Deallocate
    deallocator1(&zx, xgpts);
    deallocator1(&zy, ygpts);
    deallocator1(&wx, xgpts);
    deallocator1(&wy, ygpts);
    deallocator1(&basisx, ncoeff);
    deallocator1(&basisy, ncoeff);
    //------------------------------------------------------------------------//
    
    
}
