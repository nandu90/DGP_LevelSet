/***************************************************************************

Author: nsaini
Created: 2018-03-28

***************************************************************************/


#include "common.h"
#include "DGPFunc.h"
#include "rhs.h"
#include "polylib.h"
#include "memory.h"

void fluxes(double ***rflux, double ***tflux, struct elemsclr elem)
{
    //------------------------------------------------------------------------//
    /*This routine will construct the fluxes at the right and top faces.
      Upwinding will also be performed here. 
      Boundary integrals will be calculated in a separate routine.
     */
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Loop indexes
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
    //
    double ***recRflux, ***recTflux;      //Reconstructed fluxes at the right and top faces
    double ***recRu,  ***recTu;            //Reconstructed wall normal velocities at the right and top face

    allocator3(&recRflux, xelem, yelem, ygpts);
    allocator3(&recTflux, xelem, yelem, xgpts);
    allocator3(&recRu, xelem, yelem, ygpts);
    allocator3(&recTu, xelem, yelem, xgpts);
    //------------------------------------------------------------------------//


    //------------------------------------------------------------------------//
    //Reconstruct the face normal velocities and the fluxes at the Gauss Quadrature
    //points on the top and right face of each cell
    //Upwinding will be done later - so that its easier to change the flux scheme if i have to
    for(ielem=0; ielem<xelem; ielem++)
    {
	for(jelem=0; jelem<yelem; jelem++)
	{
	    //Loop over the Gauss Quadrature points on the right face
	    for(iygauss=0; iygauss < ygpts; iygauss++)
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

		recRflux[ielem][jelem][iygauss] = recphi;
		//Will have to have face normals here if the mesh is not cartesian - BE CAREFUL
		recRu[ielem][jelem][iygauss] = recu;		    
	    }

	    //Loop over the Gauss Quadrature points on the top face
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

		recTflux[ielem][jelem][ixgauss] = recphi;
		//Will have to have face normals here if the mesh is not cartesian - BE CAREFUL
		recTu[ielem][jelem][ixgauss] = recv;
	    }
	}
    }
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //The flux scheme is implemented now
    //On right face
    upwind(recRflux, recRu, rflux, xgpts, 1);
    //On top face
    upwind(recTflux, recTu, tflux, ygpts, 2);
    //------------------------------------------------------------------------//

    


    //------------------------------------------------------------------------//
    //Deallocators
    deallocator3(&recRflux, xelem, yelem, ygpts);
    deallocator3(&recTflux, xelem, yelem, xgpts);
    deallocator3(&recRu, xelem, yelem, ygpts);
    deallocator3(&recTu, xelem, yelem, xgpts);
    //------------------------------------------------------------------------//

    
    
    
}

void upwind(double ***recflux, double ***recu, double ***flux, int ngauss, int dircode)
{
    //------------------------------------------------------------------------//
    //Temporary variables
    //Minus values are inside the cell and plus are outside    
    double minusU, plusU;
    double minusFlux, plusFlux;

    double Fminus, Fplus;
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Loop indexes
    int ielem, jelem;
    int igauss;
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Accounting for direction code
    int xinc;
    int yinc;
    if(dircode == 1)
    {
	xinc = 1;
	yinc = 0;
    }
    else if(dircode == 2)
    {
	xinc = 0;
	yinc = 1;
    }
    else
    {
	xinc = 0;
	yinc = 0;
	if(myrank == master)
	{
	    printf("Please apply proper direction code for fluxes.\nExiting...");
	    exit(1);
	}
    }
    
    //------------------------------------------------------------------------//
    //Loop through elements
    for(ielem=0; ielem < xelem-1; ielem++)
    {
	for(jelem=0; jelem < yelem-1; jelem++)
	{
	    //Loop over Gauss Quadrature points
	    for(igauss=0; igauss<ngauss; igauss++)
	    {
		minusU = recu[ielem][jelem][igauss];
		plusU = recu[ielem+xinc][jelem+yinc][igauss];

		minusFlux = recflux[ielem][jelem][igauss];
		plusFlux = recflux[ielem+xinc][jelem+yinc][igauss];

		//Upwinding
		if(minusU >= 0.0)
		{
		    Fminus = minusFlux;
		}
		else
		{
		    Fminus = 0.0;
		}

		if(plusU > 0.0)
		{
		    Fplus = 0.0;
		}
		else
		{
		    Fplus = plusFlux;
		}

		flux[ielem][jelem][igauss] = Fminus + Fplus;
	    }
	}
    }

    
}


void boundaryIntegral(double ***integral, double ***rflux, double ***tflux)
{
    //------------------------------------------------------------------------//
    /*
      This routine will calculate the net domain integral
     */
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Loop indexes
    int ielem,jelem;
    int ixgauss, iygauss;
    int icoeff;
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Temporary Variables
    double ***rintegral, ***tintegral;
    allocator3(&rintegral, xelem, yelem, ncoeff);
    allocator3(&tintegral, xelem, yelem, ncoeff);

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
    //------------------------------------------------------------------------//


    //------------------------------------------------------------------------//
    //Loop over the elements
    for(ielem = 0; ielem<xelem; ielem++)
    {
	for(jelem = 0; jelem<yelem; jelem++)
	{
	    //Loop over the quadrature points on the right face
	    for(iygauss = 0; iygauss<ygpts; iygauss++)
	    {
		//Get the basis
		basis2D(1.0, zy[iygauss], basisy);

		//Sum to the integral
		for(icoeff = 0; icoeff<ncoeff; icoeff++)
		{
		    rintegral[ielem][jelem][icoeff] += wy[iygauss]*basisy[icoeff]*rflux[ielem][jelem][iygauss];
		}
	    }

	    //Loop over the quadrature points in the top face
	    for(ixgauss=0; ixgauss<xgpts; ixgauss++)
	    {
		//Get the basis
		basis2D(zx[ixgauss], 1.0, basisx);

		//Sum to the integral
		for(icoeff=0; icoeff<ncoeff; icoeff++)
		{
		    tintegral[ielem][jelem][icoeff] += wx[ixgauss]*basisx[icoeff]*tflux[ielem][jelem][ixgauss];
		}
	    }
	}
    }
    //------------------------------------------------------------------------//


    //------------------------------------------------------------------------//
    //Finally calculate the total integral
    for(ielem=1; ielem<xelem-1; ielem++)
    {
	for(jelem=1; jelem<yelem-1; jelem++)
	{
	    for(icoeff=0; icoeff<ncoeff; icoeff++)
	    {
		integral[ielem][jelem][icoeff] = rintegral[ielem][jelem][icoeff] - rintegral[ielem-1][jelem][icoeff];

		integral[ielem][jelem][icoeff] += tintegral[ielem][jelem][icoeff] - tintegral[ielem][jelem-1][icoeff];
	    }
	}
    }
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Deallocators
    deallocator3(&rintegral, xelem, yelem, ncoeff);
    deallocator3(&tintegral, xelem, yelem, ncoeff);
    deallocator1(&zx, xgpts);
    deallocator1(&zy, ygpts);
    deallocator1(&wx, xgpts);
    deallocator1(&wy, ygpts);
    deallocator1(&basisx, ncoeff);
    deallocator1(&basisy, ncoeff);
    //------------------------------------------------------------------------//

    
}


