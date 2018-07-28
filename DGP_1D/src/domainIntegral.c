/***************************************************************************

Author: nsaini
Created: 2018-07-25

***************************************************************************/

#include "common.h"
#include "memory.h"
#include "functions.h"

void domainIntegral(double *x, struct elemsclr elem, double **rhs)
{
    //------------------------------------------------------------------------//
    //Loop indexes
    int ielem;
    int igauss;
    int icoeff, icoeff1;
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Temporary variables
    double *dBdz;
    allocator1(&dBdz, ncoeff);

    double recphi;
    double recu;
    double *basis;
    allocator1(&basis, ncoeff);

    double *inv, *jacobian;
    allocator1(&inv, 1);
    allocator1(&jacobian, 1);

    double gradBz;
    double detJ;
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Loop over the elemeents
    for(ielem = 1; ielem<xelem-1; ielem++)
    {
	//Initialize the integral to 0
	for(icoeff=0; icoeff<ncoeff; icoeff++)
	{
	    rhs[ielem][icoeff] = 0.0;
	}

	//Loop over the test functions
	for(icoeff=0; icoeff<ncoeff; icoeff++)
	{
	    //Loop over the Gauss Quadrature points
	    for(igauss=0; igauss<tgauss; igauss++)
	    {
		//get the basis and differntial vectors
		basisdiff1D(zeta[igauss], dBdz);

		basis1D(zeta[igauss], basis);

		recphi = 0.0;
		recu = 0.0;

		for(icoeff1=0; icoeff1<ncoeff; icoeff1++)
		{
		    recphi += basis[icoeff1]*elem.phi[ielem][icoeff1];
		    recu += basis[icoeff1]*elem.u[ielem][icoeff1];
		}

		//Get the jacobian
		detJ = mappingJacobianDeterminant(ielem, zeta[igauss], x, inv, jacobian);

		gradBz = dBdz[icoeff] * inv[0];

		rhs[ielem][icoeff] += weights[igauss]*recphi*gradBz*recu*detJ;
	    }
	}
    }
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Via exact integration
    
    //------------------------------------------------------------------------//

    
    //------------------------------------------------------------------------//
    //Check
    /*printf("\nDomain integral check\n");
    printf("rhs elem.u elem.phi\n");
    ielem = xelem/2;
    for(icoeff=0; icoeff<ncoeff; icoeff++)
    {
	
	printf("%d %.4e %.4e %.4e\n",icoeff+1, rhs[ielem][icoeff], elem.u[ielem][icoeff], elem.phi[ielem][icoeff]);

    }

     printf("\nxs recphi recu\n");
    double *xs;
    allocator1(&xs, tgauss);
    naturalToCartesian(xs, x, ielem);
    for(igauss=0; igauss<tgauss; igauss++)
    {
	basis1D(zeta[igauss], basis);

	recphi = 0.0;
	recu = 0.0;
	
	for(icoeff1=0; icoeff1<ncoeff; icoeff1++)
	{
	    recphi += basis[icoeff1]*elem.phi[ielem][icoeff1];
	    recu += basis[icoeff1]*elem.u[ielem][icoeff1];
	}
	printf("%.4e %.4e %4e\n",xs[igauss], recphi, recu);
	
    }
    printf("\n");

    printf("\nElement coordinates\n");
    printf("%.4e %.4e\n",x[ielem],x[ielem+1]);
    
    exit(1);*/
    //------------------------------------------------------------------------//
    //Deallocators
    deallocator1(&dBdz, ncoeff);
    deallocator1(&basis, ncoeff);
    deallocator1(&inv, 1);
    deallocator1(&jacobian, 1);
    //------------------------------------------------------------------------//

}
