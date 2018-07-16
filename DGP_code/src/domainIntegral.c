/***************************************************************************

Author: nsaini
Created: 2018-03-26

Notes:
 - Communication should not be required here
 - The domain integrals are calculated for the interior elements only
 - x,y coordinates are not required for this routine and you can choose not to pass them
   Passed over here only to perform some checks

***************************************************************************/


#include "common.h"
#include "memory.h"
#include "polylib.h"
#include "DGPFunc.h"
#include "rhs.h"



void domainIntegral(double **x , double **y, struct elemsclr elem, double ***rhs)
{
    //------------------------------------------------------------------------//
    //Loop Indexes
    int ielem,jelem;
    int igauss;
    int icoeff, icoeff1;
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Temporary variables
    double *dBdz1, *dBdz2;
    allocator1(&dBdz1, ncoeff);
    allocator1(&dBdz2, ncoeff);

    double recphi;
    double recu;
    double recv;
    double *basis;
    allocator1(&basis, ncoeff);

    double *inv,  *jacobian;
    allocator1(&inv, 4);
    allocator1(&jacobian, 4);

    double gradBz1, gradBz2;
    double detJ;
    //------------------------------------------------------------------------//

   
    //------------------------------------------------------------------------//
    //Loop over the elements
    for(ielem = 1; ielem<xelem-1; ielem++)
    {
	for(jelem = 1; jelem<yelem-1; jelem++)
	{
	    //Initialize integral to 0
	    for(icoeff=0; icoeff<ncoeff; icoeff++)
	    {
		rhs[ielem][jelem][icoeff] = 0.0;
	    }
	    
	    //Loop over the gauss quadrature points
	    for(icoeff=0; icoeff<ncoeff; icoeff++)
	    {
		for(igauss=0; igauss<tgauss; igauss++)
		{
		    
		    //------------------------------------------------------------------------//
		    //Get the basis and differential vectors
		    //w.r.t. zeta1
		    basisDiff2D(zeta[igauss][0], zeta[igauss][1],dBdz1, 1);
		    //w.r.t zeta2
		    basisDiff2D(zeta[igauss][1], zeta[igauss][0],dBdz2, 2);
		    
		    //check
		    /*printf("The natural coordinates are: %.4f %.4f\n",zeta[igauss][0],zeta[igauss][1]);
		      for(icoeff= 0; icoeff < ncoeff; icoeff++)
		      {
		      printf("%.4f %.4f\n",dBdz1[icoeff], dBdz2[icoeff]);
		      }
		      printf("\n");*/
		    //------------------------------------------------------------------------//
		    
		    //------------------------------------------------------------------------//
		    //Get the flux vector
		    //Reconstruct the solution at the Quadrature point
		    recphi = 0.0;
		    recu = 0.0;
		    recv = 0.0;
		    basis2D(zeta[igauss][0], zeta[igauss][1], basis);
		    for(icoeff1=0; icoeff1<ncoeff; icoeff1++)
		    {
			recphi += basis[icoeff1]*elem.phi[ielem][jelem][icoeff1];
			recu += basis[icoeff1]*elem.u[ielem][jelem][icoeff1];
			recv += basis[icoeff1]*elem.v[ielem][jelem][icoeff1];
		    }
			 
		    //------------------------------------------------------------------------//
		    
		    //------------------------------------------------------------------------//
		    //get the matrix inverse
		    detJ = mappingJacobianDeterminant(ielem, jelem, zeta[igauss][0], zeta[igauss][1], x, y, inv, jacobian);
		    //------------------------------------------------------------------------//

		    //------------------------------------------------------------------------//
		    //Get the zeta1, zeta2 components of basis
		    gradBz1 = dBdz1[icoeff] * inv[0] + dBdz1[icoeff] * inv[1];
		    gradBz2 = dBdz2[icoeff] * inv[2] + dBdz2[icoeff] * inv[3];

		    rhs[ielem][jelem][icoeff] +=weights[igauss][0]*weights[igauss][1]*recphi*(gradBz1*recu + gradBz2*recv)*detJ;
		    
		    //------------------------------------------------------------------------//
		    
		}
	    }

	    
	    //------------------------------------------------------------------------//
	    //Check the integral
	    /* if(ielem == 2 && jelem == 2)
	    { 
		for(icoeff=0; icoeff<ncoeff; icoeff++)
		{
		    printf("The integral is : %.6f\n",rhs[ielem][jelem][icoeff]);
		    printf("The phi value is : %.4f\n",elem.phi[ielem][jelem][icoeff]);
		    printf("The u value is : %.4f\n",elem.u[ielem][jelem][icoeff]);
		    printf("The v value is : %.4f\n\n",elem.v[ielem][jelem][icoeff]);
		}
		printf("The coordinates are\n:");
		printf("%.4f %.4f\n",x[ielem][jelem], y[ielem][jelem]);
		printf("%.4f %.4f\n",x[ielem+1][jelem], y[ielem+1][jelem]);
		printf("%.4f %.4f\n",x[ielem][jelem+1], y[ielem][jelem+1]);
		printf("%.4f %.4f\n",x[ielem+1][jelem+1], y[ielem+1][jelem+1]);
	    exit(1);
	    }*/
	    //------------------------------------------------------------------------//

	}
    }

    /*double sum = 0.0;
    printf("Domain\n");
    for(ielem = 2; ielem<3; ielem++)
    {
	for(jelem=2; jelem<3; jelem++)
	{
	    printf("%d %d ",ielem,jelem);
	    for(icoeff=0; icoeff<ncoeff; icoeff++)
	    {
		printf("%.4f ",integral[ielem][jelem][icoeff]);
	    }
	    printf("\n");
	}
    }
    exit(1);*/
	
    //------------------------------------------------------------------------//
    //Deallocators
    deallocator1(&dBdz1, ncoeff);
    deallocator1(&dBdz2, ncoeff);
    deallocator1(&basis, ncoeff);
    deallocator1(&inv, 4);
    deallocator1(&jacobian, 4);
    //------------------------------------------------------------------------//

}
