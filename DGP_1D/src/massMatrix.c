/***************************************************************************

Author: nsaini
Created: 2018-07-24

***************************************************************************/

#include "common.h"
#include "memory.h"
#include "functions.h"

void naturalToCartesian(double *xs, double *x, int ielem)
{
    int igauss;

    double x1 = x[ielem];
    double x2 = x[ielem + 1];

    double N1;
    double N2;
    
    for(igauss=0; igauss < tgauss; igauss++)
    {
	N1 = 0.5*(1.0-zeta[igauss]);
	N2 = 0.5*(1.0+zeta[igauss]);
	xs[igauss] = N1*x1 + N2*x2;
    }
}

void basisdiff1D(double z, double *b)
{
    double barray[4] = {0.0, 1.0, 3.0*z, 0.5*(15.0*z*z - 3.0)};

    int i;

    for(i=0; i<polyorder+1; i++)
    {
	b[i] = barray[i];
    }
}

void basis1D(double z, double *basis)
{
    // Upto P = 3 at the moment
    double barray[4] = {1.0, z, 0.5*(3.0*z*z - 1.0), 0.5*(5*pow(z,3.0)-3.0*z)};
    
    int i;

    for(i=0; i<polyorder+1; i++)
    {
	basis[i] = barray[i];
    }
}

double mappingJacobianDeterminant(int ielem, double z, double *x, double *inv, double *jacobian)
{
    double x1, x2;
    double dxdz;

    double detJ;
    
    x1 = x[ielem];
    x2 = x[ielem+1];

    dxdz = 0.5*(-x1 + x2);

    jacobian[0] = dxdz;

    detJ = fabs(dxdz);

    inv[0] = 1.0/jacobian[0];

    return detJ;
}

void massmatrix(double ***mass, double *x)
{
    //------------------------------------------------------------------------//
    //Temporary variables
    int ielem;
    int igauss;
    int Bi, Bj;

    double *basis;
    allocator1(&basis, ncoeff);

    double *inv, *jacobian;
    allocator1(&inv, 1);
    allocator1(&jacobian, 1);

    double detJ;
    //------------------------------------------------------------------------//

    for(ielem=0; ielem<xelem; ielem++)
    {
	//Loop over the Gauss quadrature points
	for(igauss=0; igauss<tgauss; igauss++)
	{
	    //Get the basis vector
	    basis1D(zeta[igauss], basis);

	    //Get the determinant
	    detJ = mappingJacobianDeterminant(ielem, zeta[igauss], x, inv, jacobian);
	    
	    //Loop over the mass matrix
	    for(Bi = 0; Bi<ncoeff; Bi++)
	    {
		for(Bj = 0; Bj<ncoeff; Bj++)
		{
		    mass[ielem][Bi][Bj] += basis[Bj]*basis[Bi]*weights[igauss]*detJ;
		}
	    }
	}
    }

    for(Bi = 0; Bi<ncoeff; Bi++)
    {
	for(Bj = 0; Bj<ncoeff; Bj++)
	{
	    printf("%.6f ",mass[2][Bi][Bj]);
	}
	printf("\n");
    }
    //------------------------------------------------------------------------//
    //Dellocators
    deallocator1(&basis, ncoeff);
    deallocator1(&inv, 1);
    deallocator1(&jacobian, 1);
    //------------------------------------------------------------------------//

}
