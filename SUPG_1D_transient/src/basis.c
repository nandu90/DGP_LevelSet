/***************************************************************************

Author: nsaini
Created: 2018-07-24

***************************************************************************/

#include "common.h"
#include "memory.h"
#include "functions.h"
#include "polylib.h"

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

void weightdiff1D(double z, double *w, double *x, int ielem)
{
    //------------------------------------------------------------------------//
    //First find the derivative of the original test function
    double *inv, *jacobian;
    allocator1(&inv, 1);
    allocator1(&jacobian, 1);
    double detJ;
    detJ = mappingJacobianDeterminant(ielem, z, x, inv, jacobian);

    double *bdiff;
    allocator1(&bdiff, ncoeff);

    basisdiff1D(z, bdiff);
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    double Peh = detJ*2.0*Pe;
    double beta = 1.0/tanh(Peh) - 1.0/Peh;
    //Now construct the modified weight function
    //Get the original weight first
    basisdiff1D(z, w);
    //Now add perturbation to it
    int icoeff;
    for(icoeff=0; icoeff<ncoeff; icoeff++)
    {
	w[icoeff] += 0.0*beta*detJ;
    }
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    deallocator1(&inv, 1);
    deallocator1(&jacobian, 1);
    deallocator1(&bdiff, ncoeff);
    //------------------------------------------------------------------------//

}

void weight1D(double z, double *w, double *x, int ielem)
{
    //------------------------------------------------------------------------//
    //First find the derivative of the original test function
    double *inv, *jacobian;
    allocator1(&inv, 1);
    allocator1(&jacobian, 1);
    double detJ;
    detJ = mappingJacobianDeterminant(ielem, z, x, inv, jacobian);

    double *bdiff;
    allocator1(&bdiff, ncoeff);

    basisdiff1D(z, bdiff);
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    double Peh = detJ*2.0*Pe;
    double beta = 1.0/tanh(Peh) - 1.0/Peh;
    //Now construct the modified weight function
    //Get the original weight first
    basis1D(z, w);
    //Now add perturbation to it
    int icoeff;
    for(icoeff=0; icoeff<ncoeff; icoeff++)
    {
	w[icoeff] += beta*detJ*bdiff[icoeff]*inv[0];
    }
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    deallocator1(&inv, 1);
    deallocator1(&jacobian, 1);
    deallocator1(&bdiff, ncoeff);
    //------------------------------------------------------------------------//

}

void basisdiff1D(double z, double *b)
{
    //double barray[4] = {0.0, 1.0, 3.0*z, 0.5*(15.0*z*z - 3.0)};

    double *barray;
    allocator1(&barray, polyorder+1);

    if(polyorder == 1)
    {
	barray[0] = -0.5;
	barray[1] = 0.5;
    }
    else if(polyorder == 2)
    {
	barray[0] = 0.5*(2.0*z - 1.0);
	barray[1] = -2.0*z;
	barray[2] = 0.5*(2.0*z + 1.0);
    }
    
    int i;

    for(i=0; i<polyorder+1; i++)
    {
	b[i] = barray[i];
    }
    deallocator1(&barray, polyorder+1);
}

void basis1D(double z, double *basis)
{
    // Upto P = 3 at the moment
    //double barray[4] = {1.0, z, 0.5*(3.0*z*z - 1.0), 0.5*(5*pow(z,3.0)-3.0*z)};

    double *barray;
    allocator1(&barray, polyorder+1);
    
    if(polyorder == 1)
    {
	barray[0] = 0.5*(1.0-z);
	barray[1] = 0.5*(1.0+z);
    }
    else if(polyorder == 2)
    {
	barray[0] = z*(z-1.0)/2.0;
	barray[1] = (1.0-z)*(1.0+z);
	barray[2] = z*(z+1.0)/2.0;
    }
    
    int i;

    for(i=0; i<polyorder+1; i++)
    {
	basis[i] = barray[i];
    }

    deallocator1(&barray, polyorder+1);
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
