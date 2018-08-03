#include "common.h"
#include "DGPFunc.h"
#include "memory.h"
#include "polylib.h"

void GaussPoints1D(double *zeta, double *weights, int quadtype, int total)
{
    
    if(quadtype == 1) //Gauss-Legendre-Lobatto - Deprecated
    {
	zwgll(zeta,weights,total);
    }
    else if(quadtype == 2) //Gauss-Legendre
    {
        if(total == 1)
	{
	    zeta[0] = 0.0;
	    weights[0] = 1.0;
	}
	else
	{
	    zwgl(zeta,weights,total);
	}

    }
}

void GaussPoints2D(double **zeta, double **weights, int quadtype, int tg)
{
    int i,j,k;

    
    //------------------------------------------------------------------------//
    int xg = sqrt(tg);
    int yg = sqrt(tg);
    
    double *zeta1, *zeta2;
    double *weight1, *weight2;
    allocator1(&zeta1, xg);
    allocator1(&zeta2, yg);
    allocator1(&weight1, xg);
    allocator1(&weight2, yg);
    if(quadtype == 1) //Gauss-Legendre-Lobatto - Deprecated
    {
	zwgll(zeta1,weight1,xg);
	zwgll(zeta2,weight2,yg);
    }
    else if(quadtype == 2) //Gauss-Legendre
    {
        if(xg == 1)
	{
	    zeta1[0] = 0.0;
	    weight1[0] = 1.0;
	}
	else
	{
	    zwgl(zeta1,weight1,xg);
	}
	if(yg == 1)
	{
	    zeta2[0] = 0.0;
	    weight2[0] = 1.0;
	}
	else
	{
	    zwgl(zeta2,weight2,yg);
	}
    }


    
    //Arrange quadrature points in a easy to access array
    
    k=0;
    for(j=0; j<yg; j++)
    {
	for(i=0; i<xg; i++)
	{
	    zeta[k][0] = zeta1[i];
	    zeta[k][1] = zeta2[j];
	    weights[k][0] = weight1[i];
	    weights[k][1] = weight2[j];
	    k++;
	}
    }

    
}


void basisDiff2D(double zdiff, double z, double *basis, int order)
{
    int i,j,k;
    
    double *Bdiff, *B;
    allocator1(&Bdiff, polyorder+1);
    allocator1(&B, polyorder+1);

    basis1D(z, B);
    basisdiff1D(zdiff, Bdiff);

    //Now get the outer product of the above two to get the differential
    k=0;
    for(i=0; i<polyorder+1; i++)
    {
	for(j=0; j<polyorder+1; j++)
	{
	    if(order == 1)
	    { 
		basis[k++] = B[i]*Bdiff[j];
	    }
	    else if(order == 2)
	    {
		basis[k++] = B[j]*Bdiff[i];
	    }
	}
    }

    deallocator1(&Bdiff, polyorder+1);
    deallocator1(&B, polyorder+1);
}

void basisdiff1D(double z, double *b)
{
    //------------------------------------------------------------------------//
    //Legendre Polynomials
    //Upto p=3 at the moment
    double barray[4] = {0.0, 1.0, 3.0*z, 0.5*(15.0*z*z - 3.0)};
    int i;
    if(basistype == 1)
    {
	for(i=0; i<polyorder+1; i++)
	{
	    b[i] = barray[i];
	}
    }
    //------------------------------------------------------------------------//
    //------------------------------------------------------------------------//
    //Lagrange polynomials
    else if(basistype == 2)
    {
	if(polyorder == 1)
	{
	    b[0] = -0.5;
	    b[1] = 0.5;
	}
	else if(polyorder == 2)
	{
	    b[0] = 0.5*(2.0*z - 1.0);
	    b[1] = -2.0*z;
	    b[2] = 0.5*(2.0*z + 1.0);
	}
    }
}
void basis1D(double z, double *basis)
{
    

    //------------------------------------------------------------------------//
    //Legendre Polynomials
    // Upto P = 3 at the moment
    double barray[4] = {1.0, z, 0.5*(3.0*z*z - 1.0), 0.5*(5*pow(z,3.0)-3.0*z)};
    int i;
    if(basistype == 1)
    {
        for(i=0; i<polyorder+1; i++)
	{
	    basis[i] = barray[i];
	}
    }
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Lagrange Polynomials
    // Upto p=2 at the moment
    else if(basistype == 2)
    {
	if(polyorder == 1)
	{
	    basis[0] = 0.5*(1.0-z);
	    basis[1] = 0.5*(1.0+z);
	}
	else if(polyorder == 2)
	{
	    basis[0] = z*(z-1.0)/2.0;
	    basis[1] = (1.0-z)*(1.0+z);
	    basis[2] = z*(z+1.0)/2.0;
	}
    }
    //------------------------------------------------------------------------//


    
}

void basis2D(double z1, double z2, double *basis)
{
    int i,j,k;

    double *B1, *B2;
    allocator1 (&B1, polyorder+1);
    allocator1 (&B2, polyorder+1);

    basis1D(z1, B1);
    basis1D(z2, B2);

    //Now get the outer product of above two vectors to get the 2D basis
    
    k=0;
    for(i=0; i<polyorder+1; i++)
    {
	for(j=0; j<polyorder+1; j++)
	{
	    basis[k++] = B2[i]*B1[j];
	}
    }
    
    deallocator1(&B1, polyorder+1);
    deallocator1(&B2, polyorder+1);
}


void mappingFunc(double *scalar, double sc_vertex[], int tgauss, double **zeta, double **weights)
{
    int i;

    double N1, N2, N3, N4;
    //------------------------------------------------------------------------//
    //Loop over quadrature points
    for(i=0; i<tgauss; i++)
    {
	N1 = (1.0-zeta[i][0])*(1.0-zeta[i][1])/4.0;
	N2 = (1.0+zeta[i][0])*(1.0-zeta[i][1])/4.0;
	N3 = (1.0-zeta[i][0])*(1.0+zeta[i][1])/4.0;
	N4 = (1.0+zeta[i][0])*(1.0+zeta[i][1])/4.0;
	

	scalar[i] = N1*sc_vertex[0] + N2*sc_vertex[1] + N3*sc_vertex[2] + N4*sc_vertex[3];
    }
	    
}
