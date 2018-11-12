#include "common.h"
#include "SUPGFunc.h"
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
    allocator1(&Bdiff, supgorder+1);
    allocator1(&B, supgorder+1);

    basis1D(z, B);
    basisdiff1D(zdiff, Bdiff);

    //Now get the outer product of the above two to get the differential
    k=0;
    for(i=0; i<supgorder+1; i++)
    {
	for(j=0; j<supgorder+1; j++)
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

    deallocator1(&Bdiff, supgorder+1);
    deallocator1(&B, supgorder+1);
}

void basisdiff1D(double z, double *b)
{
    
    //------------------------------------------------------------------------//
    //Lagrange polynomials    
    if(supgorder == 1)
    {
	b[0] = -0.5;
	b[1] = 0.5;
    }
    else if(supgorder == 2)
    {
	b[0] = 0.5*(2.0*z - 1.0);
	b[1] = -2.0*z;
	b[2] = 0.5*(2.0*z + 1.0);
    }
}
void basis1D(double z, double *basis)
{
    //------------------------------------------------------------------------//
    //Lagrange Polynomials
    // Upto p=2 at the moment
   
    if(supgorder == 1)
    {
	basis[0] = 0.5*(1.0-z);
	basis[1] = 0.5*(1.0+z);
    }
    else if(supgorder == 2)
    {
	basis[0] = z*(z-1.0)/2.0;
	basis[1] = (1.0-z)*(1.0+z);
	basis[2] = z*(z+1.0)/2.0;
    }
    
    //------------------------------------------------------------------------//


    
}

void basis2D(double z1, double z2, double *basis)
{
    int i,j,k;

    double *B1, *B2;
    allocator1 (&B1, supgorder+1);
    allocator1 (&B2, supgorder+1);

    basis1D(z1, B1);
    basis1D(z2, B2);

    //Now get the outer product of above two vectors to get the 2D basis
    
    k=0;
    for(i=0; i<supgorder+1; i++)
    {
	for(j=0; j<supgorder+1; j++)
	{
	    basis[k++] = B2[i]*B1[j];
	}
    }
    
    deallocator1(&B1, supgorder+1);
    deallocator1(&B2, supgorder+1);
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


double mappingJacobianDeterminant(int i, int j, double zeta1, double zeta2, double **x, double **y, double *inv, double *jacobian)
{
    double x1, x2, x3, x4;
    double y1, y2, y3, y4;
    double dxdz1, dxdz2;
    double dydz1, dydz2;
    double detJ;
    
    x1 = x[i][j];
    x2 = x[i+1][j];
    x3 = x[i][j+1];
    x4 = x[i+1][j+1];
    
    y1 = y[i][j];
    y2 = y[i+1][j];
    y3 = y[i][j+1];
    y4 = y[i+1][j+1];

    dxdz1 = x1*(zeta2-1.0)*(1.0/4.0)-x2*(zeta2-1.0)*(1.0/4.0)-x3*(zeta2+1.0)*(1.0/4.0)+x4*(zeta2+1.0)*(1.0/4.0);    
    dxdz2 = x1*(zeta1-1.0)*(1.0/4.0)-x2*(zeta1+1.0)*(1.0/4.0)-x3*(zeta1-1.0)*(1.0/4.0)+x4*(zeta1+1.0)*(1.0/4.0);
    dydz1 = y1*(zeta2-1.0)*(1.0/4.0)-y2*(zeta2-1.0)*(1.0/4.0)-y3*(zeta2+1.0)*(1.0/4.0)+y4*(zeta2+1.0)*(1.0/4.0);
    dydz2 = y1*(zeta1-1.0)*(1.0/4.0)-y2*(zeta1+1.0)*(1.0/4.0)-y3*(zeta1-1.0)*(1.0/4.0)+y4*(zeta1+1.0)*(1.0/4.0);

    jacobian[0] = dxdz1;
    jacobian[1] = dxdz2;
    jacobian[2] = dydz1;
    jacobian[3] = dydz2;
    
    detJ = fabs(dxdz1*dydz2-dxdz2*dydz1);

    double dz1dx, dz1dy;
    double dz2dx, dz2dy;

    dz1dx = dydz2/detJ;
    dz1dy = -dxdz2/detJ;
    dz2dx = -dydz1/detJ;
    dz2dy = dxdz1/detJ;

    inv[0] = dz1dx;
    inv[1] = dz1dy;
    inv[2] = dz2dx;
    inv[3] = dz2dy;

    
    return detJ;

    
    
}
