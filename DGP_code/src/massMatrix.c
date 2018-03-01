#include "common.h"
#include "polylib.h"
#include "memory.h"
#include "DGPFunc.h"

double mappingJacobianDeterminant(int i, int j, double zeta1, double zeta2)
{
    double x1, x2, x3, x4;
    double y1, y2, y3, y4;
    double dxdz1, dxdz2;
    double dydz1, dydz2;
    double detJ;
    
    x1 = x[i][j];
    x2 = x[i+1][j];
    x3 = x[i+1][j+1];
    x4 = x[i][j+1];
    
    y1 = y[i][j];
    y2 = y[i+1][j];
    y3 = y[i+1][j+1];
    y4 = y[i][j+1];

    dxdz1 = x1*(zeta2-1.0)*(1.0/4.0)-x2*(zeta2-1.0)*(1.0/4.0)+x3*(zeta2+1.0)*(1.0/4.0)-x4*(zeta2+1.0)*(1.0/4.0);
    dxdz2 = x1*(zeta1-1.0)*(1.0/4.0)-x2*(zeta1+1.0)*(1.0/4.0)+x3*(zeta1+1.0)*(1.0/4.0)-x4*(zeta1-1.0)*(1.0/4.0);
    dydz1 = y1*(zeta2-1.0)*(1.0/4.0)-y2*(zeta2-1.0)*(1.0/4.0)+y3*(zeta2+1.0)*(1.0/4.0)-y4*(zeta2+1.0)*(1.0/4.0);
    dydz2 = y1*(zeta1-1.0)*(1.0/4.0)-y2*(zeta1+1.0)*(1.0/4.0)+y3*(zeta1+1.0)*(1.0/4.0)-y4*(zeta1-1.0)*(1.0/4.0);

    detJ = dxdz1*dydz2-dxdz2*dydz1;

    return detJ;
    
}


void massmatrix()
{
    int i,j,k,l;
    double detJ;
    
    int Bi, Bj;
    
    int npz1 = polyorder + 12;
    int npz2 = polyorder + 12;
    
    // Get the Gauss Quadrature zeros and weights
    
    double *z1g, *z2g; //zeta1 and zeta2 Gauss zeros
    double *w1g, *w2g; //Gauss weights
    
    allocator1(&z1g, npz1);
    allocator1(&z2g, npz2);
    allocator1(&w1g, npz1);
    allocator1(&w2g, npz2);
    
    //Get the Gauss-Lobatto-Legendre Quadrature points(Two points at the end + polyorder)
    zwgll(z1g,w1g,npz1);
    zwgll(z2g,w2g,npz2);
    //

    double *basis;
    allocator1(&basis, (int)pow(polyorder+1,2.0));
    
    

   
    
    for(i=0; i<xelem; i++)
    {
	for(j=0; j<yelem; j++)
	{
	    //Loop over the quadrature points
	    for(k=0; k<npz1; k++)
	    {
		for(l=0; l<npz2; l++)
		{
		    // Get the 2D Basis vector
		    basis2D(z1g[k], z2g[l], basis);
		    // Get the value of determinant
		    detJ = mappingJacobianDeterminant(i,j, z1g[k], z2g[l]);
		    
		    //Loop over the Basis matrix
		    for(Bi=0; Bi<(int)pow(polyorder+1,2.0); Bi++)
		    {
			for(Bj=0; Bj<(int)pow(polyorder+1,2.0); Bj++)
			{
			    mass[i][j][Bi][Bj] += basis[Bi]*basis[Bj]*detJ*w1g[k]*w2g[l];
			}
		    }
		    
		}
	    }
	}
    }

    for(Bi=0; Bi<(int)pow(polyorder+1,2.0); Bi++)
    {
	for(Bj=0; Bj<(int)pow(polyorder+1,2.0); Bj++)
	{
	    if(Bi == Bj)printf("%.6f       ",mass[2][2][Bi][Bj]);
	}
	printf("\n");
    }
    
    deallocator1(&basis, (int)pow(polyorder+1,2.0));
    deallocator1(&z1g, npz1);
    deallocator1(&z2g, npz2);
    deallocator1(&w1g, npz1);
    deallocator1(&w2g, npz2);
}
