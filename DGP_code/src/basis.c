#include "common.h"
#include "DGPFunc.h"
#include "memory.h"

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


void mappingFunc(double *scalar, double sc_vertex[])
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
