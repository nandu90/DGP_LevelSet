#include "common.h"
#include "polylib.h"
#include "memory.h"
#include "DGPFunc.h"

double lineJacobian(int ielem, int jelem, double zeta, double **x, double **y, int facecode)
{
    //------------------------------------------------------------------------//
    //Temporary variables
    double x1, x2;
    double y1, y2;

    double length;
    double costheta;
    double sintheta;

    double detJ;

    double dxdz;
    double dydz;
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    /*CONVENTION: Choose the coordinates such that if you go from 1-2 
     the right side points outside the cell*/    
    //------------------------------------------------------------------------//
    //Choose the coordinates based on facecode
    if(facecode == 1)
    {
	x1 = x[ielem+1][jelem];
	y1 = y[ielem+1][jelem];

	x2 = x[ielem+1][jelem+1];
	y2 = y[ielem+1][jelem+1];
    }
    else if(facecode == 2)
    {
	x1 = x[ielem][jelem+1];
	y1 = y[ielem][jelem+1];

	x2 = x[ielem+1][jelem+1];
	y2 = y[ielem+1][jelem+1];
    }
    else if(facecode == 3)
    {
	x1 = x[ielem][jelem];
	y1 = y[ielem][jelem];

	x2 = x[ielem][jelem+1];
	y2 = y[ielem][jelem+1];
    }
    else
    {
	x1 = x[ielem][jelem];
	y1 = y[ielem][jelem];

	x2 = x[ielem+1][jelem];
	y2 = y[ielem+1][jelem];
    }
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Get the length
    length = sqrt(pow(y2-y1,2) + pow(x2-x1,2));

    //Get the slope
    costheta = (x2-x1)/length;
    sintheta = (y2-y1)/length;
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Get the components of jacobian matrix
    dxdz = length*costheta/2.0;

    dydz = length*sintheta/2.0;
    //------------------------------------------------------------------------//

    detJ = dxdz + dydz;

    return detJ;
}

double mappingJacobianDeterminant(int i, int j, double zeta1, double zeta2, double **x, double **y, double *inv)
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

    detJ = (dxdz1*dydz2-dxdz2*dydz1);


    inv[0] = dydz2/detJ;
    inv[1] = -dxdz2/detJ;
    inv[2] = -dydz1/detJ;
    inv[3] = dxdz1/detJ;
    
    return detJ;

    
    
}


void massmatrix(double **x, double **y, double ****mass)
{
    //------------------------------------------------------------------------//
    //Temporary variables
    int ielem,jelem;
    int igauss;
    double detJ;
    
    int Bi, Bj;
    
    double *basis;
    allocator1(&basis, ncoeff);

    double *inv;
    allocator1(&inv, 4);

    //------------------------------------------------------------------------//

    
    for(ielem=0; ielem<xelem; ielem++)
    {
	for(jelem=0; jelem<yelem; jelem++)
	{
	    //Loop over the quadrature points
	    for(igauss=0; igauss<tgauss; igauss++)
	    {
		// Get the 2D Basis vector
		basis2D(zeta[igauss][0], zeta[igauss][1], basis);
		// Get the value of determinant
		detJ = mappingJacobianDeterminant(ielem,jelem, zeta[igauss][0], zeta[igauss][1], x, y, inv);
		
		//Loop over the Basis matrix
		for(Bi=0; Bi<ncoeff; Bi++)
		{
		    for(Bj=0; Bj<ncoeff; Bj++)
		    {
			mass[ielem][jelem][Bi][Bj] += basis[Bi]*basis[Bj]*weights[igauss][0]*weights[igauss][1]*detJ;
		    }
		}		
	    }
	}
    }

    //------------------------------------------------------------------------//
    //Check
    /*printf("The mass matrix is\n");
    for(i=0; i<1; i++)
    {
	for(j=0; j<1; j++)
	{
	    //printf("%d %d ",i,j);
	    for(Bi=0; Bi<(int)pow(polyorder+1,2.0); Bi++)
	    {
		for(Bj=0; Bj<(int)pow(polyorder+1,2.0); Bj++)
		{
		    //if(Bi == Bj)
		    // {
			printf("%.6f ",mass[2][2][Bi][Bj]);
			//}
		}
		printf("\n");
	    }
	    printf("\n");
	}
    }

    exit(1);*/
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Deallocators
    deallocator1(&basis, ncoeff);

    deallocator1(&inv, 4);
    //------------------------------------------------------------------------//

}
