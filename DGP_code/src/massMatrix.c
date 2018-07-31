#include "common.h"
#include "polylib.h"
#include "memory.h"
#include "DGPFunc.h"


void lineNormal(int ielem, int jelem, double **x, double **y, double *normz1, double *normz2, int facecode, double zeta1, double zeta2)
{
    //SUBROUTINE NO LONGER USED
    //------------------------------------------------------------------------//
    //Find the slope of line
    double xa, xb;
    double ya, yb;

    if(facecode == 1)
    {
	xa = x[ielem+1][jelem];
	ya = y[ielem+1][jelem];

	xb = x[ielem+1][jelem+1];
	yb = y[ielem+1][jelem+1];
    }
    else if(facecode == 2)
    {
	xb = x[ielem][jelem+1];
	yb = y[ielem][jelem+1];

	xa = x[ielem+1][jelem+1];
	ya = y[ielem+1][jelem+1];
    }
    else if(facecode == 3)
    {
	xa = x[ielem][jelem];
	ya = y[ielem][jelem];

	xb = x[ielem][jelem+1];
	yb = y[ielem][jelem+1];
    }
    else
    {
	xb = x[ielem][jelem];
	yb = y[ielem][jelem];

	xa = x[ielem+1][jelem];
	ya = y[ielem+1][jelem];
    }

    
    double m = 0.0;
    double a, b;
    if(fabs(xb-xa) > 0.0)
    {
	m = (yb-ya)/(xb-xa);
	a = -m;
	b = 1.0;	
    }
    else
    {
	a = 1.0;
	b = 0.0;
    }

    b = -(xb-xa);
    a = (yb-ya);
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Get the jacobian and its inverse
    double *jacobian, *inv;
    allocator1(&jacobian, 4);
    allocator1(&inv, 4);

    mappingJacobianDeterminant(ielem, jelem, zeta1, zeta2, x, y, inv, jacobian);
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Get the differential of equation of line w.r.t zeta1, zeta2
    double dfdz1, dfdz2;
    
    dfdz1 = b * jacobian[2] + a*jacobian[0];
    dfdz2 = b * jacobian[3] + a*jacobian[1];
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Get the normal components in zeta1, zeta2 space
    (*normz1) = dfdz1*(inv[0] + inv[1]);
    (*normz2) = dfdz2*(inv[2] + inv[3]);

    //(*normz1) = a*inv[0] + b*inv[1];
    //(*normz2) = a*inv[2] + b*inv[3];

    double mag = sqrt(pow((*normz1),2.0) + pow((*normz2),2.0));

    (*normz1) = (*normz1)/mag;
    (*normz2) = (*normz2)/mag;
    
    double xcomp = (*normz1);
    double ycomp = (*normz2);

    mag = sqrt(xcomp*xcomp + ycomp*ycomp);

    //printf("%.4e %.4e %.4e %.4e\n", xa, xb, ya, yb);
    //printf("Here1 %.4e %.4e %.4e\n",m, dfdz1, dfdz2);
    /*printf("facecode is %d\n",facecode);
    printf("Here %.4e %.4e\n",a,b);
    printf("Here %.4e %.4e %.4e\n\n", xcomp, ycomp, mag);*/
    //exit(1);
    
    //------------------------------------------------------------------------//


    
    //------------------------------------------------------------------------//
    //Deallocators
    deallocator1(&jacobian, 4);
    deallocator1(&inv, 4);
    //------------------------------------------------------------------------//
    
    
}

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

    detJ = fabs(dxdz + dydz);


    return detJ;
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

    double *jacobian, *inv;
    allocator1(&jacobian, 4);
    allocator1(&inv, 4);

    int i,j;
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Define quad points and weights here independent of what is in the rest of the code
    double *z1, *w1;
    double *z2, *w2;
    allocator1(&z1, xgpts);
    allocator1(&w1, xgpts);
    allocator1(&z2, ygpts);
    allocator1(&w2, ygpts);
    zwgl(z1, w1, xgpts);
    zwgl(z2, w2, ygpts);

    int ngauss = xgpts*ygpts;
    double **z, **w;
    allocator2(&z, ngauss, 2);
    allocator2(&w, ngauss, 2);
    int k=0;
    for(j=0; j<ygpts; j++)
    {
	for(i=0; i<xgpts; i++)
	{
	    z[k][0] = z1[i];
	    z[k][1] = z2[j];
	    w[k][0] = w1[i];
	    w[k][1] = w2[j];
	    k++;
	}
    }
    //------------------------------------------------------------------------//

    
    for(ielem=0; ielem<xelem; ielem++)
    {
	for(jelem=0; jelem<yelem; jelem++)
	{
	    //Loop over the quadrature points
	    for(igauss=0; igauss<ngauss; igauss++)
	    {
		// Get the 2D Basis vector
		basis2D(z[igauss][0], z[igauss][1], basis);
		// Get the value of determinant
		detJ = mappingJacobianDeterminant(ielem,jelem, z[igauss][0], z[igauss][1], x, y, inv, jacobian);
		
		//Loop over the Basis matrix
		for(Bi=0; Bi<ncoeff; Bi++)
		{
		    for(Bj=0; Bj<ncoeff; Bj++)
		    {
			mass[ielem][jelem][Bi][Bj] += basis[Bi]*basis[Bj]*w[igauss][0]*w[igauss][1]*detJ;
		    }
		}		
	    }
	}
    }

    //------------------------------------------------------------------------//
    //Check
    /*int i,j;
    printf("The mass matrix is\n");
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
    deallocator1(&jacobian, 4);

    deallocator1(&z1, xgpts);
    deallocator1(&w1, xgpts);
    deallocator1(&z2, ygpts);
    deallocator1(&w2, ygpts);
    deallocator2(&z, ngauss, 2);
    deallocator2(&w, ngauss, 2);
    //------------------------------------------------------------------------//

}
