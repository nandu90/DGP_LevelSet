/***************************************************************************

Author: nsaini
Created: 2018-03-08

***************************************************************************/


#include "common.h"
#include "icbc.h"
#include "commu.h"
#include "solvers.h"
#include "rhs.h"
#include "memory.h"
#include "generalFunc.h"
#include "DGPFunc.h"

void level_setBC(double ***scalar, int **iBC, double **x, double **y)
{
    int i,j,k;
    
    double *basis;
    allocator1(&basis, ncoeff);
    int igauss, icoeff;

     //------------------------------------------------------------------------//
    //Define quad points and weights here 
    int extra;
    if(quadtype == 1)
    {
	extra = 1;
    }
    else
    {
	extra = 0;
    }
    double **zeta, **weights;
    int tgauss = pow(polyorder + 1 + extra, 2);

    allocator2(&zeta, tgauss,2);
    allocator2(&weights, tgauss,2);
    
    GaussPoints2D(zeta, weights, quadtype, tgauss); 
    //------------------------------------------------------------------------//

    //Allocate the Vandermonde matrix
    double **vand;
    allocator2(&vand, tgauss, ncoeff);

    //Loop over the quadrature points to fill the Vandermonde Matrix
    for(igauss=0; igauss<tgauss; igauss++)
    {
	//Get the basis vector
	basis2D(zeta[igauss][0], zeta[igauss][1], basis);
	//Fill up row of the Vandermonde matrix
	for(icoeff=0; icoeff<ncoeff; icoeff++)
	{
	    vand[igauss][icoeff] = basis[icoeff];
	}
    }

    //Allocate coordinate matrix corresponding to zs - solution points
    double **xs;
    allocator2(&xs, tgauss, 2);
    double *cartphi;
    allocator1(&cartphi, tgauss);
    //------------------------------------------------------------------------//
    
    
    //------------------------------------------------------------------------//
    //Left and right walls
    if(x_bound == 1 || x_bound == 2)
    {
        for( j=0; j<yelem; j++)
        {
	    if(case_tog != 10)
	    {
		//Loop over the coefficients
		for(k=0; k<ncoeff; k++)
		{
		    if(iBC[1][j] == 2)scalar[0][j][k] = scalar[1][j][k];
		    if(iBC[xelem-2][j] == 2)scalar[xelem-1][j][k] = scalar[xelem-2][j][k];
		}
	    }
	    else
	    {
		if(iBC[1][j] == 2)
		{
		    //Convert natural coordinates at quadrature points to Cartesian
		    naturalToCartesian(xs, x, y, 1, j, zeta, tgauss);
		    //Get the phi values at Cartesian Quadrature points
		    for(igauss=0; igauss<tgauss; igauss++)
		    {
			cartphi[igauss] = getphi(xs[igauss][0], xs[igauss][1],simtime);
		    }
		    //Solve the system to get the source coefficents in natural coordinates
		    solveSystem(vand, cartphi, scalar[1][j], tgauss, ncoeff);
		}
		if(iBC[xelem-2][j] == 2)
		{
		    //Convert natural coordinates at quadrature points to Cartesian
		    naturalToCartesian(xs, x, y, xelem-2, j, zeta, tgauss);
		    //Get the phi values at Cartesian Quadrature points
		    for(igauss=0; igauss<tgauss; igauss++)
		    {
			cartphi[igauss] = getphi(xs[igauss][0], xs[igauss][1],simtime);
		    }
		    //Solve the system to get the source coefficents in natural coordinates
		    solveSystem(vand, cartphi, scalar[xelem-2][j], tgauss, ncoeff);
		}
	    }
	}
    }
    else if(x_bound == 3) //Periodic BC are mostly handled in commu2
    {
	if(bhailog[0] < 0 && per[0] == 1) //Means processor is periodic with itself
	{
	    //printf("%d went here\n",myrank);
	    for( j=0; j<yelem; j++)
	    {
		//Loop over quadrature points
		for(k=0; k<ncoeff; k++)
		{
		    scalar[1][j][k] = scalar[xelem-3][j][k];
		    scalar[xelem-2][j][k] = scalar[2][j][k];
		}
	    }
	}
    }
    else if(x_bound == 4)
    {
	for( j=0; j<yelem; j++)
        {
	    //Loop over the coefficients
	    for(k=0; k<ncoeff; k++)
	    {
		if(iBC[1][j] == 2)scalar[0][j][k] = -scalar[1][j][k];
		if(iBC[xelem-2][j] == 2)scalar[xelem-1][j][k] = -scalar[xelem-2][j][k];
	    }
	}
    }
    
    //Top and bottom walls
    if(y_bound == 1 || y_bound == 2)
    {
        for( i=0; i<xelem; i++)
        {
	    if(case_tog != 10)
	    {
		//Loop over quadrature points
		for(k=0; k<ncoeff; k++)
		{
		    if(iBC[i][1] == 2)scalar[i][0][k] = scalar[i][1][k];
		    if(iBC[i][yelem-2] == 2)scalar[i][yelem-1][k] = scalar[i][yelem-2][k];
		}
	    }
	    else
	    {
		if(iBC[i][1] == 2)
		{
		    //Convert natural coordinates at quadrature points to Cartesian
		    naturalToCartesian(xs, x, y, i, 1, zeta, tgauss);
		    //Get the phi values at Cartesian Quadrature points
		    for(igauss=0; igauss<tgauss; igauss++)
		    {
			cartphi[igauss] = getphi(xs[igauss][0], xs[igauss][1],simtime);
		    }
		    //Solve the system to get the source coefficents in natural coordinates
		    solveSystem(vand, cartphi, scalar[i][1], tgauss, ncoeff);
		}
		if(iBC[i][yelem-2] == 2)
		{
		    //Convert natural coordinates at quadrature points to Cartesian
		    naturalToCartesian(xs, x, y, i, yelem-2, zeta, tgauss);
		    //Get the phi values at Cartesian Quadrature points
		    for(igauss=0; igauss<tgauss; igauss++)
		    {
			cartphi[igauss] = getphi(xs[igauss][0], xs[igauss][1],simtime);
		    }
		    //Solve the system to get the source coefficents in natural coordinates
		    solveSystem(vand, cartphi, scalar[i][yelem-2], tgauss, ncoeff);
		}
	    }
	} 
    }
    else if(y_bound == 3) //Periodic BC are mostly handled in commu2
    {
	if(bhailog[1] < 0 && per[1] == 1) //Means processor is periodic with itself
	{
	    //printf("%d went here2\n",myrank);
	    for( i=0; i<xelem; i++)
	    {
		//Loop over quadrature points
		for(k=0; k<ncoeff; k++)
		{
		    scalar[i][0][k] = scalar[i][yelem-2][k];
		    scalar[i][yelem-1][k] = scalar[i][1][k];
		}
	    }
	}
    }
    else if(y_bound == 4)
    {
        for( i=0; i<xelem; i++)
        {
	    //Loop over quadrature points
	    for(k=0; k<ncoeff; k++)
	    {
		if(iBC[i][1] == 2)scalar[i][0][k] = -scalar[i][1][k];
		if(iBC[i][yelem-2] == 2)scalar[i][yelem-1][k] = -scalar[i][yelem-2][k];
	    }
	} 
    }
    //Finally do the communication
    commu2(scalar);

    //------------------------------------------------------------------------//
    //Deallocators
    deallocator1(&basis, ncoeff);
    deallocator2(&zeta, tgauss,2);
    deallocator2(&weights, tgauss,2);
    deallocator2(&vand, tgauss, ncoeff);
    deallocator2(&xs, tgauss, 2);
    deallocator1(&cartphi, tgauss);
    //------------------------------------------------------------------------//

}



