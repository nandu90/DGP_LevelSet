/***************************************************************************

Author: nsaini
Created: 2018-03-08

***************************************************************************/


#include "common.h"
#include "icbc.h"



void level_setBC(double ***scalar, int **iBC)
{
    int i,j,k;

  //------------------------------------------------------------------------//
  //Left and right walls
    if(x_bound == 1 || x_bound == 2)
    {
        for( j=0; j<yelem; j++)
        {
	    //Loop over the coefficients
	    for(k=0; k<ncoeff; k++)
	    {
		if(iBC[1][j] == 2)scalar[1][j][k] = scalar[2][j][k];
		if(iBC[xelem-2][j] == 2)scalar[xelem-2][j][k] = scalar[xelem-3][j][k];
	    }
	}
    }
    else if(x_bound == 3)
    {
        for( j=0; j<yelem; j++)
        {
	    //Loop over quadrature points
	    for(k=0; k<ncoeff; k++)
	    {
		if(iBC[1][j] == 2)scalar[1][j][k] = scalar[xelem-3][j][k];
		if(iBC[xelem-2][j] == 2)scalar[xelem-2][j][k] = scalar[2][j][k];
	    }
	}
    }
    
    //Top and bottom walls
    if(y_bound == 1 || y_bound == 2)
    {
        for( i=0; i<xelem; i++)
        {
	    //Loop over quadrature points
	    for(k=0; k<ncoeff; k++)
	    {
		if(iBC[i][1] == 2)scalar[i][1][k] = scalar[i][2][k];
		if(iBC[i][yelem-2] == 2)scalar[i][yelem-2][k] = scalar[i][yelem-3][k];
	    }
	} 
    }
    else if(y_bound == 3)
    {
        for( i=0; i<xelem; i++)
        {
	    //Loop over quadrature points
	    for(k=0; k<ncoeff; k++)
	    {
		if(iBC[i][0] == 2)scalar[i][0][k] = scalar[i][yelem-2][k];
		if(iBC[i][yelem-1] == 2)scalar[i][yelem-1][k] = scalar[i][1][k];
	    }
        }
    }
}



