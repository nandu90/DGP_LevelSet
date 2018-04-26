/***************************************************************************

Author: nsaini
Created: 2018-04-25

***************************************************************************/


#include "common.h"
#include "commu.h"
#include "INS.h"

void INSlevel_setBC(double **scalar, int **iBC)
{
    int i,j;
    /****left and right walls*****/
    if(x_bound == 1 || x_bound == 2)
    {
        for( j=0; j<yelem; j++)
        {
            if(iBC[1][j] == 2)scalar[1][j] = scalar[2][j];
            if(iBC[xelem-2][j] == 2)scalar[xelem-2][j] = scalar[xelem-3][j];
        }
    }
    else if(x_bound == 3)
    {
        for( j=0; j<yelem; j++)
        {
            if(iBC[1][j] == 2)scalar[1][j] = scalar[xelem-3][j];
            if(iBC[xelem-2][j] == 2)scalar[xelem-2][j] = scalar[2][j];
        }
    }
    
    /****Top and bottom walls***/
    if(y_bound == 1 || y_bound == 2)
    {
        for( i=0; i<xelem; i++)
        {
            if(iBC[i][1] == 2)scalar[i][1] = scalar[i][2];
            if(iBC[i][yelem-2] == 2)scalar[i][yelem-2] = scalar[i][yelem-3];
        }
    }
    else if(y_bound == 3)
    {
        for( i=0; i<xelem; i++)
        {
            if(iBC[i][0] == 2)scalar[i][0] = scalar[i][yelem-2];
            if(iBC[i][yelem-1] == 2)scalar[i][yelem-1] = scalar[i][1];
        }
    }

    INScommu2(scalar);
}
