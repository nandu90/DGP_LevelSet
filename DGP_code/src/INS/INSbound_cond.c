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
            if(iBC[1][j] == 2)scalar[0][j] = scalar[1][j];
            if(iBC[xelem-2][j] == 2)scalar[xelem-1][j] = scalar[xelem-2][j];
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
            if(iBC[i][1] == 2)scalar[i][0] = scalar[i][1];
            if(iBC[i][yelem-2] == 2)scalar[i][yelem-1] = scalar[i][yelem-2];
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

void vel_BC(double **u, double **v, int **iBC)
{
  int i,j;
    /***Taking care of left and right direction wall****/
    if(x_bound == 1)
    {
        for( j=0; j<yelem; j++)
        {
            if(iBC[1][j] == 2)u[1][j] = 0.0;
            if(iBC[xelem-2][j] == 2)u[xelem-2][j] = 0.0;
            if(iBC[xelem-2][j] == 2)u[xelem-3][j] = 0.0;
            
            if(iBC[1][j] == 2)v[1][j] = -v[2][j];
            if(iBC[xelem-2][j] == 2)v[xelem-2][j] = -v[xelem-3][j];
        }
    }
    else if(x_bound == 2)
    {
        for( j=0; j<yelem; j++)
        {
            if(iBC[1][j] == 2)u[1][j] = u[2][j];
            if(iBC[xelem-2][j] == 2)u[xelem-3][j] = u[xelem-4][j];
            if(iBC[xelem-2][j] == 2)u[xelem-2][j] = u[xelem-3][j];
            
            if(iBC[1][j] == 2)v[1][j] = v[2][j];
            if(iBC[xelem-2][j] == 2)v[xelem-2][j] = v[xelem-3][j];
        }
    }
    else if(x_bound == 3)
    {
        for( j=0; j<yelem; j++)
        {
            if(iBC[1][j] == 2)u[1][j] = u[xelem-3][j];
            if(iBC[xelem-2][j] == 2)u[xelem-2][j] = u[2][j];
        
            if(iBC[1][j] == 2)v[1][j] = v[xelem-3][j];
            if(iBC[xelem-2][j] == 2)v[xelem-1][j] = v[1][j];
        }
    }
    
    /****Now take care of up and down walls****/
    if(y_bound == 1)
    {
        for( i=0; i<xelem; i++)
        {
            if(iBC[i][1] == 2)u[i][1] = -u[i][2];
            if(iBC[i][yelem-2] == 2)u[i][yelem-2]= -u[i][yelem-3];

            if(iBC[i][1] == 2)v[i][1] = 0.0;
            if(iBC[i][yelem-2] == 2)v[i][yelem-2] = 0.0;
            if(iBC[i][yelem-2] == 2)v[i][yelem-3] = 0.0;
        }
    }
    else if(y_bound == 2)
    {
        for( i=0; i<xelem; i++)
        {
            if(iBC[i][1] == 2)u[i][1] = u[i][2];
            if(iBC[i][yelem-2] == 2)u[i][yelem-2]= u[i][yelem-3];

            if(iBC[i][1] == 2)v[i][1] = v[i][2];
            if(iBC[i][yelem-2] == 2)v[i][yelem-3] = v[i][yelem-4];
            if(iBC[i][yelem-2] == 2)v[i][yelem-2] = v[i][yelem-3];
            
        }
    }
    else if(y_bound == 3)
    {
        for( i=0; i<xelem; i++)
        {
            if(iBC[i][0] == 2)u[i][0] = u[i][yelem-2];
            if(iBC[i][yelem-1] == 2)u[i][yelem-1] = u[i][1];
            
            if(iBC[i][0] == 2)v[i][0] = v[i][yelem-2];
            if(iBC[i][yelem-1] == 2)v[i][yelem-1] = v[i][1];
        }
    }
    
}


void grad_level_setBC(double **scalar, int **iBC)
{
  int i,j;
    /****left and right walls*****/
    if(x_bound == 1 || x_bound == 2)
    {
        for( j=0; j<yelem; j++)
        {
            if(iBC[1][j] == 2)scalar[1][j] = 0.0;
            if(iBC[xelem-2][j] == 2)scalar[xelem-2][j] = 0.0;
        }
    }
    else if(x_bound == 3)
    {
        for( j=0; j<yelem; j++)
        {
            if(iBC[0][j] == 2)scalar[0][j] = scalar[xelem-2][j];
            if(iBC[xelem-1][j] == 2)scalar[xelem-1][j] = scalar[1][j];
        }
    }
    
    /****Top and bottom walls***/
    if(y_bound == 1 || y_bound == 2)
    {
        for( i=0; i<xelem; i++)
        {
            if(iBC[i][1] == 2)scalar[i][1] = 0.0;
            if(iBC[i][yelem-2] == 2)scalar[i][yelem-2] = 0.0;
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
}


void pressureBC(double **scalar, int **iBC)
{
  int i,j;
    //For advection case
    /*for( i=0; i<xelem; i++)
    {
        scalar[i][0][0] = scalar[i][1][0];
        scalar[i][yelem-1][0]= scalar[i][yelem-2][0];
        
    }
    
    
    for( j=0; j<yelem; j++)
    {
        scalar[0][j][0] = (-4800*xc[0][j] + 96)*nu;//scalar[1][j][0] + 4800.0*area[1][j][1][1];
        scalar[xelem-1][j][0] = scalar[xelem-2][j][0] - 4800*nu*area[xelem-2][j][1][1];
        
    }*/
    
    //For Bubble breakup and bubble rise case
    for( i=0; i<xelem; i++)
    {
        if(iBC[i][1] == 2)scalar[i][1] = scalar[i][2];
        if(iBC[i][yelem-2] == 2)scalar[i][yelem-2]= 0.0;
        
    }
    
    
    for( j=0; j<yelem; j++)
    {
        if(iBC[1][j] == 2)scalar[1][j] = scalar[2][j];//scalar[1][j][0] + 4800.0*area[1][j][1][1];
        if(iBC[xelem-2][j] == 2)scalar[xelem-2][j] = scalar[xelem-3][j];
        
    }
    
    /*double line = 2.5*2.0*rb_in;
    for( i=2; i<xelem-2; i++)
    {
        for( j=2; j<yelem-2; j++)
        {
            if(yc[i][j] > line)
            {
                scalar[i][j][0] = 0.0;
            }
        }
	}*/
}
