/***************************************************************************

Author: nsaini
Created: 2018-04-25

***************************************************************************/


#include "common.h"
#include "INS.h"
#include "memory.h"

void variable_pressure(double **ustar, double **vstar, double **p, double deltat, double **rho, double **stx, double **sty, double ****area, double **vol, int **iBC)
{
  int i,j;
 
    /****Calculate the RHS of matrix****/
  double **b;
  allocator2(&b, xelem, yelem);
    
    for(i=2; i<xelem-2; i++)
    {
        for(j=2; j<yelem-2; j++)
        {
            double hx = area[i][j][1][1];
            double hy = area[i][j][0][0];
            double volume = vol[i][j];
            b[i][j] = (volume/deltat)*((ustar[i][j] - ustar[i-1][j])/hx + (vstar[i][j]-vstar[i][j-1])/hy);            
        }
    }
    
    /**Now subtract contribution due to surface tension force***/
    /**First calculate the force at the faces***/
    double stx_face[xelem][yelem];
    double sty_face[xelem][yelem];
    for(i=1; i<xelem-2; i++)
    {
        for(j=1; j<yelem-2; j++)
        {
            stx_face[i][j] = 0.5*(stx[i][j] + stx[i+1][j]);
            sty_face[i][j] = 0.5*(sty[i][j] + sty[i][j+1]);
        }
    }
    /**Now calculate and add the contribution from force**/
    for(i=2; i<xelem-2; i++)
    {
        for(j=2; j<yelem-2; j++)
        {
            double volume = vol[i][j];
            b[i][j] -= (1.0/rho[i][j]) * volume*((stx_face[i][j]-stx_face[i-1][j])/area[i][j][1][1] + (sty_face[i][j]-sty_face[i][j-1])/area[i][j][0][0]);
        }
    }
    
    /*if(myrank == master)
      {
	for(i=2; i<xelem-2; i++)
	  {
	    for(j=2; j<yelem-2; j++)
	      {
		printf("%d %d %.6f\n",i-1,j-1,b[i][j]);
	      }
	  }
	  }*/

    /*Create a matrix to store the rho values at all 4 faces of CV*/
    double elem_rho[xelem][yelem][4];
    for(i=2; i < xelem-2; i++)
    {
        for(j=2; j< yelem-2; j++)
        {
            elem_rho[i][j][0] = 0.5*(rho[i+1][j] + rho[i][j]);
            elem_rho[i][j][1] = 0.5*(rho[i][j+1] + rho[i][j]);
            elem_rho[i][j][2] = 0.5*(rho[i][j] + rho[i-1][j]);
            elem_rho[i][j][3] = 0.5*(rho[i][j] + rho[i][j-1]);
        }
    }

    /*if(myrank == master)
      {
	for(i=2; i<xelem-2; i++)
	  {
	    for(j=2; j<yelem-2; j++)
	      {
		printf("%d %d %.6f %.6f %.6f %.6f\n",i-1,j-1,elem_rho[i][j][0],elem_rho[i][j][1],elem_rho[i][j][2],elem_rho[i][j][3]);
	      }
	  }
	  }*/

    double ***a;
    allocator3(&a, xelem, yelem, 5);
    
    for(i=2; i<xelem-2; i++)
    {
        for(j=2; j<yelem-2; j++)
        {
            double hx = area[i][j][1][1];
            double hy = area[i][j][0][0];
            double den = hx*hy*elem_rho[i][j][0]*elem_rho[i][j][1]*elem_rho[i][j][2]*elem_rho[i][j][3];
            a[i][j][0] = hy*hy*elem_rho[i][j][0]*elem_rho[i][j][1]*elem_rho[i][j][3]/den;
            a[i][j][1] = hx*hx*elem_rho[i][j][0]*elem_rho[i][j][1]*elem_rho[i][j][2]/den;
            a[i][j][2] = -(hy*hy*(elem_rho[i][j][1]*elem_rho[i][j][2]*elem_rho[i][j][3] + elem_rho[i][j][0]*elem_rho[i][j][1]*elem_rho[i][j][3]) + hx*hx*(elem_rho[i][j][0]*elem_rho[i][j][2]*elem_rho[i][j][3] + elem_rho[i][j][0]*elem_rho[i][j][1]*elem_rho[i][j][2]) )/den;
            a[i][j][3] = hx*hx*elem_rho[i][j][0]*elem_rho[i][j][2]*elem_rho[i][j][3]/den;
            a[i][j][4] = hy*hy*elem_rho[i][j][1]*elem_rho[i][j][2]*elem_rho[i][j][3]/den;
	    //printf("%d %d %.6f %.6f %.6f %.6f %.6f\n",i-1,j-1,a[i][j][0],a[i][j][1],a[i][j][2],a[i][j][3],a[i][j][4]);
	    //if(myrank==master)printf("%d %d %.6f %.6f %.6f %.6f\n",i-1,j-1,elem_rho[i][j][0],elem_rho[i][j][1],elem_rho[i][j][2],elem_rho[i][j][3]);
	}
    }
    
    /*for(i=2; i<xelem-2; i++)
      {
	for(j=2; j<yelem-2; j++)
	  {
	    if(myrank==master)printf("%d %d %.6f %.6f %.6f %.6f %.6f\n",i-1,j-1,a[i][j][0],a[i][j][1],a[i][j][2],a[i][j][3],a[i][j][4]);
	  }
	  }*/
    
    pointJacobi(a,b,p, iBC, vol);
    deallocator3(&a, xelem, yelem, 5);
    deallocator2(&b, xelem, yelem);
}

