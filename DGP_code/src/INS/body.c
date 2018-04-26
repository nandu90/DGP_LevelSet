/***************************************************************************

Author: nsaini
Created: 2018-04-25

***************************************************************************/


#include "common.h"
#include "memory.h"
#include "INS.h"
#include "commu.h"
#include "generalFunc.h"

void body(struct elemsclr sclr, double **st_forcex, double **st_forcey, double **vol)
{
  int i,j;
     /**Compute eps based on grid size*/
     double eps=epsilon*max(xlen/(gxelem), ylen/(gyelem));
    
    /**Compute Heavyside function**/
    double **H;
    allocator2(&H,xelem,yelem);
    
    heavy_func(H, sclr.phi2, eps);
    
    for(i=0; i<xelem; i++)
    {
        for(j=0; j<yelem; j++)
        {
            st_forcex[i][j] += (min(sclr.phi2[i][j],0.0)/sclr.phi2[i][j])*sclr.rho[i][j]*vol[i][j]*gx;
            st_forcey[i][j] += (min(sclr.phi2[i][j],0.0)/sclr.phi2[i][j])*sclr.rho[i][j]*vol[i][j]*gy;
            
        }
    }

    INScommu2(st_forcex);
    INScommu2(st_forcey);
    deallocator2(&H,xelem,yelem);
}
