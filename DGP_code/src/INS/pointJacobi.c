/***************************************************************************

Author: nsaini
Created: 2018-04-25

***************************************************************************/

#include "common.h"
#include "memory.h"
#include "INS.h"
#include "commu.h"

void pointJacobi(double ***a, double **b, double **p, int **iBC, double **vol)
{
  //printf("%d reached pressure solver\n",myrank);
  int i,j;
  double **tempp;
  double **delp;
  double **d;
  allocator2(&tempp, xelem, yelem);
  allocator2(&delp, xelem, yelem);
  allocator2(&d, xelem, yelem);

    double ires = 0.0;
    double resnorm=0.0;

    for(j=0; j<yelem; j++)
    {
        for(i=0; i<xelem; i++)
        {
            tempp[i][j] = p[i][j];
        }
    }

    int iter;
    for(iter=0; iter< 10000; iter++)
    {
        for(i=2; i<xelem-2; i++)
        {
            for(j=2; j<yelem-2; j++)
            {
                double res = b[i][j] - (a[i][j][0]*tempp[i-1][j] + a[i][j][1]*tempp[i][j-1] + a[i][j][3]*tempp[i][j+1] + a[i][j][4]*tempp[i+1][j]);
                delp[i][j] = res/a[i][j][2];
            }
        }

        for(i=2; i<xelem-2; i++)
        {
            for(j=2; j<yelem-2; j++)
            {
                tempp[i][j] = delp[i][j];
            }
        }
	INScommu2(tempp);
        pressureBC(tempp, iBC);

	double res1=0.0;
        
	double buf = 0.0;
        for(i=2; i<xelem-2; i++)
        {
            for(j=2; j<yelem-2; j++)
            {
                res1 = b[i][j] - (a[i][j][0]*tempp[i-1][j] + a[i][j][1]*tempp[i][j-1] + a[i][j][3]*tempp[i][j+1] + a[i][j][4]*tempp[i+1][j] + a[i][j][2]*tempp[i][j]);

                buf = buf + pow(res1,2.0)*vol[i][j];
            }
        }

	MPI_Allreduce(&buf,&resnorm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	
	
        if(iter == 0)
        {
            ires = resnorm;
            //if(myrank==master)printf("Initial residual %.6f\n",ires);
            
        }
        else
        {
	  //if(myrank == master)printf("Pressure Step: %d residual: %.6f\n",iter,resnorm/ires);
            if(resnorm / ires < ptol)
            {
	      if(myrank == master)printf("Pressure converged in %d \n",iter);
                break;
            }
        }
    }

    INScommu2(tempp);
    pressureBC(tempp,iBC);

    if(myrank==master)printf("Pressure iterations %d Residual %.2e\n",iter,resnorm/ires);
    for(i=0; i<xelem; i++)
    {
        for(j=0; j<yelem; j++)
        {
            p[i][j] = tempp[i][j];
        }
    }
    deallocator2(&tempp, xelem, yelem);
    deallocator2(&delp, xelem, yelem);
    deallocator2(&d, xelem, yelem);
}
