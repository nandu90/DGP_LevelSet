/***************************************************************************

Author: nsaini
Created: 2018-04-25

***************************************************************************/


#include "common.h"
#include "INS.h"

void find_density_visc(double **H, double **rho, double **mu)
{
  int i,j;
    for(i=0; i<xelem; i++)
    {
        for(j=0; j<yelem; j++)
        {
            rho[i][j] = rhog + (rhof -rhog)*H[i][j];
            mu[i][j] = mug + (muf -mug)*H[i][j];
        }
    }
}

void heavy_func(double **H, double **phi, double eps)
{
  int i,j;
    for(i=0; i<xelem; i++)
        {
            for(j=0; j<yelem; j++)
            {
                if(phi[i][j] < -eps)
                {
                    H[i][j] = 0.0;
                }
                else if(phi[i][j] > eps)
                {
                    H[i][j] = 1.0;
                }
                else
                {
                    //H[i][j][0] = 1.0/(1.0 + exp(-3.0*sclr.phi[i][j][0]/eps));
                    H[i][j] = 0.5 + phi[i][j]/(2.0*eps) + (1.0/(2.0*PI)*sin(PI * phi[i][j]/eps));
                }
                //<<H[i][j][0]<<" ";
            }
            //<<endl;
        }
        //exit(0);
}

/*void delta_func(double ***delta, double ***phi, double eps)
{
  int i,j;
    for(i=2; i<xelem-2; i++)
    {
        for(j=2; j<yelem-2; j++)
        {
            if(fabs(phi[i][j][0]) > eps)
            {
                delta[i][j][0] = 0.0;
            }
            else
            {
                //delta[i][j][0] = 3.0*exp(-3.0*sclr.phi[i][j][0]/eps)/(1.0 * pow(1.0 + exp(-3.0*sclr.phi[i][j][0]/eps),2.0));
                delta[i][j][0] = (1.0/2.0*eps) * (1.0 + cos(PI * phi[i][j][0]/eps));
            }
            //<<delta[i][j][0]<<" ";
        }
        //<<endl;
    }
    //exit(0);
}

void grad_func(double ***grad_phi, double ***phi)
{
  int i,j;
  double ***grad_phix, ***grad_phiy;
  allocator3(&grad_phix, xelem, yelem, zelem);
  allocator3(&grad_phiy, xelem, yelem, zelem);

  double ***phiRface, ***phiTface;
  allocator3(&phiRface, xelem, yelem, zelem);
  allocator3(&phiTface, xelem, yelem, zelem);
     for(i=1; i<xelem-2; i++)
     {
         for(j=1; j<yelem-2; j++)
         {
             phiRface[i][j][0] = 0.5*(phi[i+1][j][0] + phi[i][j][0]);
             phiTface[i][j][0] = 0.5*(phi[i][j+1][0] + phi[i][j][0]);
         }
     }
     
     for(j=2; j<yelem-2; j++)
     {
         for(i=2; i<xelem-2; i++)
         {
             grad_phix[i][j][0] = (phiRface[i][j][0] - phiRface[i-1][j][0])/area[i][j][1][1];
             grad_phiy[i][j][0] = (phiRface[i][j][0] - phiRface[i][j-1][0])/area[i][j][0][0];
             //if(delta[i][j][0] != 0.0){
             //<<grad_phix[i][j][0]<<" ";
             //}
         }
         //<<endl;
     }
     //exit(0);
     //Need to impose BC for grad_phix and grad_phiy
     grad_level_setBC(grad_phix);
     grad_level_setBC(grad_phiy);
     
     for(i=0; i<xelem; i++)
     {
         for(j=0; j<yelem; j++)
         {
             grad_phi[i][j][0] = sqrt(pow(grad_phix[i][j][0],2.0) + pow(grad_phiy[i][j][0],2.0));
         }
     }

     deallocator3(&grad_phix, xelem, yelem, zelem);
     deallocator3(&grad_phiy, xelem, yelem, zelem);
     deallocator3(&phiRface, xelem, yelem, zelem);
     deallocator3(&phiTface, xelem, yelem, zelem);
}


void vol_contraint(double ***phi2, double ***phi, double ***grad_phi, double ***delta, double deltat)
{
  int i,j;
    for(i=0; i<xelem; i++)
    {
        for(j=0; j<yelem; j++)
        {
            if(delta[i][j][0] != 0.0)
            {
                phi2[i][j][0] = phi[i][j][0];
            }
        }
    }
}*/