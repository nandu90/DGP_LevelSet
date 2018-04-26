/***************************************************************************

Author: nsaini
Created: 2018-04-25

***************************************************************************/


#include "common.h"
#include "INS.h"
#include "memory.h"

void flux(double u1R, double u1L, double u1T, double u1B, double u2R, double u2L, double u2T, double u2B, int dir, double *flux1, double *flux2)
{
    if(dir ==1)
    {
      (*flux1) = pow(u1R + u1L,2.0)/4.0;
      (*flux2) = (u2R +u2L)*0.5 * (u2T + u2B)*0.5;
    }
    else if (dir == 2)
    {
      (*flux1) = (u1R +u1L)*0.5 * (u1T + u1B)*0.5;
      (*flux2) =  pow(u2T + u2B,2.0)/4.0;
    }
}




void rhscalc(struct elemsclr sclr, double **rhsx, double **rhsy, double ****area, double **vol, int **iBC, double **u, double **v)
{


  int i,j;


    //***Calculate contribution from advection
    //Note for index i,j the CV under consideration is the CV between i,j and i+1,j
  double **advx, **advy;
  allocator2(&advx, xelem, yelem);
  allocator2(&advy, xelem, yelem);

  double flux1=0.0, flux2=0.0;
  
    for(i=1; i<xelem-2; i++)
    {
        for(j=2; j<yelem-2; j++)
        {
            //Calculate the velocities at the edge centers of CV
            double u1R=0.0, u1L=0.0;
            double u1T=0.0, u1B=0.0;
            double u2R=0.0, u2L=0.0;
            double u2T=0.0, u2B=0.0;

            double u1LL=0.0, u1RR=0.0;

            u1R = u[i+1][j];
            u1L = u[i][j];

            u1T = v[i+1][j];
            u1B = v[i+1][j-1];

            u2R = v[i+1][j];
            u2L = v[i][j];

            u2T = u[i][j+1];
            u2B = u[i][j];

            if(i==1)
            {
                u1LL = u[xelem-3][j];
                u1RR = u[i+2][j];
            }
            else if(i==xelem-3)
            {
                u1LL = u[i-1][j];
                u1RR = u[2][j];
            }
            else
            {
                u1LL = u[i-1][j];
                u1RR = u[i+2][j];
            }


            
            flux(u1R, u1L, u1T, u1B, u2R, u2L, u2T, u2B, 1, &flux1, &flux2);
            //quick(u1R, u1L, u1T, u1B, u2R, u2L, u2T, u2B, u1LL, u1RR, 1, flux1, flux2);
            advx[i][j]=(area[i][j][0][0]*flux1 + area[i][j][1][1]*flux2)/vol[i][j];

        }
    }

    for (i=2; i< xelem-2 ;i++)
    {
        for (j=1; j< yelem-2; j++)
        {
            double v1R=0.0, v1L=0.0;
            double v1T=0.0, v1B=0.0;
            double v2R=0.0, v2L=0.0;
            double v2T=0.0, v2B=0.0;

            double v2TT=0.0, v2BB=0.0;

            v1R = v[i+1][j];
            v1L = v[i][j];

            v1T = u[i][j+1];
            v1B = u[i][j];

            v2R = u[i][j+1];
            v2L = u[i-1][j+1];

            v2T = v[i][j+1];
            v2B = v[i][j];

            if(j==1)
            {
                v2TT = v[i][j+2];
                v2BB = 0.0;
            }
            else if(j==yelem-3)
            {
                v2TT = 0.0;
                v2BB = v[i][j-1];
            }
            else
            {
                v2TT = v[i][j+2];
                v2BB = v[i][j-1];
            }


            flux(v1R, v1L, v1T, v1B, v2R, v2L, v2T, v2B, 2, &flux1, &flux2);
            //quick(v1R, v1L, v1T, v1B, v2R, v2L, v2T, v2B, v2BB, v2TT, 2, flux1, flux2);
            advy[i][j]=(area[i][j][0][0]*flux1 + area[i][j][1][1]*flux2)/vol[i][j];
        }
    }


    //***Calculate contribution from diffusion
    double **diffx, **diffy;
    allocator2(&diffx, xelem, yelem);
    allocator2(&diffy, xelem, yelem);

    for(i=1; i<xelem-2; i++)
    {
        for(j=2; j<yelem-2; j++)
        {
            //Simple taylor series approximation of second diff is used
            double uc=0.0, uR=0.0, uL=0.0, uT=0.0, uB=0.0;
            if (i == 1)
            {
                uc = u[i][j];
                uR = u[i+1][j];
                if(x_bound == 1)
                {
		  if(iBC[i][j]==2)
		    {
		      uL = -u[i+1][j];
		    }
		  else
		    {
		      uL = u[i-1][j];
		    }
                }
                else if(x_bound == 2)
                {
		  if(iBC[i][j]==2)
		    {
		      uL = u[i+1][j];
		    }
		  else
		    {
		      uL = u[i-1][j];
		    }
                }
                else if(x_bound == 3)
                {
                    uL = u[xelem-3][j];
                }

                uT = u[i][j+1];
                uB = u[i][j-1];
            }
            else
            {
                uc = u[i][j];
                uR = u[i+1][j];
                uL = u[i-1][j];
                uT = u[i][j+1];
                uB = u[i][j-1];
            }

            diffx[i][j] = ((sclr.mu[i+1][j] + sclr.mu[i][j])/(sclr.rho[i+1][j] + sclr.rho[i][j]))*((uR + uL - 2.0*uc)/pow(area[i][j][1][1],2.0) + (uT + uB - 2.0*uc)/pow(area[i][j][0][0],2.0));
        }
    }

    for(i=2; i<xelem-2; i++)
    {
        for(j=1; j<yelem-2; j++)
        {
            //Simple taylor series approximation of second diff is used
            double vc=0.0, vR=0.0, vL=0.0, vT=0.0, vB=0.0;
            if (j == 1)
            {
                vc = v[i][j];
                vR = v[i+1][j];
                vL = v[i-1][j];
                vT = v[i][j+1];
                if(y_bound == 1)
                {
		  if(iBC[i][j]==2)
		    {
		      vB = -vT;
		      //vB = 0.0;
		    }
		  else
		    {
		      vB = v[i][j-1];
		    }
                }
                else if(y_bound == 2)
                {
                    if(iBC[i][j]==2)
		    {
		      vB = vT;
		      //vB = 0.0;
		    }
		  else
		    {
		      vB = v[i][j-1];
		    }
                }
                else if(y_bound == 3)
                {
                    vB = v[i][yelem-3];
                }

            }
            else
            {
                vc = v[i][j];
                vR = v[i+1][j];
                vL = v[i-1][j];
                vT = v[i][j+1];
                vB = v[i][j-1];
            }
            diffy[i][j] = ((sclr.mu[i][j+1] + sclr.mu[i][j])/(sclr.rho[i][j+1] + sclr.rho[i][j]))*((vR + vL - 2.0*vc)/pow(area[i][j][1][1],2.0) + (vT + vB - 2.0*vc)/pow(area[i][j][0][0],2.0));
        }
    }





    for(j=2; j<yelem-2; j++)
    {
        for(i=2; i<xelem-2; i++)
        {
            rhsx[i][j]=diffx[i][j]-(advx[i][j]-advx[i-1][j]);
            rhsy[i][j]=diffy[i][j]-(advy[i][j]-advy[i][j-1]);

            

        }
    }

    deallocator2(&diffx, xelem, yelem);
    deallocator2(&diffy, xelem, yelem);
    deallocator2(&advx, xelem, yelem);
    deallocator2(&advy, xelem, yelem);
}
