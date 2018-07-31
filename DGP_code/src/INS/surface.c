/***************************************************************************

Author: nsaini
Created: 2018-04-25

***************************************************************************/

#include "common.h"
#include "INS.h"
#include "memory.h"
#include "commu.h"
#include "generalFunc.h"

void surface(struct elemsclr sclr, double **st_forcex, double **st_forcey, double ****area)
{
    int i,j;
    /**Compute eps based on grid size*/
    double eps=epsilon*max(xlen/(gxelem), ylen/(gyelem));
    
    /**Compute Heavyside function**/
    double **H;
    allocator2(&H, xelem, yelem);
    
    heavy_func(H, sclr.phi2, eps);
    
    find_density_visc(H, sclr.rho, sclr.mu);
    //<<eps<<endl;
    
    
    
    if(sf_toggle == 1)
    {
      double **delta;
      allocator2(&delta, xelem, yelem);
        for(i=2; i<xelem-2; i++)
        {
            for(j=2; j<yelem-2; j++)
            {
                if(fabs(sclr.phi2[i][j]) > eps)
                {
                    delta[i][j] = 0.0;
                }
                else
                {
                    //delta[i][j][0] = 3.0*exp(-3.0*sclr.phi2[i][j][0]/eps)/(1.0 * pow(1.0 + exp(-3.0*sclr.phi2[i][j][0]/eps),2.0));
                    delta[i][j] = (1.0/2.0*eps) * (1.0 + cos(PI * sclr.phi2[i][j]/eps));
                }
                //<<delta[i][j][0]<<" ";
            }
            //<<endl;
        }
        //exit(0);
	double **del_scaling;
	allocator2(&del_scaling, xelem, yelem);

         //This is equivalent to marker function
         for(i=2 ; i<xelem-2; i++)
         {
             for(j=2; j<yelem-2; j++)
             {
                 del_scaling[i][j] = 2.0*H[i][j]*delta[i][j];
                 //<<del_scaling[i][j][0]<<" ";
             }
             //<<endl;
         }
         //exit(0);

         /***Now calculate gradient of level set for ultimately calculating curvature*/
	 double **grad_phix;
	 double **grad_phiy;
	 double **phiRface;
	 double **phiTface;
	 allocator2(&phiRface, xelem, yelem);
	 allocator2(&phiTface, xelem, yelem);
	 allocator2(&grad_phiy, xelem, yelem);
	 allocator2(&grad_phix, xelem, yelem);
         for(i=1; i<xelem-2; i++)
         {
             for(j=1; j<yelem-2; j++)
             {
                 phiRface[i][j] = 0.5*(sclr.phi2[i+1][j] + sclr.phi2[i][j]);
                 phiTface[i][j] = 0.5*(sclr.phi2[i][j+1] + sclr.phi2[i][j]);
             }
         }

	 

         for(j=2; j<yelem-2; j++)
         {
             for(i=2; i<xelem-2; i++)
             {
                 grad_phix[i][j] = (phiRface[i][j] - phiRface[i-1][j])/area[i][j][1][1];
                 grad_phiy[i][j] = (phiRface[i][j] - phiRface[i][j-1])/area[i][j][0][0];
                 //if(delta[i][j][0] != 0.0){
                 //<<grad_phix[i][j][0]<<" ";
                 //}
             }
             //<<endl;
         }
         //exit(0);
         /*Need to impose BC for grad_phix and grad_phiy*/
	 INScommu2(grad_phix);
	 INScommu2(grad_phiy);
         grad_level_setBC(grad_phix, sclr.iBC);
         grad_level_setBC(grad_phiy, sclr.iBC);

         /***Compute double and mixed derivatives***/
	 double **grad_phixx;
	 double **grad_phixy;
	 double **grad_phiyy;
	 double **phixRface;
	 double **phiyTface;
	 double **phixTface;
	 allocator2(&grad_phixx, xelem, yelem);
	 allocator2(&grad_phixy, xelem, yelem);
	 allocator2(&grad_phiyy, xelem, yelem);
	 allocator2(&phixRface, xelem, yelem);
	 allocator2(&phiyTface, xelem, yelem);
	 allocator2(&phixTface, xelem, yelem);

         for(i=1; i<xelem - 2; i++)
         {
             for(j=1; j<yelem-2; j++)
             {
                 phixRface[i][j] = 0.5*(grad_phix[i][j] + grad_phix[i+1][j]);
                 phiyTface[i][j] = 0.5*(grad_phiy[i][j] + grad_phiy[i][j+1]);
                 phixTface[i][j] = 0.5*(grad_phix[i][j] + grad_phix[i][j+1]);
             }
         }

         for(i=2; i<xelem-2; i++)
         {
             for(j=2; j<yelem-2; j++)
             {
                 grad_phixx[i][j] = (phixRface[i][j] - phixRface[i-1][j])/area[i][j][1][1];
                 grad_phiyy[i][j] = (phiyTface[i][j] - phiyTface[i][j-1])/area[i][j][0][0];
                 grad_phixy[i][j] = (phixTface[i][j] - phixTface[i][j-1])/area[i][j][0][0];
             }
         }
	 double **curvature;
	 allocator2(&curvature, xelem, yelem);


         for(j=2; j<yelem-2; j++)
         {
             for(i=2; i<xelem-2; i++)
             {
                 //<<i<<" "<<j<<endl;
                 curvature[i][j] = (pow(grad_phiy[i][j],2.0)*grad_phixx[i][j] + pow(grad_phix[i][j],2.0)*grad_phiyy[i][j] - 2.0*grad_phix[i][j]*grad_phiy[i][j]*grad_phixy[i][j]);
                 curvature[i][j] = curvature[i][j]/pow((pow(grad_phix[i][j],2.0) + pow(grad_phiy[i][j],2.0)),1.5);
                 /*if(delta[i][j][0] != 0)
                 {
                 <<curvature[i][j][0]<<" ";
                 }*/
             }
             //<<endl;

         }
         //exit(0);

         for(j=2; j<yelem-2; j++)
         {
             for(i=2; i<xelem-2; i++)
             {
                 st_forcex[i][j] = sf_coeff*curvature[i][j]*del_scaling[i][j]*grad_phix[i][j];
                 st_forcey[i][j] = sf_coeff*curvature[i][j]*del_scaling[i][j]*grad_phiy[i][j];
                 //<<st_forcey[i][j][0]<<" ";
             }
             //<<endl;
         }
         //exit(0);
	 INScommu2(st_forcex);
	 INScommu2(st_forcey);
         grad_level_setBC(st_forcex, sclr.iBC);
         grad_level_setBC(st_forcey, sclr.iBC);

	 deallocator2(&curvature, xelem, yelem);
	 deallocator2(&grad_phixx, xelem, yelem);
	 deallocator2(&grad_phixy, xelem, yelem);
	 deallocator2(&grad_phiyy, xelem, yelem);
	 deallocator2(&phixRface, xelem, yelem);
	 deallocator2(&phiyTface, xelem, yelem);
	 deallocator2(&phixTface, xelem, yelem);
	 deallocator2(&grad_phix, xelem, yelem);
	 deallocator2(&grad_phiy, xelem, yelem);
	 deallocator2(&phiRface, xelem, yelem);
	 deallocator2(&phiTface, xelem, yelem);
	 deallocator2(&del_scaling, xelem, yelem);
	 deallocator2(&delta, xelem, yelem);
	 
    }
    deallocator2(&H, xelem, yelem);
}
