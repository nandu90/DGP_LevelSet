/***************************************************************************

Author: nsaini
Created: 2018-04-25

***************************************************************************/


#include "common.h"
#include "INS.h"
#include "memory.h"
#include "DGPFunc.h"

void heavy_funcDG(double ***H, double ***phi, double eps)
{
    //------------------------------------------------------------------------//
    //Note: Heavy func will be constructed and stored at all quadrature points
    
     //------------------------------------------------------------------------//
    //Loop Indexes
    int ielem,jelem;
    int igauss;
    int icoeff;
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Temporary variables
    double recphi;
    double *basis;
    allocator1(&basis, ncoeff);
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    for(ielem =2; ielem<xelem-2; ielem++)
    {
	for(jelem=2; jelem<yelem-2; jelem++)
	{
	    for(igauss=0; igauss<tgauss; igauss++)
	    {
		recphi = 0.0;
		basis2D(zeta[igauss][0], zeta[igauss][1], basis);
		for(icoeff=0; icoeff<ncoeff; icoeff++)
		{
		    recphi += basis[icoeff]*phi[ielem][jelem][icoeff];
		}
		if(recphi < -eps)
                {
                    H[ielem][jelem][igauss] = 0.0;
                }
                else if(recphi > eps)
                {
                    H[ielem][jelem][igauss] = 1.0;
                }
                else
                {
                    H[ielem][jelem][igauss] = 0.5 + recphi/(2.0*eps) + (1.0/(2.0*PI)*sin(PI * recphi/eps));
                }
	    }
	}
    }

    //------------------------------------------------------------------------//
    //Deallocators
    deallocator1(&basis, ncoeff);
    //------------------------------------------------------------------------//

    
}

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

void delta_func(double **delta, double **phi, double eps)
{
  int i,j;
    for(i=2; i<xelem-2; i++)
    {
        for(j=2; j<yelem-2; j++)
        {
            if(fabs(phi[i][j]) > eps)
            {
                delta[i][j] = 0.0;
            }
            else
            {
                //delta[i][j][0] = 3.0*exp(-3.0*sclr.phi[i][j][0]/eps)/(1.0 * pow(1.0 + exp(-3.0*sclr.phi[i][j][0]/eps),2.0));
                delta[i][j] = (1.0/2.0*eps) * (1.0 + cos(PI * phi[i][j]/eps));
            }
            //<<delta[i][j][0]<<" ";
        }
        //<<endl;
    }
    //exit(0);
}

void grad_func(double **grad_phi, double **phi, double ****area, int **iBC)
{
  int i,j;
  double **grad_phix, **grad_phiy;
  allocator2(&grad_phix, xelem, yelem);
  allocator2(&grad_phiy, xelem, yelem);

  double **phiRface, **phiTface;
  allocator2(&phiRface, xelem, yelem);
  allocator2(&phiTface, xelem, yelem);
  
     for(i=1; i<xelem-2; i++)
     {
         for(j=1; j<yelem-2; j++)
         {
             phiRface[i][j] = 0.5*(phi[i+1][j] + phi[i][j]);
             phiTface[i][j] = 0.5*(phi[i][j+1] + phi[i][j]);
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
     //Need to impose BC for grad_phix and grad_phiy
     grad_level_setBC(grad_phix, iBC);
     grad_level_setBC(grad_phiy, iBC);
     
     for(i=0; i<xelem; i++)
     {
         for(j=0; j<yelem; j++)
         {
             grad_phi[i][j] = sqrt(pow(grad_phix[i][j],2.0) + pow(grad_phiy[i][j],2.0));
         }
     }

     deallocator2(&grad_phix, xelem, yelem);
     deallocator2(&grad_phiy, xelem, yelem);
     deallocator2(&phiRface, xelem, yelem);
     deallocator2(&phiTface, xelem, yelem);
}


/*void vol_contraint(double ***phi2, double ***phi, double ***grad_phi, double ***delta, double deltat)
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
