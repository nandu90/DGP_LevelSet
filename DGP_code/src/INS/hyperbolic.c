/***************************************************************************

Author: nsaini
Created: 2018-04-25

***************************************************************************/


#include "common.h"
#include "INS.h"
#include "memory.h"
#include "commu.h"
#include "generalFunc.h"


void hyperbolic(struct elemsclr sclr, double ****area)
{
  //printf("%d is here1\n",myrank);
  int i,j;
    /***Store phi values in a separate matrix***/
  double **phi2;
  allocator2(&phi2, xelem, yelem);
    for(i=0; i< xelem; i++)
    {
        for(j=0; j< yelem; j++)
        {
            phi2[i][j] = sclr.phi2[i][j];
        }
    }
    
   
    /**Compute eps based on grid size*/
    double eps=epsilon*max(xlen/(gxelem), ylen/(gyelem));
    
     double deltat=re_time;
    
    double ires=0.0;  
    
    bool exitflag = false;
    
    /****Heavyside and delta functions for volume constraint***/
    double **H, **delta, **grad_phi;
    allocator2(&H, xelem, yelem);
    allocator2(&delta, xelem, yelem);
    allocator2(&grad_phi, xelem, yelem);
    heavy_func(H,sclr.phi2,eps);
    delta_func(delta,sclr.phi2,eps);
    grad_func(grad_phi, sclr.phi2, area, sclr.iBC);
    
    double **lambda;
    allocator2(&lambda, xelem, yelem);
    
    double **temp_phi2;
    allocator2(&temp_phi2, xelem, yelem);

    double **rhs;
    allocator2(&rhs, xelem, yelem);

    double **phistar;
    allocator2(&phistar, xelem, yelem);

    double **rhstar;
    allocator2(&rhstar, xelem, yelem);

    
    
    int iter;
    //printf("%d is here2\n",myrank);
    for(iter=0; iter < re_loops; iter++)
    {
      if(myrank == master)printf("    Re-dist iter = %d\n",iter+1);
      
       
        for(i=2;i<xelem-2;i++)
        {
            for(j=2;j<yelem-2;j++)
            {
                temp_phi2[i][j]=phi2[i][j];
            }
        }
        
                
        /*****Now onto calculating fluxes******/
	
        
        rhs_redist2(rhs, phi2, sclr.phi2, area);
        
        
        
        
        for(i=2; i<xelem-2; i++) 
        {
            for(j=2; j<yelem-2; j++)
            {
                phistar[i][j] = phi2[i][j] + deltat * (rhs[i][j]);
            }
        }
        
        /*if(vf_control == 1)
        {
            vol_contraint(phistar, sclr.phi, grad_phi, delta, deltat);
	    }*/
        
	//printf("%d is here3\n",myrank);
        //bothscalarBC(phistar);
	INScommu2(phistar);
        INSlevel_setBC(phistar, sclr.iBC);

	//printf("%d id here4\n",myrank);
	
	
        //Calculate the star fluxes
        rhs_redist2(rhstar, phistar, sclr.phi2, area);
        
        for(i=2; i<xelem-2; i++)
        {
            for(j=2; j<yelem-2; j++)
            {
                phi2[i][j] = phistar[i][j] + 0.5*deltat *(rhs[i][j] + rhstar[i][j]);
                //<<" "<<phi[i][j][0]<<endl;
            }
        }
        //printf("%d is here5\n",myrank);
        //level_setBC(phi2);
        /*if(vf_control == 1)
        {
            vol_contraint(phi2, sclr.phi, grad_phi, delta, deltat);
	    }*/
	INScommu2(phi2);
        INSlevel_setBC(phi2, sclr.iBC);
        //printf("%d is here6\n",myrank);
        //bothscalarBC(phistar);
        /****Apply volume constraint****/
        
        /*if(exitflag == false)
        {
            monitor_res_redist(&ires, &exitflag, iter,  phi2,  temp_phi2);
        }
        else
        {
            break;
	    }*/
        
	//printf("%d is here7\n",myrank);
        
    }
    
    INScommu2(phi2);
    INSlevel_setBC(phi2, sclr.iBC);
    
    /*Reassign values*/
    for(i=0; i< xelem; i++)
    {
        for(j=0; j< yelem; j++)
        {
            sclr.phi2[i][j] = phi2[i][j];
        }
    }

    
    deallocator2(&H, xelem, yelem);
    deallocator2(&delta, xelem, yelem);
    deallocator2(&grad_phi, xelem, yelem); 
    deallocator2(&phi2,xelem,yelem);

    deallocator2(&lambda, xelem, yelem);
    deallocator2(&temp_phi2, xelem, yelem);
    deallocator2(&phistar, xelem, yelem);
    deallocator2(&rhs, xelem, yelem);
    deallocator2(&rhstar, xelem, yelem);
}
