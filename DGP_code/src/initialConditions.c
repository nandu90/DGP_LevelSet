/***************************************************************************

Author: nsaini
Created: 2018-03-04

***************************************************************************/

#include "common.h"
#include "icbc.h"

void initializeLS(struct elemsclr elem)
{
    int i,j,k,l;

  
    
    
    for (i=0; i<xelem; i++)
    {
        for (j=0; j<yelem; j++)
        {
	    //Loop over Gauss Quadrature points
	    for (k=0; k<tgauss; k++)
	    {
		elem.phi[i][j][k] = sqrt(pow(xb_in - elem.xgauss[i][j][k],2) + pow(yb_in - elem.ygauss[i][j][k],2)) - rb_in;
	    }
	}
    }
    
    /*//if(myrank == master)printf("here");
    if(case_tog == 4)
      {
        for (i=0; i<xelem; i++)
	  {
            for (j=0; j<yelem; j++)
	      {
                double a,b,c,d;
                a=(0.5-0.03)*100.0;
                b=(0.5+0.03)*100.0;
                c=(0.75-0.15)*100.0;
                d=(0.75+0.15-0.06)*100.0;
                if(xc[i][j]-a == 0.0 || xc[i][j]-b == 0.0 || yc[i][j]-c == 0.0 || yc[i][j]-d == 0.0)
		  {
                    phi[i][j][0] =0.0;
		  }
                else
		  {
		    
                    double dist[4];
                    dist[0] = sqrt(pow(xc[i][j]-a,2.0) + pow(yc[i][j]-c,2.0));
                    dist[1] = sqrt(pow(xc[i][j]-b,2.0) + pow(yc[i][j]-c,2.0));
                    dist[2] = sqrt(pow(xc[i][j]-b,2.0) + pow(yc[i][j]-d,2.0));
                    dist[3] = sqrt(pow(xc[i][j]-a,2.0) + pow(yc[i][j]-d,2.0));
		    
                    int min = 0;
                    for (k=1; k<4; k++)
		      {
                        if(dist[min] > dist[k])
			  {
                            min=k;
			  }
		      }
		    
                    double base, height;
                    if(min == 0)
		      {
                        base = a - xc[i][j];
                        height = c - yc[i][j];
		      }
                    else if(min == 1)
		      {
                        base = b - xc[i][j];
                        height = c - yc[i][j];
		      }
                    else if(min == 2)
		      {
                        base = b - xc[i][j];
                        height = d - yc[i][j];
		      }
                    else if(min == 3)
		      {
                        base = a - xc[i][j];
                        height = d - yc[i][j];
		      }

                    double angle;
		    angle = atan2(height,base);
		    if(min == 0)
		      {
			if(angle < PI/4.0 && angle > -3.0*PI/4.0)
			  {
			    phi[i][j][0] = base;//(height);
			  }
			else
			  {
			    phi[i][j][0] = height;//(base);
			  }
		      }
		    
		    else if(min == 2)
		      {
			if(angle < PI/4.0 && angle > -3.0*PI/4.0)
			  {
			    phi[i][j][0] = -height;//(base);
			  }
			else
			  {
			    phi[i][j][0] = -base;//(height);
			  }
		      }
		    
		    else if(min == 1)
		      {
			if(angle < 3.0*PI/4.0 && angle > -PI/4.0)
			  {
			    phi[i][j][0] = height;//(base);
			  }
			else
			  {
			    phi[i][j][0] = -base;//(height);
			  }
		      }

		    else if (min == 3)
		      {
			if(angle <= 3.0*PI/4.0 && angle >= -PI/4.0)
			  {
			    phi[i][j][0] = base;//(height);
			  }
			else
			  {
			    phi[i][j][0] = -height;//(base);
			  }
		      }
		  }
		
		
	    
	      }
	  }

	//if(myrank == master)printf("here");
        for(i=0; i<xelem; i++)
        {
            for(j=0; j<yelem; j++)
            {
                phi[i][j][0] = -phi[i][j][0];
            }
        }


        double temp_phi[xelem][yelem];
        for(i=0; i<xelem; i++)
        {
            for(j=0; j<yelem; j++)
            {
                temp_phi[i][j] = phi[i][j][0];
            }
        }

        for(i=0; i<xelem; i++)
        {
            for(j=0; j<yelem; j++)
            {
                double cir = sqrt(pow(xc[i][j]-50.0,2.0) + pow(yc[i][j]-75.0,2.0)) - 15.0;
                if(phi[i][j][0] < 0.0)
                {
                    phi[i][j][0] = max(cir,phi[i][j][0]);
                }
            }
        }
	}*/
}

void initialize(struct elemsclr elem)
{
    int i,j;
    //calcp(sclr); //For the simple advection case
    //pressureBC(sclr.p);

    /*This initializes a single bubble*/
      
    initializeLS(elem);
    /*level_setBC(elem.phi);
    double eps=epsilon*max(xlen/(gxelem), ylen/(gyelem));

    double ***H;
    allocator3(&H, xelem, yelem, zelem);

    heavy_func(H, elem.phi, eps);


    /*****Bubble breakup case****/
    /*if (case_tog == 3)
    {
      double line = 2.5*2.0*rb_in;
      double temp_phi[xelem][yelem];
      for(i=0; i<xelem; i++)
	{
	  for(j=0; j<yelem; j++)
	    {
	      temp_phi[i][j] = line - yc[i][j];
	    }
	}

      for(i=0; i<xelem; i++)
	{
	  for(j=0; j<yelem; j++)
	    {
	      if(fabs(temp_phi[i][j]) < fabs(elem.phi[i][j][0]))
		{
		  elem.phi[i][j][0] = temp_phi[i][j];
		}
	    }
	}
      level_setBC(elem.phi);
    }

    /****Pressure*****/
    /*for(i=1; i<xelem-1; i++)
    {
        for(j=1; j<yelem-1; j++)
        {
            if(yc[i][j] < line)
            {
                elem.p[i][j][0] = rhof * fabs(gy)*fabs(line-yc[i][j]);
            }
            else
            {
                elem.p[i][j][0] = 0.0;
            }
        }
    }
    double p_at_cen = rhof * fabs(gy)*fabs(line-yb_in);*/



    /*for(i=1; i<xelem-1; i++)
    {
        for(j=1; j<yelem-1; j++)
        {
            if(H[i][j][0] < 1.0 && yc[i][j] < line)
            {
                elem.p[i][j][0] = (1.0 - H[i][j][0])*p_at_cen;
            }
        }
    }
    /*********End of Bubble breakup initial conditions************/

    /***********************************************/

    /****Initializes density and viscosity based on phi****/
    /*find_density_visc(H, elem.rho, elem.mu);
    level_setBC(elem.rho);
    level_setBC(elem.mu);

    /*if(myrank == 1)
      {
	for(j=0; j<yelem; j++)
	  {
	    printf("%d %d %.6f\n",xelem-3,j,H[xelem-2][j][0]);
	  }
	  }*/


    //deallocator3(&H, xelem, yelem, zelem);

}

