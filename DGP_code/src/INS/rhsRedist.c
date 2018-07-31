/***************************************************************************

Author: nsaini
Created: 2018-04-25

***************************************************************************/

#include "common.h"
#include "INS.h"
#include "memory.h"
#include "generalFunc.h"

double minmod(double a, double b)
{
    double result;

    if(fabs(a) < fabs(b))
    {
        result = a;
    }
    else
    {
        result = b;
    }

    return result;
}

double dplus(double phi[5], int index)
{
    double result = phi[index+1] - phi[index];
    return result;
}


double dminus(double phi[5], int index)
{
    double result = phi[index] - phi[index-1];
    return result;
}

double signof(double a)
{
    double result = a/fabs(a);
    return result;
}

void rhs_redist2(double **rhs, double **phi2, double **phi, double ****area, int **iBC)
{
    double **dxbarplus;
    double **dxbarminus;
    double **dxbar;

    allocator2(&dxbarplus, xelem, yelem);
    allocator2(&dxbarminus, xelem, yelem);
    allocator2(&dxbar, xelem, yelem);

    int i,j;

    double temp_phi[5];
    for(i=0; i<5; i++)
    {
	temp_phi[i] = 0.0;
    }
    
    for(i=2; i<xelem-2; i++)
    {
        for(j=2; j<yelem-2; j++)
        {
            
            if(i == 2)
            {
                temp_phi[2] = phi2[i][j];
                temp_phi[3] = phi2[i+1][j];
                temp_phi[4] = phi2[i+2][j];
                temp_phi[1] = phi2[i-1][j];
                if(x_bound == 1 || x_bound ==2)
                {
		  if(iBC[i-1][j]==2)
		    {
		      temp_phi[0] = temp_phi[1];
		    }
		  else
		    {
		      temp_phi[0] = phi2[i-2][j];
		    }
                }
                else if(x_bound ==3)
                {
                    temp_phi[0] = phi2[xelem-3][j];
                }

            }
            else if (i == xelem-3)
            {
                temp_phi[2] = phi2[i][j];
                temp_phi[3] = phi2[i+1][j];
                if(x_bound == 1 || x_bound ==2)
                {
		  if(iBC[i+1][j]==2)
		    {
		      temp_phi[4] = temp_phi[3];
		    }
		  else
		    {
		      temp_phi[4] = phi2[i+2][j];
		    }
                }
                else if(x_bound ==3)
                {
                    temp_phi[4] = phi2[2][j];
                }

                temp_phi[1] = phi2[i-1][j];
                temp_phi[0] = phi2[i-2][j];
            }
            else
            {
                temp_phi[2] = phi2[i][j];
                temp_phi[3] = phi2[i+1][j];
                temp_phi[4] = phi2[i+2][j];
                temp_phi[1] = phi2[i-1][j];
                temp_phi[0] = phi2[i-2][j];
            }

            dxbarplus[i][j] = dplus(temp_phi, 2) -0.5*minmod(dplus(temp_phi,2) - dplus(temp_phi,1) , dplus(temp_phi,3) - dplus(temp_phi,2));
            dxbarminus[i][j] = dminus(temp_phi, 2) +0.5*minmod(dplus(temp_phi,2) - dplus(temp_phi,1) , dplus(temp_phi,1) - dplus(temp_phi,0));

           if(signof(phi[i][j])*dplus(temp_phi,2) < 0.0  && signof(phi[i][j])*dminus(temp_phi,2) < -signof(phi[i][j])*dplus(temp_phi,2))
           {
               dxbar[i][j] = dxbarplus[i][j];
           }
           else if(signof(phi[i][j])*dminus(temp_phi,2) > 0.0  && signof(phi[i][j])*dplus(temp_phi,2) > -signof(phi[i][j])*dminus(temp_phi,2))
           {
               dxbar[i][j] = dxbarminus[i][j];
           }
           else
           {
               dxbar[i][j] = 0.5*(dxbarplus[i][j] + dxbarminus[i][j]);
           }
        }

    }

    double **dybarplus;
    double **dybarminus;
    double **dybar;

    allocator2(&dybarplus, xelem, yelem);
    allocator2(&dybarminus, xelem, yelem);
    allocator2(&dybar, xelem, yelem);

    for(i=0; i<5; i++)
    {
	temp_phi[i] = 0.0;
    }
    
    for(i=2; i<xelem-2; i++)
    {
        for(j=2; j<yelem-2; j++)
        {
            if(j == 2)
            {
                temp_phi[2] = phi2[i][j];
                temp_phi[3] = phi2[i][j+1];
                temp_phi[4] = phi2[i][j+2];
                temp_phi[1] = phi2[i][j-1];
                if(y_bound == 1 || y_bound ==2)
                {
		  if(iBC[i][j-1]==2)
		    {
		      temp_phi[0] = phi2[i][j-1];
		    }
		  else
		    {
		      temp_phi[0] = phi2[i][j-2];
		    }
                }
                else if(y_bound ==3)
                {
                    temp_phi[0] = phi2[i][yelem-3];
                }

            }
            else if (j == yelem-3)
            {
                temp_phi[2] = phi2[i][j];
                temp_phi[3] = phi2[i][j+1];
                if(y_bound == 1 || y_bound ==2)
                {
		  if(iBC[i][j+1]==2)
		    {
                    temp_phi[4] = phi2[i][j+1];
		    }
		  else
		    {
		      temp_phi[4] = phi2[i][j+2];
		    }
                }
                else if(y_bound ==3)
                {
                    temp_phi[4] = phi2[i][2];
                }

                temp_phi[1] = phi2[i][j-1];
                temp_phi[0] = phi2[i][j-2];
            }
            else
            {
                temp_phi[2] = phi2[i][j];
                temp_phi[3] = phi2[i][j+1];
                temp_phi[4] = phi2[i][j+2];
                temp_phi[1] = phi2[i][j-1];
                temp_phi[0] = phi2[i][j-2];
            }

            dybarplus[i][j] = dplus(temp_phi, 2) -0.5*minmod(dplus(temp_phi,2) - dplus(temp_phi,1) , dplus(temp_phi,3) - dplus(temp_phi,2));
            dybarminus[i][j] = dminus(temp_phi, 2) +0.5*minmod(dplus(temp_phi,2) - dplus(temp_phi,1) , dplus(temp_phi,1) - dplus(temp_phi,0));

           if(signof(phi[i][j])*dplus(temp_phi,2) < 0.0  && signof(phi[i][j])*dminus(temp_phi,2) < -signof(phi[i][j])*dplus(temp_phi,2))
           {
               dybar[i][j] = dybarplus[i][j];
           }
           else if(signof(phi[i][j])*dminus(temp_phi,2) > 0.0  && signof(phi[i][j])*dplus(temp_phi,2) > -signof(phi[i][j])*dminus(temp_phi,2))
           {
               dybar[i][j] = dybarminus[i][j];
           }
           else
           {
               dybar[i][j] = 0.5*(dybarplus[i][j] + dybarminus[i][j]);
           }
        }

    }

    

    double eps=epsilon*max(xlen/(gxelem), ylen/(gyelem));
    for(i=2; i < xelem-2; i++)
    {
        for(j=2; j < yelem-2; j++)
        {
            double sign_phi;
            if(phi[i][j] >= eps)
            {
                sign_phi = 1.0;
            }
            else if(phi[i][j] <= -eps)
            {
                sign_phi =-1.0;
            }
            else
            {
                sign_phi = (phi[i][j]/ eps) - (1.0/PI)*sin(PI*phi[i][j]/eps);
            }

            rhs[i][j] = sign_phi*(1-sqrt(pow(dxbar[i][j]/area[i][j][1][1],2.0) + pow(dybar[i][j]/area[i][j][0][0],2.0)));
        }
    }

    deallocator2(&dxbarplus, xelem, yelem);
    deallocator2(&dxbarminus, xelem, yelem);
    deallocator2(&dxbar, xelem, yelem);
    deallocator2(&dybarplus, xelem, yelem);
    deallocator2(&dybarminus, xelem, yelem);
    deallocator2(&dybar, xelem, yelem);
}
