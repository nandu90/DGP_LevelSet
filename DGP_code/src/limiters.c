/***************************************************************************

Author: nsaini
Created: 2018-08-14

***************************************************************************/


#include "common.h"
#include "generalFunc.h"
#include "memory.h"
#include "DGPFunc.h"
#include "memory.h"


void momentLimiter(double ***phi)
{
    
    //------------------------------------------------------------------------//
    int ielem, jelem;
    //int icoeff;
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //For traversing the coefficients
    int p;
    int m, n;

    double tempphi;
    double tempphi1, tempphi2;
    
    int z1jump = 1;
    int z2jump = polyorder+1;
    
    int current;
    int current1, current2;
    int otrack;
    int porder;

    int flag;

    double *basis;
    allocator1(&basis, ncoeff);

    //double recphi;

    for(ielem=2; ielem<xelem-2; ielem++)
    {
	for(jelem = 2; jelem<yelem-2; jelem++)
	{
	    flag = 0;
	    p = ncoeff-1;
	    porder = polyorder;
	    
	    while(p > 0)
	    {
		//Limit p,p
		current = p;
		tempphi = phi[ielem][jelem][current];
		limitCoeff(phi, ielem, jelem, porder, porder, current);
		
		if(fabs(tempphi-phi[ielem][jelem][current]) < 1e-15 && fabs(tempphi) > 1e-16)
		{
		    break;
		}
		else
		{
		    //Check for mixed derivatives
		    //Limit z1
		    otrack = porder;
		    m = porder;
		    n = porder;
		    for(otrack=porder; otrack > 0; otrack--)
		    {
			m -= 1;
			current1 = p-z1jump*(porder+1-otrack);
			tempphi1 = phi[ielem][jelem][current1];
			limitCoeff(phi, ielem, jelem, m ,porder, current1);

			if(fabs(tempphi1-phi[ielem][jelem][current1]) < 1e-15 && fabs(tempphi1) > 1e-16)
			{
			    flag = 1;
			    break;
			}
		    }

		    //Limit z2
		    otrack = porder;
		    m = porder;
		    n = porder;
		    for(otrack=porder; otrack > 0; otrack--)
		    {
			n -= 1;
			current2 = p-z2jump*(porder+1-otrack);
			tempphi2 = phi[ielem][jelem][current2];
			limitCoeff(phi, ielem, jelem, porder, n, current2);

			if(fabs(tempphi2-phi[ielem][jelem][current2]) < 1e-15 && fabs(tempphi2) > 1e-16)
			{
			    flag = 1;
			    break;
			}
		    }
		    
		    if(flag == 1)
		    {
			break;
		    }
		    else
		    {
			porder--;
			p -= (polyorder+2);
		    }
		}
		
	    }
	}
    }
    deallocator1(&basis, ncoeff);
}

void limitCoeff(double ***phi, int i, int j, int m, int n, int current)
{
    double a, b, c, d, e;
    a = 0.0;
    b = 0.0;
    c = 0.0;
    d = 0.0;
    e = 0.0;
    
    int z1jump = 1;
    int z2jump = polyorder+1;
        
    double alpham; 
    double alphan;

    
    a = phi[i][j][current];
    if(n != 0)
    {
	alphan = 1.0;
	
	b = alphan*(phi[i][j+1][current-z2jump] - phi[i][j][current-z2jump]);
	c = alphan*(phi[i][j][current-z2jump] - phi[i][j-1][current-z2jump]);
    }
    if(m != 0)
    {
	alpham = 1.0;
	
	d = alpham*(phi[i+1][j][current-z1jump] - phi[i][j][current-z1jump]);
	e = alpham*(phi[i][j][current-z1jump] - phi[i-1][j][current-z1jump]);
    }

    double *array;
    if(m == 0 || n==0)
    {
	allocator1(&array, 3);
	array[0] = a;
	if(m==0)
	{
	    array[1] = b;
	    array[2] = c;
	}
	else if(n==0)
	{
	    array[1] = d;
	    array[2] = e;
	}
	phi[i][j][current] = minmod(array,3);
    }
    else
    {
	allocator1(&array, 5);
	array[0] = a;
	array[1] = b;
	array[2] = c;
	array[3] = d;
	array[4] = e;
	phi[i][j][current] = minmod(array,5);
        
    }

    //------------------------------------------------------------------------//
    //Deallocate
    if(m == 0 || n==0)
    {
	deallocator1(&array, 3);
    }
    else
    {
	deallocator1(&array, 5);
    }
    //------------------------------------------------------------------------//

}
