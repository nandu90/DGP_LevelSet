/***************************************************************************

Author: nsaini
Created: 2018-08-13

***************************************************************************/


#include "common.h"
#include "functions.h"
#include "memory.h"



void cockburn(double **phi)
{
    //------------------------------------------------------------------------//
    int ielem;
    int icoeff;

    double recphiL;
    double recphiR;
    double sR, sL;

    double *basis;
    allocator1(&basis, ncoeff);
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
       
    double a, b, c;
    for(ielem = 1; ielem<xelem-1; ielem++)
    {	
	//REconstruct solution at the right face
	recphiR = 0.0;
	basis1D(1.0, basis);
	for(icoeff = 0; icoeff<ncoeff; icoeff++)
	{
	    recphiR += basis[icoeff] * phi[ielem][icoeff]; 
	}

	sR = recphiR - phi[ielem][0];

	//Reconstruct solution at the left face
	recphiL = 0.0;
	basis1D(-1.0, basis);
	for(icoeff = 0; icoeff<ncoeff; icoeff++)
	{
	    recphiL += basis[icoeff] * phi[ielem][icoeff]; 
	}
	sL = recphiL - phi[ielem][0];

	a = sR;
	b = phi[ielem][0] - phi[ielem-1][0];
	c = phi[ielem+1][0] - phi[ielem][0];
	sR = minmod(a, b, c);


	a = -sL;
	b = phi[ielem][0] - phi[ielem-1][0];
	c = phi[ielem+1][0] - phi[ielem][0];
	sL = minmod(a, b, c);



	phi[ielem][1] = (sR + sL)*0.5;

    }

    deallocator1(&basis, ncoeff);
}

void momentLimiter(double **phi)
{
    //------------------------------------------------------------------------//
    int ielem;
    int icoeff;

    double *basis;
    allocator1(&basis, ncoeff);

    //int k;
    double a, b, c;
    double alpha;

    double tempphi;

    //------------------------------------------------------------------------//

       
    for(ielem = 1; ielem<xelem-1; ielem++)
    {
	//Limit from higher to lower
	for(icoeff=ncoeff-1; icoeff>0; icoeff--)
	{
	    //k = icoeff - 1;
	    //alpha =2.0*k+1.0;
	    alpha = 1.0;
	    a = phi[ielem][icoeff];
	    b = alpha*(phi[ielem+1][icoeff-1] - phi[ielem][icoeff-1]);
	    c = alpha*(phi[ielem][icoeff-1] - phi[ielem-1][icoeff-1]);

	    tempphi = minmod(a, b, c);
	    
	    if(fabs(tempphi - phi[ielem][icoeff]) < 1e-15)
	    {
		break;
	    }
	    else
	    {
		phi[ielem][icoeff] = tempphi;
	    }   
	}
	
    }
    
    deallocator1(&basis, ncoeff);
}
