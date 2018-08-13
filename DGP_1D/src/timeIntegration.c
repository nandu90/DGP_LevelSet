/***************************************************************************

Author: nsaini
Created: 2018-07-25

***************************************************************************/

#include "common.h"
#include "functions.h"
#include "memory.h"


void eulerIncrement(double **inc, double ***mass, double **rhs, double deltat)
{
    //------------------------------------------------------------------------//
    //Loop indexes
    int ielem;
    int icoeff;
    //------------------------------------------------------------------------//

    for(ielem = 1; ielem<xelem-1; ielem++)
    {
	    for(icoeff = 0; icoeff< ncoeff; icoeff++)
	    {
		inc[ielem][icoeff] = deltat*rhs[ielem][icoeff]/mass[ielem][icoeff][icoeff];
	    }
    }
}


void Runge_Kutta(struct elemsclr elem, double *x, double deltat, double **rhs)
{
    //------------------------------------------------------------------------//
    int ielem, icoeff;
    //------------------------------------------------------------------------//

    double **k1, **k2, **k3, **k4;
    allocator2(&k1, xelem, ncoeff);
    allocator2(&k2, xelem, ncoeff);
    allocator2(&k3, xelem, ncoeff);
    allocator2(&k4, xelem, ncoeff);
    
    //------------------------------------------------------------------------//
    //Store the primary value in a temp array
    double **tempphi;
    allocator2(&tempphi, xelem, ncoeff);
    for(ielem =0; ielem<xelem; ielem++)
    {
	for(icoeff=0; icoeff<ncoeff; icoeff++)
	{
	    tempphi[ielem][icoeff] = elem.phi[ielem][icoeff];
	}
    }
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Implement Limiter
    if(limiter == 1)
    {
	cockburn(elem.phi);
    }
    else if(limiter == 2)
    {
	momentLimiter(elem.phi);
    }
    //------------------------------------------------------------------------//


    //------------------------------------------------------------------------//
    //1st increment
    getRHS(elem, x, rhs);
    eulerIncrement(k1, elem.mass, rhs, deltat);
    for(ielem=1; ielem<xelem-1; ielem++)
    {
	for(icoeff=0; icoeff<ncoeff; icoeff++)
	{
	    elem.phi[ielem][icoeff] = tempphi[ielem][icoeff] + 0.5*k1[ielem][icoeff];
	}
    }
    //Apply BC
    level_setBC(elem.phi);
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //2nd increment
    getRHS(elem, x, rhs);
    eulerIncrement(k2, elem.mass, rhs, deltat);
    for(ielem=1; ielem<xelem-1; ielem++)
    {
	for(icoeff=0; icoeff<ncoeff; icoeff++)
	{
	    elem.phi[ielem][icoeff] = tempphi[ielem][icoeff] + 0.5*k2[ielem][icoeff];
	}
    }

    //Apply BC
    level_setBC(elem.phi);
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //3rd Increment
    getRHS(elem, x, rhs);
    eulerIncrement(k3, elem.mass, rhs, deltat);
    for(ielem=1; ielem<xelem-1; ielem++)
    {
	for(icoeff=0; icoeff<ncoeff; icoeff++)
	{
	    elem.phi[ielem][icoeff] = tempphi[ielem][icoeff] + k3[ielem][icoeff];
	}
    }

    //Apply BC
    level_setBC(elem.phi);
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //4th Increment
    getRHS(elem, x, rhs);
    eulerIncrement(k4, elem.mass, rhs, deltat);
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Now get the next time step value
    for(ielem=1; ielem<xelem-1; ielem++)
    {
	for(icoeff=0; icoeff<ncoeff; icoeff++)
	{
	    elem.phi[ielem][icoeff] = tempphi[ielem][icoeff] + (1.0/6.0)*(k1[ielem][icoeff] + 2.0*k2[ielem][icoeff] + 2.0*k3[ielem][icoeff] + k4[ielem][icoeff]);
	}
    }

    //Apply BC
    level_setBC(elem.phi);
    //------------------------------------------------------------------------//

	    
    deallocator2(&tempphi, xelem, ncoeff);
    deallocator2(&k1, xelem, ncoeff);
    deallocator2(&k2, xelem, ncoeff);
    deallocator2(&k3, xelem, ncoeff);
    deallocator2(&k4, xelem, ncoeff);
}
