/***************************************************************************

Author: nsaini
Created: 2018-04-19

***************************************************************************/


#include "common.h"
#include "solvers.h"
#include "rhs.h"
#include "icbc.h"
#include "commu.h"
#include "memory.h"
#include "DGPFunc.h"

void euler(double ***scalar, double ****mass, double ***rhs, double deltat)
{
    //------------------------------------------------------------------------//
    //Loop Indexes
    int ielem, jelem;
    int icoeff;
    //------------------------------------------------------------------------//

    double *inc;
    allocator1(&inc, ncoeff);
    
    for(ielem = 1; ielem < xelem-1; ielem++)
    {
	for(jelem = 1; jelem < yelem-1; jelem++)
	{
	    solveSystem(mass[ielem][jelem], rhs[ielem][jelem], inc, ncoeff, ncoeff);
	    for(icoeff=0; icoeff < ncoeff; icoeff++)
	    {
		scalar[ielem][jelem][icoeff] += deltat*inc[icoeff];//rhs[ielem][jelem][icoeff]/mass[ielem][jelem][icoeff][icoeff];
	    }
	}
    }
}

void eulerIncrement(double ***inc, double ****mass, double ***rhs, double deltat)
{
    //------------------------------------------------------------------------//
    //Loop indexes
    int ielem, jelem;
    int icoeff;
    //------------------------------------------------------------------------//

    for(ielem = 1; ielem<xelem-1; ielem++)
    {
	for(jelem = 1; jelem<yelem-1; jelem++)
	{
	    solveSystem(mass[ielem][jelem], rhs[ielem][jelem], inc[ielem][jelem], ncoeff, ncoeff);
	    for(icoeff = 0; icoeff< ncoeff; icoeff++)
	    {
		inc[ielem][jelem][icoeff] *= deltat;//*rhs[ielem][jelem][icoeff]/mass[ielem][jelem][icoeff][icoeff];
	    }
	}
    }
}

void Runge_Kutta(struct elemsclr elem, double **x, double **y, double deltat, double ***rhs, double ****area)
{
    
    //------------------------------------------------------------------------//
    //Loop Indexes
    int ielem, jelem, icoeff;
    //------------------------------------------------------------------------//

    
    
    if(RKstages == 1)
    {
	//Get the Right hand side
	getRHS(elem, x, y, rhs, area);
	//Forward Euler
	euler(elem.phi, elem.mass, rhs, deltat);
	//Apply slope limiter
	if(limit == 1)
	{
	    momentLimiter(elem.phi);
	}
	//Apply boundary conditions
	level_setBC(elem.phi, elem.iBC);
    }

    else if(RKstages == 2)
    {
	double ***temp;
	allocator3(&temp, xelem, yelem, ncoeff);

	
	for(ielem=1; ielem<xelem-1; ielem++)
	{
	    for(jelem=1; jelem<yelem-1; jelem++)
	    {
		for(icoeff=0; icoeff<ncoeff; icoeff++)
		{
		    temp[ielem][jelem][icoeff]= elem.phi[ielem][jelem][icoeff];
		}
	    }
	}
	//Get the Right hand side
	getRHS(elem, x, y, rhs, area);
	//Forward Euler
	euler(elem.phi, elem.mass, rhs, deltat);
	//Apply slope limiter
	if(limit == 1)
	{
	    momentLimiter(elem.phi);
	}
	//Apply boundary conditions
	level_setBC(elem.phi, elem.iBC);

	//Get the Right hand side
	getRHS(elem, x, y, rhs, area);
	//Forward Euler
	euler(elem.phi, elem.mass, rhs, deltat);
	
	for(ielem=1; ielem<xelem-1; ielem++)
	{
	    for(jelem=1; jelem<yelem-1; jelem++)
	    {
		for(icoeff=0; icoeff<ncoeff; icoeff++)
		{
		    elem.phi[ielem][jelem][icoeff] = 0.5*temp[ielem][jelem][icoeff] + 0.5*elem.phi[ielem][jelem][icoeff];
		}
	    }
	}
	//Apply slope limiter
	if(limit == 1)
	{
	    momentLimiter(elem.phi);
	}
	//Apply boundary conditions
	level_setBC(elem.phi, elem.iBC);
	
	deallocator3(&temp, xelem, yelem, ncoeff);
    }

    else if(RKstages == 3)
    {
	double ***temp;
	allocator3(&temp, xelem, yelem, ncoeff);

	for(ielem=1; ielem<xelem-1; ielem++)
	{
	    for(jelem=1; jelem<yelem-1; jelem++)
	    {
		for(icoeff=0; icoeff<ncoeff; icoeff++)
		{
		    temp[ielem][jelem][icoeff] = elem.phi[ielem][jelem][icoeff];
		}
	    }
	}

	//Get the Right hand side
	getRHS(elem, x, y, rhs, area);
	//Forward Euler
	euler(elem.phi, elem.mass, rhs, deltat);
	//Apply slope limiter
	if(limit == 1)
	{
	    momentLimiter(elem.phi);
	}
	//Apply boundary conditions
	level_setBC(elem.phi, elem.iBC);

	//Get the Right hand side
	getRHS(elem, x, y, rhs, area);
	//Forward Euler
	euler(elem.phi, elem.mass, rhs, deltat);
	//Apply slope limiter
	if(limit == 1)
	{
	    momentLimiter(elem.phi);
	}
	//Apply boundary conditions
	level_setBC(elem.phi, elem.iBC);

	for(ielem=1; ielem<xelem-1; ielem++)
	{
	    for(jelem=1; jelem<yelem-1; jelem++)
	    {
		for(icoeff=0; icoeff<ncoeff; icoeff++)
		{
		    elem.phi[ielem][jelem][icoeff] = (3.0/4.0)*temp[ielem][jelem][icoeff] + (1.0/4.0)*elem.phi[ielem][jelem][icoeff];
		}
	    }
	}
	//Apply slope limiter
	if(limit == 1)
	{
	    momentLimiter(elem.phi);
	}
	//Apply boundary conditions
	level_setBC(elem.phi, elem.iBC);

	//Get the Right hand side
	getRHS(elem, x, y, rhs, area);
	//Forward Euler
	euler(elem.phi, elem.mass, rhs, deltat);
	for(ielem=1; ielem<xelem-1; ielem++)
	{
	    for(jelem=1; jelem<yelem-1; jelem++)
	    {
		for(icoeff=0; icoeff<ncoeff; icoeff++)
		{
		    elem.phi[ielem][jelem][icoeff] = (1.0/3.0)*temp[ielem][jelem][icoeff] + (2.0/3.0)*elem.phi[ielem][jelem][icoeff];
		}
	    }
	}
	//Apply slope limiter
	if(limit == 1)
	{
	    momentLimiter(elem.phi);
	}
	//Apply boundary conditions
	level_setBC(elem.phi, elem.iBC);
	
	deallocator3(&temp, xelem, yelem, ncoeff);
    }
    else if(RKstages == 4)
    {
	double ***k1, ***k2, ***k3, ***k4;
	allocator3(&k1, xelem, yelem, ncoeff);
	allocator3(&k2, xelem, yelem, ncoeff);
	allocator3(&k3, xelem, yelem, ncoeff);
	allocator3(&k4, xelem, yelem, ncoeff);

	//------------------------------------------------------------------------//
	//Store the primary value in a temp array
	double ***tempphi;
	allocator3(&tempphi, xelem, yelem, ncoeff);
	for(ielem =0; ielem<xelem; ielem++)
	{
	    for(jelem=0; jelem<yelem; jelem++)
	    {
		for(icoeff=0; icoeff<ncoeff; icoeff++)
		{
		    tempphi[ielem][jelem][icoeff] = elem.phi[ielem][jelem][icoeff];
		}
	    }
	}
	//------------------------------------------------------------------------//

	/*for(icoeff =0; icoeff<ncoeff; icoeff++)
	{
	    printf("%.4e ",elem.phi[2][2][icoeff]);
	}
	printf("\n\n");*/
	
	//------------------------------------------------------------------------//
	//1st increment
	getRHS(elem, x, y, rhs, area);
	eulerIncrement(k1, elem.mass, rhs, deltat);
	for(ielem=1; ielem<xelem-1; ielem++)
	{
	    for(jelem=1; jelem<yelem-1; jelem++)
	    {
		for(icoeff=0; icoeff<ncoeff; icoeff++)
		{
		    elem.phi[ielem][jelem][icoeff] = tempphi[ielem][jelem][icoeff] + 0.5*k1[ielem][jelem][icoeff];
		}
	    }
	}
	//Apply slope limiter
	if(limit == 1)
	{
	    momentLimiter(elem.phi);
	}
	//Apply BC
	level_setBC(elem.phi, elem.iBC);
	//------------------------------------------------------------------------//

	/*for(icoeff =0; icoeff<ncoeff; icoeff++)
	{
	    printf("%.4e ",elem.phi[2][2][icoeff]);
	}
	printf("\n\n");*/
	
	//------------------------------------------------------------------------//
	//2nd increment
	getRHS(elem, x, y, rhs, area);
	eulerIncrement(k2, elem.mass, rhs, deltat);
	for(ielem=1; ielem<xelem-1; ielem++)
	{
	    for(jelem=1; jelem<yelem-1; jelem++)
	    {
		for(icoeff=0; icoeff<ncoeff; icoeff++)
		{
		    elem.phi[ielem][jelem][icoeff] = tempphi[ielem][jelem][icoeff] + 0.5*k2[ielem][jelem][icoeff];
		}
	    }
	}
	//Apply slope limiter
	if(limit == 1)
	{
	    momentLimiter(elem.phi);
	}
	//Apply BC
	level_setBC(elem.phi, elem.iBC);
	//------------------------------------------------------------------------//

	/*for(icoeff =0; icoeff<ncoeff; icoeff++)
	{
	    printf("%.4e ",elem.phi[2][2][icoeff]);
	}
	printf("\n\n");*/
	
	//------------------------------------------------------------------------//
	//3rd Increment
	getRHS(elem, x, y, rhs, area);
	eulerIncrement(k3, elem.mass, rhs, deltat);
	for(ielem=1; ielem<xelem-1; ielem++)
	{
	    for(jelem=1; jelem<yelem-1; jelem++)
	    {
		for(icoeff=0; icoeff<ncoeff; icoeff++)
		{
		    elem.phi[ielem][jelem][icoeff] = tempphi[ielem][jelem][icoeff] + k3[ielem][jelem][icoeff];
		}
	    }
	}
	//Apply slope limiter
	if(limit == 1)
	{
	    momentLimiter(elem.phi);
	}
	//Apply BC
	level_setBC(elem.phi, elem.iBC);
	//------------------------------------------------------------------------//

	/*for(icoeff =0; icoeff<ncoeff; icoeff++)
	{
	    printf("%.4e ",elem.phi[2][2][icoeff]);
	}
	printf("\n\n");*/
	
	//------------------------------------------------------------------------//
	//4th Increment
	getRHS(elem, x, y, rhs, area);
	eulerIncrement(k4, elem.mass, rhs, deltat);
	//------------------------------------------------------------------------//

	//------------------------------------------------------------------------//
	//Now get the next time step value
	for(ielem=1; ielem<xelem-1; ielem++)
	{
	    for(jelem=1; jelem<yelem-1; jelem++)
	    {
		for(icoeff=0; icoeff<ncoeff; icoeff++)
		{
		    elem.phi[ielem][jelem][icoeff] = tempphi[ielem][jelem][icoeff] + (1.0/6.0)*(k1[ielem][jelem][icoeff] + 2.0*k2[ielem][jelem][icoeff] + 2.0*k3[ielem][jelem][icoeff] + k4[ielem][jelem][icoeff]);
		}
	    }
	}
	//Apply slope limiter
	if(limit == 1)
	{
	    momentLimiter(elem.phi);
	}
	//Apply BC
	level_setBC(elem.phi, elem.iBC);
	//------------------------------------------------------------------------//

	/*for(icoeff =0; icoeff<ncoeff; icoeff++)
	{
	    printf("%.4e ",elem.phi[2][2][icoeff]);
	}
	printf("\n\n");*/

	//exit(1);
	deallocator3(&tempphi, xelem, yelem, ncoeff);
	deallocator3(&k1, xelem, yelem, ncoeff);
	deallocator3(&k2, xelem, yelem, ncoeff);
	deallocator3(&k3, xelem, yelem, ncoeff);
	deallocator3(&k4, xelem, yelem, ncoeff);
    }
    else
    {
	if(myrank==master)
	{
	    printf("%d Stage Runge-Kutta not available.\nExiting...",RKstages);
	    exit(1);
	}
    }


}
