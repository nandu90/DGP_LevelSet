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

void euler(double ***scalar, double ****mass, double ***rhs, double deltat)
{
    //------------------------------------------------------------------------//
    //Loop Indexes
    int ielem, jelem;
    int icoeff;
    //------------------------------------------------------------------------//

    for(ielem = 1; ielem < xelem-1; ielem++)
    {
	for(jelem = 1; jelem < yelem-1; jelem++)
	{
	    for(icoeff=0; icoeff < ncoeff; icoeff++)
	    {
		scalar[ielem][jelem][icoeff] += deltat*rhs[ielem][jelem][icoeff]/mass[ielem][jelem][icoeff][icoeff];
	    }
	}
    }
}


void Runge_Kutta(struct elemsclr elem, double **x, double **y, double deltat, double ***rhs)
{
    
    //------------------------------------------------------------------------//
    //Loop Indexes
    int ielem, jelem, icoeff;
    //------------------------------------------------------------------------//

    
    if(RKstages == 1)
    {
	//Get the Right hand side
	getRHS(elem, x, y, rhs);
	//Forward Euler
	euler(elem.phi, elem.mass, rhs, deltat);	
	//Apply boundary conditions
	level_setBC(elem.phi, elem.iBC);
    }

    else if(RKstages == 2)
    {
	double ***temp;
	allocator3(&temp, xelem, yelem, ncoeff);

	
	for(ielem=2; ielem<xelem-2; ielem++)
	{
	    for(jelem=2; jelem<yelem-2; jelem++)
	    {
		for(icoeff=0; icoeff<ncoeff; icoeff++)
		{
		    temp[ielem][jelem][icoeff]= elem.phi[ielem][jelem][icoeff];
		}
	    }
	}
	//Get the Right hand side
	getRHS(elem, x, y, rhs);
	//Forward Euler
	euler(elem.phi, elem.mass, rhs, deltat);	
	//Apply boundary conditions
	level_setBC(elem.phi, elem.iBC);

	//Get the Right hand side
	getRHS(elem, x, y, rhs);
	//Forward Euler
	euler(elem.phi, elem.mass, rhs, deltat);
	
	for(ielem=2; ielem<xelem-2; ielem++)
	{
	    for(jelem=2; jelem<yelem-2; jelem++)
	    {
		for(icoeff=0; icoeff<ncoeff; icoeff++)
		{
		    elem.phi[ielem][jelem][icoeff] = 0.5*temp[ielem][jelem][icoeff] + 0.5*elem.phi[ielem][jelem][icoeff];
		}
	    }
	}
	
	//Apply boundary conditions
	level_setBC(elem.phi, elem.iBC);
	
	deallocator3(&temp, xelem, yelem, ncoeff);
    }

    else if(RKstages == 3)
    {
	double ***temp;
	allocator3(&temp, xelem, yelem, ncoeff);

	for(ielem=2; ielem<xelem-2; ielem++)
	{
	    for(jelem=2; jelem<yelem-2; jelem++)
	    {
		for(icoeff=0; icoeff<ncoeff; icoeff++)
		{
		    temp[ielem][jelem][icoeff] = elem.phi[ielem][jelem][icoeff];
		}
	    }
	}

	//Get the Right hand side
	getRHS(elem, x, y, rhs);
	//Forward Euler
	euler(elem.phi, elem.mass, rhs, deltat);	
	//Apply boundary conditions
	level_setBC(elem.phi, elem.iBC);

	//Get the Right hand side
	getRHS(elem, x, y, rhs);
	//Forward Euler
	euler(elem.phi, elem.mass, rhs, deltat);	
	//Apply boundary conditions
	level_setBC(elem.phi, elem.iBC);

	for(ielem=2; ielem<xelem-2; ielem++)
	{
	    for(jelem=2; jelem<yelem-2; jelem++)
	    {
		for(icoeff=0; icoeff<ncoeff; icoeff++)
		{
		    elem.phi[ielem][jelem][icoeff] = (3.0/4.0)*temp[ielem][jelem][icoeff] + (1.0/4.0)*elem.phi[ielem][jelem][icoeff];
		}
	    }
	}

	//Apply boundary conditions
	level_setBC(elem.phi, elem.iBC);

	//Get the Right hand side
	getRHS(elem, x, y, rhs);
	//Forward Euler
	euler(elem.phi, elem.mass, rhs, deltat);
	for(ielem=2; ielem<xelem-2; ielem++)
	{
	    for(jelem=2; jelem<yelem-2; jelem++)
	    {
		for(icoeff=0; icoeff<ncoeff; icoeff++)
		{
		    elem.phi[ielem][jelem][icoeff] = (1.0/3.0)*temp[ielem][jelem][icoeff] + (2.0/3.0)*elem.phi[ielem][jelem][icoeff];
		}
	    }
	}

	//Apply boundary conditions
	level_setBC(elem.phi, elem.iBC);
	
	deallocator3(&temp, xelem, yelem, ncoeff);
    }
    else
    {
	if(myrank==master)
	{
	    printf("%d Stage Runge-Kutta not available.\nExiting...",RKstages);
	    exit(1);
	}
    }


    //------------------------------------------------------------------------//
    //Deallocators
    
    //------------------------------------------------------------------------//

}
