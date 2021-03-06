/***************************************************************************

Author: nsaini
Created: 2018-07-24

***************************************************************************/


#include "common.h"
#include "memory.h"
#include "polylib.h"
#include "functions.h"

void initializeVel(struct elemsclr elem, double *x)
{
    //------------------------------------------------------------------------//
    //Temporary Variables
    int ielem;
    int igauss;
    int icoeff;

    double *us;
    allocator1(&us, tgauss);

    double *basis;
    allocator1(&basis, ncoeff);

    double **vand;
    allocator2(&vand, tgauss, ncoeff);
    //Fill up the Vandermonde matrix
    for(igauss=0; igauss<tgauss; igauss++)
    {
	//Get the basis vector
	basis1D(zeta[igauss], basis);
	for(icoeff=0; icoeff< ncoeff; icoeff++)
	{
	    vand[igauss][icoeff] = basis[icoeff];
	}
    }

    double *xs;
    allocator1(&xs, tgauss);
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Loop over the elements
    for(ielem=0; ielem <xelem; ielem++)
    {
	//Convert natural coordinates to physical space
	naturalToCartesian(xs,x,ielem);

	for(igauss=0; igauss<tgauss; igauss++)
	{
	    if(case_tog == 1)
	    {
		us[igauss] = 2.0*PI;
	    }
	    else if(case_tog == 2 || case_tog == 3)
	    {
		us[igauss] = 1.0;
	    }
	    
	}

	solveSystem(vand, us, elem.u[ielem]);
    }
    
    
    //------------------------------------------------------------------------//
    //Deallocators
    deallocator1(&us, tgauss);
    deallocator1(&basis, ncoeff);
    deallocator2(&vand, ncoeff, ncoeff);
    deallocator1(&xs, tgauss);
    //------------------------------------------------------------------------//

}


void initializeLS(struct elemsclr elem, double *x)
{
    //------------------------------------------------------------------------//
    //Temporary Variables
    int ielem;
    int igauss;
    int icoeff;

    double *phis;
    allocator1(&phis, tgauss);

    double *basis;
    allocator1(&basis, ncoeff);

    double **vand;
    allocator2(&vand, tgauss, ncoeff);
    //Fill up the Vandermonde matrix
    for(igauss=0; igauss<tgauss; igauss++)
    {
	//Get the basis vector
	basis1D(zeta[igauss], basis);
	for(icoeff=0; icoeff< ncoeff; icoeff++)
	{
	    vand[igauss][icoeff] = basis[icoeff];
	}
    }

    double *xs;
    allocator1(&xs, tgauss);
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Loop over the elements
    for(ielem=0; ielem <xelem; ielem++)
    {
	//Convert natural coordinates to physical space
	naturalToCartesian(xs,x,ielem);

	for(igauss=0; igauss<tgauss; igauss++)
	{
	    if(case_tog == 1)
	    {
		phis[igauss] = sin(xs[igauss]);
	    }
	    else if(case_tog == 2)
	    {
		phis[igauss] = exp(-200.0*pow(xs[igauss]-xb_in,2.0));
	    }
	    else if(case_tog == 3)
	    {
		if(xs[igauss] <= 0.6)
		{
		    phis[igauss] = exp(-200.0*pow(xs[igauss]-xb_in,2.0));
		}
		else if(xs[igauss] >= 0.6 && xs[igauss] <= 0.8)
		{
		    phis[igauss] = 1.0;
		}
		else
		{
		    phis[igauss] = 0.0;
		}
	    }
	}

	solveSystem(vand, phis, elem.phi[ielem]);
    }
    
    
    //------------------------------------------------------------------------//
    //Deallocators
    deallocator1(&phis, tgauss);
    deallocator1(&basis, ncoeff);
    deallocator2(&vand, ncoeff, ncoeff);
    deallocator1(&xs, tgauss);
    //------------------------------------------------------------------------//

}


void initialize(struct elemsclr elem, double *x)
{
    //------------------------------------------------------------------------//
    //Initialize Velocity
    initializeVel(elem, x);
    level_setBC(elem.u);
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Initialize scalar
    initializeLS(elem, x);
    //Apply BC
    level_setBC(elem.phi);

   
    //------------------------------------------------------------------------//

}
