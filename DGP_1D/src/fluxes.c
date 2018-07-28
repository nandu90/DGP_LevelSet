/***************************************************************************

Author: nsaini
Created: 2018-07-25

***************************************************************************/

#include "common.h"
#include "memory.h"
#include "functions.h"

void fluxes(double *flux, double *x, struct elemsclr elem)
{
    //------------------------------------------------------------------------//
    //Loop Indexes
    int ielem, igauss, icoeff;
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Temporary Variables
    double recphi, recu;

    double Rflux, Lflux;
    double Ru, Lu;

    double *basis;
    allocator1(&basis, ncoeff);
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    for(ielem=0; ielem<xelem-1; ielem++)
    {
	//Reconstruct solution at the left side of face
	basis1D(1.0, basis);
	recphi = 0.0;
	recu = 0.0;
	
	for(icoeff=0; icoeff<ncoeff; icoeff++)
	{
	    recphi += basis[icoeff]*elem.phi[ielem][icoeff];
	    recu += basis[icoeff]*elem.u[ielem][icoeff];
	}
	Lflux = recphi*recu;
	Lu = recu;

	//Recontruct the solution at the right side of face
	basis1D(-1.0, basis);
	recphi = 0.0;
	recu = 0.0;
	
	for(icoeff=0; icoeff<ncoeff; icoeff++)
	{
	    recphi += basis[icoeff]*elem.phi[ielem+1][icoeff];
	    recu += basis[icoeff]*elem.u[ielem+1][icoeff];
	}
	Rflux = recphi*recu;
	Ru = recu;

	//Do upwinding
	flux[ielem] = upwind(Lflux, Rflux, Lu, Ru);

	/*if(ielem == xelem/2)
	{
	    printf("%.4e %.4e %.4e %.4e\n",Lflux, Lu, Rflux, Ru);
	    printf("%.4e\n",flux);
	    exit(1);
	    }*/
    }

    //------------------------------------------------------------------------//
    //Deallocators
    deallocator1(&basis, ncoeff);
    //------------------------------------------------------------------------//

}

double upwind(double minusFlux, double plusFlux, double minusU, double plusU)
{
    double Fminus;
    double Fplus;
    double flux;
    
    //Upwinding
    if(minusU >= 0.0)
    {
	Fminus = minusFlux;
    }
    else
    {
	Fminus = 0.0;
    }
    
    if(plusU > 0.0)
    {
	Fplus = 0.0;
    }
    else
    {
	Fplus = plusFlux;
    }

    flux = Fminus + Fplus;

    return flux;
}


void boundaryIntegral(double **rhs, double *flux, double *x)
{
    int ielem;
    int icoeff;

    double *basis;
    allocator1(&basis, ncoeff);

    double *rightInt, *leftInt;
    allocator1(&rightInt, ncoeff);
    allocator1(&leftInt, ncoeff);

    /*ielem = xelem/2;
    if(ielem==xelem/2)
    {
	printf("boundary\n");
	for(icoeff=0; icoeff<ncoeff; icoeff++)
	{
	    printf("%.4e %.4e %.4e\n", rightInt[icoeff], leftInt[icoeff], rhs[ielem][icoeff]);
	}
	}*/
    for(ielem=1; ielem<xelem-1; ielem++)
    {
	basis1D(1.0, basis);
	for(icoeff=0; icoeff<ncoeff; icoeff++)
	{	    
	    rightInt[icoeff] = basis[icoeff]*flux[ielem];
	}

	basis1D(-1.0, basis);
	for(icoeff=0; icoeff<ncoeff; icoeff++)
	{
	    leftInt[icoeff] = basis[icoeff]*flux[ielem-1];
	}

	for(icoeff=0; icoeff<ncoeff; icoeff++)
	{
	    rhs[ielem][icoeff] += -(rightInt[icoeff] - leftInt[icoeff]);
	}

	/*if(ielem==xelem/2)
	{
	    printf("boundary\n");
	    for(icoeff=0; icoeff<ncoeff; icoeff++)
	    {
		printf("%.4e %.4e %.4e\n", rightInt[icoeff], leftInt[icoeff], rhs[ielem][icoeff]);
	    }
	    exit(1);
	    }*/
    }

    deallocator1(&basis, ncoeff);
    deallocator1(&rightInt, ncoeff);
    deallocator1(&leftInt, ncoeff);
}
