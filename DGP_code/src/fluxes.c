/***************************************************************************

Author: nsaini
Created: 2018-03-28

***************************************************************************/


#include "common.h"
#include "DGPFunc.h"
#include "rhs.h"
#include "polylib.h"
#include "memory.h"

void fluxes(double ***rflux, double ***tflux, double **x, double **y,  struct elemsclr elem)
{
    //------------------------------------------------------------------------//
    /*This routine will construct the fluxes at the right and top faces.
      Upwinding will also be performed here. 
      Boundary integrals will be calculated in a separate routine.
     */
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Loop indexes
    int ielem, jelem;
    int ixgauss, iygauss;
    int icoeff;
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Temporary variables
    double *zx, *zy;
    double *wx, *wy;
    allocator1(&zx, xgpts);
    allocator1(&zy, ygpts);
    allocator1(&wx, xgpts);
    allocator1(&wy, ygpts);
    //Get the quadrature points and weights
    if(xgpts == 1)
    {
	zx[0] = 0.0;
	wx[0] = 1.0;
    }
    else
    {
	zwgl(zx,wx,xgpts);
    }
    if(ygpts == 1)
    {
	zy[0] = 0.0;
	wy[0] = 1.0;
    }
    else
    {
	zwgl(zy,wy,ygpts);
    }

    double *basisx, *basisy;
    allocator1(&basisx, ncoeff);
    allocator1(&basisy, ncoeff);

    double recphi, recu, recv;
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    double Rflux, Lflux;
    double Tflux, Bflux;
    double Ru, Lu;
    double Tv, Bv;
    double normz1, normz2;
    double normalVel;
    //------------------------------------------------------------------------//


    //------------------------------------------------------------------------//
    //Reconstruct the face normal velocities and the fluxes at the Gauss Quadrature
    //points on the top and right face of each cell
    //Upwinding will be done later - so that its easier to change the flux scheme if i have to
    for(ielem=0; ielem<xelem-1; ielem++)
    {
	for(jelem=0; jelem<yelem-1; jelem++)
	{
	    //Loop over the Gauss Quadrature points on the right face of cell
	    for(iygauss=0; iygauss < ygpts; iygauss++)
	    {
		//Reconstruct the solution at the left side
		//Get the basis
		basis2D(1.0, zy[iygauss], basisy);		

		//get the face normal - deprecated
		//lineNormal(ielem, jelem, x, y, &normz1, &normz2, 1, 1.0 ,zy[iygauss]);
		
		recphi = 0.0;
		recu = 0.0;
		recv = 0.0;
		for(icoeff = 0; icoeff<ncoeff; icoeff++)
		{
		    recphi += basisy[icoeff]*elem.phi[ielem][jelem][icoeff];
		    recu += basisy[icoeff]*elem.u[ielem][jelem][icoeff];
		    recv += basisy[icoeff]*elem.v[ielem][jelem][icoeff];
		}
				
		//normalVel = recu*normz1 + recv*normz2;// - deprecated
		Lflux = recphi*recu;
		Lu = recu;
		//Lflux = recphi*normalVel;// - deprecated
		//Lu = normalVel;// - deprecated

		//Reconstruct the solution at the right side
		//Get the basis
		basis2D(-1.0, zy[iygauss], basisy);		

		//get the face normal - deprecated
		//lineNormal(ielem, jelem, x, y, &normz1, &normz2, 3, 1.0 ,zy[iygauss]);
		
		recphi = 0.0;
		recu = 0.0;
		recv = 0.0;
		for(icoeff = 0; icoeff<ncoeff; icoeff++)
		{
		    recphi += basisy[icoeff]*elem.phi[ielem+1][jelem][icoeff];
		    recu += basisy[icoeff]*elem.u[ielem+1][jelem][icoeff];
		    recv += basisy[icoeff]*elem.v[ielem+1][jelem][icoeff];
		}
		
	        //normalVel = recu*normz1 + recv*normz2;// - deprecated
		Rflux = recphi * recu;
		Ru = recu;
		//Rflux = recphi*normalVel;// - deprecated
		//Ru = normalVel;// - deprecated


		rflux[ielem][jelem][iygauss] = upwind(Lflux, Rflux, Lu, Ru);
	    }	    

	    
	    //Loop over the Gauss Quadrature points on the top face of the cell
	    for(ixgauss=0; ixgauss<xgpts; ixgauss++)
	    {
		//Recontruct the solution at the bottom
		//Get the basis
		basis2D(zx[ixgauss], 1.0, basisx);

		//get the face normal - deprecated
		//lineNormal(ielem, jelem, x, y, &normz1, &normz2, 4, zx[ixgauss] ,1.0);
		
		recphi = 0.0;
		recu = 0.0;
		recv = 0.0;
		for(icoeff = 0; icoeff<ncoeff; icoeff++)
		{
		    recphi += basisx[icoeff]*elem.phi[ielem][jelem][icoeff];
		    recu += basisx[icoeff]*elem.u[ielem][jelem][icoeff];
		    recv += basisx[icoeff]*elem.v[ielem][jelem][icoeff];
		}
		//normalVel = recu*normz1 + recv*normz2;// - deprecated
		Bflux = recphi * recv;
		Bv = recv;
		//Bflux = recphi*normalVel;// - deprecated
		//Bv = normalVel;// - deprecated


		//Recontruct the solution at the top
		//Get the basis
		basis2D(zx[ixgauss], -1.0, basisx);

		//get the face normal - deprecated
		//lineNormal(ielem, jelem, x, y, &normz1, &normz2, 2, zx[ixgauss] ,1.0);
		
		recphi = 0.0;
		recv = 0.0;
		recu = 0.0;
		for(icoeff = 0; icoeff<ncoeff; icoeff++)
		{
		    recphi += basisx[icoeff]*elem.phi[ielem][jelem+1][icoeff];
		    recu += basisx[icoeff]*elem.u[ielem][jelem+1][icoeff];
		    recv += basisx[icoeff]*elem.v[ielem][jelem+1][icoeff];
		}
		//normalVel = recu*normz1 + recv*normz2;// - deprecated
		Tflux = recphi * recv;
		Tv = recv;
		//Tflux = recphi*normalVel;// - deprecated
		//Tv = normalVel;// - dprecated

		tflux[ielem][jelem][ixgauss] = upwind(Bflux, Tflux, Bv, Tv);
	    }

	    //exit(1);
	}
    }
    //------------------------------------------------------------------------//


    //------------------------------------------------------------------------//
    //Deallocators
    deallocator1(&zx, xgpts);
    deallocator1(&zy, ygpts);
    deallocator1(&wx, xgpts);
    deallocator1(&wy, ygpts);
    deallocator1(&basisx, ncoeff);
    deallocator1(&basisy, ncoeff);
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



void boundaryIntegral(double ***rhs, double ***rflux, double ***tflux, double **x, double **y, double ****area)
{
    //------------------------------------------------------------------------//
    /*
      This routine will calculate the net boundary integral
     */
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Loop indexes
    int ielem,jelem;
    int ixgauss, iygauss;
    int icoeff;
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Temporary Variables
    
    double *zx, *zy;
    double *wx, *wy;
    allocator1(&zx, xgpts);
    allocator1(&zy, ygpts);
    allocator1(&wx, xgpts);
    allocator1(&wy, ygpts);
    //Get the quadrature points and weights
    if(xgpts == 1)
    {
	zx[0] = 0.0;
	wx[0] = 1.0;
    }
    else
    {
	zwgl(zx,wx,xgpts);
    }
    if(ygpts == 1)
    {
	zy[0] = 0.0;
	wy[0] = 1.0;
    }
    else
    {
	zwgl(zy,wy,ygpts);
    }

    /*printf("%d %d\n", xgpts, ygpts);
    printf("%.4e %.4e %.4e\n",zx[0], zx[1], zx[2]);
    printf("%.4e %.4e %.4e\n\n",wx[0], wx[1], wx[2]);
    printf("%.4e %.4e %.4e\n",zy[0], zy[1], zy[2]);
    printf("%.4e %.4e %.4e\n",wy[0], wy[1], wy[2]);
    exit(1);*/
    
    double *basisx, *basisy;
    allocator1(&basisx, ncoeff);
    allocator1(&basisy, ncoeff);

    //------------------------------------------------------------------------//
    double rintegral;
    double tintegral;
    double *rightInt;
    double *leftInt;
    double *topInt;
    double *bottomInt;

    allocator1(&rightInt, ncoeff);
    allocator1(&leftInt, ncoeff);
    allocator1(&topInt, ncoeff);
    allocator1(&bottomInt, ncoeff);


    //------------------------------------------------------------------------//
    //Loop over the elements
    for(ielem = 1; ielem<xelem-1; ielem++)
    {
	for(jelem = 1; jelem<yelem-1; jelem++)
	{
	    for(icoeff=0; icoeff < ncoeff; icoeff++)
	    {
		rightInt[icoeff] = 0.0;
		leftInt[icoeff] = 0.0;
		topInt[icoeff] = 0.0;
		bottomInt[icoeff] = 0.0;
	    }

	    for(icoeff=0; icoeff<ncoeff; icoeff++)
	    {
		//Loop over the quadrature points on the right face
		for(iygauss = 0; iygauss<ygpts; iygauss++)
		{
		    //Get the basis
		    basis2D(1.0, zy[iygauss], basisy);

		    if(polyorder == 0)
		    {
			rightInt[icoeff] += wy[iygauss]*basisy[icoeff]*rflux[ielem][jelem][iygauss] / area[ielem][jelem][0][0];
		    }
		    else
		    {
			rightInt[icoeff] += wy[iygauss]*basisy[icoeff]*rflux[ielem][jelem][iygauss] * lineJacobian(ielem, jelem, basisy[icoeff], x, y, 1);
		    }
		    
		}

		//Loop over the quadrature points on the left face
		for(iygauss = 0; iygauss<ygpts; iygauss++)
		{
		    //Get the basis
		    basis2D(-1.0, zy[iygauss], basisy);

		    if(polyorder == 0)
		    {
			leftInt[icoeff] += wy[iygauss]*basisy[icoeff]*rflux[ielem-1][jelem][iygauss] / area[ielem][jelem][0][0];
		    }
		    else
		    {
			leftInt[icoeff] += wy[iygauss]*basisy[icoeff]*rflux[ielem-1][jelem][iygauss] * lineJacobian(ielem, jelem, basisy[icoeff], x, y, 3);
		    }
		}
		
		
		
		//Loop over the quadrature points on the top face
		for(ixgauss=0; ixgauss<xgpts; ixgauss++)
		{
		    //Get the basis
		    basis2D(zx[ixgauss], 1.0, basisx);

		    if(polyorder == 0)
		    {
		    //Sum to the integral		    
			topInt[icoeff] += wx[ixgauss]*basisx[icoeff]*tflux[ielem][jelem][ixgauss] / area[ielem][jelem][1][1];
		    }
		    else
		    {
			topInt[icoeff] += wx[ixgauss]*basisx[icoeff]*tflux[ielem][jelem][ixgauss] * lineJacobian(ielem, jelem, basisx[icoeff], x, y, 2);
		    }
		    
		}

		//Loop over the quadrature points on the bottom face
		for(ixgauss=0; ixgauss<xgpts; ixgauss++)
		{
		    //Get the basis
		    basis2D(zx[ixgauss], -1.0, basisx);
		    
		    //Sum to the integral
		    if(polyorder == 0)
		    {
			bottomInt[icoeff] += wx[ixgauss]*basisx[icoeff]*tflux[ielem][jelem-1][ixgauss]  / area[ielem][jelem][1][1];
		    }
		    else
		    {
			bottomInt[icoeff] += wx[ixgauss]*basisx[icoeff]*tflux[ielem][jelem-1][ixgauss] * lineJacobian(ielem, jelem, basisx[icoeff], x, y, 4);
		    }
		    
		}
	    }
	    
	    for(icoeff=0; icoeff<ncoeff; icoeff++)
	    {
		rintegral =  rightInt[icoeff] - leftInt[icoeff];
		tintegral = topInt[icoeff] - bottomInt[icoeff];

		rhs[ielem][jelem][icoeff] += -(rintegral + tintegral);
	    }
	}
    }
    //------------------------------------------------------------------------//


    //------------------------------------------------------------------------//
    //Check
    /*ielem = 2;
    jelem = 2;
    printf("Boundary \n");
    for(icoeff=0; icoeff<ncoeff; icoeff++)
    {
	printf("%.4e ",integral[ielem][jelem][icoeff]);
    }
    printf("\n");*/
    

    
    /*for(ielem = 2; ielem<xelem-2; ielem++)
    {
	for(jelem =2; jelem<yelem-2; jelem++)
	{
	    printf("%d %d ", ielem, jelem);
	    for(iygauss=0; iygauss<ygpts; iygauss++)
	    {
		//if((integral[ielem][jelem][icoeff]) < 0.0)exit(1);
		printf("%.4f %.4f %.4f ", rflux[ielem][jelem][iygauss], rflux[ielem-1][jelem][iygauss], rflux[ielem][jelem][iygauss]- rflux[ielem-1][jelem][iygauss]);
	    }
	    printf("\n");
	}
    }

    for(ielem = 2; ielem<xelem-2; ielem++)
    {
	for(jelem =2; jelem<yelem-2; jelem++)
	{
	    printf("%d %d ", ielem, jelem);
	    for(icoeff=0; icoeff<ncoeff; icoeff++)
	    {
		//if((integral[ielem][jelem][icoeff]) < 0.0)exit(1);
		printf("%.4f ",integral[ielem][jelem][icoeff]);
	    }
	    printf("\n");
	}
    }
    exit(1);*/
    //------------------------------------------------------------------------//
    
    //------------------------------------------------------------------------//
    //Deallocators
    deallocator1(&zx, xgpts);
    deallocator1(&zy, ygpts);
    deallocator1(&wx, xgpts);
    deallocator1(&wy, ygpts);
    deallocator1(&basisx, ncoeff);
    deallocator1(&basisy, ncoeff);
    deallocator1(&rightInt, ncoeff);
    deallocator1(&leftInt, ncoeff);
    deallocator1(&topInt, ncoeff);
    deallocator1(&bottomInt, ncoeff);
    //------------------------------------------------------------------------//

    
}


