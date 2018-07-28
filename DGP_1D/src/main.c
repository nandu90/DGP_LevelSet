/***************************************************************************

Author: nsaini
Created: 2018-07-24

***************************************************************************/

#include "memory.h"
#include "common.h"
#include "functions.h"
#include "polylib.h"

int main(int argc, char **argv)
{
    feenableexcept(FE_INVALID   | 
                   FE_DIVBYZERO | 
                   FE_OVERFLOW  | 
                   FE_UNDERFLOW);
    
    //------------------------------------------------------------------------//
    //INput Section
    polyorder = 3;
    xelem = 25;
    deltat = 0.0005;
    
    xlen = 2.0*PI;
    //xb_in = 75.0;
    int maxiter = 2000;
    int print_gap = 10;
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Temporary Variables
    int ielem;
    int iter;
    int icoeff;
    int igauss;
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Variables as result of input
    tgauss = polyorder + 1;
    ncoeff = polyorder + 1;
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //COnstruct grid
    double deltax = xlen/((double)xelem);
    xelem += 2;
    double *x;
    allocator1(&x,xelem+1);

    for(ielem = 1; ielem < xelem; ielem++)
    {
	x[ielem] = (ielem-1)*deltax;
    }
    //Construct ghost cells
    x[0] = x[1] - deltax;
    x[xelem] = x[xelem-1] + deltax;

    
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Get the Gauss quadrature points and weights
    allocator1(&zeta, tgauss);
    allocator1(&weights, tgauss);

    zwgll(zeta, weights, tgauss);

    printf("The Gauss points and weights are:\n");
    for(igauss=0; igauss<tgauss; igauss++)
    {
	printf("%.8e %.8e\n",zeta[igauss],weights[igauss]);
    }
    printf("\n");
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Contruct elemental variables      
    struct elemsclr elem;
    allocator2(&elem.u, xelem, ncoeff);
    allocator2(&elem.phi, xelem, ncoeff);
    allocator3(&elem.mass, xelem, ncoeff, ncoeff);
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Initialize Coefficients
    initialize(elem, x);

    double **iniphi;
    allocator2(&iniphi, xelem, ncoeff);
    for(ielem=0; ielem<xelem; ielem++)
    {
	for(icoeff=0; icoeff<ncoeff; icoeff++)
	{
	    iniphi[ielem][icoeff] = elem.phi[ielem][icoeff];
	}
    }
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Get the Gauss quadrature points and weights
    deallocator1(&zeta, tgauss);
    deallocator1(&weights, tgauss);
    
    tgauss = polyorder+10;
    allocator1(&zeta, tgauss);
    allocator1(&weights, tgauss);

    zwgll(zeta, weights, tgauss);

    printf("The Gauss points and weights are:\n");
    for(igauss=0; igauss<tgauss; igauss++)
    {
	printf("%.8e %.8e\n",zeta[igauss],weights[igauss]);
    }
    printf("\n");
    //------------------------------------------------------------------------//
   
    //------------------------------------------------------------------------//
    //Get the mass Matrix
    
    massmatrix(elem.mass, x);
    //------------------------------------------------------------------------//
    
    //------------------------------------------------------------------------//
    //Time loop
    double **rhs;
    allocator2(&rhs, xelem, ncoeff);

    int print_count = 0;

    double time = 0.0;
    for(iter=0; iter<maxiter; iter++)
    {
	
	Runge_Kutta(elem, x, deltat, rhs);
	
	time += deltat;
	//------------------------------------------------------------------------//
	//Print out the paraview output
	print_count++;
        if(print_count == 1)
        {
           output(elem, x,iter);
        }
        if(print_count == print_gap)
        {
            print_count = 0;
        }
	//------------------------------------------------------------------------//

	//errExact(elem.phi, x, time ,iter);
    }

    output(elem, x,iter);

    //------------------------------------------------------------------------//
    //Calculate Error norms
    double err1, lerr1;
    errorNormL1(iniphi, elem.phi, &err1, &lerr1, x);
    
    double err2, lerr2;
    errorNormL2(iniphi, elem.phi, &err2, &lerr2, x);
    
    printf("The L1 norm of error is %.4e and the Log of norm is %.4e\n", err1, lerr1);
    printf("The L2 norm of error is %.4e and the Log of norm is %.4e\n", err2, lerr2);
    
    
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Deallocators
    deallocator1(&x, xelem);
    deallocator2(&elem.u, xelem, ncoeff);
    deallocator2(&elem.phi, xelem, ncoeff);
    deallocator3(&elem.mass, xelem, ncoeff, ncoeff);
    deallocator2(&rhs, xelem, ncoeff);
    deallocator2(&iniphi, xelem, ncoeff);
    //------------------------------------------------------------------------//

}
