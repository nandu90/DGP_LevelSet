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
    int maxiter = 1000;
    int print_gap = 10;
    case_tog = 4;        //1-Sine; 2-Gaussian

    polyorder = 1;
    xelem = 10;
    deltat = 0.001;

    if(case_tog == 1)
    {
	xlen = 2*PI;
	xb_in = 75.0;
    }
    else if(case_tog == 2)
    {
	xlen = 1.0;
	xb_in = 0.5;
    }
    else if(case_tog == 3)
    {
	xlen = 1.0;
	xb_in = 0.3;
    }
    else if(case_tog == 4)
    {
	xlen = 1.0;
    }

    double Pe = 100.0;
    double f = 0.0;
    
    double u = 1.0;
    double deltax = xlen/((double)xelem);
    double k = deltax*u/Pe;

    double leftbc = 0.0;
    double rightbc = 1.0;
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Temporary Variables
    int ielem;
    //int iter;
    //int icoeff;
    int igauss;
    int inode;
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Variables as result of input
    tgauss = polyorder + 1;
    ncoeff = polyorder + 1;
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //COnstruct grid    
    //xelem += 2;
    xnode = xelem + 1;
    double *x;
    allocator1(&x,xelem+1);

    for(inode = 0; inode < xnode; inode++)
    {
	x[inode] = (inode)*deltax;
    }
    //Construct ghost cells
    //x[0] = x[1] - deltax;
    //x[xelem] = x[xelem-1] + deltax;   
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Get the Gauss quadrature points and weights
    allocator1(&zeta, tgauss);
    allocator1(&weights, tgauss);

    zwgl(zeta, weights, tgauss);

    printf("The Gauss points and weights are:\n");
    for(igauss=0; igauss<tgauss; igauss++)
    {
	printf("%.8e %.8e\n",zeta[igauss],weights[igauss]);
    }
    printf("\n");
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Contruct nodal variables      
    double *phi;
    allocator1(&phi, xnode);

    //Allocate Global Matrices
    double **M, **C, **S, *F;
    allocator2(&M, xnode, xnode);
    allocator2(&C, xnode, xnode);
    allocator2(&S, xnode, xnode);
    allocator1(&F, xnode);
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Initialize Coefficients
    //Initialize to 0
    //initialize(phi, x);

    double *iniphi;
    allocator1(&iniphi, xnode);
    for(inode=0; inode<xnode; inode++)
    {
	iniphi[inode] = phi[inode];
    }
    //------------------------------------------------------------------------//
    
    //------------------------------------------------------------------------//
    //Get the mass Matrix
    massmatrix(M, x);

    //Get the Convection Matrix
    convection(C, x);

    //Get the Stiffness Matrix
    stiffness(S, x);

    //Get the force vector
    forceVector(F, x);
    //------------------------------------------------------------------------//
   
    //------------------------------------------------------------------------//
    output(phi, x,0);
    //Steady state solution
    int inode2;

    
    
    double **A;
    allocator2(&A, xnode, xnode);
    for(inode=0; inode<xnode; inode++)
    {
	for(inode2=0; inode2<xnode; inode2++)
	{
	    A[inode][inode2] += (1.0/Pe)*S[inode][inode2] + C[inode][inode2];
	}
	F[inode] *= f;
    }
    //------------------------------------------------------------------------//
    //Apply BC
    //left
    A[0][0] = 1.0;
    for(inode=1; inode<xnode; inode++)
    {
	A[0][inode] = 0.0;
    }
    F[0] = leftbc;

    //right
    A[xnode-1][xnode-1] = 1.0;
    for(inode=0; inode<xnode-1; inode++)
    {
	A[xnode-1][inode] = 0.0;
    }
    F[xnode-1] = rightbc;
    
    //------------------------------------------------------------------------//

    for(inode=0; inode<xnode; inode++)
    {
	for(inode2=0; inode2<xnode; inode2++)
	{
	    printf("%.4f ",A[inode][inode2]);
	}
	printf("\n");
    }

    for(inode=0; inode<xnode; inode++)
    {
	printf("%.4f\n",F[inode]);
    }

    solveSystem(A, F, phi, xnode, xnode);
    
    printf("Result:\n");
    for(inode=0; inode<xnode; inode++)
    {
	printf("%d %.4f\n",inode+1,phi[inode]);
    }
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Print output
    output(phi, x,1);
    //------------------------------------------------------------------------//


    /*//------------------------------------------------------------------------//
    //Calculate Error norms
    double err1, lerr1;
    errorNormL1(iniphi, elem.phi, &err1, &lerr1, x);
    
    double err2, lerr2;
    errorNormL2(iniphi, elem.phi, &err2, &lerr2, x);
    
    printf("The L1 norm of error is %.4e and the Log of norm is %.4e\n", err1, lerr1);
    printf("The L2 norm of error is %.4e and the Log of norm is %.4e\n", err2, lerr2);
    //------------------------------------------------------------------------//*/

    //------------------------------------------------------------------------//
    //Deallocators
    deallocator1(&x, xnode);
    deallocator1(&phi, xnode);
    deallocator2(&M, xnode, xnode);
    deallocator2(&C, xnode, xnode);
    deallocator2(&S, xnode, xnode);
    deallocator1(&F, xnode);
    deallocator1(&iniphi, xnode);
    deallocator2(&A, xnode, xnode);
    //------------------------------------------------------------------------//

}
