/***************************************************************************

Author: nsaini
Created: 2018-03-11

***************************************************************************/


#include "solvers.h"
#include "common.h"
#include "memory.h"

void solveSystem(double **vand, double *rhs, double *soln)
{
    int i,j,k;
    //------------------------------------------------------------------------//
    //Set up variables for LAPACK
    double *A;
    allocator1(&A, tgauss*ncoeff);

    k=0;
    for(i=0; i<ncoeff; i++)
    {
	for(j=0; j<tgauss; j++)
	{
	    A[k] = vand[j][i];  //Be carefule over here. Store Column-wise for correct result
	    k++;
	}
    }

    int M = tgauss; //Rows of matrix
    int N = ncoeff;        //Columns of matrix
    int NRHS = 1;
    int LDA = M;
    int *IPIV;
    iallocator1(&IPIV, M);

    double *B;
    allocator1(&B, M);

    for(i=0; i<tgauss; i++)
    {
	B[i] = rhs[i];
    }

    int LDB = M;
    int INFO;

    int LWORK = -1;
    double WORK[1];
    char TRANS = 'N';
    
    //------------------------------------------------------------------------//
    //Solve the system to get the coefficients
    if(quadtype == 2)
    {
	dgesv_(&N, &NRHS, A, &LDA, IPIV, B, &LDB, &INFO);
    }
    else if(quadtype == 1)
    {
	dgels_(&TRANS,&M,&N,&NRHS,A,&LDA,B,&LDB,WORK,&LWORK,&INFO);
    }

    //Assign to the soln array
    for(i=0; i<ncoeff; i++)
    {
	soln[i] = B[i];
    }
    //------------------------------------------------------------------------//


    /*//------------------------------------------------------------------------//
    //check
    printf("The Vandermonde matrix is\n");
    for(i=0; i<M*N; i++)
    {
	printf("%d %.4f\n",i, temp[i]);
    }
    printf("\n\n");

    printf("The solution vector is\n");
    for(i=0; i<N; i++)
    {
	printf("%.4f\n", B[i]);
    }

    printf("The rhs vector was\n");
    for(i=0; i<M; i++)
    {
	printf("%.4f\n", rhs[i]);
    }
    //------------------------------------------------------------------------//*/

    


    deallocator1(&A, ncoeff*tgauss);
    ideallocator1(&IPIV, M);
    deallocator1(&B, M);
}
