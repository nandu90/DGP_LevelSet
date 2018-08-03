/***************************************************************************

Author: nsaini
Created: 2018-03-11

***************************************************************************/


#include "solvers.h"
#include "common.h"
#include "memory.h"

void solveSystem(double **vand, double *rhs, double *soln, int M, int N)
{
    int i,j,k;
    //------------------------------------------------------------------------//
    //Set up variables for LAPACK
    double *A;
    allocator1(&A, M*N);

    k=0;
    for(i=0; i<N; i++)
    {
	for(j=0; j<M; j++)
	{
	    A[k] = vand[j][i];  //Be carefule over here. Store Column-wise for correct result
	    k++;
	}
    }

    
    int NRHS = 1;
    int LDA = M;
    int *IPIV;
    iallocator1(&IPIV, M);

    double *B;
    allocator1(&B, M);

    for(i=0; i<M; i++)
    {
	B[i] = rhs[i];
    }

    int LDB = M;
    int INFO;


    //------------------------------------------------------------------------//
    //Solve the system to get the coefficients
 
    dgesv_(&N, &NRHS, A, &LDA, IPIV, B, &LDB, &INFO);
	

    //Assign to the soln array
    for(i=0; i<M; i++)
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

    


    deallocator1(&A, M*N);
    ideallocator1(&IPIV, M);
    deallocator1(&B, M);
}
