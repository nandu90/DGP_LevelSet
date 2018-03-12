/***************************************************************************

Author: nsaini
Created: 2018-03-11

***************************************************************************/


#include "solvers.h"
#include "common.h"
#include "memory.h"

void solveSystem(double **vand, double *rhs)
{
    int i,j,k;
    //------------------------------------------------------------------------//
    //Set up variables for LAPACK
    double *A;
    allocator1(&A, pow(pow(polyorder+1,2),2));

    double *temp;
    allocator1(&temp, pow(pow(polyorder+1,2),2));

    k=0;
    for(i=0; i<pow(polyorder+1,2); i++)
    {
	for(j=0; j<pow(polyorder+1,2); j++)
	{
	    A[k] = vand[j][i];  //Be carefule over here. Store Column-wise for correct result
	    temp[k] = A[k];
	    k++;
	}
    }

    int N = pow(polyorder+1,2);
    int NRHS = N;
    int LDA = N;
    int *IPIV;
    iallocator1(&IPIV, N);

    double *B;
    allocator1(&B, N);

    for(i=0; i<pow(polyorder+1,2); i++)
    {
	B[i] = rhs[i];
    }

    int LDB = N;
    int INFO;
    
    //------------------------------------------------------------------------//

    //Solve the system to get the coefficients
    dgesv_(&N, &NRHS, A, &LDA, IPIV, B, &LDB, &INFO);

    

    //------------------------------------------------------------------------//
    //check
    printf("The Vandermonde matrix is\n");
    for(i=0; i<pow(pow(polyorder+1,2),2); i++)
    {
	printf("%d %.4f\n",i, temp[i]);
    }
    printf("\n\n");

    printf("The solution vector is\n");
    for(i=0; i<pow(polyorder+1,2); i++)
    {
	printf("%.4f\n", B[i]);
    }

    printf("The rhs vector was\n");
    for(i=0; i<pow(polyorder+1,2); i++)
    {
	printf("%.4f\n", rhs[i]);
    }
    //------------------------------------------------------------------------//

    for(i=0; i<pow(polyorder+1,2); i++)
    {
	rhs[i] = B[i];
    }

    exit(1);

    deallocator1(&A, pow(pow(polyorder+1,2),2));
    deallocator1(&temp, pow(pow(polyorder+1,2),2));
    ideallocator1(&IPIV, N);
    deallocator1(&B, N);
}
