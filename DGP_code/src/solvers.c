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
    allocator1(&A, tgauss*pow(polyorder+1,2));

    double *temp;
    allocator1(&temp, tgauss*pow(polyorder+1,2));

    k=0;
    for(i=0; i<pow(polyorder+1,2); i++)
    {
	for(j=0; j<tgauss; j++)
	{
	    A[k] = vand[j][i];  //Be carefule over here. Store Column-wise for correct result
	    temp[k] = A[k];
	    k++;
	}
    }

    int M = tgauss; //Rows of matrix
    int N = pow(polyorder+1,2);        //Columns of matrix
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
    printf("here\n");
    //Solve the system to get the coefficients
    //if(quadtype == 2)
    //{
	dgesv_(&N, &NRHS, A, &LDA, IPIV, B, &LDB, &INFO);
	//}
	//else if(quadtype == 1)
	//{
	//dgels_(&TRANS,&M,&N,&NRHS,A,&LDA,B,&LDB,WORK,&LWORK,&INFO);
	//}

    

    //------------------------------------------------------------------------//
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
    //------------------------------------------------------------------------//

    for(i=0; i<pow(polyorder+1,2); i++)
    {
	soln[i] = B[i];
    }

    //exit(1);

    deallocator1(&A, pow(polyorder+1,2)*tgauss);
    deallocator1(&temp, pow(polyorder+1,2)*tgauss);
    ideallocator1(&IPIV, M);
    deallocator1(&B, M);
}
