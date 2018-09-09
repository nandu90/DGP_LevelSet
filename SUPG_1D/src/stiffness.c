/***************************************************************************

Author: nsaini
Created: 2018-09-06

***************************************************************************/

#include "common.h"
#include "memory.h"
#include "functions.h"
#include "polylib.h"

void stiffness(double **S, double *x)
{
    //------------------------------------------------------------------------//
    int ielem, icoeff, igauss;
    int icoeff1;
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Define temp array to store element contibutions
    double ***stiff;
    allocator3(&stiff, xelem, ncoeff, ncoeff);
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Define array to store differential of basis functions
    double *bdiff;
    allocator1(&bdiff, ncoeff);
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Define array to store differential of weight functions
    double *wdiff;
    allocator1(&wdiff, ncoeff);
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Other temporary variables
    double *inv, *jacobian;
    allocator1(&inv, 1);
    allocator1(&jacobian, 1);

    double detJ;
    //------------------------------------------------------------------------//

    
    //------------------------------------------------------------------------//
    //Loop over the elements
    for(ielem=0; ielem<xelem; ielem++)
    {
	//Loop over the Gauss Quadrature points
	for(igauss=0; igauss<tgauss; igauss++)
	{
	    //Get the basis vector
	    basisdiff1D(zeta[igauss], bdiff);
	    
	    //Get the weight vector
	    basisdiff1D(zeta[igauss], wdiff);

	    //Get the determinant
	    detJ = mappingJacobianDeterminant(ielem, zeta[igauss], x, inv, jacobian);
	    
	    //Loop over the weight vector
	    for(icoeff=0; icoeff<ncoeff; icoeff++)
	    {
		//Loop over the basis vector
		for(icoeff1=0; icoeff1<ncoeff; icoeff1++)
		{
		    stiff[ielem][icoeff][icoeff1] += weights[igauss]*wdiff[icoeff]*bdiff[icoeff1]*inv[0]*inv[0]*detJ;
		}
	    }
	    /*printf("%.4f %.4f\n",wdiff[0],wdiff[1]);
	    printf("%.4f %.4f\n",bdiff[0],bdiff[1]);
	    printf("%.4f %.4f\n",inv[0],jacobian[0]);
	    printf("%.4f\n",detJ);
	    printf("\n");*/
	}
	//exit(1);
    }
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Assemble the global stiffness matrix
    for(ielem=0; ielem<xelem; ielem++)
    {
	//Define mapping to global matrix
	int row[2] = {ielem, ielem+1};
	int col[2] = {ielem, ielem+1};

	for(icoeff=0; icoeff<ncoeff; icoeff++)
	{
	    for(icoeff1=0; icoeff1<ncoeff; icoeff1++)
	    {
		S[row[icoeff]][col[icoeff1]] += stiff[ielem][icoeff][icoeff1];
	    }
	}
    }
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    /*int inode, inode2;
    for(inode=0; inode<xnode; inode++)
    {
	for(inode2=0; inode2<xnode; inode2++)
	{
	    printf("%.4f ",S[inode][inode2]);
	}
	printf("\n");
	}*/
    //------------------------------------------------------------------------//

    
    //------------------------------------------------------------------------//
    //Deallocators
    deallocator3(&stiff, xelem, ncoeff, ncoeff);
    deallocator1(&bdiff, ncoeff);
    deallocator1(&wdiff, ncoeff);
    deallocator1(&inv, 1);
    deallocator1(&jacobian, 1);
    //------------------------------------------------------------------------//

}
