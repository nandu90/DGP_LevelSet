/***************************************************************************

Author: nsaini
Created: 2018-03-04


Notes:
-  3d groundwater flow problem using 2d decomposition for reference - /home/nsaini3/PGREM3D/flow
-  Conjugate gradient solver - /home/nsaini3/blas_tar/solverc/source
-  BiCGstab or GMRES solver - /home/nsaini3/PGREM3D/trans/

***************************************************************************/




#include "common.h"
#include "polylib.h"
#include "fileIO.h"
#include "generalFunc.h"
#include "memory.h"
#include "mesh.h"
#include "DGPFunc.h"
#include "icbc.h"
#include "commu.h"
#include "rhs.h"
/*#include "partition.h"
#include "grid.h"
#include "commu.h"
#include "output.h"*/



int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    master = 0;

    if(myrank == master)
    {
	printf("---------------------------------------\n");
	printf("------------And so it begins-----------\n");
	printf("---------------------------------------\n");
    }
    int i,j,k;
    double time1 = MPI_Wtime();

    //------------------------------------------------------------------------//
    //Create necessary output directories//
    //Let only the master create directory
    if(myrank == master)
    {
	char* path;
	path = concat(getexepath(), "/output");
	mkdir(path, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	memset(path,0,strlen(path));
	path = concat(getexepath(),"/laststep");
	mkdir(path, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	free(path);
    }
    //------------------------------------------------------------------------//

  
    //------------------------------------------------------------------------//
    //Read control file
    //All processors may read the control file
    control();
    if(myrank == master)printf("Read the input file\n\n");
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Routine to partition the mesh
    iallocator2(&io_info,nprocs,4);  //needed inside partition. Takes care of output IO
    partition();
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Initialize mesh arrays. Done here now instead of common block
    double **x;
    double **y;
    
    allocator2(&x,xnode,ynode);
    allocator2(&y,xnode,ynode);

    //Use the following if required
    /*allocator2(&xc,xelem,yelem);
    allocator2(&yc,xelem,yelem);
    allocator2(&vol,xelem,yelem);
    allocator4(&area,xelem,yelem,2,2);*/
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Read grid
    gridread(x,y);
    //------------------------------------------------------------------------//

      
    //------------------------------------------------------------------------//
    //Initialize element data
    struct elemsclr elem;

    ncoeff = (int)pow(polyorder+1,2.0); //Number of coefficients = number of basis
    
    if(quadtype == 1)
    {
	xgpts=polyorder+2;
	ygpts=polyorder+2;
    }
    else if(quadtype == 2)
    {
	xgpts = polyorder+1;
	ygpts = polyorder+1;
    }
    tgauss = xgpts*ygpts;
    
    allocator3(&elem.u,xelem,yelem,ncoeff);
    allocator3(&elem.v,xelem,yelem,ncoeff);
    allocator3(&elem.phi,xelem,yelem,ncoeff);
    allocator4(&elem.mass,xelem,yelem,tgauss,tgauss);
    iallocator2(&elem.iBC,xelem,yelem);

    //Uncomment when you need them
    /*allocator3(&elem.p,xelem,yelem,zelem);
    allocator3(&elem.rho,xelem,yelem,zelem);
    allocator3(&elem.mu,xelem,yelem,zelem);*/
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //For Gauss-Lobatto-Legendre Quadrature we need two points at the end + polyorder
    double *zeta1, *zeta2;
    double *weight1, *weight2;
    allocator1(&zeta1, xgpts);
    allocator1(&zeta2, ygpts);
    allocator1(&weight1, xgpts);
    allocator1(&weight2, ygpts);
    if(quadtype == 1) //Gauss-Legendre-Lobatto
    {
	zwgll(zeta1,weight1,xgpts);
	zwgll(zeta2,weight2,ygpts);
    }
    else if(quadtype == 2) //Gauss-Legendre
    {
        zwgl(zeta1,weight1,xgpts);
	zwgl(zeta2,weight2,ygpts);
    }

    //Arrange quadrature points in a easy to access array
    allocator2(&zeta, tgauss,2);
    allocator2(&weights, tgauss,2);
    k=0;
    for(j=0; j<ygpts; j++)
    {
	for(i=0; i<xgpts; i++)
	{
	    zeta[k][0] = zeta1[i];
	    zeta[k][1] = zeta2[j];
	    weights[k][0] = weight1[i];
	    weights[k][1] = weight2[j];
	    k++;
	}
    }

    if(myrank == master)
    {
	printf("The quadrature points and weights are:\n");
	for(i=0; i<tgauss; i++)
	{
	    printf("%d  %.6f  %.6f  %.6f  %.6f\n",i,zeta[i][0],zeta[i][1],weights[i][0],weights[i][1]);
	}
	printf("\n");
    }
    
    //Get the mass matrix
    massmatrix(x,y,elem.mass); 
    //------------------------------------------------------------------------//

   

    //------------------------------------------------------------------------//
    //Set up Communicator arrays
    genibc(elem.iBC);
    sendptr = (double **) malloc(4 * sizeof(double *));
    recvptr = (double **) malloc(4 * sizeof(double *));
    setupcommu();
    /*for(i=0; i<yelem; i++)
    {
      if(myrank == master)printf("%d %.4f\n",i,bhai.sendrbuf[i]);
      }*/
    //------------------------------------------------------------------------//

    
    //------------------------------------------------------------------------//
    //Initialize solution vectors
    //Arrange quadratures in an array for easy access
    
    initialize(elem, x, y); 
    //------------------------------------------------------------------------//

    
    
    //------------------------------------------------------------------------//
    //Time loop will start here
    double ***rhs;
    allocator3(&rhs, xelem, yelem, ncoeff);
    getRHS(elem, x, y, rhs);
    //------------------------------------------------------------------------//

    
    
    
    //------------------------------------------------------------------------//
    //Print out the paraview outputs
    output_xml(elem, 0, x, y);
    //------------------------------------------------------------------------//



    //------------------------------------------------------------------------//
    //deallocate all arrays - in the reverse order they were allocated
    //Gauss quadrature
    deallocator2(&zeta, tgauss,2);
    deallocator2(&weights, tgauss,2);
    deallocator1(&zeta1, xgpts);
    deallocator1(&zeta2, ygpts);
    deallocator1(&weight1, xgpts);
    deallocator1(&weight2, ygpts);
    //Element
    deallocator3(&elem.u,xelem,yelem,ncoeff);
    deallocator3(&elem.v,xelem,yelem,ncoeff);
    deallocator3(&elem.phi,xelem,yelem,ncoeff);
    deallocator4(&elem.mass,xelem,yelem,tgauss,tgauss); //Should this be ncoeff? - Yes it should be
    ideallocator2(&elem.iBC,xelem,yelem);

    //Uncomment when you need them
    /*deallocator3(&elem.p,xelem,yelem,zelem);
    deallocator3(&elem.rho,xelem,yelem,zelem);
    deallocator3(&elem.mu,xelem,yelem,zelem);*/

    //Mesh
    deallocator2(&x,xnode,ynode);
    deallocator2(&y,xnode,ynode);

    //Use the following if required
    /*deallocator2(&xc,xelem,yelem);
    deallocator2(&yc,xelem,yelem);
    deallocator2(&vol,xelem,yelem);
    deallocator4(&area,xelem,yelem,2,2);*/

    //IO array
    ideallocator2(&io_info,nprocs,4);

    //Communication Arrays
    destroycommu();
    free(sendptr);
    free(recvptr);

    //Time loop arrays
    deallocator3(&rhs, xelem, yelem, ncoeff);
    //------------------------------------------------------------------------//

    
    
    
    double time2 = MPI_Wtime();
    double secs = time2-time1;
    if(myrank == master)
    {
	printf("Total run time: %.6f secs\n",secs);
/*fprintf(out,"Total run time: %.6f secs\n",secs);
  fclose(out);*/
    }
    
    if(myrank == master)
    {
      printf("---------------------------------------\n");
      printf("-----------That's all folks!-----------\n");
      printf("---------------------------------------\n");
    }
    MPI_Finalize();
}
      
