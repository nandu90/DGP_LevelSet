/***************************************************************************

Author: nsaini
Created: 2018-03-04


Notes:
-  3d groundwater flow problem using 2d decomposition for reference - /home/nsaini3/PGREM3D/flow
-  Conjugate gradient solver - /home/nsaini3/blas_tar/solverc/source
-  BiCGstab or GMRES solver - /home/nsaini3/PGREM3D/trans/

Plan:

Make an integer array which stores the global mapping for each node


***************************************************************************/





#include "common.h"
#include "polylib.h"
#include "generalFunc.h"
#include "memory.h"
#include "mesh.h"
#include "fileIO.h"

int main(int argc, char **argv)
{
    feenableexcept(FE_INVALID   | 
                   FE_DIVBYZERO | 
                   FE_OVERFLOW);
    
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

	char* dirpath;
	dirpath = concat(getexepath(), "/laststep/");
	mkdir(dirpath, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	free(dirpath);
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
    //printf("%d %d\n",gxelem,gyelem);
    iallocator2(&io_info,nprocs,4);  //needed inside partition. Takes care of output IO
    partition();
    //printf("%d %d %d %d\n",xelem,yelem,xnode,ynode);
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Initialize mesh arrays. Done here now instead of common block
    double **x;
    double **y;
    
    allocator2(&x,xnode,ynode);
    allocator2(&y,xnode,ynode);
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Initialize element array
    struct elemdata **edata = NULL;
    //Initialize dof array
    struct dofdata *dof = NULL;
    
    supgcoeff = (int)pow(supgorder+1,2.0);
    //------------------------------------------------------------------------//


    //------------------------------------------------------------------------//
    //Read grid
    if(meshread == 0)
    {
	gridgen(x, y, edata, dof);
    }
    /*else
    {
	gridread(x,y,area,vol,xc,yc);
	}*/
    //------------------------------------------------------------------------//

    printf("herhe\n");
    exit(1);

    /*
      
    //------------------------------------------------------------------------//
    //Initialize element data
    struct elemsclr elem;

    ncoeff = (int)pow(polyorder+1,2.0); //Number of coefficients = number of basis
    
    
    //------------------------------------------------------------------------//

   

    //------------------------------------------------------------------------//
    //Set up Communicator arrays
    genibc(elem.iBC);
    sendptr = (double **) malloc(4 * sizeof(double *));
    recvptr = (double **) malloc(4 * sizeof(double *));
    setupcommu();

    
    //------------------------------------------------------------------------//

    
    //------------------------------------------------------------------------//
    //Initialize solution vectors
    //Arrange quadratures in an array for easy access
    
    
    initialize(elem, x, y);
    
    //------------------------------------------------------------------------//
    //deallocate all arrays - in the reverse order they were allocated
    
    //Element
    deallocator3(&elem.u,xelem,yelem,ncoeff);
    deallocator3(&elem.v,xelem,yelem,ncoeff);
    deallocator3(&elem.phi,xelem,yelem,ncoeff);
    deallocator4(&elem.mass,xelem,yelem,ncoeff,ncoeff);//Should this be ncoeff?- Yes it should be
    ideallocator2(&elem.iBC,xelem,yelem);

    deallocator2(&elem.p,xelem,yelem);
    deallocator2(&elem.rho,xelem,yelem);
    deallocator2(&elem.mu,xelem,yelem);
    deallocator2(&elem.phi2,xelem,yelem);

    //Mesh
    deallocator2(&x,xnode,ynode);
    deallocator2(&y,xnode,ynode);U

    //Use the following if required
    deallocator2(&xc,xelem,yelem);
    deallocator2(&yc,xelem,yelem);
    deallocator2(&vol,xelem,yelem);
    deallocator4(&area,xelem,yelem,2,2);

    //IO array
    ideallocator2(&io_info,nprocs,4);

    //Communication Arrays
    destroycommu();
    free(sendptr);
    free(recvptr);

   
    //------------------------------------------------------------------------//

    
    */
    
    
    
    if(myrank == master)
    {
      printf("---------------------------------------\n");
      printf("-----------That's all folks!-----------\n");
      printf("---------------------------------------\n");
    }
    MPI_Finalize();
}
      
