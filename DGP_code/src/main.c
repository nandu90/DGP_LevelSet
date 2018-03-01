/*
Notes:
-  3d groundwater flow problem using 2d decomposition for reference - /home/nsaini3/PGREM3D/flow
-  Conjugate gradient solver - /home/nsaini3/blas_tar/solverc/source
-  BiCGstab or GMRES solver - /home/nsaini3/PGREM3D/trans/
 */
/*
 * Author: nsaini3
 *
 * Created on February 27, 2018
 */

#include "common.h"
#include "polylib.h"
#include "fileIO.h"
#include "generalFunc.h"
#include "memory.h"
#include "mesh.h"
#include "DGPFunc.h"
/*#include "partition.h"
#include "grid.h"
#include "commu.h"
#include "output.h"*/



int main(int argc, char **argv)
{
    ///Catches mathematical exceptions
    //feenableexcept(FE_INVALID | FE_OVERFLOW |FE_DIVBYZERO);

  //MPI relevant variables/

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
  double time1 = MPI_Wtime();
  //int i,j,k;

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
/////////////////////////////////////
  
    
    //Read control file/
  //Lets not worry about openmp at the moment
  //    omp_set_num_threads(1);

  //All processors may read the control file
   control();

    
    //Routine to partition the mesh
   iallocator2(&io_info,nprocs,4);
    partition();
    
    
    xelem=xelem+4; //Include 2 ghost cells on each side
    yelem=yelem+4; //Include 2 ghost cells on each side
    xnode=xelem+1; //
    ynode=yelem+1; //
    
    ///Resize the vectors and initialize data structures//
    allocator2(&x,xnode,ynode);
    allocator2(&y,xnode,ynode);

    allocator2(&xc,xelem,yelem);
    allocator2(&yc,xelem,yelem);
    allocator2(&vol,xelem,yelem);
    iallocator2(&iBC,xelem,yelem);
    
    allocator4(&area,xelem,yelem,2,2);

    //Initialize the mass matrix
    allocator4(&mass,xelem,yelem,(int)pow(polyorder+1,2.0),(int)pow(polyorder+1,2.0));
    
    //Read grid and populate element and node vectors/
    gridread();
    massmatrix();
    /*genibc();
    sendptr = (double **) malloc(4 * sizeof(double *));
    recvptr = (double **) malloc(4 * sizeof(double *));
    setupcommu();
    for(i=0; i<yelem; i++)
    {
      if(myrank == master)printf("%d %.4f\n",i,bhai.sendrbuf[i]);
      }*/
    

    //Initialize solution vectors/
    struct elemsclr sclr;
    
    allocator3(&sclr.p,xelem,yelem,zelem);
    allocator3(&sclr.u,xelem,yelem,zelem);
    allocator3(&sclr.v,xelem,yelem,zelem);
    allocator3(&sclr.phi,xelem,yelem,zelem);
    allocator3(&sclr.rho,xelem,yelem,zelem);
    allocator3(&sclr.mu,xelem,yelem,zelem);
    
    



    deallocator3(&sclr.p,xelem,yelem,zelem);
    deallocator3(&sclr.u,xelem,yelem,zelem);
    deallocator3(&sclr.v,xelem,yelem,zelem);
    deallocator3(&sclr.phi,xelem,yelem,zelem);
    deallocator3(&sclr.rho,xelem,yelem,zelem);
    deallocator3(&sclr.mu,xelem,yelem,zelem);
    
    deallocator2(&x,xnode,ynode);
    deallocator2(&y,xnode,ynode);
    ideallocator2(&iBC,xelem,yelem);
    deallocator2(&xc,xelem,yelem);
    deallocator2(&yc,xelem,yelem);
    deallocator2(&vol,xelem,yelem);

    deallocator4(&area,xelem,yelem,2,2);
    deallocator4(&mass,xelem,yelem,(int)pow(polyorder+1,2.0),(int)pow(polyorder+1,2.0));
    ideallocator2(&io_info,nprocs,4);
/*destroy commu();
    free(sendptr);
    free(recvptr);
    */
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
      
