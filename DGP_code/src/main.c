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
#include "solvers.h"
#include "INS.h"

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
    double **xc, **yc;
    double **vol;
    double ****area;
    allocator2(&xc,xelem,yelem);
    allocator2(&yc,xelem,yelem);
    allocator2(&vol,xelem,yelem);
    allocator4(&area,xelem,yelem,2,2);
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Read grid
    gridread(x,y,area, vol, xc, yc);
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
    allocator2(&elem.p,xelem,yelem);
    allocator2(&elem.rho,xelem,yelem);
    allocator2(&elem.mu,xelem,yelem);
    allocator2(&elem.phi2,xelem,yelem);
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
        if(xgpts == 1)
	{
	    zeta1[0] = 0.0;
	    weight1[0] = 1.0;
	}
	else
	{
	    zwgl(zeta1,weight1,xgpts);
	}
	if(ygpts == 1)
	{
	    zeta2[0] = 0.0;
	    weight2[0] = 1.0;
	}
	else
	{
	    zwgl(zeta2,weight2,ygpts);
	}
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

    INSsendptr = (double **) malloc(4 * sizeof(double *));
    INSrecvptr = (double **) malloc(4 * sizeof(double *));
    INSsetupcommu();
    /*for(i=0; i<yelem; i++)
    {
      if(myrank == master)printf("%d %.4f\n",i,bhai.sendrbuf[i]);
      }*/
    //------------------------------------------------------------------------//

    
    //------------------------------------------------------------------------//
    //Initialize solution vectors
    //Arrange quadratures in an array for easy access
    
    initialize(elem, x, y);
    INSinitialize(elem);
    //------------------------------------------------------------------------//

    
    
    //------------------------------------------------------------------------//
    //Preliminaries before time loop
    double deltat = 0.0;
    double time = deltat*startstep;
    if(time_control == 2)
    {
	deltat = advect_deltat;
    }
    else
    {
	if(myrank == master)
	{
	    printf("Time control is not constant time.\nExiting...");
	    exit(1);
	}
    }
   

    int iter = startstep;
    int print_count = 1;

    FILE *out;
    if(myrank == master)
    {
	out = fopen("sim_out.txt","w");
	if(out == NULL)
	{
	    printf("Error opening sim_out.txt!\n");
	    exit(0);
	}
    }

    
    //Print out the paraview output if startstep == 0 i.e., initial conditions
    if(startstep == 0)
    {
	output_xml(elem,startstep,x,y);
    }

    double ***rhs;
    allocator3(&rhs, xelem, yelem, ncoeff);

    //------------------------------------------------------------------------//
    //Temporary variables for INS
    double **uedge, **vedge;
    allocator2(&uedge, xelem, yelem);
    allocator2(&vedge, xelem, yelem);

    double **ustar, **vstar;
    allocator2(&ustar, xelem, yelem);
    allocator2(&vstar, xelem, yelem);
    
    int ielem, jelem;

    double **rhsx, **rhsy;
    allocator2(&rhsx, xelem, yelem);
    allocator2(&rhsy, xelem, yelem);

    double **st_forcex, **st_forcey;
    allocator2(&st_forcex, xelem, yelem);
    allocator2(&st_forcey, xelem, yelem);

    double uL, uR, vT, vB;
    //------------------------------------------------------------------------//

    //Time loop
    for(iter = startstep; iter < itermax; iter++)
    {
	time += deltat;
	
	if(myrank == master)
	{
	    printf("Step: %d Time: %.4f\n",iter+1, time);
	    fprintf(out,"Step: %d Time: %.4f\n",iter+1, time);
	}

	if(flow_solve == 1)
	{
	    if(myrank == master)printf("Solving Flow\n");

	    //------------------------------------------------------------------------//
	    //Note that the first coefficient gives the velocity at the centroid
	    //Upwind and assign velocities to uedge, vedge
	    for(ielem=1; ielem<xelem-1; ielem++)
	    {
		for(jelem=1; jelem<yelem-1; jelem++)
		{
		    uL = elem.u[ielem][jelem][0];
		    uR = elem.u[ielem+1][jelem][0];
		    if(uL > uR)
		    {
			uedge[ielem][jelem] = uL;
		    }
		    else
		    {
			uedge[ielem][jelem] = uR;
		    }
		    vB = elem.v[ielem][jelem][0];
		    vT = elem.v[ielem][jelem+1][0];
		    if(vB > vT)
		    {
			vedge[ielem][jelem] = vB;
		    }
		    else
		    {
			vedge[ielem][jelem] = vT;
		    }
		}
	    }
	    //------------------------------------------------------------------------//

	    //------------------------------------------------------------------------//
	    //Get  RHS
	    rhscalc(elem, rhsx, rhsy, area, vol, elem.iBC, uedge, vedge);
	    //------------------------------------------------------------------------//

	    //------------------------------------------------------------------------//
	    //Predictor Step
	    for(ielem=1; ielem<xelem-1; ielem++)
	    {
		for(jelem=1; jelem<yelem-1; jelem++)
		{
		    ustar[ielem][jelem] = uedge[ielem][jelem] + deltat*rhsx[ielem][jelem];
		    vstar[ielem][jelem] = vedge[ielem][jelem] + deltat*rhsy[ielem][jelem];
		}
	    }

	    INScommu2(ustar);
	    INScommu2(vstar);
	    vel_BC(ustar, vstar, elem.iBC);
	    //------------------------------------------------------------------------//

	    //------------------------------------------------------------------------//
	    //Add in contributions from body and surface tension force
	    surface(elem, st_forcex, st_forcey, area);
	    body(elem, st_forcex, st_forcey, vol);	    
	    //------------------------------------------------------------------------//

	    //------------------------------------------------------------------------//
	    //Solve the pressure Poisson equation
	    variable_pressure(ustar, vstar, elem.p, deltat, elem.rho, st_forcex, st_forcey, area, vol, elem.iBC);
	    //------------------------------------------------------------------------//

	    //------------------------------------------------------------------------//
	    //Projection Step
	    for(ielem=2; ielem<xelem-2; ielem++)
	    {
		for(jelem=2; jelem<yelem-2; jelem++)
		{
		    ustar[ielem][jelem] = ustar[ielem][jelem] - deltat*((2.0/(elem.rho[ielem][jelem]+elem.rho[ielem+1][jelem]))*((elem.p[ielem+1][jelem]-elem.p[ielem][jelem])/area[ielem][jelem][1][1] + 0.5*(st_forcex[ielem+1][jelem]+st_forcex[ielem][jelem])));
                    vstar[ielem][jelem] = vstar[ielem][jelem] - deltat*((2.0/(elem.rho[ielem][jelem]+elem.rho[ielem][jelem+1]))*((elem.p[ielem][jelem+1]-elem.p[ielem][jelem])/area[ielem][jelem][0][0] + 0.5*(st_forcey[ielem][jelem+1]+st_forcey[ielem][jelem])));
		}
	    }

	    INScommu2(ustar);
	    INScommu2(vstar);
	    vel_BC(ustar, vstar, elem.iBC);

	    //Assign to elem vectors
	    for(ielem=1; ielem<xelem; ielem++)
	    {
		for(jelem=1; jelem<yelem; jelem++)
		{
		    uedge[ielem][jelem] = ustar[ielem][jelem];
		    vedge[ielem][jelem] = vstar[ielem][jelem];

		    elem.u[ielem][jelem][0] = 0.5*(uedge[ielem][jelem] + uedge[ielem-1][jelem]);
		    elem.v[ielem][jelem][0] = 0.5*(uedge[ielem][jelem] + uedge[ielem][jelem-1]);
		}
	    }
	    //------------------------------------------------------------------------//

	}

	//------------------------------------------------------------------------//
	//Level-Set advection
	Runge_Kutta(elem, x, y,deltat,rhs);
	//------------------------------------------------------------------------//

	//------------------------------------------------------------------------//
	//Redistancing
	if(redist_method == 1)
	{
	    //Copy the high order LS field onto the lower order one for re-distancing
	    for(ielem=0; ielem<xelem; ielem++)
	    {
		for(jelem=0; jelem<yelem; jelem++)
		{
		    elem.phi2[ielem][jelem] = elem.phi[ielem][jelem][0];
		}
	    }

	    hyperbolic(elem, area);
	}
	//------------------------------------------------------------------------//

	

	
	//Print out the paraview output
	print_count++;
        if(print_count == 1)
        {
            output_xml(elem,iter+1,x,y);
        }
        if(print_count == print_gap)
        {
            print_count = 0;
        }
    }
    
    
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

    deallocator2(&elem.p,xelem,yelem);
    deallocator2(&elem.rho,xelem,yelem);
    deallocator2(&elem.mu,xelem,yelem);
    deallocator2(&elem.phi2,xelem,yelem);

    //Mesh
    deallocator2(&x,xnode,ynode);
    deallocator2(&y,xnode,ynode);

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

    INSdestroycommu();
    free(INSsendptr);
    free(INSrecvptr);

    //Time arrays
    deallocator3(&rhs, xelem, yelem, ncoeff);

    //INS related
    deallocator2(&uedge, xelem, yelem);
    deallocator2(&vedge, xelem, yelem);
    deallocator2(&ustar, xelem, yelem);
    deallocator2(&vstar, xelem, yelem);
    deallocator2(&rhsx, xelem, yelem);
    deallocator2(&rhsy, xelem, yelem);
    deallocator2(&st_forcex, xelem, yelem);
    deallocator2(&st_forcey, xelem, yelem);
    //------------------------------------------------------------------------//

    
    
    
    double time2 = MPI_Wtime();
    double secs = time2-time1;
    if(myrank == master)
    {
	printf("Total run time: %.6f secs\n",secs);
	fprintf(out,"Total run time: %.6f secs\n",secs);
	fclose(out);
    }
    
    if(myrank == master)
    {
      printf("---------------------------------------\n");
      printf("-----------That's all folks!-----------\n");
      printf("---------------------------------------\n");
    }
    MPI_Finalize();
}
      
