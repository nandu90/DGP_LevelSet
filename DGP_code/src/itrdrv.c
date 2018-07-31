/***************************************************************************

Author: nsaini
Created: 2018-07-02

***************************************************************************/

#define SOLCHECK 1

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

void itrdrv(struct elemsclr elem ,double **x, double **y, double **xc, double **yc, double **vol ,double ****area)
{

    double time1 = MPI_Wtime();
    
    //------------------------------------------------------------------------//
    //Temporary variables
    int ielem, jelem, icoeff;
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Copy initial LS values for comparison - Can be removed later
    double ***iniphi;
    allocator3(&iniphi, xelem, yelem, ncoeff);
    for(ielem=0; ielem<xelem; ielem++)
    {
	for(jelem=0; jelem<yelem; jelem++)
	{
	    for(icoeff=0; icoeff<ncoeff; icoeff++)
	    {
		iniphi[ielem][jelem][icoeff] = elem.phi[ielem][jelem][icoeff];
	    }
	}
    }

    
    
    //------------------------------------------------------------------------//
    //Preliminaries before time loop
    double deltat = 0.0;
    double time = deltat*startstep;
    if(time_control == 2)
    {
	deltat = advect_deltat;
	if(case_tog == 6 || case_tog == 3)
	{
	    deltat = 2.0*PI*25.0/(PI*25.0/314.0);
	    deltat = deltat * advect_deltat;
	}
	if(case_tog == 7)
	{
	    deltat = 2.0*PI*advect_deltat;
	}
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

    FILE *vfout;
    if(myrank == master)
    {
	vfout = fopen("vf.txt","w");
	if(vfout == NULL)
	{
	    printf("Error opening vf.txt!\n");
	    exit(0);
	}
    }

    
    //Print out the paraview output if startstep == 0 i.e., initial conditions
    if(startstep == 0)
    {
	output_xml(elem,startstep,x,y);
    }

   

    //------------------------------------------------------------------------//
    //Temporary variables for INS
    double ***rhs;
    allocator3(&rhs, xelem, yelem, ncoeff);
    
    double **uedge, **vedge;
    allocator2(&uedge, xelem, yelem);
    allocator2(&vedge, xelem, yelem);

    double **ustar, **vstar;
    allocator2(&ustar, xelem, yelem);
    allocator2(&vstar, xelem, yelem);
    
    double **rhsx, **rhsy;
    allocator2(&rhsx, xelem, yelem);
    allocator2(&rhsy, xelem, yelem);

    double **st_forcex, **st_forcey;
    allocator2(&st_forcex, xelem, yelem);
    allocator2(&st_forcey, xelem, yelem);

    double uL, uR, vT, vB;
    
    double cellSize = max(xlen/gxelem, ylen/gyelem);
    double *basisL, *basisR;
    allocator1(&basisL, ncoeff);
    allocator1(&basisR, ncoeff);
    //------------------------------------------------------------------------//

    double inivf;

    
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
		    uL = 0.0;
		    uR = 0.0;
		    vT = 0.0;
		    vB = 0.0;
		    basis2D(1.0,0.0,basisL);
		    basis2D(-1.0,0.0,basisR);
		    for(icoeff=0; icoeff<ncoeff; icoeff++)
		    {
			uL += elem.u[ielem][jelem][icoeff]*basisL[icoeff];
			uR += elem.u[ielem+1][jelem][icoeff]*basisR[icoeff];
		    }
		    basis2D(0.0,1.0,basisL);
		    basis2D(0.0,-1.0,basisR);
		    for(icoeff=0; icoeff<ncoeff; icoeff++)
		    {
			vB += elem.v[ielem][jelem][icoeff]*basisL[icoeff];
			vT += elem.v[ielem][jelem+1][icoeff]*basisR[icoeff];
			
		    }
		    uedge[ielem][jelem] = 0.5*(uL + uR);
		    vedge[ielem][jelem] = 0.5*(vB + vT);

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
		    elem.v[ielem][jelem][0] = 0.5*(vedge[ielem][jelem] + vedge[ielem][jelem-1]);
		    for(icoeff=1; icoeff<ncoeff; icoeff++)
		    {
			elem.u[ielem][jelem][icoeff] = 0.0;
			elem.v[ielem][jelem][icoeff] = 0.0;
		    }
		}
	    }
	    commu2(elem.u);
	    commu2(elem.v);
	    //------------------------------------------------------------------------//

	}

	//------------------------------------------------------------------------//
	//Level-Set advection
	Runge_Kutta(elem, x, y,deltat,rhs, area);
	
	//------------------------------------------------------------------------//

	//------------------------------------------------------------------------//
	//Redistancing
	if(flow_solve == 1 && redist_method == 1)
	{
	    //Copy the high order LS field onto the lower order one for re-distancing
	    for(ielem=0; ielem<xelem; ielem++)
	    {
		for(jelem=0; jelem<yelem; jelem++)
		{
		    if(fabs(elem.phi2[ielem][jelem]) <  10.0*cellSize)
		    {
			elem.phi2[ielem][jelem] = elem.phi[ielem][jelem][0];
		    }
		}
	    }
	    INScommu2(elem.phi2);
	    hyperbolic(elem, area);

	    for(ielem=0; ielem<xelem; ielem++)
	    {
		for(jelem=0; jelem<yelem; jelem++)
		{
		    if(fabs(elem.phi2[ielem][jelem]) > 10.0*cellSize)
		    {
			elem.phi[ielem][jelem][0] = elem.phi2[ielem][jelem];
			for(icoeff=1; icoeff<ncoeff; icoeff++)
			{
			    elem.phi[ielem][jelem][icoeff] = 0.0;
			}
		    }
		}
	    }

	    level_setBC(elem.phi, elem.iBC);
	    
	}
	//------------------------------------------------------------------------//

	
	

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
	//------------------------------------------------------------------------//


	//------------------------------------------------------------------------//
	//Mass Conservation monitoring
	/*if(case_tog == 3 || case_tog == 4)
	{
	    double ***H;
	    allocator3(&H, xelem, yelem, tgauss);

	    double vf;
	    
	    double eps=1.0*max(xlen/(gxelem), ylen/(gyelem));
	    
	    heavy_funcDG(H, elem.phi, eps);

	    //double **inv;
	    //allocator2(&inv, 2, 2);
	    
	    //double jacobian = mappingJacobianDeterminant(2, 2, 0.0, 0.0, x, y, inv);
	    
	    //printf("Jacobian is %.4e and area is %.4e\n", jacobian, area[2][2][0][0]*area[2][2][1][1]);

	    calc_vf(H, jacobian, &vf);

	    if(iter == 0)
	    {
		inivf = vf;
	    }
	    
	    if(myrank == master)
	    {
		printf("Void fraction evolution is %.4e\n",vf/inivf);
		fprintf(vfout,"%d %.4e %.4e\n",iter, time, vf/inivf);
	    }
	    
	    //exit(1);
	    //deallocator2(&inv, 2, 2);
	    deallocator3(&H, xelem, yelem, tgauss);
	}*/
	//------------------------------------------------------------------------//

	//------------------------------------------------------------------------//
	//Active error norm calc
	/*if(case_tog == 1 || case_tog == 6)
	{
	    errorGaussian(elem.phi, time, x, y);
	    }*/
	//------------------------------------------------------------------------//

	if(myrank == master)printf("\n\n");
    }
    
    
    //------------------------------------------------------------------------//
    
    //------------------------------------------------------------------------//
    //Calculate Error Norms
    if(case_tog == 3 || case_tog == 6 || case_tog == 1 || case_tog == 7)
    {
	double err1, lerr1;
	errorNormL1(iniphi, elem.phi, &err1, &lerr1, x, y);
	
	double err2, lerr2;
	errorNormL2(iniphi, elem.phi, &err2, &lerr2, x, y);
	
	if(myrank == master)
	{
	    printf("The L1 norm of error is %.4e and the Log of norm is %.4e\n", err1, lerr1);
	    printf("The L2 norm of error is %.4e and the Log of norm is %.4e\n", err2, lerr2);

	    fprintf(out,"The L1 norm of error is %.4e and the Log of norm is %.4e\n", err1, lerr1);
	    fprintf(out,"The L2 norm of error is %.4e and the Log of norm is %.4e\n", err2, lerr2);
	}
    }
    
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Deallocators
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
    deallocator1(&basisL, ncoeff);
    deallocator1(&basisR, ncoeff);
    
    

    deallocator3(&iniphi, xelem, yelem, ncoeff);

    //------------------------------------------------------------------------//

    double time2 = MPI_Wtime();
    double secs = time2-time1;
    if(myrank == master)
    {
	printf("Total run time: %.6f secs\n",secs);
	fprintf(out,"Total run time: %.6f secs\n",secs);
	fclose(out);
	fclose(vfout);
    }
    
}
