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
	    deltat = 2.0*PI*0.25/(PI*0.25/3.14);
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
    //Temporary variables
    double ***rhs;
    allocator3(&rhs, xelem, yelem, ncoeff);
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
	}

	//------------------------------------------------------------------------//
	//Level-Set advection
	Runge_Kutta(elem, x, y,deltat,rhs, area);
	
	//------------------------------------------------------------------------//

	//------------------------------------------------------------------------//
	//Redistancing
	if(flow_solve == 1 && redist_method == 1)
	{
	    
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
    deallocator3(&iniphi, xelem, yelem, ncoeff);
    deallocator3(&rhs, xelem, yelem, ncoeff);
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
