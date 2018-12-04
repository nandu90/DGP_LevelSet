#include "generalFunc.h"
#include "common.h"
#include "DGPFunc.h"
#include "memory.h"
#include "polylib.h"
#include "rhs.h"

double minmod(double *array, int size)
{
    double result;

    int i;
    //Check whether all components have same sign
    int flag = 0;
    for(i=0; i<size-1; i++)
    {
	if(array[i]*array[i+1] < 0.0)
	{
	    flag = 1;
	    break; //Different sign detected
	}
    }

    double mini = fabs(array[0]);
    double a = array[0];
    if(flag == 0)
    {
	for(i=1; i<size; i++)
	{
	    mini = min(mini, fabs(array[i]));
	}
	
	result = mini;
	
	if(fabs(a) > 0.0)
	{
	    result = a*result/fabs(a);
	}
	else
	{
	    result = 0.0;
	}
    }
    else
    {
	result = 0.0;
    }

    return result;
}

double remainder(double a, double b)
{
    int i;
    int n = 100000000;

    a = fabs(a);
    b = fabs(b);
    if(b == 0.0)
    {
	printf("Division by 0 in remainder() function");
	exit(1);
    }
    for(i=0; i<n; i++)
    {
	if(a >= b*i && a < b*(i+1))
	{
	    break;
	}
    }

    double rem = a - b*i;

    return rem;
}
double max(double a, double b)
{
  if(a > b)
    {
      return a;
    }
  else
    {
      return b;
    }
}

double min(double a, double b)
{
  if(a > b)
    {
      return b;
    }
  else
    {
      return a;
    }
}

char* getexepath()
{
  static char cwd[1024];
  char *err = getcwd(cwd,sizeof(cwd));
  if(err == NULL)
    {
      printf("Error getting the current working directory\n.Exiting...");
      exit(1);
    }
  else
    {
      return cwd;
    }
}

char* concat(char s1[], char s2[])
{
    char* result = malloc(strlen(s1)+strlen(s2)+1);//+1 for the null-terminator
    //in real code you would check for errors in malloc here
    strcpy(result, s1);
    strcat(result, s2);
    return result;
}


void naturalToCartesian(double **xs, double **x, double **y, int i, int j, double **zeta, int tgauss)
{
    int k;
    double xvertex[4];
    double yvertex[4];
    
    xvertex[0] = x[i][j];
    yvertex[0] = y[i][j];

    xvertex[1] = x[i+1][j];
    yvertex[1] = y[i+1][j];

    xvertex[2] = x[i][j+1];
    yvertex[2] = y[i][j+1];

    xvertex[3] = x[i+1][j+1];
    yvertex[3] = y[i+1][j+1];

       
    double N1, N2, N3, N4;

    
    //Populate the coordinate vector
    for(k=0; k<tgauss; k++)
    {
	N1 = (1.0-zeta[k][0])*(1.0-zeta[k][1])/4.0;
	N2 = (1.0+zeta[k][0])*(1.0-zeta[k][1])/4.0;
	N3 = (1.0-zeta[k][0])*(1.0+zeta[k][1])/4.0;
	N4 = (1.0+zeta[k][0])*(1.0+zeta[k][1])/4.0;
	
	
	xs[k][0] = N1*xvertex[0] + N2*xvertex[1] + N3*xvertex[2] + N4*xvertex[3];
	xs[k][1] = N1*yvertex[0] + N2*yvertex[1] + N3*yvertex[2] + N4*yvertex[3];
    }


}


void errorNormL2(double ***iniphi, double ***phi, double *err, double *lerr, double **x, double **y)
{
    int ielem, jelem, icoeff;
    int igauss;

    double *basis;
    allocator1(&basis, ncoeff);

    double recini;
    double rec;
	
    (*err) = 0.0;

    double sum = 0.0;
    double sum1 = 0.0;

    double elemsum = 0.0;
    double elemsum1 = 0.0;

    double *inv,  *jacobian;
    allocator1(&inv, 4);
    allocator1(&jacobian, 4);

    double detJ;

   

    //------------------------------------------------------------------------//
    //Define quad points and weights here independent of what is in the rest of the code
   
    int extra = 4;
    
    double **z, **w;
    int ngauss = pow(polyorder + 1 + extra, 2);

    allocator2(&z, ngauss,2);
    allocator2(&w, ngauss,2);
    
    GaussPoints2D(z, w, quadtype, ngauss);

    double **xs;
    allocator2(&xs, ngauss, 2);
    //------------------------------------------------------------------------//
    
    for(ielem =2; ielem<xelem-2; ielem++)
    {
	for(jelem=2; jelem<yelem-2; jelem++)
	{
	    elemsum = 0.0;
	    elemsum1 = 0.0;
	    
	    naturalToCartesian(xs, x, y, ielem, jelem, z, ngauss);
	    
	    for(igauss=0; igauss<ngauss; igauss++)
	    {
		recini = 0.0;
		rec = 0.0;
		basis2D(z[igauss][0], z[igauss][1], basis);
		for(icoeff=0; icoeff<ncoeff; icoeff++)
		{
		    recini += basis[icoeff]*iniphi[ielem][jelem][icoeff];
		    rec += basis[icoeff]*phi[ielem][jelem][icoeff];
		}
		detJ = mappingJacobianDeterminant(ielem, jelem, z[igauss][0], z[igauss][1], x, y, inv, jacobian);
		if(case_tog == 3)
		{
		    if((xs[igauss][0] >= 0.3) && (xs[igauss][0] <= 0.4) && (xs[igauss][1] >= 0.6) && (xs[igauss][1] <= 0.9))
		    {
			elemsum += pow(recini-rec,2.0)*w[igauss][0]*w[igauss][1]*detJ;
			elemsum1 += pow(recini,2.0)*w[igauss][0]*w[igauss][1]*detJ;
		    }
		    else
		    {
			elemsum = 0.0;
			elemsum1 = 0.0;
			break;
		    }
		}
		else if(case_tog == 6)
		{
		    if((xs[igauss][0] >= 0.3) && (xs[igauss][0] <= 0.7) && (xs[igauss][1] >= 0.5) && (xs[igauss][1] <= 0.7))
		    {
			elemsum += pow(recini-rec,2.0)*w[igauss][0]*w[igauss][1]*detJ;
			elemsum1 += pow(recini,2.0)*w[igauss][0]*w[igauss][1]*detJ;
		    }
		    else
		    {
			elemsum = 0.0;
			elemsum1 = 0.0;
			break;
		    }
		}
		else
		{
		    elemsum += pow(recini-rec,2.0)*w[igauss][0]*w[igauss][1]*detJ;
		    elemsum1 += pow(recini,2.0)*w[igauss][0]*w[igauss][1]*detJ;
		}
		
	    }
	    sum += elemsum;
	    sum1 += elemsum1;
	}
    }

    MPI_Allreduce(&sum, err, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    double exact;

    MPI_Allreduce(&sum1, &exact, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		  
    (*err) = sqrt((*err));///sqrt(exact);

    

    deallocator1(&basis, ncoeff);
    deallocator1(&inv, 4);
    deallocator1(&jacobian, 4);
    deallocator2(&z, ngauss, 2);
    deallocator2(&w, ngauss, 2);
    deallocator2(&xs, ngauss, 2);
}


void errorNormL1(double ***iniphi, double ***phi, double *err, double *lerr, double **x, double **y)
{
    int ielem, jelem, icoeff;
    int igauss;

    double *basis;
    allocator1(&basis, ncoeff);

    double recini;
    double rec;
	
    (*err) = 0.0;

    double sum = 0.0;
    double sum1 = 0.0;
    double elemsum = 0.0;
    double elemsum1 = 0.0;
    
    double *inv,  *jacobian;
    allocator1(&inv, 4);
    allocator1(&jacobian, 4);

    double detJ;

   
    
    //------------------------------------------------------------------------//
    //Define quad points and weights here independent of what is in the rest of the code
   
    int extra = 4;
    
    double **z, **w;
    int ngauss = pow(polyorder + 1 + extra, 2);

    allocator2(&z, ngauss,2);
    allocator2(&w, ngauss,2);
    
    GaussPoints2D(z, w, quadtype, ngauss);

    double **xs;
    allocator2(&xs, ngauss, 2);
    //------------------------------------------------------------------------//
    
    for(ielem =2; ielem<xelem-2; ielem++)
    {
	for(jelem=2; jelem<yelem-2; jelem++)
	{
	    elemsum = 0.0;
	    elemsum1 = 0.0;
	    
	    naturalToCartesian(xs, x, y, ielem, jelem, z, ngauss);
	    
	    for(igauss=0; igauss<ngauss; igauss++)
	    {
		recini = 0.0;
		rec = 0.0;
		basis2D(z[igauss][0], z[igauss][1], basis);
		for(icoeff=0; icoeff<ncoeff; icoeff++)
		{
		    recini += basis[icoeff]*iniphi[ielem][jelem][icoeff];
		    rec += basis[icoeff]*phi[ielem][jelem][icoeff];
		}
		detJ = mappingJacobianDeterminant(ielem, jelem, z[igauss][0], z[igauss][1], x, y, inv, jacobian);
		if(case_tog == 3)
		{
		    if((xs[igauss][0] >= 0.3) && (xs[igauss][0] <= 0.4) && (xs[igauss][1] >= 0.6) && (xs[igauss][1] <= 0.9))
		    {
			elemsum += fabs(recini-rec)*w[igauss][0]*w[igauss][1]*detJ;
			elemsum1 += fabs(recini)*w[igauss][0]*w[igauss][1]*detJ;
		    }
		    else
		    {
			elemsum = 0.0;
			elemsum1 = 0.0;
			break;
		    }
		}
		else if(case_tog == 6)
		{
		    if((xs[igauss][0] >= 0.3) && (xs[igauss][0] <= 0.7) && (xs[igauss][1] >= 0.5) && (xs[igauss][1] <= 0.7))
		    {
			elemsum += fabs(recini-rec)*w[igauss][0]*w[igauss][1]*detJ;
			elemsum1 += fabs(recini)*w[igauss][0]*w[igauss][1]*detJ;
		    }
		    else
		    {
			elemsum = 0.0;
			elemsum1 = 0.0;
			break;
		    }
		}
		else
		{
		    elemsum += fabs(recini-rec)*w[igauss][0]*w[igauss][1]*detJ;
		    elemsum1 += fabs(recini)*w[igauss][0]*w[igauss][1]*detJ;
		}
		
	    }
	    sum += elemsum;
	    sum1 += elemsum1;
	}
    }

    MPI_Allreduce(&sum, err, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    double exact;
     
    MPI_Allreduce(&sum1, &exact, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    (*err) = ((*err));///(exact);

    deallocator1(&basis, ncoeff);
    deallocator1(&inv, 4);
    deallocator1(&jacobian, 4);
    deallocator2(&z, ngauss, 2);
    deallocator2(&w, ngauss, 2);
    deallocator2(&xs, ngauss, 2);
}


void calc_vf(double ***phi, double **x, double **y, double *inivf)
{
    //------------------------------------------------------------------------//
    //Loop indexes
    int ielem, jelem;
    int igauss;
    int icoeff;
    //------------------------------------------------------------------------//

    double *basis;
    allocator1(&basis, ncoeff);
    
    double *inv,  *jacobian;
    allocator1(&inv, 4);
    allocator1(&jacobian, 4);

    double detJ;
    
    //------------------------------------------------------------------------//
    //Define quad points and weights here independent of what is in the rest of the code

    int extra = 4;
        
    double **z, **w;
    int ngauss = pow(polyorder + 1 + extra, 2);

    allocator2(&z, ngauss,2);
    allocator2(&w, ngauss,2);
    
    GaussPoints2D(z, w, quadtype, ngauss);
    
    
	
    //------------------------------------------------------------------------//
    
    double vf;
    double totalvf;
    double recphi;

    vf = 0.0;
    
    for(ielem = 2; ielem<xelem-2; ielem++)
    {
	for(jelem = 2; jelem<yelem-2; jelem++)
	{
	    for(igauss=0; igauss<ngauss; igauss++)
	    {
		basis2D(z[igauss][0], z[igauss][1], basis);

		recphi = 0.0;
		for(icoeff=0; icoeff<ncoeff; icoeff++)
		{
		    recphi += phi[ielem][jelem][icoeff]*basis[icoeff];
		}
		detJ = mappingJacobianDeterminant(ielem, jelem, z[igauss][0], z[igauss][1], x, y, inv, jacobian);
		//Interchange sign and normalize
		if(recphi != 0.0)
		{
		    recphi = -recphi/fabs(recphi);
		}
		vf += w[igauss][0]*w[igauss][1]*max(0.0,recphi)*detJ;
	    }
	}
    }

    MPI_Allreduce(&vf, &totalvf, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    if((*inivf) == 0.0)
    {
	(*inivf) = totalvf;
    }
    if(myrank == master)
    {
	printf("The internal volume is %.6e and error is %.6e\n",totalvf, fabs(totalvf - (*inivf)));
    }
    //------------------------------------------------------------------------//
    //Deallocators
    deallocator1(&basis, ncoeff);
    deallocator1(&inv, 4);
    deallocator1(&jacobian, 4);
    deallocator2(&z, ngauss, 2);
    deallocator2(&w, ngauss, 2);
    //------------------------------------------------------------------------//

}


/*void errorGaussian(double ***phi, double time, double **x, double **y)
{
    int ielem, jelem, icoeff;
    int igauss;

    double *basis;
    allocator1(&basis, ncoeff);

    double rec;

    //Allocate coordinate matrix corresponding to zs - solution points
    double **xs;
    allocator2(&xs, tgauss, 2);
    double sigmax;
    double sigmay;
    double term1;
    double term2;
    double exact = 0.0;

    double exactx = 0.0;

    double sum1 = 0.0;
    double sum2 = 0.0;
    double suminf = 0.0;
    double error1, error2, errorinf;

    //Position of center on the circumference
    double speed, dist,theta,xc = 0.0,yc = 0.0;
    if(case_tog == 6)
    {
	 speed = PI*25.0/314.0;
	 dist = speed * time;
	 theta = PI/2.0 + dist/25.0;
	 yc = 25.0*sin(theta)+50.0;
	 xc = 25.0*cos(theta)+50.0;
    
	if(myrank == master)
	{
	    printf("Theta %.4e, xc %.4e, yc %.4e\n",theta*180.0/PI, xc, yc);
	}
    }
    //Lets compare at the Gauss Quadrature points
    for(ielem =2; ielem<xelem-2; ielem++)
    {
	for(jelem=2; jelem<yelem-2; jelem++)
	{
	    //Convert natural coordinates at quadrature points to Cartesian
	   
	    naturalToCartesian(xs, x, y, ielem, jelem);
	    
	    for(igauss=0; igauss<tgauss; igauss++)
	    {
		rec = 0.0;
		basis2D(zeta[igauss][0], zeta[igauss][1], basis);
		for(icoeff=0; icoeff<ncoeff; icoeff++)
		{
		    rec += basis[icoeff]*phi[ielem][jelem][icoeff];
		}

		if(case_tog == 1)
		{
		    sigmax = 25.0;
		    sigmay = 25.0;
		    exactx = xs[igauss][0] - 1.0*remainder(time, 150.0);
		    if(exactx < 0.0)
		    {
			exactx = 150.0 + exactx;
		    }
		    else if(exactx > 150.0)
		    {
			exactx = exactx - 150.0;
		    }
		    term1 = 0.5*pow((exactx - xb_in)/sigmax,2.0);
		    term2 = 0.5*pow((xs[igauss][1] - yb_in)/sigmay,2.0);
		    exact = 1.0*exp(-(term1 + term2));
		}
		else if(case_tog == 6)
		{
		    
		    exact = sqrt(pow(xc-xs[igauss][0],2.0) + pow(yc-xs[igauss][1],2.0)) - 15.0;
		    if(exact > 1.0)
		    {
			exact = rec;
		    }
		}
		sum1 += fabs(exact - rec);
		sum2 += pow(exact - rec, 2.0);

		suminf = max(suminf, fabs(exact-rec));
	    }
	    
	}
    }

    MPI_Allreduce(&sum1, &error1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&sum2, &error2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    error2 = sqrt(error2);

    MPI_Allreduce(&suminf, &errorinf, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    
    if(myrank == master)
    {
	printf("The errors for time %.4e are:\n L1: %.4e and Log of error is %.4e\n L2: %.4e and Log of error is %.4e\n Linf: %.4e and Log of error is %.4e\n",time,error1, -log(error1), error2, -log(error2), errorinf, -log(errorinf));
    }
    
    deallocator1(&basis, ncoeff);
    deallocator2(&xs, tgauss, 2);
    
    //if(myrank == master) exit(1);
    }*/


void errorMMS(double ***phi, double **x, double **y, double *L1norm, double *L2norm)
{
    //------------------------------------------------------------------------//
    //Temp
    int ielem, jelem, icoeff;
    int igauss;

    double *basis;
    allocator1(&basis, ncoeff);
    
    double *inv,  *jacobian;
    allocator1(&inv, 4);
    allocator1(&jacobian, 4);

    double detJ;
    double xmax,ymax;
    double recphi;
    double exact;
    //------------------------------------------------------------------------//

   
    
    //------------------------------------------------------------------------//
    //Define quad points and weights here independent of what is in the rest of the code
   
    int extra = 1;
    
    double **z, **w;
    int ngauss = pow(polyorder + 1 + extra, 2);

    allocator2(&z, ngauss,2);
    allocator2(&w, ngauss,2);
    
    GaussPoints2D(z, w, quadtype, ngauss);

    double **xs;
    allocator2(&xs, ngauss, 2);
    
    //------------------------------------------------------------------------//

    double elemsum1;
    double elemsum2;
    double totalerr1 = 0.0;
    double totalerr2 = 0.0;
    double err1, err2;
    
    //------------------------------------------------------------------------//
    //Loop over elements
    for(ielem=2; ielem<xelem-2; ielem++)
    {
	for(jelem=2; jelem<yelem-2; jelem++)
	{
	    //INitialize
	    elemsum1 = 0.0;
	    elemsum2 = 0.0;
	    
	    //Get cartesian coordinates of Quad points
	    naturalToCartesian(xs, x, y, ielem, jelem, z, ngauss);

	    //Decide whether to calculate error in cell
	    xmax = 0.0;
	    ymax = 0.0;
	    for(igauss=0; igauss<ngauss; igauss++)
	    {
		xmax = max(xmax, fabs(xs[igauss][0]));
		ymax = max(ymax, fabs(xs[igauss][1]));
	    }

	    if(xmax < 1.0+(simtime*2.0/PI) && ymax < 1.0)
	    {
		//Loop over the Gauss quad points
		for(igauss=0; igauss<ngauss; igauss++)
		{
		    //Get the exact soln at this Gauss point
		    exact = getphi(xs[igauss][0],xs[igauss][1],simtime);

		    //Get the basis at this point
		    basis2D(z[igauss][0], z[igauss][1], basis);

		    //Get the Jacobian
		    detJ = mappingJacobianDeterminant(ielem, jelem, z[igauss][0], z[igauss][1], x, y, inv, jacobian);
		    
		    //REconstruct the soln at this gauss point
		    recphi = 0.0;
		    for(icoeff=0; icoeff<ncoeff; icoeff++)
		    {
			recphi += basis[icoeff]*phi[ielem][jelem][icoeff];
		    }

		    err1 = fabs(recphi - exact);
		    err2 = pow(recphi-exact, 2.0);

		    elemsum1 += w[igauss][0]*w[igauss][1]*detJ*err1;
		    elemsum2 += w[igauss][0]*w[igauss][1]*detJ*err2;
		}
	    }

	    //Add contribution to total error
	    totalerr1 += elemsum1;
	    totalerr2 += elemsum2;
	}
    }
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Collect error from all procs
    double l1sum, l2sum;
    MPI_Allreduce(&totalerr1, &l1sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&totalerr2, &l2sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    l2sum = sqrt(l2sum);	
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //return error norms
    (*L1norm) = l1sum;
    (*L2norm) = l2sum;
    //------------------------------------------------------------------------//


    //------------------------------------------------------------------------//
    //Deallocators
    deallocator1(&basis, ncoeff);
    deallocator1(&inv, 4);
    deallocator1(&jacobian, 4);
    deallocator2(&z, ngauss,2);
    deallocator2(&w, ngauss,2);
    deallocator2(&xs, ngauss, 2);
    //------------------------------------------------------------------------//


}
