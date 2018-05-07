#include "generalFunc.h"
#include "common.h"
#include "DGPFunc.h"
#include "memory.h"

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


void naturalToCartesian(double **xs, double **x, double **y, int i, int j)
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


void errorNormL2(double ***iniphi, double ***phi, double *err, double *lerr)
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
    
    for(ielem =2; ielem<xelem-2; ielem++)
    {
	for(jelem=2; jelem<yelem-2; jelem++)
	{
	    elemsum = 0.0;
	    for(igauss=0; igauss<tgauss; igauss++)
	    {
		recini = 0.0;
		rec = 0.0;
		basis2D(zeta[igauss][0], zeta[igauss][1], basis);
		for(icoeff=0; icoeff<ncoeff; icoeff++)
		{
		    recini += basis[icoeff]*iniphi[ielem][jelem][icoeff];
		    rec += basis[icoeff]*phi[ielem][jelem][icoeff];
		}
		elemsum += pow(recini-rec,2.0)*weights[igauss][0]*weights[igauss][1];
		sum1 += pow(recini,2.0);
	    }
	    sum += elemsum;
	}
    }

    MPI_Allreduce(&sum, err, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    double result;

    MPI_Allreduce(&sum1, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		  
    (*err) = sqrt((*err));

    (*lerr) = log(*err);
    

    deallocator1(&basis, ncoeff);
    
}


void errorNormL1(double ***iniphi, double ***phi, double *err, double *lerr)
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
    
    for(ielem =2; ielem<xelem-2; ielem++)
    {
	for(jelem=2; jelem<yelem-2; jelem++)
	{
	    for(igauss=0; igauss<tgauss; igauss++)
	    {
		recini = 0.0;
		rec = 0.0;
		basis2D(zeta[igauss][0], zeta[igauss][1], basis);
		for(icoeff=0; icoeff<ncoeff; icoeff++)
		{
		    recini += basis[icoeff]*iniphi[ielem][jelem][icoeff];
		    rec += basis[icoeff]*phi[ielem][jelem][icoeff];
		}
		sum += fabs(recini-rec);
		sum1 += fabs(recini);
	    }
	    
	}
    }

    MPI_Allreduce(&sum, err, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    double result;

    MPI_Allreduce(&sum1, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		  
    (*err) = sqrt((*err));

    (*lerr) = log(*err);
    

    deallocator1(&basis, ncoeff);
    
}


void calc_vf(double ***H, double detJ, double *vf)
{
    //------------------------------------------------------------------------//
    //Loop indexes
    int ielem, jelem;
    int igauss;
    //------------------------------------------------------------------------//
    *vf = 0.0;
    double sum = 0.0;
    
    for(ielem =2; ielem<xelem-2; ielem++)
    {
	for(jelem=2; jelem<yelem-2; jelem++)
	{
	    sum = 0.0;
	    for(igauss=0; igauss<tgauss; igauss++)
	    {
		sum += weights[igauss][0]*weights[igauss][1]*(1.0-H[ielem][jelem][igauss]);
	    }
	    *vf += (sum) * detJ;
	}
    }

    double totalvf;

    MPI_Allreduce(vf, &totalvf, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    *vf = totalvf;
    
    if(myrank == master)
    {
	printf("Total area is %.4e\n",totalvf);
    }
}
