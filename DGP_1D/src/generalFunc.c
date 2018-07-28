/***************************************************************************

Author: nsaini
Created: 2018-07-24

***************************************************************************/

#include "common.h"
#include "memory.h"
#include "functions.h"

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


void errorNormL2(double **iniphi, double **phi, double *err, double *lerr, double *x)
{
    int ielem, icoeff;
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
    allocator1(&inv, 1);
    allocator1(&jacobian, 1);

    double detJ;
    
    for(ielem =1; ielem<xelem-1; ielem++)
    {
	elemsum = 0.0;
	elemsum1 = 0.0;
	for(igauss=0; igauss<tgauss; igauss++)
	{
	    recini = 0.0;
	    rec = 0.0;
	    basis1D(zeta[igauss], basis);
	    for(icoeff=0; icoeff<ncoeff; icoeff++)
	    {
		recini += basis[icoeff]*iniphi[ielem][icoeff];
		rec += basis[icoeff]*phi[ielem][icoeff];
	    }
	    detJ = mappingJacobianDeterminant(ielem, zeta[igauss], x, inv, jacobian);
	    
	    elemsum += pow(recini-rec,2.0)*weights[igauss]*detJ;
	    elemsum1 += pow(recini,2.0)*weights[igauss]*detJ;
		
	    
	}
	sum += elemsum;
	sum1 += elemsum1;
    }						

    
    
    (*err) = sqrt(sum);///sqrt(exact);

    (*lerr) = log(*err);
    

    deallocator1(&basis, ncoeff);
    deallocator1(&inv, 1);
    deallocator1(&jacobian, 1);
    
    
}


void errorNormL1(double **iniphi, double **phi, double *err, double *lerr, double *x)
{
    int ielem, icoeff;
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
    allocator1(&inv, 1);
    allocator1(&jacobian, 1);

    double detJ;
    
    for(ielem =1; ielem<xelem-1; ielem++)
    {
	elemsum = 0.0;
	elemsum1 = 0.0;
	for(igauss=0; igauss<tgauss; igauss++)
	{
	    recini = 0.0;
	    rec = 0.0;
	    basis1D(zeta[igauss], basis);
	    for(icoeff=0; icoeff<ncoeff; icoeff++)
	    {
		recini += basis[icoeff]*iniphi[ielem][icoeff];
		rec += basis[icoeff]*phi[ielem][icoeff];
	    }
	    detJ = mappingJacobianDeterminant(ielem, zeta[igauss], x, inv, jacobian);
	    //if(recini <= 2.0)
	    //{
	    elemsum += fabs(recini-rec)*weights[igauss]*detJ;
	    elemsum1 += fabs(recini)*weights[igauss]*detJ;
	    //}
	    
	}
	sum += elemsum;
	sum1 += elemsum1;
	
    }
    

    (*err) = (sum);///(exact);

    deallocator1(&basis, ncoeff);
    deallocator1(&inv, 1);
    deallocator1(&jacobian, 1);
}


void errExact(double **phi, double *x, double time, int iter)
{
    int ielem, icoeff;
    int igauss;

    double *basis;
    allocator1(&basis, ncoeff);

    double rec;
        
    double sum = 0.0;
    double elemsum = 0.0;
    
    double *inv,  *jacobian;
    allocator1(&inv, 1);
    allocator1(&jacobian, 1);

    //double detJ;

    double exact;

    double *xs;
    allocator1(&xs, tgauss);

    double xtemp;

    for(ielem = 1; ielem<xelem-1; ielem++)
    {
	elemsum = 0.0;
	
	naturalToCartesian(xs,x,ielem);
	
	for(igauss =0; igauss<tgauss; igauss++)
	{
	    rec = 0.0;
	    basis1D(zeta[igauss], basis);
	    for(icoeff=0; icoeff<ncoeff; icoeff++)
	    {
		rec += basis[icoeff]*phi[ielem][icoeff];
	    }

	    xtemp = xs[igauss] - time;
	    if(xtemp < 0.0)
	    {
		xtemp = xlen + xtemp;
	    }
	    else if(xtemp > xlen)
	    {
		xtemp = xtemp - xlen;
	    }
	    double sigmax = 25.0;
	    double term1 = 0.5*pow((xtemp - xb_in)/sigmax,2.0);
	    exact = 1.0*exp(-(term1));
	    
	    elemsum += fabs(rec - exact);
	}

	sum += elemsum;
    }

    printf("The L1 norm at time %d is %.4e\n\n",iter,sum);
    
    deallocator1(&basis, ncoeff);
    deallocator1(&inv, 1);
    deallocator1(&jacobian, 1);
}
