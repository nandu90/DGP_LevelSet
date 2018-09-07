/***************************************************************************

Author: nsaini
Created: 2018-07-24

***************************************************************************/

#include "common.h"
#include "memory.h"
#include "functions.h"
#include "polylib.h"

double minmod(double a, double b, double c)
{
    double result;
    
    if(a*b > 0.0 && b*c > 0.0)
    {
	if(a == 0.0)
	{
	    result = 0.0;
	}
	else
	{
	    result = (a/fabs(a))*min(fabs(a),min(fabs(b),fabs(c)));
	}
    }
    else
    {
	result = 0.0;
    }

    return result;
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

    int ngauss = polyorder+10;
    double *z, *w;
    allocator1(&z, ngauss);
    allocator1(&w, ngauss);
    zwgll(z, w, ngauss);
    
    double *inv,  *jacobian;
    allocator1(&inv, 1);
    allocator1(&jacobian, 1);

    double detJ;
    
    for(ielem =1; ielem<xelem-1; ielem++)
    {
	elemsum = 0.0;
	elemsum1 = 0.0;
	for(igauss=0; igauss<ngauss; igauss++)
	{
	    recini = 0.0;
	    rec = 0.0;
	    basis1D(z[igauss], basis);
	    for(icoeff=0; icoeff<ncoeff; icoeff++)
	    {
		recini += basis[icoeff]*iniphi[ielem][icoeff];
		rec += basis[icoeff]*phi[ielem][icoeff];
	    }
	    detJ = mappingJacobianDeterminant(ielem, z[igauss], x, inv, jacobian);
	    
	    elemsum += pow(recini-rec,2.0)*w[igauss]*detJ;
	    elemsum1 += pow(recini,2.0)*w[igauss]*detJ;
		
	    
	}
	sum += elemsum;
	sum1 += elemsum1;
    }						

    
    
    (*err) = sqrt(sum);///sqrt(exact);

    (*lerr) = log(*err);
    

    deallocator1(&basis, ncoeff);
    deallocator1(&inv, 1);
    deallocator1(&jacobian, 1);
    deallocator1(&z, ngauss);
    deallocator1(&w, ngauss);
    
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

    int ngauss = polyorder+10;
    double *z, *w;
    allocator1(&z, ngauss);
    allocator1(&w, ngauss);
    zwgll(z, w, ngauss);
    
    double *inv,  *jacobian;
    allocator1(&inv, 1);
    allocator1(&jacobian, 1);

    double detJ;
    
    for(ielem =1; ielem<xelem-1; ielem++)
    {
	elemsum = 0.0;
	elemsum1 = 0.0;
	for(igauss=0; igauss<ngauss; igauss++)
	{
	    recini = 0.0;
	    rec = 0.0;
	    basis1D(z[igauss], basis);
	    for(icoeff=0; icoeff<ncoeff; icoeff++)
	    {
		recini += basis[icoeff]*iniphi[ielem][icoeff];
		rec += basis[icoeff]*phi[ielem][icoeff];
	    }
	    detJ = mappingJacobianDeterminant(ielem, z[igauss], x, inv, jacobian);
	    //if(recini <= 2.0)
	    //{
	    elemsum += fabs(recini-rec)*w[igauss]*detJ;
	    elemsum1 += fabs(recini)*w[igauss]*detJ;
	    //}
	    
	}
	sum += elemsum;
	sum1 += elemsum1;
	
    }
    

    (*err) = (sum);///(exact);

    deallocator1(&basis, ncoeff);
    deallocator1(&inv, 1);
    deallocator1(&jacobian, 1);
    deallocator1(&z, ngauss);
    deallocator1(&w, ngauss);
}
