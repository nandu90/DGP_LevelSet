/***************************************************************************

Author: nsaini
Created: 2018-03-26
Exercise on numerical differentiation in a standard quadrilateral element
Exercise 5.6.1.a
***************************************************************************/


#include <iostream>
#include <stdio.h>
#include <math.h>
#include "../polylib.h"
#include <stdlib.h>




using namespace std;
using namespace polylib;

double *dvector(int np);


main()
{

    int np = 2;
    
  

  double *z1,*w1,*p1;
  double *z2,*w2,*p2;
  
  z1 = dvector(np);
  w1 = dvector(np);
  p1 = dvector(np);
  
  z2 = dvector(np);
  w2 = dvector(np);
  p2 = dvector(np);

  double sum = 0.0;
  
      // get zeros and weights
  zwgl(z1,w1, np);
  zwgl(z2,w2, np);

  cout<<"Gauss Weights and zeros in zeta 1 direction are:"<<endl;
  for(int i=0; i<np; i++)
    {
      cout<<w1[i]<<" "<<z1[i]<<endl;
    }

   cout<<"Gauss Weights and zeros in zeta 2 direction are:"<<endl;
  for(int i=0; i<np; i++)
    {
      cout<<w2[i]<<" "<<z2[i]<<endl;
    }

  cout<<endl;


  //For Legendre polynomials
  double alpha = 0.0;
  double beta = 0.0;

  //------------------------------------------------------------------------//
  // get the differentiation matrix
  double *D;
  double *Dt;
  

  D = dvector(pow(np,2.0));
  Dt = dvector(pow(np,2.0));
  double Dmat[np][np];
  double n1,n2;
  double exact;

  int i,j,k,l,m;

  int ig, jg;
  double maxerror = 100.0;

  //Evaluating diff w.r.t zeta1
  Dgj(D, Dt, z1, np, alpha, beta);
  m=0;
  for(k=0; k<np; k++)
  {
      for(l=0; l<np; l++)
      {
	  Dmat[k][l] = D[m++];
      }
  }
  
  for(ig = 0; ig<np; ig++)
  {
      for(jg = 0; jg<np; jg++)
      {
	  n1 = z1[ig];
	  n2 = z2[jg];
	  
	  
	  /*cout<<"The mtrix elements are: "<<endl;
	    for(i=0; i<np*np;i++)
	    {
	    cout<<D[i]<<" ";
	    }
	    cout<<endl;*/
	  
	  
	  sum = 0.0;
	  
	  
	  
	  for(j=0; j<np; j++)
	  {	  
	    sum += Dmat[ig][j]*z1[j];//pow(z1[j],7.0);
	  }
	  sum = sum * n2;// pow(n2,6.0);
	  
	  cout<<"The differential for point "<<n1<<", "<<n2<<" is "<<sum<<endl;
	  
	  exact = 7.0*pow(n1,6.0)*pow(n2,6.0);
	  cout<<"The error in differential is "<<exact-sum<<endl<<endl;

	  if(abs(exact-sum) < abs(maxerror))
	  {
	      maxerror = exact-sum;
	  }
	  
      }
      
  
  }

  cout<<"The maximum error is "<<maxerror<<endl;
  

  

}

double *dvector(int n)
{
  double *v;

  v = (double *)malloc(n*sizeof(double));
  return v;
}


