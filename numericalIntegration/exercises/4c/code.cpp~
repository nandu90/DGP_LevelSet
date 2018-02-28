#include <iostream>
#include <stdio.h>
#include <math.h>
#include "../../polylib.h"
#include <stdlib.h>



using namespace std;
using namespace polylib;

double *dvector(int np);
void *func(int np, double *z, double *p, double y);
double integr(int np, double *w, double *p);


main()
{

  int np[2] = {6, 6};

  double *z1,*w1,*p1;
  double *z2,*w2,*p2;
  
  z1 = dvector(np[0]);
  w1 = dvector(np[0]);
  p1 = dvector(np[0]);
  
  z2 = dvector(np[1]);
  w2 = dvector(np[1]);
  p2 = dvector(np[1]);

  double sum = 0.0;
  
      // get zeros and weights
  zwgll(z1,w1, np[0]);
  zwgll(z2,w2, np[1]);
    

  double sum1 = 0.0;
  double sum2 = 0.0;
      
  for(int j=0; j<np[0]; j++)
    {
      sum1 = 0.0;
      func(np[1], z2, p2, z1[j]);
      
      // integrate p over [-1,1] 
      sum1 = integr(np[1], w2, p2);

      sum2 = sum2 + sum1*w1[j];
    }
  
  cout << "\nIntegral = " << sum2  << "\n";
   

}

double *dvector(int n)
{
  double *v;

  v = (double *)malloc(n*sizeof(double));
  return v;
}

void *func(int np, double *z, double *p, double y)
{  
  
   for(int i = 0;i<np;i++)
     {
       p[i] = pow(z[i],6)*pow(y,6);
     } 
}

double integr(int np, double *w, double *p)
{
  register double sum = 0.;

  for(int i=0;i<np;i++){
    sum = sum + p[i]*w[i]; 
  }
  return sum;
}
