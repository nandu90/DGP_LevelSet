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

void *func(int np, double *z2, double *p, double zeta1)
{  
  double x1 = 0.0;
  double x2 = 1.0;
  double x3 = 0.0;
  double x4 = 2.0;

  double y1 = 0.0;
  double y2 = 0.0;
  double y3 = 1.0;
  double y4 = 1.0;

  
  
   for(int i = 0;i<np;i++)
     {
       double zeta2 = z2[i];

       double detJ = x1*y2*(1.0/8.0)-x2*y1*(1.0/8.0)-x1*y3*(1.0/8.0)+x3*y1*(1.0/8.0)+x2*y4*(1.0/8.0)-x4*y2*(1.0/8.0)-x3*y4*(1.0/8.0)+x4*y3*(1.0/8.0)-x1*y2*zeta2*(1.0/8.0)+x1*y3*zeta1*(1.0/8.0)+x2*y1*zeta2*(1.0/8.0)-x3*y1*zeta1*(1.0/8.0)-x1*y4*zeta1*(1.0/8.0)-x2*y3*zeta1*(1.0/8.0)+x3*y2*zeta1*(1.0/8.0)+x4*y1*zeta1*(1.0/8.0)+x1*y4*zeta2*(1.0/8.0)-x2*y3*zeta2*(1.0/8.0)+x2*y4*zeta1*(1.0/8.0)+x3*y2*zeta2*(1.0/8.0)-x4*y1*zeta2*(1.0/8.0)-x4*y2*zeta1*(1.0/8.0)-x3*y4*zeta2*(1.0/8.0)+x4*y3*zeta2*(1.0/8.0);

       cout<<detJ<<endl;
       p[i] = pow(z2[i],6)*pow(zeta1,6)*detJ;
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
