#include <iostream>
#include <stdio.h>
#include <math.h>
#include "polylib.h"
#include <stdlib.h>



using namespace std;
using namespace polylib;

double *dvector(int np);
void *func(int np, double *z, double *p, double y, int r, int c);
double integr(int np, double *w, double *p);


main()
{

  int np[2] = {2, 2};

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
  zwgl(z1,w1, np[0]);
  zwgl(z2,w2, np[1]);

  cout<<"Gauss Weights and zeros in zeta 1 direction are:"<<endl;
  for(int i=0; i<np[0]; i++)
    {
      cout<<w1[i]<<" "<<z1[i]<<endl;
    }

   cout<<"Gauss Weights and zeros in zeta 2 direction are:"<<endl;
  for(int i=0; i<np[1]; i++)
    {
      cout<<w2[i]<<" "<<z2[i]<<endl;
    }

  cout<<endl;
    

  double sum1 = 0.0;
  double sum2 = 0.0;

  
  for(int r=0; r<4; r++)
    {
      for(int c=0; c<4; c++)
	{
	  sum2 = 0.0;
	  for(int j=0; j<np[0]; j++)
	    {
	      sum1 = 0.0;
	      func(np[1], z2, p2, z1[j],r+1,c+1);
	      
	      // integrate p over [-1,1] 
	      sum1 = integr(np[1], w2, p2);
	      
	      sum2 = sum2 + sum1*w1[j];
	    }
	  
	  //cout << "\nIntegral for Row "<<r+1<<" and column "<<c+1<<" = "<< sum2  << "\n";
	  if(fabs(sum2)<1e-10)sum2 = 0.0;
	  cout<<sum2<<",           ";
	}
      cout<<"\n";
    }

}

double *dvector(int n)
{
  double *v;

  v = (double *)malloc(n*sizeof(double));
  return v;
}

void *func(int np, double *z2, double *p, double zeta1, int r, int c)
{  
  double x1 = 0.0;
  double x2 = 1.0;
  double x3 = 0.0;
  double x4 = 1.0;

  double y1 = 0.0;
  double y2 = 0.0;
  double y3 = 2.0;
  double y4 = 2.0;

  
  
  
   for(int i = 0;i<np;i++)
     {
       double zeta2 = z2[i];

       double detJ = x1*y2*(1.0/8.0)-x2*y1*(1.0/8.0)-x1*y3*(1.0/8.0)+x3*y1*(1.0/8.0)+x2*y4*(1.0/8.0)-x4*y2*(1.0/8.0)-x3*y4*(1.0/8.0)+x4*y3*(1.0/8.0)-x1*y2*zeta2*(1.0/8.0)+x1*y3*zeta1*(1.0/8.0)+x2*y1*zeta2*(1.0/8.0)-x3*y1*zeta1*(1.0/8.0)-x1*y4*zeta1*(1.0/8.0)-x2*y3*zeta1*(1.0/8.0)+x3*y2*zeta1*(1.0/8.0)+x4*y1*zeta1*(1.0/8.0)+x1*y4*zeta2*(1.0/8.0)-x2*y3*zeta2*(1.0/8.0)+x2*y4*zeta1*(1.0/8.0)+x3*y2*zeta2*(1.0/8.0)-x4*y1*zeta2*(1.0/8.0)-x4*y2*zeta1*(1.0/8.0)-x3*y4*zeta2*(1.0/8.0)+x4*y3*zeta2*(1.0/8.0);
       
       if(r==1 && c==1)
	 {
	   p[i] = 1.0;
	 }
       
       else if(r == 1 && c==2)
	 {
	   p[i] = zeta1;
	 }

       else if(r==1 && c==3)
	 {
	   p[i] = zeta2;
	 }

       else if(r==1 && c==4)
	 {
	   p[i] = zeta1*zeta2;
	 }

       else if(r==2 && c==1)
	 {
	   p[i] = zeta1;
	 }

       else if(r==2 && c==2)
	 {
	   p[i] = zeta1*zeta1;
	 }

       else if(r==2 && c==3)
	 {
	   p[i] = zeta1*zeta2;
	 }

       else if(r==2 && c==4)
	 {
	   p[i] = zeta1*zeta1*zeta2;
	 }

       else if(r==3 && c==1)
	 {
	   p[i] = zeta2;
	 }

       else if(r==3 && c==2)
	 {
	   p[i] = zeta1*zeta2;
	 }

       else if(r==3 && c==3)
	 {
	   p[i] = zeta2*zeta2;
	 }

       else if(r==3 && c==4)
	 {
	   p[i] = zeta1*zeta2*zeta2;
	 }

       else if(r==4 && c==1)
	 {
	   p[i] = zeta1*zeta2;
	 }

       else if(r==4 && c==2)
	 {
	   p[i] = zeta1*zeta1*zeta2;
	 }

       else if(r==4 && c==3)
	 {
	   p[i] = zeta1*zeta2*zeta2;
	 }

       else if(r==4 && c==4)
	 {
	   p[i] = zeta1*zeta1*zeta2*zeta2;
	 }

       p[i] = p[i]*detJ;

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
