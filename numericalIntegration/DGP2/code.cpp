#include <iostream>
#include <stdio.h>
#include <math.h>
#include "../polylib.h"
#include <stdlib.h>



using namespace std;
using namespace polylib;

double *dvector(int np);
void *func(int np, double *z, double *p, double y, int r, int c);
double integr(int np, double *w, double *p);


main()
{

  int np[2] = {4, 4};

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

  
  for(int r=0; r<9; r++)
    {
      for(int c=0; c<9; c++)
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
	  if(c==r)cout<<sum2<<",           ";
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
  double x3 = 1.0;
  double x4 = 0.0;

  double y1 = 0.0;
  double y2 = 0.0;
  double y3 = 1.0;
  double y4 = 1.0;

  double M[9][9];
  
  
   for(int i = 0;i<np;i++)
     {
       double zeta2 = z2[i];

       double detJ = x1*y2*(1.0/8.0)-x2*y1*(1.0/8.0)-x1*y4*(1.0/8.0)+x2*y3*(1.0/8.0)-x3*y2*(1.0/8.0)+x4*y1*(1.0/8.0)+x3*y4*(1.0/8.0)-x4*y3*(1.0/8.0)-x1*y2*zeta2*(1.0/8.0)-x1*y3*zeta1*(1.0/8.0)+x2*y1*zeta2*(1.0/8.0)+x3*y1*zeta1*(1.0/8.0)+x1*y3*zeta2*(1.0/8.0)+x1*y4*zeta1*(1.0/8.0)+x2*y3*zeta1*(1.0/8.0)-x3*y1*zeta2*(1.0/8.0)-x3*y2*zeta1*(1.0/8.0)-x4*y1*zeta1*(1.0/8.0)-x2*y4*zeta1*(1.0/8.0)+x4*y2*zeta1*(1.0/8.0)-x2*y4*zeta2*(1.0/8.0)+x4*y2*zeta2*(1.0/8.0)+x3*y4*zeta2*(1.0/8.0)-x4*y3*zeta2*(1.0/8.0);

       M[0][0] = 1.0;
       M[0][1] = zeta1;
       M[0][2] = zeta2;
       M[0][3] = zeta1*zeta2;
       M[0][4] = (zeta1*zeta1)*(3.0/2.0)-1.0/2.0;
       M[0][5] = (zeta2*zeta2)*(3.0/2.0)-1.0/2.0;
       M[0][6] = zeta2*((zeta1*zeta1)*3.0-1.0)*(1.0/2.0);
       M[0][7] = zeta1*((zeta2*zeta2)*3.0-1.0)*(1.0/2.0);
       M[0][8] = ((zeta2*zeta2)*3.0-1.0)*((zeta1*zeta1)*(3.0/2.0)-1.0/2.0)*(1.0/2.0);
       M[1][0] = zeta1;
       M[1][1] = zeta1*zeta1;
       M[1][2] = zeta1*zeta2;
       M[1][3] = (zeta1*zeta1)*zeta2;
       M[1][4] = zeta1*((zeta1*zeta1)*(3.0/2.0)-1.0/2.0);
       M[1][5] = zeta1*((zeta2*zeta2)*(3.0/2.0)-1.0/2.0);
       M[1][6] = zeta1*zeta2*((zeta1*zeta1)*3.0-1.0)*(1.0/2.0);
       M[1][7] = (zeta1*zeta1)*((zeta2*zeta2)*3.0-1.0)*(1.0/2.0);
       M[1][8] = zeta1*((zeta2*zeta2)*3.0-1.0)*((zeta1*zeta1)*(3.0/2.0)-1.0/2.0)*(1.0/2.0);
       M[2][0] = zeta2;
       M[2][1] = zeta1*zeta2;
       M[2][2] = zeta2*zeta2;
       M[2][3] = zeta1*(zeta2*zeta2);
       M[2][4] = zeta2*((zeta1*zeta1)*(3.0/2.0)-1.0/2.0);
       M[2][5] = zeta2*((zeta2*zeta2)*(3.0/2.0)-1.0/2.0);
       M[2][6] = (zeta2*zeta2)*((zeta1*zeta1)*3.0-1.0)*(1.0/2.0);
       M[2][7] = zeta1*zeta2*((zeta2*zeta2)*3.0-1.0)*(1.0/2.0);
       M[2][8] = zeta2*((zeta2*zeta2)*3.0-1.0)*((zeta1*zeta1)*(3.0/2.0)-1.0/2.0)*(1.0/2.0);
       M[3][0] = zeta1*zeta2;
       M[3][1] = (zeta1*zeta1)*zeta2;
       M[3][2] = zeta1*(zeta2*zeta2);
       M[3][3] = (zeta1*zeta1)*(zeta2*zeta2);
       M[3][4] = zeta1*zeta2*((zeta1*zeta1)*(3.0/2.0)-1.0/2.0);
       M[3][5] = zeta1*zeta2*((zeta2*zeta2)*(3.0/2.0)-1.0/2.0);
       M[3][6] = zeta1*(zeta2*zeta2)*((zeta1*zeta1)*3.0-1.0)*(1.0/2.0);
       M[3][7] = (zeta1*zeta1)*zeta2*((zeta2*zeta2)*3.0-1.0)*(1.0/2.0);
       M[3][8] = zeta1*zeta2*((zeta2*zeta2)*3.0-1.0)*((zeta1*zeta1)*(3.0/2.0)-1.0/2.0)*(1.0/2.0);
       M[4][0] = (zeta1*zeta1)*(3.0/2.0)-1.0/2.0;
       M[4][1] = zeta1*((zeta1*zeta1)*(3.0/2.0)-1.0/2.0);
       M[4][2] = zeta2*((zeta1*zeta1)*(3.0/2.0)-1.0/2.0);
       M[4][3] = zeta1*zeta2*((zeta1*zeta1)*(3.0/2.0)-1.0/2.0);
       M[4][4] = pow((zeta1*zeta1)*(3.0/2.0)-1.0/2.0,2.0);
       M[4][5] = ((zeta1*zeta1)*(3.0/2.0)-1.0/2.0)*((zeta2*zeta2)*(3.0/2.0)-1.0/2.0);
       M[4][6] = zeta2*((zeta1*zeta1)*3.0-1.0)*((zeta1*zeta1)*(3.0/2.0)-1.0/2.0)*(1.0/2.0);
       M[4][7] = zeta1*((zeta2*zeta2)*3.0-1.0)*((zeta1*zeta1)*(3.0/2.0)-1.0/2.0)*(1.0/2.0);
       M[4][8] = ((zeta2*zeta2)*3.0-1.0)*pow((zeta1*zeta1)*(3.0/2.0)-1.0/2.0,2.0)*(1.0/2.0);
       M[5][0] = (zeta2*zeta2)*(3.0/2.0)-1.0/2.0;
       M[5][1] = zeta1*((zeta2*zeta2)*(3.0/2.0)-1.0/2.0);
       M[5][2] = zeta2*((zeta2*zeta2)*(3.0/2.0)-1.0/2.0);
       M[5][3] = zeta1*zeta2*((zeta2*zeta2)*(3.0/2.0)-1.0/2.0);
       M[5][4] = ((zeta1*zeta1)*(3.0/2.0)-1.0/2.0)*((zeta2*zeta2)*(3.0/2.0)-1.0/2.0);
       M[5][5] = pow((zeta2*zeta2)*(3.0/2.0)-1.0/2.0,2.0);
       M[5][6] = zeta2*((zeta1*zeta1)*3.0-1.0)*((zeta2*zeta2)*(3.0/2.0)-1.0/2.0)*(1.0/2.0);
       M[5][7] = zeta1*((zeta2*zeta2)*3.0-1.0)*((zeta2*zeta2)*(3.0/2.0)-1.0/2.0)*(1.0/2.0);
       M[5][8] = ((zeta2*zeta2)*3.0-1.0)*((zeta1*zeta1)*(3.0/2.0)-1.0/2.0)*((zeta2*zeta2)*(3.0/2.0)-1.0/2.0)*(1.0/2.0);
       M[6][0] = zeta2*((zeta1*zeta1)*3.0-1.0)*(1.0/2.0);
       M[6][1] = zeta1*zeta2*((zeta1*zeta1)*3.0-1.0)*(1.0/2.0);
       M[6][2] = (zeta2*zeta2)*((zeta1*zeta1)*3.0-1.0)*(1.0/2.0);
       M[6][3] = zeta1*(zeta2*zeta2)*((zeta1*zeta1)*3.0-1.0)*(1.0/2.0);
       M[6][4] = zeta2*((zeta1*zeta1)*3.0-1.0)*((zeta1*zeta1)*(3.0/2.0)-1.0/2.0)*(1.0/2.0);
       M[6][5] = zeta2*((zeta1*zeta1)*3.0-1.0)*((zeta2*zeta2)*(3.0/2.0)-1.0/2.0)*(1.0/2.0);
       M[6][6] = (zeta2*zeta2)*pow((zeta1*zeta1)*3.0-1.0,2.0)*(1.0/4.0);
       M[6][7] = zeta1*zeta2*((zeta1*zeta1)*3.0-1.0)*((zeta2*zeta2)*3.0-1.0)*(1.0/4.0);
       M[6][8] = zeta2*((zeta1*zeta1)*3.0-1.0)*((zeta2*zeta2)*3.0-1.0)*((zeta1*zeta1)*(3.0/2.0)-1.0/2.0)*(1.0/4.0);
       M[7][0] = zeta1*((zeta2*zeta2)*3.0-1.0)*(1.0/2.0);
       M[7][1] = (zeta1*zeta1)*((zeta2*zeta2)*3.0-1.0)*(1.0/2.0);
       M[7][2] = zeta1*zeta2*((zeta2*zeta2)*3.0-1.0)*(1.0/2.0);
       M[7][3] = (zeta1*zeta1)*zeta2*((zeta2*zeta2)*3.0-1.0)*(1.0/2.0);
       M[7][4] = zeta1*((zeta2*zeta2)*3.0-1.0)*((zeta1*zeta1)*(3.0/2.0)-1.0/2.0)*(1.0/2.0);
       M[7][5] = zeta1*((zeta2*zeta2)*3.0-1.0)*((zeta2*zeta2)*(3.0/2.0)-1.0/2.0)*(1.0/2.0);
       M[7][6] = zeta1*zeta2*((zeta1*zeta1)*3.0-1.0)*((zeta2*zeta2)*3.0-1.0)*(1.0/4.0);
       M[7][7] = (zeta1*zeta1)*pow((zeta2*zeta2)*3.0-1.0,2.0)*(1.0/4.0);
       M[7][8] = zeta1*pow((zeta2*zeta2)*3.0-1.0,2.0)*((zeta1*zeta1)*(3.0/2.0)-1.0/2.0)*(1.0/4.0);
       M[8][0] = ((zeta2*zeta2)*3.0-1.0)*((zeta1*zeta1)*(3.0/2.0)-1.0/2.0)*(1.0/2.0);
       M[8][1] = zeta1*((zeta2*zeta2)*3.0-1.0)*((zeta1*zeta1)*(3.0/2.0)-1.0/2.0)*(1.0/2.0);
       M[8][2] = zeta2*((zeta2*zeta2)*3.0-1.0)*((zeta1*zeta1)*(3.0/2.0)-1.0/2.0)*(1.0/2.0);
       M[8][3] = zeta1*zeta2*((zeta2*zeta2)*3.0-1.0)*((zeta1*zeta1)*(3.0/2.0)-1.0/2.0)*(1.0/2.0);
       M[8][4] = ((zeta2*zeta2)*3.0-1.0)*pow((zeta1*zeta1)*(3.0/2.0)-1.0/2.0,2.0)*(1.0/2.0);
       M[8][5] = ((zeta2*zeta2)*3.0-1.0)*((zeta1*zeta1)*(3.0/2.0)-1.0/2.0)*((zeta2*zeta2)*(3.0/2.0)-1.0/2.0)*(1.0/2.0);
       M[8][6] = zeta2*((zeta1*zeta1)*3.0-1.0)*((zeta2*zeta2)*3.0-1.0)*((zeta1*zeta1)*(3.0/2.0)-1.0/2.0)*(1.0/4.0);
       M[8][7] = zeta1*pow((zeta2*zeta2)*3.0-1.0,2.0)*((zeta1*zeta1)*(3.0/2.0)-1.0/2.0)*(1.0/4.0);
       M[8][8] = pow((zeta2*zeta2)*3.0-1.0,2.0)*pow((zeta1*zeta1)*(3.0/2.0)-1.0/2.0,2.0)*(1.0/4.0);


       p[i] = M[r-1][c-1];
	

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
