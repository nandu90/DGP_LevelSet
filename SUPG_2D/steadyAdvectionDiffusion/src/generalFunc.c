#include "generalFunc.h"
#include "common.h"
#include "memory.h"
#include "polylib.h"

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
