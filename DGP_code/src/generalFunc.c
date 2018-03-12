#include "generalFunc.h"
#include "common.h"

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


void getvertices(double **xs,double **zs, double **x, double **y, int i, int j)
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

    printf("Vertex points are:\n");
    for(k=0; k<4; k++)
    {
	printf("%.2f %.2f\n",xvertex[k],yvertex[k]);
    }

    printf("\n\n");
    
    double N1, N2, N3, N4;
    
    //Populate the coordinate vector
    for(k=0; k<pow(polyorder+1,2); k++)
    {
	N1 = (1.0-zs[k][0])*(1.0-zs[k][1])/4.0;
	N2 = (1.0+zs[k][0])*(1.0-zs[k][1])/4.0;
	N3 = (1.0-zs[k][0])*(1.0+zs[k][1])/4.0;
	N4 = (1.0+zs[k][0])*(1.0+zs[k][1])/4.0;
	
	
	xs[k][0] = N1*xvertex[0] + N2*xvertex[1] + N3*xvertex[2] + N4*xvertex[3];
	xs[k][1] = N1*yvertex[0] + N2*yvertex[1] + N3*yvertex[2] + N4*yvertex[3];
    }
}

void getSolnNaturalCoord(double **zs)
{
    int i;
    //------------------------------------------------------------------------//
    //Gauss-Lobatto quadrature
    if(quadtype == 1)
    {
	if(polyorder == 1) //Only vertices are required
	{
	    zs[0][0] = zeta[0][0];
	    zs[0][1] = zeta[0][1];
	    
	    zs[1][0] = zeta[2][0];
	    zs[1][1] = zeta[2][1];
	    
	    zs[2][0] = zeta[6][0];
	    zs[2][1] = zeta[6][1];
	    
	    zs[3][0] = zeta[8][0];
	    zs[3][1] = zeta[8][1];
	}
	else if(polyorder == 2) //This is messed up. Be careful about points selected
	{
	    //Vertices - anticlockwise
	    zs[0][0] = zeta[0][0];
	    zs[0][1] = zeta[0][1];
	    
	    zs[1][0] = zeta[3][0];
	    zs[1][1] = zeta[3][1];
	    
	    zs[2][0] = zeta[4][0];
	    zs[2][1] = zeta[4][1];
	    
	    zs[3][0] = zeta[7][0];
	    zs[3][1] = zeta[7][1];
	    
	    zs[4][0] = zeta[8][0];
	    zs[4][1] = zeta[8][1];
	    
	    zs[5][0] = zeta[11][0];
	    zs[5][1] = zeta[11][1];
	    
	    zs[6][0] = zeta[12][0];
	    zs[6][1] = zeta[12][1];
	    
	    zs[7][0] = zeta[13][0];
	    zs[7][1] = zeta[13][1];
	    
	    zs[8][0] = zeta[15][0];
	    zs[8][1] = zeta[15][1];
	}
    }
    else if(quadtype == 2) // Gauss-Legendre
    {
	for(i=0; i<tgauss; i++)
	{
	    zs[i][0] = zeta[i][0];
	    zs[i][1] = zeta[i][1];
	}
    }
}


