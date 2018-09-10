#include "memory.h"
#include "common.h"

void allocator1(double **p, int x)
{
  *p = (double *)malloc(x * sizeof(double));
}

void deallocator1(double **p, int x)
{
  free(*p);
}

void allocator2(double ***p, int x, int y)
{
  int i,j;
  *p = (double **)malloc(x * sizeof(double *));
  for(i=0; i<x; i++)
  {
    (*p)[i] =  (double *)malloc(y * sizeof(double));
  }

    for(i=0; i<x; i++)
    {
      for(j=0; j<y; j++)
	{
	  (*p)[i][j] = 0.0;
	}
    }	
}

void allocator3(double ****p, int x, int y, int z)
{
  int i,j,k;
  *p = (double ***)malloc(x * sizeof(double **));
  for(i=0; i<x; i++)
  {
    (*p)[i] =  (double **)malloc(y * sizeof(double *));
    for(j=0; j<y; j++)
      {
	(*p)[i][j] = (double *)malloc(z * sizeof(double));
      }
  }

  for(i=0; i<x; i++)
    {
      for(j=0; j<y; j++)
	{
	  for(k=0; k<z; k++)
	    {
	      (*p)[i][j][k] = 0.0;
	    }
	}
    }

}

void allocator4(double *****p, int x, int y, int z, int w)
{
  int i,j,k,l;
  *p = (double ****)malloc(x * sizeof(double ***));
  for(i=0; i<x; i++)
  {
    (*p)[i] =  (double ***)malloc(y * sizeof(double **));
    for(j=0; j<y; j++)
      {
	(*p)[i][j] = (double **)malloc(z * sizeof(double *));
	for(k=0; k<z; k++)
	  {
	    (*p)[i][j][k] = (double *) malloc(w * sizeof(double));
	  }
      }
  }

  for(i=0; i<x; i++)
    {
      for(j=0; j<y; j++)
	{
	  for(k=0; k<z; k++)
	    {
	      for(l=0; l<w; l++)
		{
		  (*p)[i][j][k][l] = 0.0;
		}
	    }
	}
    }

}


void iallocator1(int **p, int x)
{
  *p = (int *)malloc(x * sizeof(int));
}

void ideallocator1(int **p, int x)
{
  free(*p);
}


void iallocator2(int ***p, int x, int y)
{
  int i,j;
  (*p) = (int **)malloc(x * sizeof(int *));
  for(i=0; i<x; i++)
  {
    (*p)[i] =  (int *)malloc(y * sizeof(int));
  }

  for(i=0; i<x; i++)
    {
      for(j=0; j<y; j++)
	{
	  (*p)[i][j] = 0;
	}
    }
}

void deallocator2(double ***p, int x, int y)
{
  int i;
  for (i=0; i<x; i++)
    {
      free((*p)[i]);
    }
  free(*p);
}

void deallocator3(double ****p, int x, int y, int z)
{
  int i,j;
  for(i=0; i<x; i++)
    {
      for(j=0; j<y; j++)
	{
	  free((*p)[i][j]);
	}
      free((*p)[i]);
    }
  free(*p);
}

void deallocator4(double *****p, int x, int y, int z, int w)
{
  int i,j,k;
  for(i=0; i<x; i++)
    {
      for(j=0; j<y; j++)
	{
	  for(k=0; k<z; k++)
	    {
	      free((*p)[i][j][k]);
	    }
	  free((*p)[i][j]);
	}
      free((*p)[i]);
    }
  free(*p);
}

void ideallocator2(int ***p, int x, int y)
{
  int i;
  for (i=0; i<x; i++)
    {
      free((*p)[i]);
    }
  free(*p);
}
