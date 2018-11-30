/***************************************************************************

Author: nsaini
Created: 2018-03-04

***************************************************************************/

#include "common.h"
#include "mesh.h"
#include "memory.h"

void gridgen(double **x, double **y, double ****area, double **vol, double **xc, double **yc)
{
    int i,j;
    double deltax=xlen/(gxelem);
    double deltay=ylen/(gyelem);
    double firstx = (myrank%procm)*elemm*deltax;
    double firsty = (double)floor(myrank/procm)*elemn*deltay;
    //Assign coordinates
    for(j=2; j<ynode-2; j++)
    {
        for(i=2; i<xnode-2; i++)
        {
            x[i][j] = firstx + (i-2)*deltax;
            y[i][j] = firsty + (j-2)*deltay;
	    //printf("%d %d %.2f %.2f\n",i,j,x[i][j],y[i][j]); 
        }
    }

    meshOperations(x, y, area, vol, xc, yc);

    if(shiftorigin == 1)
    {
	meshShift(x,y,xc,yc);
    }
    //printf("X and y limits in %d with m x n = %d %d are: %.5f %.5f %.5f %.5f\n",myrank+1,xelem-2, yelem-2,x[0][0], x[xnode-1][0], y[0][0],y[0][ynode-1]);
}

void gridread(double **x, double **y, double ****area, double **vol, double **xc, double **yc)
{
    //------------------------------------------------------------------------//
    //Let only the master processor read the mesh
    double *xread, *yread;
    
    allocator1(&xread, (gxelem+1)*(gyelem+1));
    allocator1(&yread, (gxelem+1)*(gyelem+1));

    int i,j,k;
    
    if(myrank == master)
    {
	/*printf("%d %d\n",gxelem,gyelem);
	  exit(1);*/
	
	FILE *gridfile;
	gridfile = fopen("../../mesh.dat","r");
	if(gridfile == NULL)
	{
	    printf("Error opening mesh.dat!\n");
	    exit(0);
	}

	char* line = NULL;
	ssize_t size;
	size_t len = 0;
	
	char delim [1] = " ";
	char* word = NULL;

	i = 0;
	while((size = getline(&line, &len, gridfile)) != -1)
	{
	    word = strtok(line, delim);
	    
	    xread[i] = atof(word);
	    word = strtok(NULL,delim);
	    yread[i] = atof(word);
	    word = strtok(NULL,delim); 
	    
	    i++;
	}
	
	fclose(gridfile);
	free(word);
	free(line);
    }

    //Pass all info to all processors
    MPI_Bcast(xread, (gxelem+1)*(gyelem+1), MPI_DOUBLE, master, MPI_COMM_WORLD);
    MPI_Bcast(yread, (gxelem+1)*(gyelem+1), MPI_DOUBLE, master, MPI_COMM_WORLD);
    
    //Rearrange in 2D array
    double **xtemp, **ytemp;
    allocator2(&xtemp, gxelem+1, gyelem+1);
    allocator2(&ytemp, gxelem+1, gyelem+1);

    i=0;
    j=0;
    k=0;
    
    for(j=0; j<gyelem+1; j++)
    {
	for(i=0; i<gxelem+1; i++)
	{
	    xtemp[i][j] = xread[k];
	    ytemp[i][j] = yread[k];
	    k++;
	}
    }

    //Now split the mesh
    int firstx = (myrank%procm)*elemm;
    int firsty = (double)floor(myrank/procm)*elemn;

    for(j=2; j<ynode-2; j++)
    {
        for(i=2; i<xnode-2; i++)
        {
            x[i][j] = xtemp[firstx+i-2][firsty+j-2];
            y[i][j] = ytemp[firstx+i-2][firsty+j-2];
        }
    }

    //Now do the rest of mesh operations
    meshOperations(x, y, area, vol, xc, yc);

    /*FILE *meshcheck;
    meshcheck = fopen("../../meshcheck.dat","w");

    for(j=0; j<ynode; j++)
    {
	for(i=0; i<xnode; i++)
	{
	    fprintf(meshcheck,"%12.8f %12.8f\n",x[i][j],y[i][j]);
	}
    }
    
    fclose(meshcheck);*/
    //------------------------------------------------------------------------//
    //Deallocators
    deallocator1(&xread, (gxelem+1)*(gyelem+1));
    deallocator1(&yread, (gxelem+1)*(gyelem+1));
    deallocator2(&xtemp, gxelem+1, gyelem+1);
    deallocator2(&ytemp, gxelem+1, gyelem+1);
    //------------------------------------------------------------------------//

}

void meshOperations(double **x, double **y, double ****area, double **vol, double **xc, double **yc)
{
    int i,j;
    
    //Generate ghost nodes
    //Add cells on both sides of x
    for(j=2; j < ynode-2; j++)
    {
	for(i=1; i>=0;i--)
	{
	  x[i][j]=2.0*x[i+1][j]-x[i+2][j];
	  y[i][j]=2.0*y[i+1][j]-y[i+2][j];
	  x[xnode-1-i][j] = 2.0*x[xnode-2-i][j]-x[xnode-3-i][j];
	  y[xnode-1-i][j] = 2.0*y[xnode-2-i][j]-y[xnode-3-i][j];
	  //printf("%d %d %.2f %.2f\n",i,j,x[i][j],y[i][j]);
	}
    }
    
    for(i=2; i < xnode-2; i++)
    {
      for(j=1; j>=0; j--)
	{
	  x[i][j]=2.0*x[i][1+j]-x[i][2+j];
	  y[i][j]=2.0*y[i][1+j]-y[i][2+j];
	  x[i][ynode-1-j] = 2*x[i][ynode-2-j]-x[i][ynode-3-j];
	  y[i][ynode-1-j] = 2*y[i][ynode-2-j]-y[i][ynode-3-j];
	}
    }
    //Not required but assign proper coordinates to corner ghost nodes
    for (i=1; i>=0; i--)
    {
	for(j=1; j>=0; j--)
	{
	    x[i][j]=2.0*x[i][1+j]-x[i][2+j];
	    y[i][j]=2.0*y[i][1+j]-y[i][2+j];
	    x[i][ynode-1-j] = 2*x[i][ynode-2-j]-x[i][ynode-3-j];
	    y[i][ynode-1-j] = 2*y[i][ynode-2-j]-y[i][ynode-3-j];
	}
    }

    for(i=xnode-2; i<xnode; i++)
    {
	for(j=1; j>=0; j--)
	{
	    x[i][j]=2.0*x[i][1+j]-x[i][2+j];
	    y[i][j]=2.0*y[i][1+j]-y[i][2+j];
	    x[i][ynode-1-j] = 2*x[i][ynode-2-j]-x[i][ynode-3-j];
	    y[i][ynode-1-j] = 2*y[i][ynode-2-j]-y[i][ynode-3-j];
	}
    }
    
    x[0][0] = 2.0*x[0][1]-x[0][2];
    y[0][0] = 2.0*y[0][1]-y[0][2];
    x[0][ynode-1] = 2.0*x[0][ynode-2]-x[0][ynode-3];
    y[0][ynode-1] = 2.0*y[0][ynode-2]-y[0][ynode-3];
    x[xnode-1][0] = 2.0*x[xnode-2][0]-x[xnode-3][0];
    y[xnode-1][0] = 2.0*y[xnode-2][0]-y[xnode-3][0];
    x[xnode-1][ynode-1] =  2.0*x[xnode-2][ynode-1]-x[xnode-3][ynode-1];
    y[xnode-1][ynode-1] =  2.0*y[xnode-2][ynode-1]-y[xnode-3][ynode-1];
    
    /*for(i=0;i<xnode;i++)
    {
        for(j=0;j<ynode;j++)
        {
            printf("%d %d %.2f %.2f\n",i,j,x[i][j],y[i][j]);
        }
    }
    exit(1);*/

    //------------------------------------------------------------------------//
    //I don't think i require the following anymore. If needed just uncomment
    //Area components of faces of parallelogram
    for(i=0; i<xelem; i++)
    {
        for(j=0; j<yelem; j++)
        {
            area[i][j][0][0] = y[i+1][j+1] - y[i+1][j];
            area[i][j][0][1] = -(x[i+1][j+1] - x[i+1][j]);
            area[i][j][1][0] = y[i][j+1] - y[i+1][j+1];
            area[i][j][1][1] = -(x[i][j+1]-x[i+1][j+1]);
	    //printf("%.4e %.4e %.4e %.4e\n",area[i][j][0][0],area[i][j][0][1],area[i][j][1][0],area[i][j][1][1]);
	    
        }
    }
    //exit(1);
    
    for(i=0; i<xelem; i++)
    {
        for(j=0; j<yelem; j++)
        {
            //Heron's Formula for area of triangle (area of cell)
            double a=sqrt(pow(x[i][j]-x[i+1][j],2)+pow(y[i][j]-y[i+1][j],2));
            double b=sqrt(pow(x[i+1][j+1]-x[i+1][j],2)+pow(y[i+1][j+1]-y[i+1][j],2));
            double c=sqrt(pow(x[i][j]-x[i+1][j+1],2)+pow(y[i][j]-y[i+1][j+1],2));
            double s=(a+b+c)/2.0;
            vol[i][j]=sqrt(s*(s-a)*(s-b)*(s-c));

            double a1=sqrt(pow(x[i][j]-x[i][j+1],2)+pow(y[i][j]-y[i][j+1],2));
            double b1=sqrt(pow(x[i][j+1]-x[i+1][j+1],2)+pow(y[i][j+1]-y[i+1][j+1],2));
            double s1=(a1+b1+c)/2.0;
            vol[i][j]=vol[i][j]+sqrt(s1*(s1-a1)*(s1-b1)*(s1-c));

            //////Centroid of cell///////
            xc[i][j]=0.25*(x[i][j]+x[i+1][j]+x[i+1][j+1]+x[i][j+1]);
            yc[i][j]=0.25*(y[i][j]+y[i+1][j]+y[i+1][j+1]+y[i][j+1]);
        }
    }
}



void meshShift(double **x, double **y, double **xc, double **yc)
{
    double xshift = xlen/2.0;
    double yshift = ylen/2.0;

    int i,j;
    
    for(i=0; i<xnode; i++)
    {
	for(j=0; j<ynode; j++)
	{
	    x[i][j] -= xshift;
	    y[i][j] -= yshift;
	}
    }

    for(i=0; i<xelem; i++)
    {
	for(j=0; j<yelem; j++)
	{
	    xc[i][j] -= xshift;
	    yc[i][j] -= yshift;
	}
    }
}
