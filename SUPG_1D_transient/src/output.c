/***************************************************************************

Author: nsaini
Created: 2018-07-24

***************************************************************************/


#include "common.h"
#include "polylib.h"
#include "functions.h"
#include "memory.h"

void output(double *phi, double *x, int iter)
{
    char* dirpath;
    dirpath = concat(getexepath(), "/output/");
    //printf("%s\n",dirpath);
    DIR* dir = opendir(dirpath);
    if(dir)
    {
	closedir(dir);
    }
    else if(ENOENT == errno)
    {
	mkdir(dirpath, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    }
    free(dirpath);

    //------------------------------------------------------------------------//
    //Temporary Variables
    int i,j;
    //------------------------------------------------------------------------//


    
    FILE *out;
    char* path1;
    path1 = concat(getexepath(),"/output/out_00");
    char buffer0[12];
    snprintf(buffer0,12,"%d",iter);
    path1 = concat(path1,buffer0);
    path1 = concat(path1,".vts");
    out = fopen(path1,"w");
    free(path1);

    
    //Write Headers
    fprintf(out,"# vtk DataFile Version 3.0\n");
    fprintf(out,"vtk output\n");
    fprintf(out,"ASCII\n");
    fprintf(out,"DATASET STRUCTURED_GRID\n");
    fprintf(out,"DIMENSIONS %d 2 1\n",xnode);
    fprintf(out,"POINTS %d double\n",(xnode)*2);

    
    double dummy ;
    for(j=0; j<2; j++)
    {
	if(j==0)
	{
	    dummy = 0.0;
	}
	else
	{
	    dummy = xlen/2.0;
	}
	for(i=0; i<xnode; i++)
	{
	    fprintf(out,"%.6f %.6f 0.0\n",x[i],dummy);
	}
    }
    
    fprintf(out,"POINT_DATA %d\n",(xnode)*2);
    fprintf(out,"SCALARS phi double\n");
    fprintf(out,"LOOKUP_TABLE default\n");

    for(j=0; j<2; j++)
    {
	for(i=0; i<xnode; i++)
	{
	    double phinode = phi[i];
	    fprintf(out,"%.6f\n",phinode);
	}
    }
    fprintf(out,"CELL_DATA %d\n",(xnode));
    fclose(out);
    

}
