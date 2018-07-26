/***************************************************************************

Author: nsaini
Created: 2018-07-24

***************************************************************************/


#include "common.h"
#include "polylib.h"
#include "functions.h"
#include "memory.h"

void output(struct elemsclr elem, double *x, int iter)
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
    int ielem, icoeff;
    int i,j;
    double *basis;
    allocator1(&basis,ncoeff);
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Reconstruct the solution at the nodes
    double **recphi, **recu;
    allocator2(&recphi, xelem, 2);
    allocator2(&recu, xelem, 2);

    for(ielem=0; ielem < xelem; ielem++)
    {
	for(i=0; i<2; i++)
	{
	    if(i==0)
	    {
		basis1D(-1.0,basis);
	    }
	    else
	    {
		basis1D(1.0,basis);
	    }
	    for(icoeff=0; icoeff<ncoeff; icoeff++)
	    {
		recphi[ielem][i] += basis[icoeff]*elem.phi[ielem][icoeff];
		recu[ielem][i] += basis[icoeff]*elem.u[ielem][icoeff];
	    }
	}
    }
   
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
    fprintf(out,"DIMENSIONS %d 2 1\n",xelem-1);
    fprintf(out,"POINTS %d double\n",(xelem-1)*2);

    
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
	for(i=1; i<xelem; i++)
	{
	    fprintf(out,"%.6f %.6f 0.0\n",x[i],dummy);
	}
    }
    
    fprintf(out,"POINT_DATA %d\n",(xelem-1)*2);
    fprintf(out,"SCALARS u double\n");
    fprintf(out,"LOOKUP_TABLE default\n");

    for(j=0; j<2; j++)
    {
	for(i=1; i<xelem-1; i++)
	{
	    double unode = 0.5*(recu[i-1][1] + recu[i][0]);
	    fprintf(out,"%.6f\n",unode);
	}
    }

    fprintf(out,"SCALARS phi double\n");
    fprintf(out,"LOOKUP_TABLE default\n");

    for(j=0; j<2; j++)
    {
	for(i=1; i<xelem-1; i++)
	{
	    double phinode = 0.5*(recphi[i-1][1] + recphi[i][0]);
	    fprintf(out,"%.6f\n",phinode);
	}
    }
    fprintf(out,"CELL_DATA %d\n",(xelem-2));
    fclose(out);
    
    //------------------------------------------------------------------------//
    //Deallocators
    deallocator2(&recphi, xelem, 2);
    deallocator2(&recu, xelem, 2);
    deallocator1(&basis,ncoeff);
    //------------------------------------------------------------------------//

}
