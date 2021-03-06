/***************************************************************************

Author: nsaini
Created: 2018-03-10

***************************************************************************/
#include "common.h"
#include "fileIO.h"
#include "generalFunc.h"
#include "DGPFunc.h"
#include "memory.h"

void output_xml(struct elemsclr elem, int iter , double **x, double **y)
{
    int i,j,k,l;
    
    char* dirpath;
    dirpath = concat(getexepath(), "/output/");
    char buf_dir[12];
    snprintf(buf_dir,12,"%d",iter);
    dirpath = concat(dirpath,buf_dir);
    DIR* dir = opendir("dirpath");
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
    //Reconstruct solution for all elements at the nodes
    double ***recphi, ***recu, ***recv;
    allocator3(&recphi, xelem, yelem, 4);
    allocator3(&recu, xelem, yelem, 4);
    allocator3(&recv, xelem, yelem, 4);

    double **zs;
    allocator2(&zs, 4, 2);

    double *basis;
    allocator1(&basis, ncoeff);

    zs[0][0] = -1.0;
    zs[0][1] = -1.0;

    zs[1][0] = 1.0;
    zs[1][1] = -1.0;

    zs[2][0] = -1.0;
    zs[2][1] = 1.0;

    zs[3][0] = 1.0;
    zs[3][1] = 1.0;

    
    for(i=0; i<xelem; i++)
    {
	for(j=0; j<yelem; j++)
	{
	    for(k=0; k<4; k++)
	    {
		recphi[i][j][k] = 0.0;
		recu[i][j][k] = 0.0;
		recv[i][j][k] = 0.0;

		basis2D(zs[k][0], zs[k][1], basis);
		
		for(l=0; l<ncoeff; l++)
		{
		    recphi[i][j][k] += basis[l]*elem.phi[i][j][l];
		    recu[i][j][k] += basis[l]*elem.u[i][j][l];
		    recv[i][j][k] += basis[l]*elem.v[i][j][l];
		}
	    }
	}
    }
    //------------------------------------------------------------------------//


    FILE *out;
    FILE *out1;
    //if(myrank == master)
    //{
    if(myrank == master)
    {
	char* path1;
	path1 = concat(getexepath(),"/output/out_00");
	char buffer0[12];
	snprintf(buffer0,12,"%d",iter);
	path1 = concat(path1,buffer0);
	path1 = concat(path1,".pvts");
	out1 = fopen(path1,"w");
	free(path1);
	
	fprintf(out1,"<VTKFile type=\"PStructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
	fprintf(out1,"<PStructuredGrid WholeExtent=\"%d %d %d %d %d %d\" GhostLevel=\"1\">\n",1,io_info[nprocs-1][1]+1,1,io_info[nprocs-1][3]+1,0,0);
	fprintf(out1,"<PPoints>\n");
	fprintf(out1,"<PDataArray NumberOfComponents=\"3\" format=\"ascii\" type =\"Float32\" Name=\"mesh\"/>\n");
	fprintf(out1,"</PPoints>\n");
	fprintf(out1,"<PPointData>\n");
	fprintf(out1,"<PDataArray NumberOfComponents=\"1\" format=\"ascii\" type =\"Float32\" Name=\"phi\"/>\n");
	fprintf(out1,"<PDataArray NumberOfComponents=\"1\" format=\"ascii\" type =\"Float32\" Name=\"u\"/>\n");
	fprintf(out1,"<PDataArray NumberOfComponents=\"1\" format=\"ascii\" type =\"Float32\" Name=\"v\"/>\n");
	/*fprintf(out1,"<PDataArray NumberOfComponents=\"1\" format=\"ascii\" type =\"Float32\" Name=\"mu\"/>\n");
	fprintf(out1,"<PDataArray NumberOfComponents=\"1\" format=\"ascii\" type =\"Float32\" Name=\"rho\"/>\n");
	fprintf(out1,"<PDataArray NumberOfComponents=\"1\" format=\"ascii\" type =\"Float32\" Name=\"p\"/>\n");
	fprintf(out1,"<PDataArray NumberOfComponents=\"1\" format=\"ascii\" type =\"Float32\" Name=\"phi2\"/>\n");*/
	fprintf(out1,"</PPointData>\n");
	for(i=0;i<nprocs; i++)
	{
	    fprintf(out1,"<Piece Extent=\"%d %d %d %d %d %d\" Source=\"%d/out.%d.vts\"/>\n",io_info[i][0],io_info[i][1]+1,io_info[i][2],io_info[i][3]+1,0,0,iter,i);
	}
	fprintf(out1,"</PStructuredGrid>\n");
	fprintf(out1,"</VTKFile>\n");
	fclose(out1);
    }
    
    char* path;
    char buffer1[12];
    char buffer2[12];
    path = concat(getexepath(), "/output/");
    snprintf(buffer1,12,"%d",iter);
    path = concat(path,buffer1);
    path = concat(path,"/out.");
    snprintf(buffer2,12,"%d",myrank);
    path = concat(path,buffer2);
    path = concat(path,".vts");
    
  
    out = fopen(path,"w");
    free(path);
    if(out == NULL)
    {
	printf("Error opening output.dat!\n");
	exit(0);
    }
    //Write Headers
    fprintf(out,"<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
    fprintf(out,"<StructuredGrid WholeExtent=\"%d %d %d %d %d %d\">\n",1,io_info[nprocs-1][1]+1,1,io_info[nprocs-1][3]+1,0,0);
    fprintf(out,"<Piece Extent=\"%d %d %d %d %d %d\">\n",io_info[myrank][0],io_info[myrank][1]+1,io_info[myrank][2],io_info[myrank][3]+1,0,0);
    fprintf(out,"<PointData></PointData>\n");
    fprintf(out,"<CellData></CellData>\n");
    fprintf(out,"<Points>\n");
    fprintf(out,"<DataArray NumberOfComponents=\"3\" format=\"ascii\" type =\"Float32\" Name=\"mesh\">\n");

    double zlen = x[2][2] - x[1][2];
    for(j=2; j<ynode-2; j++)
    {
	for(i=2; i<xnode-2; i++)
	{
	    fprintf(out,"%.6f %.6f %.6f\n",x[i][j],y[i][j],zlen);
	}
    }
    fprintf(out,"</DataArray>");
    fprintf(out,"</Points>\n");
    fprintf(out,"<PointData>\n");
    
    fprintf(out,"<DataArray NumberOfComponents=\"1\" format=\"ascii\" type =\"Float32\" Name=\"phi\">\n");
    for(j=2; j<ynode-2; j++)
    {
	for(i=2; i<xnode-2; i++)
	{
	    	    
	    double phinode=0.25*(recphi[i][j][0] + recphi[i-1][j][1] + recphi[i-1][j-1][3] + recphi[i][j-1][2]);

	    if((iter == 0 && case_tog == 3) || (iter == 0 && case_tog == 8))
	    {
		if(phinode > 1.0)phinode = 1.0;

		if(phinode < 0.0)phinode = 0.0;
	    }
	    
	    fprintf(out,"%.6f\n",phinode);
	}
    }

    
    fprintf(out,"</DataArray>\n");
    fprintf(out,"<DataArray NumberOfComponents=\"1\" format=\"ascii\" type =\"Float32\" Name=\"u\">\n");
    for(j=2; j<ynode-2; j++)
    {
	for(i=2; i<xnode-2; i++)
	{
	    
	    double unode=0.25*(recu[i][j][0] + recu[i-1][j][1] + recu[i-1][j-1][3] + recu[i][j-1][2]);
	    fprintf(out,"%.6f\n",unode);
	}
    }
    
    fprintf(out,"</DataArray>\n");
    fprintf(out,"<DataArray NumberOfComponents=\"1\" format=\"ascii\" type =\"Float32\" Name=\"v\">\n");
    for(j=2; j<ynode-2; j++)
    {
	for(i=2; i<xnode-2; i++)
	{
	    
	    double vnode=0.25*(recv[i][j][0] + recv[i-1][j][1] + recv[i-1][j-1][3] + recv[i][j-1][2]);
	    fprintf(out,"%.6f\n",vnode);
	}
    }
    fprintf(out,"</DataArray>\n");
    
    /*fprintf(out,"<DataArray NumberOfComponents=\"1\" format=\"ascii\" type =\"Float32\" Name=\"mu\">\n");
    for(j=2; j<ynode-2; j++)
    {
	for(i=2; i<xnode-2; i++)
	{
	    double munode=0.25*(elem.mu[i][j]+elem.mu[i-1][j]+elem.mu[i-1][j-1]+elem.mu[i][j-1]);
	    fprintf(out,"%.6f\n",munode);
	    
	}
    }
    fprintf(out,"</DataArray>\n");
    
    fprintf(out,"<DataArray NumberOfComponents=\"1\" format=\"ascii\" type =\"Float32\" Name=\"rho\">\n");
    for(j=2; j<ynode-2; j++)
    {
	for(i=2; i<xnode-2; i++)
	{
	    
	    double rhonode=0.25*(elem.rho[i][j]+elem.rho[i-1][j]+elem.rho[i-1][j-1]+elem.rho[i][j-1]);
	    fprintf(out,"%.6f\n",rhonode);
	}
    }
    fprintf(out,"</DataArray>\n");
    
    fprintf(out,"<DataArray NumberOfComponents=\"1\" format=\"ascii\" type =\"Float32\" Name=\"p\">\n");
    for(j=2; j<ynode-2; j++)
    {
	for(i=2; i<xnode-2; i++)
	{
	    
	    double pnode=0.25*(elem.p[i][j]+elem.p[i-1][j]+elem.p[i-1][j-1]+elem.p[i][j-1]);
	    fprintf(out,"%.6f\n",pnode);
	}
    }
    fprintf(out,"</DataArray>\n");

    fprintf(out,"<DataArray NumberOfComponents=\"1\" format=\"ascii\" type =\"Float32\" Name=\"phi2\">\n");
    for(j=2; j<ynode-2; j++)
    {
	for(i=2; i<xnode-2; i++)
	{
	    
	    double phi2node=0.25*(elem.phi2[i][j]+elem.phi2[i-1][j]+elem.phi2[i-1][j-1]+elem.phi2[i][j-1]);
	    fprintf(out,"%.6f\n",phi2node);
	}
	}
	fprintf(out,"</DataArray>\n");*/
    
    fprintf(out,"</PointData>\n");
    fprintf(out,"</Piece>\n");
    fprintf(out,"</StructuredGrid>\n");
    fprintf(out,"</VTKFile>\n");
    fclose(out);
    //}

    deallocator3(&recphi, xelem, yelem, 4);
    deallocator3(&recu, xelem, yelem, 4);
    deallocator3(&recv, xelem, yelem, 4);
    deallocator2(&zs, 4, 2);
    deallocator1(&basis, ncoeff);
}

