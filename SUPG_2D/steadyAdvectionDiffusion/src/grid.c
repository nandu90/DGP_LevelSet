/***************************************************************************

Author: nsaini
Created: 2018-03-04

***************************************************************************/

#include "common.h"
#include "mesh.h"
#include "memory.h"

void assignGlobalDof(struct dofdata *dof)
{
    if(myrank == 0)
    {
	printf("In assignGlobalDof\n");
    }
    //------------------------------------------------------------------------//
    //Temporary Varibles
    int idof;
    int iproc;
    int *procdof;
    iallocator1(&procdof, nprocs);
    MPI_Request request;
    MPI_Status status;
    int index;
    int increment;
    int rows = ynode-4; //Number of rows of dofs
    int cols = xnode-4;
    int split;
    int ne1, ne2;
    int exclude1, exclude2;
    int recvproc;
    int tag;
    int *recvbuf;
    int *sendbuf;
    int count;
    int lim;
    int inc;
    int start;
    int sendsize;
    int recvsize;
    int i;
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Initialize dofs to -1
    for(idof=0; idof<tdof; idof++)
    {
	dof[idof].gindex = -1;
    }
    //------------------------------------------------------------------------//

    
    //------------------------------------------------------------------------//
    //Assign numbers to interior dofs
    index=0;
    for(idof=0; idof<tdof; idof++)
    {
	if(dof[idof].BC == 0)
	{
	    dof[idof].gindex = index++;
	    dof[idof].controlproc = myrank;
	}
    }
    //Now assign numbers to exclusive boundary dofs
    for(idof=0; idof<tdof; idof++)
    {
	if(dof[idof].BC == 1)
	{
	    dof[idof].gindex = index++;
	    dof[idof].controlproc = myrank;
	}
    }
    
    //------------------------------------------------------------------------//

    if(myrank == 0)
    {
	printf("Total procs=%d\n",nprocs);
    }
    //------------------------------------------------------------------------//
    //Pass around the info for starting index
    MPI_Allgather(&index, 1, MPI_INT, procdof, 1, MPI_INT, MPI_COMM_WORLD);
    increment = 0;    
    for(iproc=0; iproc<myrank; iproc++)
    {
	increment += procdof[iproc];
    }
    
    //Increment starting index on all but master proc
    for(idof=0; idof<tdof; idof++)
    {
	if(dof[idof].BC ==0 || dof[idof].BC == 1)
	{
	    dof[idof].gindex += increment;
	    index = dof[idof].gindex + 1;
	}
    }

    //Bcast the last index from last processor to everyone
    MPI_Bcast(&index, 1, MPI_INT, nprocs-1, MPI_COMM_WORLD);
    
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Assign shared dof
    /*Method:
      1. Loop over the faces starting from the right
      2. The lower ranked processor does all the work.
      :- Split the number of shared nodes evenly.
      :- if there is a node on this face shared by more than two proc, exclude it
      :- Pass info on who is the master processor for the nodes to the higher ranked proc
      :- Assign node numbers starting from 1
      :- Later, when master proc is decided for all shared nodes, raise the assigned node value       */
    
    //------------------------------------------------------------------------//
    //Assign control processor to the shared nodes
    for(i=0; i<4; i++)
    {
	
	if(i==0)recvproc = 2;
	if(i==1)recvproc = 3;
	if(i==2)recvproc = 0;
	if(i==3)recvproc = 1;

	if(bhailog[recvproc] >= 0 && myrank > bhailog[recvproc])
	{
	    tag = bhailog[recvproc];
	    if(recvproc == 0 || recvproc == 2)
	    {
		recvsize = rows;
	    }
	    else
	    {
		recvsize = cols;
	    }
	    
	    iallocator1(&recvbuf, recvsize);
	    
	    MPI_Recv(recvbuf, recvsize, MPI_INT, bhailog[recvproc],tag,MPI_COMM_WORLD, &status);


	    //Now just assign the proc info to the dofs    
	    if(recvproc == 0) //Right face
	    {
		start = cols-1;
		lim = tdof-1;
		inc = cols;
	    }
	    else if(recvproc==1) //Upper face
	    {
		start = cols*(rows-1);
		lim = tdof-1;
		inc = 1;
	    }
	    else if(recvproc==2) //left face
	    {
		start = 0;
		lim = cols*(rows-1);
		inc = cols;
	    }
	    else
	    {
		start = 0;
		lim = cols-1;
		inc = 1;
	    }

	    count=0;
	    for(idof=start; idof<=lim; idof += inc)
	    {
		dof[idof].controlproc = recvbuf[count++];
	    }
	    ideallocator1(&recvbuf, recvsize);
	}

	ne2 = i+1;
	ne1 = i-1;
	if(ne2 > 3)
	{
	    ne2 -= 4;
	}
	if(ne1 < 0)
	{
	    ne1 += 4;
	}
	
	if(bhailog[i] >= 0 && myrank < bhailog[i])
	{
	    exclude1 = 0;
	    exclude2 = 0;
	    if(i == 0 || i==2)
	    {
		split = rows/2;
	    }
	    else
	    {
		split =  cols/2;
	    }
	    if(bhailog[ne1] != -1)
	    {
		if(i == 0 || i==3)
		{
		    exclude1 = 1;
		}
		else
		{
		    exclude2 = 1;
		}
	    }
	    if(bhailog[ne2] != -1)
	    {
		if(i==0 || i==3)
		{
		    exclude2 = 1;
		}
		else
		{
		    exclude1 = 1;
		}
	    }   

	   
	    //Now assign ID's to the nodes you have decided to keep on this proc	    
	    if(i == 0) //Right face
	    {
		start = cols-1;
		lim = tdof-1;
		inc = cols;
	    }
	    else if(i==1) //Upper face
	    {
		start = cols*(rows-1);
		lim = tdof-1;
		inc = 1;
	    }
	    else if(i==2) //left face
	    {
		start = 0;
		lim = cols*(rows-1);
		inc = cols;
	    }
	    else
	    {
		start = 0;
		lim = cols-1;
		inc = 1;
	    }
	    
	    count = 1;
	    for(idof=start; idof<=lim; idof += inc)
	    {
		if(idof == start && exclude1 != 1)
		{
		    dof[idof].controlproc = bhailog[i];
		}
		else if(idof == lim && exclude2 !=1)
		{
		    dof[idof].controlproc = myrank;
		}
		else if(idof != start && idof != lim)
		{
		    if(count <= split)
		    {
			dof[idof].controlproc = bhailog[i];
		    }
		    else
		    {
			dof[idof].controlproc = myrank;
		    }
		}
		else
		{
		    dof[idof].controlproc = -1;
		}
		count++;
	    }

	    if(i == 0 || i == 2)
	    {
		sendsize = rows;
	    }
	    else
	    {
		sendsize = cols;
	    }

	    iallocator1(&sendbuf, sendsize);
	    count=0;
	    for(idof=start; idof<=lim; idof += inc)
	    {
		sendbuf[count++] = dof[idof].controlproc;
	    }
	    
	    tag = myrank;
	    MPI_Isend(sendbuf, sendsize, MPI_INT, bhailog[i], tag, MPI_COMM_WORLD, &request);    
	    MPI_Wait(&request, &status);
	    ideallocator1(&sendbuf, sendsize);
	}
	
    }
    //------------------------------------------------------------------------//

        
    MPI_Barrier(MPI_COMM_WORLD);
    
    //Update the global Id's on each proc, so that they are unique
    count = 0;
    for(idof =0; idof < tdof; idof++)
    {
	if(dof[idof].controlproc == myrank && (dof[idof].BC == 2 || dof[idof].BC == 3))
	{
	    count++;
	}
    }
    
    MPI_Allgather(&count, 1, MPI_INT, procdof, 1, MPI_INT, MPI_COMM_WORLD);

    increment = 0;
    for(iproc=0; iproc<myrank; iproc++)
    {
	increment += procdof[iproc];
    }

    index = index + increment;
    for(idof=0; idof < tdof; idof++)
    {
	if(dof[idof].controlproc == myrank && (dof[idof].BC == 2 || dof[idof].BC == 3))
	{
	    dof[idof].gindex = index++;
	}
    }
    
        
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Broadcast the next global index
    MPI_Bcast(&index, 1, MPI_INT, nprocs-1, MPI_COMM_WORLD);
    //Finally assign global IDs to nodes shared by more than 2 proc
    //Assign the master processor
    //Convention: diagonally southwest partition is the master
    //Loop over the corner nodes
    for(i=0; i<4; i++)
    {
	ne2 = i;
	ne1 = i-1; //Processors that share this node
	if(i == 0)
	{
	    ne1 = 3;
	}

	if(bhailog[ne1] != -1 && bhailog[ne2] != -1) //Then shared node
	{
	    if(i == 0)
	    {
		idof = cols-1;
		dof[idof].controlproc = bhailog[ne1];
	    }
	    else if(i == 1)
	    {
		idof = tdof-1;
		dof[idof].controlproc = myrank;
	    }
	    else if(i == 2)
	    {
		idof = cols*(rows - 1);
		dof[idof].controlproc = bhailog[ne2];
	    }
	    else
	    {
		idof = 0;		
		if(procm-1 < 0)
		{
		    dof[idof].controlproc = (procm*procn)-1;
		}
		else
		{
		    dof[idof].controlproc = procm - 1;
		}
	    }
	}
    }
    //Now assign global id to all shared nodes on each processor
    for(iproc=0; iproc<nprocs; iproc++)
    {
	if(iproc == myrank)
	{
	    for(idof=0; idof<tdof; idof++)
	    {
		if(dof[idof].BC == 4 && dof[idof].controlproc == myrank)
		{
		    dof[idof].gindex = index++;
		}
	    }
	}
	MPI_Bcast(&index, 1, MPI_INT, iproc, MPI_COMM_WORLD);
    }
    //------------------------------------------------------------------------//

    if(myrank == 7)
    {
	for(idof=0; idof<tdof; idof++)
	{
	    printf("%d %d %d %d\n",idof, dof[idof].BC, dof[idof].controlproc, dof[idof].gindex);
	}
    }

    MPI_Barrier(MPI_COMM_WORLD);
    exit(1);
    ideallocator1(&procdof, nprocs);
    
}

void markSharedDof(struct dofdata *dof)
{
    if(myrank == 0)
    {
	printf("In markSharedDof\n");
    }
    
    //------------------------------------------------------------------------//
    //Temporary variables
    int idof;
    int i;
    //------------------------------------------------------------------------//

    int rows = ynode-4; //Number of rows of dofs
    int cols = xnode-4;
    
    for(i=0; i<4; i++)
    {
	if(bhailog[i] != -1)
	{
	    if(i == 0) //Right face
	    {
		for(idof=cols-1; idof<tdof; idof+=cols)
		{
		    dof[idof].BC += 2;
		}
	    }
	    else if(i==1) //Upper face
	    {
		for(idof=cols*(rows-1); idof<tdof; idof++)
		{
		    dof[idof].BC += 2;
		}
	    }
	    else if(i==2) //left face
	    {
		for(idof=0; idof<tdof; idof+=cols)
		{
		    dof[idof].BC += 2;
		}
	    }
	    else
	    {
		for(idof=0; idof<cols; idof++) //bottom face
		{
		    dof[idof].BC += 2;
		}
	    }
	}
    }

    if(myrank == 0)
    {
	for(idof=0; idof<tdof; idof++)
	{
	    printf("%d %d\n",idof,dof[idof].BC);
	}
    }
    
    
}

void markBoundaryDof(struct dofdata *dof)
{
    if(myrank == 0)
    {
	printf("In markBoundaryDof\n");
    }
    /*Convention: 
      if BC = 0 ==> interior dof
      if BC = 1 ==> shared dof
      if BC = 2 ==> boundary dof
    */
    
    //------------------------------------------------------------------------//
    //Temporary variables
    int idof;
    int i;
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Initialize BC to assume this is an interior dof
    for(idof = 0; idof<tdof; idof++)
    {
	dof[idof].BC = 0;
    }
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Initialize periodic array to zero
    for(i=0; i<4; i++)
    {
	per[i] = 0;
    }
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    ///Assume all processors have 4 neighbours to start with
    //Neighbour processors will have rank as defined below
    int l_bhai = myrank - 1;
    int r_bhai = myrank + 1;
    int u_bhai = myrank + procm;
    int d_bhai = myrank - procm;

    //int rows = ynode-4; //Number of rows of dofs
    int cols = xnode-4;
    
    //Now take care of dofs of boundary processors
    if((int)floor(myrank/procm) == 0) //Means that this processor is aligned with lower boundary
    {
	//down = 1;
	d_bhai = -1000;           //Set rank of down proc to -ve value
	if(y_bound != 3)           //Dont declare it a boundary node for periodic BC
	{
	    for(idof=0; idof<cols; idof++)
	    {
		dof[idof].BC = 1;
	    }
	}
	else //If BC is periodic
	{
	    d_bhai = myrank + (procn-1)*procm;
	    per[3] = 1;
	}
    }

    if((int)floor(myrank/procm) == procn-1) //Means processor is on upper boundary
    {
	//up = 1; 
      u_bhai = -1000;
      if(y_bound != 3)                //Dont declare it a boundary node for periodic BC
	{
	  for(idof=tdof-cols; idof<tdof; idof++)
	    {
		dof[idof].BC = 1;
	    }
	}
      else
	{
	  u_bhai = myrank % procm;
	  per[1] = 1;
	}
    }

    if((myrank % procm) == 0) //Means processor is on left boundary
    {
	//left = 1;
	l_bhai = -1000;
	if(x_bound != 3)           //Dont declare it a boundary node for periodic BC
	{
	    for(idof=0; idof<tdof-cols+1; idof+=cols)
	    {
		dof[idof].BC = 1;
	    }
	}
	else
	{
	    l_bhai = myrank + procm - 1;
	    per[2] = 1;
	}
    }

    if((myrank % procm) == procm-1) //Means processor is on right boundary
    {
	//right = 1;
      r_bhai = -1000;
      if(x_bound != 3)           //Dont declare it a boundary node for periodic BC
	{
	  for(idof=cols-1; idof<tdof; idof+=cols)
	    {
		dof[idof].BC = 1;
	    }
	}
      else
	{
	  r_bhai = myrank - procm + 1;
	  per[0] = 1;
	}
    }
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Fill out the bhailog array. 
    //Convention: Starting from right  go anticlockwise
    bhailog[0] = r_bhai;
    bhailog[1] = u_bhai;
    bhailog[2] = l_bhai;
    bhailog[3] = d_bhai;
    
    for(i=0; i<4; i++)
    {
	if(bhailog[i] < 0 || bhailog[i] == myrank)
	{
	    bhailog[i] = -1;  //If no neighbour mark that placeholder as -1
	}
    }
    if(myrank == 0)
    {
	printf("%d %d %d %d\n",bhailog[0], bhailog[1], bhailog[2], bhailog[3]);
    }
    //------------------------------------------------------------------------//

    if(myrank == 0)
    {
	for(idof=0; idof<tdof; idof++)
	{
	    printf("%d %d\n",idof,dof[idof].BC);
	}
    }
}

void gridgen(double **x, double **y, struct elemdata **edata, struct dofdata *dof)
{
    //------------------------------------------------------------------------//
    //Temporary Variables
    int i,j;
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    double deltax=xlen/(gxelem);
    double deltay=ylen/(gyelem);
    double firstx = (myrank%procm)*elemm*deltax;
    double firsty = (double)floor(myrank/procm)*elemn*deltay;
    //Assign coordinates to mesh nodes
    for(j=2; j<ynode-2; j++)
    {
        for(i=2; i<xnode-2; i++)
        {
            x[i][j] = firstx + (i-2)*deltax;
            y[i][j] = firsty + (j-2)*deltay;
	    //printf("%d %d %.2f %.2f\n",i,j,x[i][j],y[i][j]); 
        }
    }
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Initialize element array
    elemallocator(&edata, xelem, yelem);
    //Initialize DOF array
    tdof = (xnode-4)*(ynode-4);  //Number of dofs on this processor
    dofallocator(&dof, tdof);
    //------------------------------------------------------------------------//

    
    
    //------------------------------------------------------------------------//
    //Assign local dofs to elements
    //int row = (ynode-4);
    int col = xnode - 4;
    
    int ielem, jelem;
    for(ielem=2; ielem<xelem-2; ielem++)
    {
	for(jelem=2; jelem<yelem-2; jelem++)
	{
	    edata[ielem][jelem].edofn = 4;
	    iallocator1(&(edata[ielem][jelem].edofs),edata[ielem][jelem].edofn);
	    i = (ielem-2) + (jelem-2)*col;
	    
	    edata[ielem][jelem].edofs[0] = i;
	    edata[ielem][jelem].edofs[1] = i+1;
	    edata[ielem][jelem].edofs[2] = i+col;
	    edata[ielem][jelem].edofs[3] = i+col+1;
	    /*for(j=0;j<4; j++)
	    {
		if(myrank == master)printf("%d ",edata[ielem][jelem].edofs[j]);
	    }
	    if(myrank == master)printf("\n");*/
	    
	}
	//if(myrank == master)printf("\n");
    }    
    //------------------------------------------------------------------------//
    
    
    markBoundaryDof(dof);

    markSharedDof(dof);

    assignGlobalDof(dof);

    
}

/*void gridread(double **x, double **y, double ****area, double **vol, double **xc, double **yc)
{
    //------------------------------------------------------------------------//
    //Let only the master processor read the mesh
    double *xread, *yread;
    
    allocator1(&xread, (gxelem+1)*(gyelem+1));
    allocator1(&yread, (gxelem+1)*(gyelem+1));

    int i,j,k;
    
    if(myrank == master)
    {
	
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

//------------------------------------------------------------------------//
    //Deallocators
    deallocator1(&xread, (gxelem+1)*(gyelem+1));
    deallocator1(&yread, (gxelem+1)*(gyelem+1));
    deallocator2(&xtemp, gxelem+1, gyelem+1);
    deallocator2(&ytemp, gxelem+1, gyelem+1);
    //------------------------------------------------------------------------//

}*/

/*void meshOperations(double **x, double **y, double ****area, double **vol, double **xc, double **yc)
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
}*/
