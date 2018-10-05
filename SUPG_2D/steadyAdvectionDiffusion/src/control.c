/***************************************************************************

Author: nsaini
Created: 2018-03-04

***************************************************************************/

#include "fileIO.h"
#include "common.h"

void control()
{
  if(myrank == master)
    {
	FILE *controlfile;
	controlfile = fopen("control.txt","r");
	if(controlfile == NULL)
	{
	    printf("Error opening control.txt!\n");
	    exit(0);
	}
  
	char* line = NULL;
	ssize_t size;
	size_t len = 0;
	
	char delim [1] = " ";
	char* word = NULL;
    
	while((size = getline(&line, &len, controlfile)) != -1)
	{
      
	    word = strtok(line, delim);
	    while(word != NULL)
	    {
		//printf("%s\n",word);
		if(strcmp(word,"Basis_Order") == 0)
		{
		    word = strtok(NULL,delim);
		    polyorder = atoi(word);
		}
		if(strcmp(word,"xlen") == 0)
		{
		    word = strtok(NULL,delim);
		    xlen = atof(word);
		}
		if(strcmp(word,"ylen") == 0)
		{
		    word = strtok(NULL,delim);
		    ylen = atof(word);
		}
		if(strcmp(word,"xelem") == 0)
		{
		    word = strtok(NULL,delim);
		    gxelem = atoi(word);
		}
		if(strcmp(word,"yelem") == 0)
		{
		    word = strtok(NULL,delim);
		    gyelem = atoi(word);
		}
		if(strcmp(word,"Read_mesh") == 0)
		{
		    word = strtok(NULL,delim);
		    meshread= atoi(word);
		}
	       
		if(strcmp(word,"x-boundary") == 0)
		{
		    word = strtok(NULL, delim);
		    if(strcmp(word,"no-slip\n") == 0)
		    {
			x_bound=1;
		    }
		    else if(strcmp(word,"slip\n") == 0)
		    {
			x_bound=2;
		    }
		    else if(strcmp(word,"periodic\n") == 0)
		    {
			x_bound=3;
		    }
		}
		
		if(strcmp(word,"y-boundary") == 0)
		{
		    word = strtok(NULL, delim);
		    if(strcmp(word,"no-slip\n") == 0)
		    {
			y_bound=1;
		    }
		    else if(strcmp(word,"slip\n") == 0)
		    {
			y_bound=2;
		    }
		    else if(strcmp(word,"periodic\n") == 0)
		    {
			y_bound=3;
		    }
		}
	        
		
		word = strtok(NULL,delim);
	    }
	    
	}
	fclose(controlfile);
	
	free(word);
	free(line);
    }

  MPI_Bcast(&polyorder,1,MPI_INT,master,MPI_COMM_WORLD);
  
  MPI_Bcast(&xlen,1,MPI_DOUBLE,master,MPI_COMM_WORLD);
  MPI_Bcast(&ylen,1,MPI_DOUBLE,master,MPI_COMM_WORLD);
  
  MPI_Bcast(&gxelem,1,MPI_INT,master,MPI_COMM_WORLD);
  MPI_Bcast(&gyelem,1,MPI_INT,master,MPI_COMM_WORLD);
  MPI_Bcast(&meshread,1,MPI_INT,master,MPI_COMM_WORLD);

  MPI_Bcast(&x_bound,1,MPI_INT,master,MPI_COMM_WORLD);
  MPI_Bcast(&y_bound,1,MPI_INT,master,MPI_COMM_WORLD);
  

    
}
