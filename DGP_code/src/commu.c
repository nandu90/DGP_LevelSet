/***************************************************************************

Author: nsaini
Created: 2018-03-07

***************************************************************************/
#include "commu.h"
#include "common.h"
#include "memory.h"

void genibc(int **iBC)
{
    
  /*This routine works only on ghost elements on each processor
    Convention: 
    if iBC = 0 ==> interior element
    if iBC = 1 ==> interior boundary, only responsible for communication
    if iBC = 2 ==> domain boundary. BC needs to be implemented here
  */

    int i,j;
    
    //------------------------------------------------------------------------//
    //Initialize periodic array to zero
    for(i=0; i<4; i++)
    {
	per[i] = 0;
    }
    //------------------------------------------------------------------------//

  //First consider all ghost elements to be 1
  for(j=0; j<yelem; j++)
    {
      for(i=1; i>=0; i--)
	{
	  iBC[i][j] = 1; //left band
	  iBC[xelem-1-i][j] = 1; //right band
	}
    }
  for(i=0; i<xelem; i++)
    {
      for(j=1; j>=0; j--)
	{
	  iBC[i][j] = 1; //lower band
	  iBC[i][yelem-1-j] = 1; //upper band
	}
    }

  //------------------------------------------------------------------------//
  //Not required anymore
  /*int left = 0;
  int right = 0;
  int up = 0;
  int down = 0;*/
  //------------------------------------------------------------------------//

  
  ///Assume all processors have 4 neighbours to start with*/
  //Neighbour processors will have rank as defined below
  int l_bhai = myrank - 1;
  int r_bhai = myrank + 1;
  int u_bhai = myrank + procm;
  int d_bhai = myrank - procm;

  /*Now take care of ghosts of boundary processors*/
  if((int)floor(myrank/procm) == 0) //Means that this processor is aligned with lower boundary
    {
	//down = 1;
      d_bhai = -1000;           //Set rank of down proc to -ve value
      if(y_bound != 3)           //Dont declare it a boundary node for periodic BC
	{
	  for(i=0; i<xelem; i++)
	    {
	      for(j=1; j>=0; j--)
		{
		  iBC[i][j] = 2;
		}
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
	  for(i=0; i<xelem; i++)
	    {
	      for(j=1; j>=0; j--)
		{
		  iBC[i][yelem-1-j] = 2;
		}
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
	  for(j=0; j<yelem; j++)
	    {
	      for(i=1; i>=0; i--)
		{
		  iBC[i][j] = 2;
		}
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
	  for(j=0; j<yelem; j++)
	    {
	      for(i=1; i>=0; i--)
		{
		  iBC[xelem-1-i][j] = 2;
		}
	    }
	}
      else
	{
	  r_bhai = myrank - procm + 1;
	  per[0] = 1;
	}
    }

  //  printf("%d processor is on left = %d; right = %d; up = %d; down = %d\n",myrank+1,left,right,up,down);

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
	  bhailog[i] = -1;  //If no neighbour mark thatt placeholder as -1
	}
    }

  
  //printf("Bhai log of proc %d are: %d %d %d %d\nPeriodic array is %d %d %d %d\n\n",myrank,bhailog[0],bhailog[1],bhailog[2],bhailog[3],per[0],per[1],per[2],per[3]);

  
}

/*void commu(double ***var)
{
  //if(myrank == master)printf("\nhere\n");
  //    Starting from the right side of mesh, ghost cell strips are numbered in an anticlockwise number
   
  int i,j,k,strip;
  int index;
  int imin[4][2], imax[4][2];  //Here 4 refers to the number of ghost cell strips
  int jmin[4][2], jmax[4][2];
  int size[4];
  imin[0][0] = xelem-3;
  imax[0][0] = xelem-3;
  jmin[0][0] = 0;
  jmax[0][0] = yelem-1;
  imin[0][1] = xelem-1;
  imax[0][1] = xelem-1;
  jmin[0][1] = 0;
  jmax[0][1] = yelem-1;
  size[0] = 2*yelem;

  imin[1][0] = 0;
  imax[1][0] = xelem-1;
  jmin[1][0] = yelem-3;
  jmax[1][0] = yelem-3;
  imin[1][1] = 0;
  imax[1][1] = xelem-1;
  jmin[1][1] = yelem-1;
  jmax[1][1] = yelem-1;
  size[1] = 2*xelem;

  imin[2][0] = 2;
  imax[2][0] = 2;
  jmin[2][0] = 0;
  jmax[2][0] = yelem-1;
  imin[2][1] = 0;
  imax[2][1] = 0;
  jmin[2][1] = 0;
  jmax[2][1] = yelem-1;
  size[2] = 2*yelem;

  imin[3][0] = 0;
  imax[3][0] = xelem-1;
  jmin[3][0] = 2;
  jmax[3][0] = 2;
  imin[3][1] = 0;
  imax[3][1] = xelem-1;
  jmin[3][1] = 0;
  jmax[3][1] = 0;
  size[3] = 2*xelem;

  //if(myrank == master)printf("here starting to send\n");
  //for(i=0; i<yelem; i++)
    //{
      //if(myrank == master)printf("%.f\n",bhai.sendrbuf[i]);
      //}
  //Send operation
    for(k=0; k<4; k++) //Loop over ghost cell strips starting from right and then antic
    {
      if(bhailog[k] >= 0)
	{
	  //Package contents to send
	  int mulx, muly;
	  if(k == 0)
	    {
	      mulx = -1;
	      muly = 0;
	    }
	  else if(k==2)
	    {
	      mulx = 1;
	      muly =0;
	    }
	  else if(k==1)
	    {
	      mulx = 0;
	      muly = -1;
	    }
	  else
	    {
	      mulx = 0;
	      muly = 1;
	    }
	  index=0;
	  for(strip=1; strip>=0; strip--)
	    {
	      for(i=imin[k][0]+(strip*mulx); i<=imax[k][0]+(strip*mulx); i++)  //loop over i index of ghost node. e.g. for right strip imin = imax = xelem-1
		{
		  for(j=jmin[k][0]+(strip*muly); j<=jmax[k][0]+(strip*muly); j++) //loop over j index of ghost node. e.g. for right strip jmin = 0 and jmax = yelem-1
		    {
		      sendptr[k][index] = var[i][j][0]; //fill up the send ptr array
		      //if(myrank== master && k==1)printf("send %d %d %d %d\n",i,j,k,index);
		      index++;
		    }
		}
	    }
	  //sendptr[0] array is sent to the processor on the right... and so on
	  //tag is specified to be myrank
	   MPI_Send(sendptr[k], size[k], MPI_DOUBLE, bhailog[k], myrank, MPI_COMM_WORLD);
	   //if(bhailog[k]==0)printf("here sent from %d to %d\n",myrank,bhailog[k]);
	}
      
    }

    //if(myrank == master)printf("here sent\n");
    
  //Recv OPeration
  for(k=3; k>=0; k--)
    {
      if(bhailog[k] >= 0)
	{
	  //The tag = bhailog[k] ensures that this receive is from the correct processor
	  //eg if bhailog[0] = 1. This means that thr processor to the right of this processor is 1.
	  //if(myrank == master)printf("here about to recv from %d\n",bhailog[k]);
	  MPI_Recv(recvptr[k],size[k], MPI_DOUBLE,bhailog[k],bhailog[k],MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	  //if(myrank == master)printf("here recv from %d by %d\n",bhailog[k],myrank);
	  int mulx, muly;
	  if(k == 0)
	    {
	      mulx = -1;
	      muly = 0;
	    }
	  else if(k==2)
	    {
	      mulx = 1;
	      muly =0;
	    }
	  else if(k==1)
	    {
	      mulx = 0;
	      muly = -1;
	    }
	  else
	    {
	      mulx = 0;
	      muly = 1;
	    }
	  
	  //unpack the contents of recvptr and place in the correct location. e.g. if k=0. this means the contents go on the right ghost strip
	  index=0;
	  for(strip = 0; strip<=1; strip++)
	    {
	      for(i=imin[k][1]+(strip*mulx); i<=imax[k][1]+(strip*mulx); i++)
		{
		  for(j=jmin[k][1]+(strip*muly); j<=jmax[k][1]+(strip*muly); j++)
		    {
		      var[i][j][0] = recvptr[k][index];
		      //if(myrank== 2 && k==3)printf("recv %d %d %d %d\n",i,j,k,index);
		      index++;
		    }
		}
	    }
	}
      //if(myrank == master)printf("received %d\n",k);
    }
  //if(myrank == master)printf("here final\n");
}*/


void commu2(double ***var)
{
  //if(debug == 1)printf("%d reached the commu \n", myrank);
  //    Starting from the right side of mesh, ghost cell strips are numbered in an anticlockwise number
   
    int i,j,k,strip, ng;
    int index;
    int imin[4][2], imax[4][2];  //Here 4 refers to the number of ghost cell strips
    int jmin[4][2], jmax[4][2];
    int size[4];
    imin[0][0] = xelem-3;
    imax[0][0] = xelem-3;
    jmin[0][0] = 0;
    jmax[0][0] = yelem-1;
    imin[0][1] = xelem-1;
    imax[0][1] = xelem-1;
    jmin[0][1] = 0;
    jmax[0][1] = yelem-1;
    size[0] = 2*yelem*ncoeff;
    
    imin[1][0] = 0;
    imax[1][0] = xelem-1;
    jmin[1][0] = yelem-3;
    jmax[1][0] = yelem-3;
    imin[1][1] = 0;
    imax[1][1] = xelem-1;
    jmin[1][1] = yelem-1;
    jmax[1][1] = yelem-1;
    size[1] = 2*xelem*ncoeff;
    
    imin[2][0] = 2;
    imax[2][0] = 2;
    jmin[2][0] = 0;
    jmax[2][0] = yelem-1;
    imin[2][1] = 0;
    imax[2][1] = 0;
    jmin[2][1] = 0;
    jmax[2][1] = yelem-1;
    size[2] = 2*yelem*ncoeff;
    
    imin[3][0] = 0;
    imax[3][0] = xelem-1;
    jmin[3][0] = 2;
    jmax[3][0] = 2;
    imin[3][1] = 0;
    imax[3][1] = xelem-1;
    jmin[3][1] = 0;
    jmax[3][1] = 0;
    size[3] = 2*xelem*ncoeff;

      

    int recvk;
    
  //Package contents
    for(k=0; k<4; k++) //Loop over ghost cell strips starting from right and then anticlockwise
    {
      if(bhailog[k] >= 0)
	{
	  //Package contents to send
	  int smulx, smuly;
	  if(k == 0)
	    {
	      smulx = -1;
	      smuly = 0;
	    }
	  else if(k==2)
	    {
	      smulx = 1;
	      smuly =0;
	    }
	  else if(k==1)
	    {
	      smulx = 0;
	      smuly = -1;
	    }
	  else
	    {
	      smulx = 0;
	      smuly = 1;
	    }
	  index=0;
	  for(strip=1; strip>=0; strip--)
	    {
	      for(i=imin[k][0]+(strip*smulx); i<=imax[k][0]+(strip*smulx); i++)  //loop over i index of ghost node. e.g. for right strip imin = imax = xelem-1
		{
		  for(j=jmin[k][0]+(strip*smuly); j<=jmax[k][0]+(strip*smuly); j++) //loop over j index of ghost node. e.g. for right strip jmin = 0 and jmax = yelem-1
		    {
			//Loop over the quadrature points
			for(ng=0; ng<ncoeff; ng++)
			{
			    sendptr[k][index] = var[i][j][ng]; //fill up the send ptr array
			    index++;
			}
		    }
		}
	    }

	  //if(debug == 1)printf("%d packd for %d\n",myrank, bhailog[k]);
	}

       
    }

    
    
   
     

    //if(myrank == 1 && debug == 1)printf("contents are pakeaged\n");

    MPI_Status status;
    MPI_Request request;
    int tag;
    //if(myrank == master && debug == 1)printf("here starting sendrecv operations\n");
    //------------------------------------------------------------------------//
    //First only let the interior boundaries be communicated
    //If you try to do everything at once, contention occurs
      for(k=0; k<4; k++)
	{
	  if(k==0)recvk = 2;
	  if(k==1)recvk = 3;
	  if(k==2)recvk = 0;
	  if(k==3)recvk = 1;
	   
	  
	  if(bhailog[recvk] >= 0 && per[recvk] == 0)// && bhailog[k] != myrank)
	    {
		
		tag = bhailog[recvk];
		
		MPI_Recv(recvptr[recvk],size[recvk], MPI_DOUBLE,bhailog[recvk],tag,MPI_COMM_WORLD,&status);	      
	    }

	    if(bhailog[k] >= 0 && per[k] == 0)// && bhailog[recvk] != myrank)
	    {
		tag = myrank;
		
	      MPI_Isend(sendptr[k], size[k], MPI_DOUBLE, bhailog[k], tag, MPI_COMM_WORLD,&request);	      
	      MPI_Wait(&request, &status);
	    }
	}
      //Now Communicate the periodic boundaries (if any)
      for(k=0; k<4; k++)
      {
	  if(k==0)recvk = 2;
	  if(k==1)recvk = 3;
	  if(k==2)recvk = 0;
	  if(k==3)recvk = 1;
	   
	  
	  if(bhailog[recvk] >= 0 && per[recvk] == 1)// && bhailog[k] != myrank)
	    {
		
		tag = bhailog[recvk];
		
		MPI_Recv(recvptr[recvk],size[recvk], MPI_DOUBLE,bhailog[recvk],tag,MPI_COMM_WORLD,&status);	      
	    }

	    if(bhailog[k] >= 0 && per[k] == 1)// && bhailog[recvk] != myrank)
	    {
		tag = myrank;
		
	      MPI_Isend(sendptr[k], size[k], MPI_DOUBLE, bhailog[k], tag, MPI_COMM_WORLD,&request);	      
	      MPI_Wait(&request, &status);
	    }
      }
      //------------------------------------------------------------------------//


      //Unpack contents
      for(recvk=0; recvk<4; recvk++)
	{
	  if(bhailog[recvk] >= 0)
	    {
	  
	      int rmulx, rmuly;
	      if(recvk == 0)
		{
		  rmulx = -1;
		  rmuly = 0;
		}
	      else if(recvk==2)
		{
		  rmulx = 1;
		  rmuly =0;
		}
	      else if(recvk==1)
		{
		  rmulx = 0;
		  rmuly = -1;
		}
	      else
		{
		  rmulx = 0;
		  rmuly = 1;
		}
	      
	      //unpack the contents of recvptr and place in the correct location. e.g. if k=0. this means the contents go on the right ghost strip
	      index=0;
	      for(strip = 0; strip<=1; strip++)
		{
		  for(i=imin[recvk][1]+(strip*rmulx); i<=imax[recvk][1]+(strip*rmulx); i++)
		    {
		      for(j=jmin[recvk][1]+(strip*rmuly); j<=jmax[recvk][1]+(strip*rmuly); j++)
		      {
			  //Loop over quadrature points
			  for(ng=0; ng<ncoeff; ng++)
			  {
			      var[i][j][ng] = recvptr[recvk][index];
			      index++;
			  }
		      }
		    }
		}
	    }
	  
	}

    
  
}

void setupcommu()
{
  /*Allocate send and receive arrays on each processor
    These arrays stay uniform since the mesh is uniform.
    Therefore allocate at the start of the program and dellocate at the end
   */

  /*For each ghost strip there is a corresponding sendptr and recvptr array.
    Thus dimension of sendptr is [4][size of that ghost strip]
   */
    if(bhailog[0] >= 0)
    {
	allocator1(&bhai.sendrbuf,2*yelem*ncoeff);
	allocator1(&bhai.recvrbuf,2*yelem*ncoeff);
	sendptr[0] = bhai.sendrbuf;
	recvptr[0] = bhai.recvrbuf;
    }
    if(bhailog[1] >= 0)
    {
	allocator1(&bhai.sendubuf,2*xelem*ncoeff);
	allocator1(&bhai.recvubuf,2*xelem*ncoeff);
	sendptr[1] = bhai.sendubuf;
	recvptr[1] = bhai.recvubuf;
    }
    if(bhailog[2] >= 0)
    {
	allocator1(&bhai.sendlbuf,2*yelem*ncoeff);
	allocator1(&bhai.recvlbuf,2*yelem*ncoeff);
	sendptr[2] = bhai.sendlbuf;
	recvptr[2] = bhai.recvlbuf;
    }
    if(bhailog[3] >= 0)
    {
	allocator1(&bhai.senddbuf,2*xelem*ncoeff);
	allocator1(&bhai.recvdbuf,2*xelem*ncoeff);
	sendptr[3] = bhai.senddbuf;
	recvptr[3] = bhai.recvdbuf;
    }
  
  /*for(i=0; i<yelem; i++)
    {
       if(myrank == master)printf("%.4f",sendptr[0][i]);
       }*/
}

void destroycommu()
{
  
  if(bhailog[0] >= 0)
    {
      deallocator1(&bhai.sendrbuf,2*yelem*ncoeff);
      deallocator1(&bhai.recvrbuf,2*yelem*ncoeff);
    }
  if(bhailog[1] >= 0)
    {
      deallocator1(&bhai.sendubuf,2*xelem*ncoeff);
      deallocator1(&bhai.recvubuf,2*xelem*ncoeff);
    }
  if(bhailog[2] >= 0)
    {
      deallocator1(&bhai.sendlbuf,2*yelem*ncoeff);
      deallocator1(&bhai.recvlbuf,2*yelem*ncoeff);
    }
  if(bhailog[3] >= 0)
    {
      deallocator1(&bhai.senddbuf,2*xelem*ncoeff);
      deallocator1(&bhai.recvdbuf,2*xelem*ncoeff);
    }
}



