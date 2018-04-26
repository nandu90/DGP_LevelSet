/***************************************************************************

Author: nsaini
Created: 2018-04-25

***************************************************************************/

#include "commu.h"
#include "common.h"
#include "memory.h"

void INScommu2(double **var)
{
  //if(debug == 1)printf("%d reached the commu \n", myrank);
  /*
    Starting from the right side of mesh, ghost cell strips are numbered in an anticlockwise number
   */
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

  
 

  int recvk;

  //Package contents
    for(k=0; k<4; k++) //Loop over ghost cell strips starting from right and then antic
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
		      INSsendptr[k][index] = var[i][j]; //fill up the send ptr array
		      index++;
		    }
		}
	    }

	  //if(debug == 1)printf("%d packd for %d\n",myrank, bhailog[k]);
	}

       
    }

    //if(myrank == 1 && debug == 1)printf("contents are pakeaged\n");

    MPI_Status status;
    
    MPI_Request request;
    //if(myrank == master && debug == 1)printf("here starting sendrecv operations\n");
      for(k=0; k<4; k++)
	{
	  if(k==0)recvk = 2;
	  if(k==1)recvk = 3;
	  if(k==2)recvk = 0;
	  if(k==3)recvk = 1;
	   
	  
	  if(bhailog[recvk] >= 0)
	    {
	      //if(myrank==0 && debug ==1)printf("here recving from %d by %d\n",bhailog[recvk],myrank);
	      MPI_Recv(INSrecvptr[recvk],size[recvk], MPI_DOUBLE,bhailog[recvk],bhailog[recvk],MPI_COMM_WORLD,&status);
	      
	    }

	  if(bhailog[k] >= 0)
	    {
	      //if(bhailog[k]==0 && debug == 1)printf("here sending from %d to %d\n",myrank,bhailog[k]);
	      MPI_Isend(INSsendptr[k], size[k], MPI_DOUBLE, bhailog[k], myrank, MPI_COMM_WORLD,&request);
	      
	      MPI_Wait(&request, &status);
	    }
	}


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
			  var[i][j] = INSrecvptr[recvk][index];
			  index++;
		      }
		    }
		}
	    }
	  
	}

    
  
}

void INSsetupcommu()
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
      allocator1(&INSbhai.sendrbuf,2*yelem);
      allocator1(&INSbhai.recvrbuf,2*yelem);
      INSsendptr[0] = INSbhai.sendrbuf;
      INSrecvptr[0] = INSbhai.recvrbuf;
    }
  if(bhailog[1] >= 0)
    {
      allocator1(&INSbhai.sendubuf,2*xelem);
      allocator1(&INSbhai.recvubuf,2*xelem);
      INSsendptr[1] = INSbhai.sendubuf;
      INSrecvptr[1] = INSbhai.recvubuf;
    }
  if(bhailog[2] >= 0)
    {
      allocator1(&INSbhai.sendlbuf,2*yelem);
      allocator1(&INSbhai.recvlbuf,2*yelem);
      INSsendptr[2] = INSbhai.sendlbuf;
      INSrecvptr[2] = INSbhai.recvlbuf;
    }
  if(bhailog[3] >= 0)
    {
      allocator1(&INSbhai.senddbuf,2*xelem);
      allocator1(&INSbhai.recvdbuf,2*xelem);
      INSsendptr[3] = INSbhai.senddbuf;
      INSrecvptr[3] = INSbhai.recvdbuf;
      }
  /*for(i=0; i<yelem; i++)
    {
       if(myrank == master)printf("%.4f",sendptr[0][i]);
       }*/
}

void INSdestroycommu()
{
  
  if(bhailog[0] >= 0)
    {
      deallocator1(&INSbhai.sendrbuf,2*yelem);
      deallocator1(&INSbhai.recvrbuf,2*yelem);
    }
  if(bhailog[1] >= 0)
    {
      deallocator1(&INSbhai.sendubuf,2*xelem);
      deallocator1(&INSbhai.recvubuf,2*xelem);
    }
  if(bhailog[2] >= 0)
    {
      deallocator1(&INSbhai.sendlbuf,2*yelem);
      deallocator1(&INSbhai.recvlbuf,2*yelem);
    }
  if(bhailog[3] >= 0)
    {
      deallocator1(&INSbhai.senddbuf,2*xelem);
      deallocator1(&INSbhai.recvdbuf,2*xelem);
    }
}
