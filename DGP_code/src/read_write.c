/***************************************************************************

Author: nsaini
Created: 2018-08-04

***************************************************************************/


#include "common.h"
#include "memory.h"
#include "commu.h"
#include "fileIO.h"
#include "generalFunc.h"

void prevfileread(struct elemsclr elem, double *time, double **x, double **y)
{
    //------------------------------------------------------------------------//
    //Path to timestep directory
    char* dirpath1;
    dirpath1 = concat(getexepath(), "/laststep/");
    char buf_dir1[12];
    snprintf(buf_dir1,12,"%d",startstep);
    dirpath1 = concat(dirpath1,buf_dir1);
    //------------------------------------------------------------------------//

    
    //------------------------------------------------------------------------//
    //Directory operations
    char* dirpath;
    dirpath = concat(dirpath1, "/");
    char buf_dir[12];
    snprintf(buf_dir,12,"%d",myrank);
    dirpath = concat(dirpath,buf_dir);

    char *path;
    path = concat(dirpath,"/u.dat");
    FILE *startu = fopen(path,"r");
    if(startu == NULL)
    {
	printf("Error opening file laststep/u.dat/n");
	exit(1);
    }
    memset(path,0,strlen(path));

    path = concat(dirpath,"/v.dat");
    FILE *startv = fopen(path,"r");
    if(startv == NULL)
    {
	printf("Error opening file laststep/v.dat/n");
	exit(1);
    }
    memset(path,0,strlen(path));

    path = concat(dirpath,"/phi.dat");
    FILE *startphi = fopen(path,"r");
    if(startv == NULL)
    {
	printf("Error opening file laststep/phi.dat/n");
	exit(1);
    }
    memset(path,0,strlen(path));

    path = concat(dirpath,"/time.dat");
    FILE *starttime = fopen(path,"r");
    if(starttime == NULL)
    {
	printf("Error opening file laststep/time.dat/n");
	exit(1);
    }
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Now the fun part
    char* line = NULL;
    ssize_t size;
    size_t len = 0;

    while((size = getline(&line, &len, starttime)) != -1)
    {
	(*time) = atof(line);
	break;
    }

    int ielem;
    int jelem;
    int icoeff;

    //Read u
    /*ielem = 0;
    jelem = 0;
    icoeff = 0;
    while((size = getline(&line, &len, startu)) != -1)
    {
	elem.u[ielem][jelem][icoeff] = atof(line);
	icoeff++;
	if(icoeff == ncoeff)
	{
	    ielem++;
	    icoeff=0;
	    if(ielem == xelem)
	    {
		jelem++;
		ielem=0;
		if(jelem == yelem)break;
	    }
	}
    }

    //Read v
    ielem = 0;
    jelem = 0;
    icoeff = 0;
    while((size = getline(&line, &len, startv)) != -1)
    {
	elem.v[ielem][jelem][icoeff] = atof(line);
	icoeff++;
	if(icoeff == ncoeff)
	{
	    ielem++;
	    icoeff=0;
	    if(ielem == xelem)
	    {
		jelem++;
		ielem=0;
		if(jelem == yelem)break;
	    }
	}
	}*/

    //Read phi
    ielem = 0;
    jelem = 0;
    icoeff = 0;
    while((size = getline(&line, &len, startphi)) != -1)
    {
	elem.phi[ielem][jelem][icoeff] = atof(line);
	icoeff++;
	if(icoeff == ncoeff)
	{
	    ielem++;
	    icoeff=0;
	    if(ielem == xelem)
	    {
		jelem++;
		ielem=0;
		if(jelem == yelem)break;
	    }
	}
    }
    //------------------------------------------------------------------------//


    //------------------------------------------------------------------------//
    free(dirpath);
    free(path);
    fclose(startu);
    fclose(startv);
    fclose(startphi);
    fclose(starttime);
    free(dirpath1);
    //------------------------------------------------------------------------//

}

void filewrite(struct elemsclr elem, double time, int iter)
{
    //------------------------------------------------------------------------//
    //Let the master processor create the timestep directory
    
    char* dirpath1;
    dirpath1 = concat(getexepath(), "/laststep/");
    char buf_dir1[12];
    snprintf(buf_dir1,12,"%d",iter);
    dirpath1 = concat(dirpath1,buf_dir1);
    
    if(myrank == master)
    {
	DIR* dir = opendir("dirpath1");
	if(dir)
	{
	    closedir(dir);
	}
	else if(ENOENT == errno)
	{
	    mkdir(dirpath1, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	}
    }
    MPI_Barrier(MPI_COMM_WORLD);
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Now all processors create sub-directories
    char* dirpath;
    dirpath = concat(dirpath1, "/");
    char buf_dir[12];
    snprintf(buf_dir,12,"%d",myrank);
    dirpath = concat(dirpath,buf_dir);
    DIR* dir1 = opendir("dirpath");
    if(dir1)
    {
	closedir(dir1);
    }
    else if(ENOENT == errno)
    {
	mkdir(dirpath, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    }
    
    
    char *path;
    path = concat(dirpath,"/u.dat");
    FILE *startu = fopen(path,"w");
    if(startu == NULL)
    {
	printf("Error opening file laststep/u.dat/n");
	exit(1);
    }
    memset(path,0,strlen(path));

    path = concat(dirpath,"/v.dat");
    FILE *startv = fopen(path,"w");
    if(startv == NULL)
    {
	printf("Error opening file laststep/v.dat/n");
	exit(1);
    }
    memset(path,0,strlen(path));

    path = concat(dirpath,"/phi.dat");
    FILE *startphi = fopen(path,"w");
    if(startv == NULL)
    {
	printf("Error opening file laststep/phi.dat/n");
	exit(1);
    }
    memset(path,0,strlen(path));

    path = concat(dirpath,"/time.dat");
    FILE *starttime = fopen(path,"w");
    if(starttime == NULL)
    {
	printf("Error opening file laststep/time.dat/n");
	exit(1);
    }
    //------------------------------------------------------------------------//

    //------------------------------------------------------------------------//
    //Print out now
    int ielem, jelem, icoeff;

    for(jelem=0; jelem<yelem; jelem++)
    {
	for(ielem=0; ielem<xelem; ielem++)
	{
	    for(icoeff=0; icoeff<ncoeff; icoeff++)
	    {
		fprintf(startu,"%.12f\n", elem.u[ielem][jelem][icoeff]);
		fprintf(startv,"%.12f\n", elem.v[ielem][jelem][icoeff]);
		fprintf(startphi,"%.12f\n", elem.phi[ielem][jelem][icoeff]);
	    }
	}
    }
    //------------------------------------------------------------------------//

    fprintf(starttime,"%.12f\n",time);

    
    free(dirpath);
    free(path);
    fclose(starttime);
    fclose(startphi);
    fclose(startu);
    fclose(startv);
    free(dirpath1);
}
