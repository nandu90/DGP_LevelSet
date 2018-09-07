/***************************************************************************

Author: nsaini
Created: 2018-07-24

***************************************************************************/


#include "common.h"
#include "memory.h"
#include "functions.h"

void level_setBC(double **scalar)
{
    int icoeff;
    
    for(icoeff=0; icoeff<ncoeff; icoeff++)
    {
	scalar[xelem-1][icoeff] = scalar[1][icoeff];
	scalar[0][icoeff] = scalar[xelem-2][icoeff];
    }
}
