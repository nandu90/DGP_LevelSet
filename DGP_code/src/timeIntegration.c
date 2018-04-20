/***************************************************************************

Author: nsaini
Created: 2018-04-19

***************************************************************************/


#include "common.h"
#include "solvers.h"

void euler(double ***scalar, double ****mass, double ***rhs, double deltat)
{
    //------------------------------------------------------------------------//
    //Loop Indexes
    int ielem, jelem;
    int icoeff;
    //------------------------------------------------------------------------//

    for(ielem = 1; ielem < xelem-1; ielem++)
    {
	for(jelem = 1; jelem < yelem-1; jelem++)
	{
	    for(icoeff=0; icoeff < ncoeff; icoeff++)
	    {
		scalar[ielem][jelem][icoeff] += deltat*rhs[ielem][jelem][icoeff]/mass[ielem][jelem][icoeff][icoeff];
	    }
	}
    }
}
