/***************************************************************************

Author: nsaini
Created: 2018-03-07

***************************************************************************/


#ifndef COMMU_H
#define COMMU_H

void genibc(int **);

void setupcommu();

void commu2 (double ***);

void destroycommu();

//------------------------------------------------------------------------//
//INS Comm
void INSsetupcommu();

void INScommu2 (double **);

void INSdestroycommu();
//------------------------------------------------------------------------//


#endif
