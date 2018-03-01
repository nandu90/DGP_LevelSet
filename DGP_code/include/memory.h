
#ifndef MEMORY_H
#define MEMORY_H

/*********Double arrays**********/
//1d array
void allocator1 (double **, int);
void deallocator1 (double **, int);

//2d array
void allocator2(double ***, int, int);
void deallocator2(double ***, int, int);

//3d array
void allocator3(double ****, int, int, int);
void deallocator3(double ****, int, int, int);

//4d array
void allocator4(double *****, int, int, int, int);
void deallocator4(double *****, int, int, int, int);
/*************************************/

/*********Integer arrays***********/
//2d arrays
void iallocator2(int ***, int, int);
void ideallocator2(int ***, int, int);
/**********************************/
#endif
