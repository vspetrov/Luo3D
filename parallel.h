#ifndef _PARALLEL_H
#define _PARALLEL_H

#include "mpi.h"

extern int ProcNum;
extern int ProcRank;
extern int N;
extern double *null_p;
extern int OwnCoords[2];
void CreateGrid(int GridSize,MPI_Comm *Comm, MPI_Comm* old_comm);
void Exchange(MPI_Comm *Comm, int GridSize, double* V,int N, MPI_Datatype *plane_x, MPI_Datatype *plane_y);
extern void save(short* V,int N,int fd);

#endif