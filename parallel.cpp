#include "parallel.h"
#include "LR_lattice.h"
#include "io.h"

int N;
int ProcNum  = 0;
int ProcRank = 0;
int OwnCoords[2];
double *null_p = new double;
void CreateGrid(int GridSize,MPI_Comm *Comm, MPI_Comm* old_comm)
{
	int dims[2], periods[2], reorder=0;
	dims[0] = dims[1] = GridSize;
	periods[0] = periods[1] = 0;
	//MPI_Dims_create(ProcNum, 2, dims);
	MPI_Cart_create(*old_comm,2,dims,periods,reorder, Comm);
	MPI_Cart_coords(*Comm,ProcRank,2,OwnCoords);

}

void Exchange(MPI_Comm *Comm, int GridSize, double* V,int N, MPI_Datatype *plane_x, MPI_Datatype *plane_y)
{
	int NextCoords[2];
	int NextRank;
	MPI_Status status;
	if (OwnCoords[1] % 2 == 0)
	{
		if (OwnCoords[1]!=GridSize-1)
		{	//sending right column and recieving left column
			NextCoords[0] = OwnCoords[0];
			NextCoords[1] = OwnCoords[1]+1;
			MPI_Cart_rank(*Comm,NextCoords,&NextRank);
			MPI_Sendrecv(&(V[N-2]),1,*plane_x,NextRank,0,&(V[N-1]),1,*plane_x,NextRank,1,*Comm,&status);
		}
		else
		{	//boundary condition of whole lattice(right column)
			for (int k=0; k<Height; k++)
				for(int i=0;i<N;i++)
				{
					V[(i+k*N)*N+N-1]=V[(i+k*N)*N+N-2];
				}
		}
		if (OwnCoords[1]!=0)
		{	//sending left column and recieving right column
			NextCoords[0] = OwnCoords[0];
			NextCoords[1] = OwnCoords[1]-1;
			MPI_Cart_rank(*Comm,NextCoords,&NextRank);
			MPI_Sendrecv(&(V[1]),1,*plane_x,NextRank,1,&(V[0]),1,*plane_x,NextRank,0,*Comm,&status);
		}
		else
		{	//boundary condition of whole lattice(left column)
			for (int k=0; k<Height; k++)
				for(int i=0;i<N;i++)
				{
					V[(i+k*N)*N+0]=V[(i+k*N)*N+1];
				}
		}
	}
	else
	{
		if (OwnCoords[1]!=0)
		{	//sending left column and receiving right column
			NextCoords[0] = OwnCoords[0];
			NextCoords[1] = OwnCoords[1]-1;
			MPI_Cart_rank(*Comm,NextCoords,&NextRank);
			MPI_Sendrecv(&(V[1]),1,*plane_x,NextRank,1,&(V[0]),1,*plane_x,NextRank,0,*Comm,&status);
		}
		else
		{	//boundary condition of whole lattice(left column)
			for (int k=0; k<Height; k++)
				for(int i=0;i<N;i++)
				{
					V[(i+k*N)*N+0]=V[(i+k*N)*N+1];
				}
		}
		if (OwnCoords[1]!=GridSize-1)
		{	//sending right column and recieving left column
			NextCoords[0] = OwnCoords[0];
			NextCoords[1] = OwnCoords[1]+1;
			MPI_Cart_rank(*Comm,NextCoords,&NextRank);
			MPI_Sendrecv(&(V[N-2]),1,*plane_x,NextRank,0,&(V[N-1]),1,*plane_x,NextRank,1,*Comm,&status);
		}
		else
		{	//boundary condition of whole lattice(right column)
			for (int k=0; k<Height; k++)
				for(int i=0;i<N;i++)
				{
					V[(i+k*N)*N+N-1]=V[(i+k*N)*N+N-2];
				}
		}
	}
	if(OwnCoords[0] % 2 == 0)
	{
		if (OwnCoords[0]!=GridSize-1)
		{	//sending bottom row and recieving top row
			NextCoords[0] = OwnCoords[0]+1;
			NextCoords[1] = OwnCoords[1];
			MPI_Cart_rank(*Comm,NextCoords,&NextRank);
			MPI_Sendrecv(&(V[(N-2)*N]),1,*plane_y,NextRank,2,&(V[(N-1)*N]),1,*plane_y,NextRank,3,*Comm,&status);
		}
		else
		{	//boundary condition of whole lattice(bottom row)
			for (int k=0; k<Height; k++)
				for(int j=0;j<N;j++)
				{
					V[(N-1+k*N)*N+j]=V[(N-2+k*N)*N+j];
				}
		}
		if (OwnCoords[0]!=0)
		{	//sending top row and recieving bottom row
			NextCoords[0] = OwnCoords[0]-1;
			NextCoords[1] = OwnCoords[1];
			MPI_Cart_rank(*Comm,NextCoords,&NextRank);
			MPI_Sendrecv(&(V[N]),1,*plane_y,NextRank,3,&(V[0]),1,*plane_y,NextRank,2,*Comm,&status);
		}
		else
		{	//boundary condition of whole lattice(top row)
			for (int k=0; k<Height; k++)
				for(int j=0;j<N;j++)
				{
					V[k*N*N+j]=V[(1+k*N)*N+j];
				}
		}
	}
	else
	{
		if (OwnCoords[0]!=0)
		{	//sending top row and recieving bottom row
			NextCoords[0] = OwnCoords[0]-1;
			NextCoords[1] = OwnCoords[1];
			MPI_Cart_rank(*Comm,NextCoords,&NextRank);
			MPI_Sendrecv(&(V[N]),1,*plane_y,NextRank,3,&(V[0]),1,*plane_y,NextRank,2,*Comm,&status);
		}
		else
		{	//boundary condition of whole lattice(top row)
			for (int k=0; k<Height; k++)
				for(int j=0;j<N;j++)
				{
					V[k*N*N+j]=V[(1+k*N)*N+j];
				}
		}
		if (OwnCoords[0]!=GridSize-1)
		{	//sending bottom row and recieving top row
			NextCoords[0] = OwnCoords[0]+1;
			NextCoords[1] = OwnCoords[1];
			MPI_Cart_rank(*Comm,NextCoords,&NextRank);
			MPI_Sendrecv(&(V[(N-2)*N]),1,*plane_y,NextRank,2,&(V[(N-1)*N]),1,*plane_y,NextRank,3,*Comm,&status);
		}
		else
		{	//boundary condition of whole lattice(bottom row)
			for (int k=0; k<Height; k++)
				for(int j=0;j<N;j++)
				{
					V[(N-1+k*N)*N+j]=V[(N-2+k*N)*N+j];
				}
		}
	}
}



void save(short* Vv,int Nn,int fd)
{
	char* buffer = (char *)Vv;
	write(fd,buffer,(Nn)*sizeof(short));
}