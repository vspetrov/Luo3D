#include <time.h>
#include "math.h"
#include <stdio.h>
#include "LR_cell.h"
#include "LR_lattice.h"
#include "parallel.h"

#include "stdlib.h"
#include "fcntl.h"

#include <stdio.h>

int main(int argc, char *argv[])
{
    double t1,t2,duration;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);//Number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);//Rank of process
    t1 = MPI_Wtime();

    MPI_Group WorldGroup, CalculatorGroup;
    MPI_Comm Calculators;
    int ranks[1];
    ranks[0] = ProcNum-1;
    MPI_Comm_group(MPI_COMM_WORLD, &WorldGroup);
    MPI_Group_excl(WorldGroup, 1, ranks, &CalculatorGroup);
    MPI_Comm_create(MPI_COMM_WORLD,CalculatorGroup,&Calculators);


    int GridSize= sqrt((double)(ProcNum-1));//size of virtual topology(Grid)
    N=Size/GridSize+2;//N-size of subsystem; +2 - for boundary condition


    if (ProcRank!=ProcNum-1)
    {

        MPI_Datatype column;
        MPI_Type_vector(N,1,N,MPI_DOUBLE,&column);
        MPI_Type_commit(&column);

        MPI_Datatype plane_x;
        MPI_Type_hvector(Height,1,N*N*8,column,&plane_x);
        MPI_Type_commit(&plane_x);

        MPI_Datatype plane_y;
        MPI_Type_vector(Height,N,N*N,MPI_DOUBLE,&plane_y);
        MPI_Type_commit(&plane_y);


        //creating datatype for square
        MPI_Datatype square;
        MPI_Type_vector(N-2,N-2,N,MPI_DOUBLE,&square);
        MPI_Type_commit(&square);

        MPI_Comm GridComm;
        CreateGrid(GridSize, &GridComm, &Calculators);
        LattInit( N, Height, "sp_init.bin");
        SolveEquations(1000, GridComm, GridSize, plane_x, plane_y, square);
        MPI_Type_free(&square);
        MPI_Type_free(&column);
        MPI_Type_free(&plane_x);
        MPI_Type_free(&plane_y);
        MPI_Group_free(&CalculatorGroup);
        MPI_Comm_free(&Calculators);
        delete[] V;
        delete[] Cell;
        delete[] I_ext;
    }
    else
    {
        double* V_all = new double[(N-2)*(N-2)*ProcNum];
        double* V_temp = new double[(N-2)*(N-2)];
        short* V_save = new short[(N-2)*(N-2)*ProcNum];
        int ii;
        //printf("Total number of frames: %i\n", DrawNum);
#ifdef OS_WINDOWS
        int fd = open("spiral.bin",O_RDWR|O_CREAT | O_BINARY,S_IREAD|S_IWRITE);
#elif defined(OS_LINUX)
        int fd = open("spiral.bin",O_RDWR|O_CREAT,S_IREAD|S_IWRITE);
#endif
        for (int i=0;i<DrawNum;i++)
        {
            for (int k=0;k<Height; k++)
            {
            MPI_Gather(V_temp,(N-2)*(N-2),MPI_DOUBLE,V_all,(N-2)*(N-2),MPI_DOUBLE,ProcNum-1,MPI_COMM_WORLD);
            for (ii=0;ii<(N-2)*(N-2)*(ProcNum-1);ii++)
            {
                V_save[ii]=short(V_all[ii]*250.);
            }
            save(V_save,(N-2)*(N-2)*(ProcNum-1),fd);
            }
            //printf("Process %i: I've saved the frame number %i\n",ProcRank,i+1);
        }
        close(fd);
//		MPI_Gather(V_temp,(N-2)*(N-2),MPI_DOUBLE,V_all,(N-2)*(N-2),MPI_DOUBLE,ProcNum-1,MPI_COMM_WORLD);
//		for (ii=0;ii<(N-2)*(N-2)*(ProcNum-1);ii++)
//		{
//			V_save[ii]=short(V_all[ii]*250.);
//		}
//int		fd = open("spiral_freqs.bin",O_RDWR|O_CREAT | O_BINARY,S_IREAD|S_IWRITE);
//		save(V_save,(N-2)*(N-2)*(ProcNum-1),fd);
//		close(fd);
        delete[] V_all;
        delete[] V_temp;
        delete[] V_save;
    }





    t2 = MPI_Wtime();
    if (ProcRank==0)
    {
        printf("Experiment duration: %f seconds\n",t2-t1);
    }
    MPI_Finalize();

    return 0;
}
