#include "LR_lattice.h"
#include "LR_cell.h"
#include "parallel.h"
#include "math.h"
#include <time.h>

#include "stdlib.h"


#ifdef OS_LINUX
#define _lseek lseek
#include <unistd.h>
#endif

#include <stdio.h>

const int Size=200;
const int Height=50;

int DrawNum = 1;
int resolution = 1;

double **Max;
int *MaxCount;
int *flag;
const double threshould = 0.;
const double threshould_1 = -50.;
double *freq;

CellVariables *Cell;
double *I_ext;
double *V;

int ln, rn, un, dn; //variables for defining coupling neighbours

const double dt=0.1;
double D=0.05;

void LattInit(const int &size, const int &height)
{
    srand((ProcRank+1)*10);
    Cell = new CellVariables[size*size*height];
    I_ext = new double[size*size*height];
    V = new double[size*size*height];

    //define all lattice variables
    int i; //int counters
    for (i=0;i<size*size*height;i++)
    {
        V[i]=-70.;

        Cell[i].m=0.0143576;
        Cell[i].h=0.745286;
        Cell[i].j=0.76125;
        Cell[i].d=0.00934794;
        Cell[i].f=0.999772;
        Cell[i].X=0.018801218;
        Cell[i].Cai=0.00025;
        Cell[i].tension = 0;
        Cell[i].s_tension = 0;
    }
    //define external currents
    for (i=0;i<size*size*height;i++)
        I_ext[i]=-2.4+(double)rand()/double(RAND_MAX)*(2.4-3.2);

    //if (ProcRank == 12)
    //{
    //	for (int k=0; k<Height; k++)
    //		for (i=0; i<size; i++)
    //			for (int j=0; j<size; j++)
    //				I_ext[(i+k*size)*size+j] = -2.8;
    //}
    //printf("Process %i finished\n\n", ProcRank);
}


void SolveEquations(double MaxTime, MPI_Comm &GridComm, int &GridSize, MPI_Datatype &plane_x, MPI_Datatype &plane_y, MPI_Datatype &square)
{
    int i,j,k; //counting variables
    int time, MT;//time itarator
    int counter=1;
    MT=int(MaxTime/dt);

    double *Vs = new double[N*N];
    int drawSkip = int(999/dt);
    int drawEvent = int((MT-drawSkip)/DrawNum);
    int drawCounter=0;
    int freqskip = int(2000/dt);
    drawSkip--;
    //integrating on the interval from time to time+dt/2
    Exchange(&GridComm,GridSize,V,N,&plane_x, &plane_y);

    for (k=0; k<Height; k++)
        for (i=1; i<N-1; i++)
            for (j=1; j<N-1; j++)
                V[N*(i+k*N)+j]+=dt/2.*Coupling(i,j,k);

    int GS=sqrt((double)(ProcNum-1));


    for (time=0; time<MT; time++)
    {
        //D=D+1e-5;
        //solve system with NO coupling on the interval from time to time+dt
        for (k=0; k<Height; k++)
            for (i=1; i<N-1; i++)
                for (j=1; j<N-1; j++)
                    OdeSolve(i, j, k);




        Exchange(&GridComm,GridSize,V,N,&plane_x, &plane_y);



        //save data to the file
        if (time > drawSkip)	drawCounter++;
        if (drawCounter == drawEvent)
        {
            for (k=0;k<Height;k++)
            {
            for (i=0; i<N; i++)
                for (j=0; j<N; j++)
                    Vs[i*N+j]=V[(i+k*N)*N+j];

            MPI_Gather(&Vs[1*N+1],1,square,null_p,(N-2)*(N-2),MPI_DOUBLE,ProcNum-1,MPI_COMM_WORLD);
            }
            drawCounter=0;
        }




            //if (time > freqskip)
            //{
            //for (i=1; i<N-1; i++)
            //	for (j=1; j<N-1; j++)
            //	{
            //		if ((flag[i*N+j]  == -1)&&(V[i*N+j]>threshould))
            //		{
            //			if (MaxCount[i*N+j] == 0)
            //				Max[i*N+j][0] = time*dt;
            //			else
            //				Max[i*N+j][1]=time*dt;

            //				MaxCount[i*N+j]++;
            //				flag[i*N+j] = 1;
            //		}
            //		if ((flag[i*N+j] == 1)&&(V[i*N+j] < threshould_1))
            //			flag[i*N+j] = -1;
            //	}
            //}


        //again calculate coupling
        for (k=0; k<Height; k++)
            for (i=1; i<N-1; i++)
                for (j=1; j<N-1; j++)
                    V[(i+k*N)*N+j]+=dt*Coupling(i,j,k);


    }

    //left calculation of coupling
    for (k=0; k<Height; k++)
        for (i=1; i<N-1; i++)
            for (j=1; j<N-1; j++)
                V[(i+k*N)*N+j]+=dt/2.*Coupling(i,j,k);

            //for (i=1; i<N-1; i++)
            //for (j=1; j<N-1; j++)
            //{
            //	if (MaxCount[i*N+j] > 1)
            //		freq[i*N+j] = 1000.*(MaxCount[i*N+j]-1)/(Max[i*N+j][1]-Max[i*N+j][0]);
            //	else
            //		freq[i*N+j] = 0.;
            //}





        for (i=0; i<N; i++)
            for (j=0; j<N; j++)
                Vs[i*N+j]=V[(i+2*N)*N+j];

//	MPI_Gather(&Vs[1*N+1],1,square,null_p,(N-2)*(N-2),MPI_DOUBLE,ProcNum-1,MPI_COMM_WORLD);

}

void OdeSolve(int &ii, int &jj, int &kk)
{
    double vd;
    int kstep=1;
    double delta_t=dt;
    for(int i=0; i<kstep; i++)
    {
        vd=VFunction(ii,jj, kk);

        if(i==0) //decide on the time substep
        {
            kstep=Substeps(vd);
            delta_t=dt/(double)kstep;
        }

        V[(ii+kk*N)*N+jj]+=delta_t*vd;
        Cell[(ii+kk*N)*N+jj].m=mFunction(delta_t);
        Cell[(ii+kk*N)*N+jj].h=hFunction(delta_t);
        Cell[(ii+kk*N)*N+jj].j=jFunction(delta_t);
        Cell[(ii+kk*N)*N+jj].d=dFunction(delta_t);
        Cell[(ii+kk*N)*N+jj].f=fFunction(delta_t);
        Cell[(ii+kk*N)*N+jj].X=XFunction(delta_t);
        Cell[(ii+kk*N)*N+jj].Cai+=delta_t*CaiFunction();
        Cell[(ii+kk*N)*N+jj].s_tension += delta_t * s_tensionFunction();
    }
}

inline int Substeps(double &vd)
{
    const int kmax=100;
    int k;

    const int k0=vd>0. ? 5 : 1;
    k=k0+(int)fabs(vd);

    return k<kmax ? k : kmax;
}

double Coupling(int &ii, int &jj, int &kk)
{
    if ((kk+1) < Height) un = kk+1;
    else un = kk;

    if ((kk-1) > -1) dn = kk-1;
    else dn = kk;

    return D*(V[(ii+kk*N)*N+jj+1]+V[(ii+kk*N)*N+jj-1]+V[(ii+1+kk*N)*N+jj]+V[(ii-1+kk*N)*N+jj]+
                V[(ii+un*N)*N+jj]+V[(ii+dn*N)*N+jj]-6.*V[(ii+kk*N)*N+jj]);
}

void LattInit(const int &size,  const int &height, const char *a)
{
    //srand((ProcRank+1)*10);
    Cell = new CellVariables[size*size*height];
    I_ext = new double[size*size*height];
    V = new double[size*size*height];
    int i,j,k; //int counters
    for (i=0;i<size*size*height;i++)
    {
        V[i]=-70.;

        Cell[i].m=0.0143576;
        Cell[i].h=0.745286;
        Cell[i].j=0.76125;
        Cell[i].d=0.00934794;
        Cell[i].f=0.999772;
        Cell[i].X=0.018801218;
        Cell[i].Cai=0.00025;
        Cell[i].tension = 0;
        Cell[i].s_tension = 0;
    }
//	V_av = new double[size*size];

#ifdef OS_WINDOWS
    int file_id=open(a,O_RDONLY | O_BINARY, S_IREAD );
#elif defined(OS_LINUX)
    int file_id=open(a,O_RDONLY);
#endif
    if (ProcRank != ProcNum - 1)
    {
        //read V
        _lseek(file_id,((N-2)*Size*OwnCoords[0]+(N-2)*OwnCoords[1])*sizeof(double),SEEK_CUR);
        for (int i=0; i < N-2; i++)
        {
            if (-1 == read(file_id,&(V[N+i*N+1]),(N-2)*sizeof(double))){
                goto Error;
            }
            _lseek(file_id, (Size-N+2)*sizeof(double), SEEK_CUR);
        }

        //read m __1 - skip
        _lseek(file_id,Size*Size*sizeof(double),SEEK_SET);
        _lseek(file_id,((N-2)*Size*OwnCoords[0]+(N-2)*OwnCoords[1])*sizeof(double),SEEK_CUR);
        for (int i=1; i < N-1; i++)
        {
            for (int j=1;j<N-1;j++)
                if (-1 == read(file_id,&(Cell[i*N+j].m),sizeof(double))){
                    goto Error;
                }
            _lseek(file_id, (Size-N+2)*sizeof(double), SEEK_CUR);
        }

        //read h __2 - skip
        _lseek(file_id,2*Size*Size*sizeof(double),SEEK_SET);
        _lseek(file_id,((N-2)*Size*OwnCoords[0]+(N-2)*OwnCoords[1])*sizeof(double),SEEK_CUR);
        for (int i=1; i < N-1; i++)
        {
            for (int j=1;j<N-1;j++)
                if (-1 == read(file_id,&(Cell[i*N+j].h),sizeof(double))){
                    goto Error;
                }
            _lseek(file_id, (Size-N+2)*sizeof(double), SEEK_CUR);
        }

        //read j __3 - skip
        _lseek(file_id,3*Size*Size*sizeof(double),SEEK_SET);
        _lseek(file_id,((N-2)*Size*OwnCoords[0]+(N-2)*OwnCoords[1])*sizeof(double),SEEK_CUR);
        for (int i=1; i < N-1; i++)
        {
            for (int j=1;j<N-1;j++)
                if (-1 == read(file_id,&(Cell[i*N+j].j),sizeof(double))){
                    goto Error;
                }
            _lseek(file_id, (Size-N+2)*sizeof(double), SEEK_CUR);
        }

        //read d __4 - skip
        _lseek(file_id,4*Size*Size*sizeof(double),SEEK_SET);
        _lseek(file_id,((N-2)*Size*OwnCoords[0]+(N-2)*OwnCoords[1])*sizeof(double),SEEK_CUR);
        for (int i=1; i < N-1; i++)
        {
            for (int j=1;j<N-1;j++)
                if (-1 == read(file_id,&(Cell[i*N+j].d),sizeof(double))){
                    goto Error;
                }
            _lseek(file_id, (Size-N+2)*sizeof(double), SEEK_CUR);
        }

        //read f __5 - skip
        _lseek(file_id,5*Size*Size*sizeof(double),SEEK_SET);
        _lseek(file_id,((N-2)*Size*OwnCoords[0]+(N-2)*OwnCoords[1])*sizeof(double),SEEK_CUR);
        for (int i=1; i < N-1; i++)
        {
            for (int j=1;j<N-1;j++)
                if (-1 == read(file_id,&(Cell[i*N+j].f),sizeof(double))){
                    goto Error;
                }
            _lseek(file_id, (Size-N+2)*sizeof(double), SEEK_CUR);
        }

        //read X __6 - skip
        _lseek(file_id,6*Size*Size*sizeof(double),SEEK_SET);
        _lseek(file_id,((N-2)*Size*OwnCoords[0]+(N-2)*OwnCoords[1])*sizeof(double),SEEK_CUR);
        for (int i=1; i < N-1; i++)
        {
            for (int j=1;j<N-1;j++)
                if (-1 == read(file_id,&(Cell[i*N+j].X),sizeof(double))){
                    goto Error;
                }
            _lseek(file_id, (Size-N+2)*sizeof(double), SEEK_CUR);
        }

        //read Cai __7 - skip
        _lseek(file_id,7*Size*Size*sizeof(double),SEEK_SET);
        _lseek(file_id,((N-2)*Size*OwnCoords[0]+(N-2)*OwnCoords[1])*sizeof(double),SEEK_CUR);
        for (int i=1; i < N-1; i++)
        {
            for (int j=1;j<N-1;j++)
                if (-1 == read(file_id,&(Cell[i*N+j].Cai),sizeof(double))){
                    goto Error;
                }
            _lseek(file_id, (Size-N+2)*sizeof(double), SEEK_CUR);
        }
    }
    close(file_id);





        for (k=1; k<20; k++)
            for (i=0; i<N; i++)
                for (j=0; j<N; j++)
                {
                    V[(i+k*N)*N+j]=V[i*N+j];
                    Cell[(i+k*N)*N+j].m=Cell[i*N+j].m;
                    Cell[(i+k*N)*N+j].h=Cell[i*N+j].h;
                    Cell[(i+k*N)*N+j].j=Cell[i*N+j].j;
                    Cell[(i+k*N)*N+j].d=Cell[i*N+j].d;
                    Cell[(i+k*N)*N+j].f=Cell[i*N+j].f;
                    Cell[(i+k*N)*N+j].X=Cell[i*N+j].X;
                    Cell[(i+k*N)*N+j].Cai=Cell[i*N+j].Cai;
                }

        for (k=20; k<Height; k++)
            for (i=0; i<N; i++)
                for (j=0; j<N; j++)
                {
                    V[(i+k*N)*N+j]=V[N+1];
                    Cell[(i+k*N)*N+j].m=Cell[N+1].m;
                    Cell[(i+k*N)*N+j].h=Cell[N+1].h;
                    Cell[(i+k*N)*N+j].j=Cell[N+1].j;
                    Cell[(i+k*N)*N+j].d=Cell[N+1].d;
                    Cell[(i+k*N)*N+j].f=Cell[N+1].f;
                    Cell[(i+k*N)*N+j].X=Cell[N+1].X;
                    Cell[(i+k*N)*N+j].Cai=Cell[N+1].Cai;
                }

    for (i=0;i<size*size*height;i++)
        I_ext[i]=-2;
Error:
    printf("Error reading init file\n");
    MPI_Abort(MPI_COMM_WORLD,-1);



}
