#ifndef LR_LATTICE_H
#define LR_LATTICE_H
#include "parallel.h"
extern const int Size; //Lattice dimensions
extern const int Height; //Lattice dimensions


struct CellVariables //Phase space variables of one cell
{
	double m;
	double h;
	double j;
	double d;
	double f;
	double X;
	double Cai;
};

extern double *V;

extern CellVariables *Cell; 
extern double *I_ext;//External current
extern const double dt;//Integrating step (variable)
extern double D;//coupling parameter 
extern int DrawNum;
double Coupling(int &I, int &J, int &K); // This function returns the value of coupling D*(...) of the cell [i,j]
extern void LattInit(const int &size, const int &height);//initializing function
void OdeSolve(int &ii, int &jj, int &kk);//this function solves the system with NO coupling on the interval 'dt'
inline int Substeps(double &vd);//devides step length due to value of Voltage (V)
extern void SolveEquations(double MaxTime /*time of calculations*/, MPI_Comm &GridComm, int &GridSize, MPI_Datatype &plane_x, MPI_Datatype &plane_y, MPI_Datatype &square);//Solves the task
extern void LattInit(const int &size , const int &height, char a[]);
#endif 