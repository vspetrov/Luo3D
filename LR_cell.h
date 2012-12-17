#ifndef LR_cell_H
#define LR_cell_H

//right-hand parts of the Luo-Rudy System
extern double VFunction(int &ii /*column cell number*/,int &jj /*Row Cell number*/, int &kk);
extern double mFunction(double &delta_t);
extern double hFunction(double &delta_t);
extern double jFunction(double &delta_t);
extern double dFunction(double &delta_t);
extern double fFunction(double &delta_t);
extern double XFunction(double &delta_t);
extern double CaiFunction();

#endif