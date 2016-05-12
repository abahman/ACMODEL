#ifndef _NUM2D
#define _NUM2D

double FindZero2DTest(double*,double(*)(double*,double*,void*),double,void*);
int FindZero2DGolden(double*,int(*)(double*,double*,void*),double,void*);
int FindZero2DGolden_InLoop(double*,int(*)(double*,double*,void*),double,void*);//B.S.
double FindZero2D(double*,double(*)(double*,double*,void*),double,void*);
double FindZero2DConst(double*,double(*)(double*,double*,void*),int(*)(double*,double*,void*),double,void*);
double Bisection2D(double*,double*,double*,double(*)(double*,double*,void*),double,void*);

#endif
