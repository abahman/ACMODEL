#ifndef _NUMERIC
#define _NUMERIC

struct dataSpline
{
double Data_points[58][20];
};


double FindZero(double,double(*)(double,void*),double,void*);
int FindZeroConst(double*,double,double,double(*)(double,void*),double,void*);
double GoldSecMin(double,double,double(*)(double,void*),double,void*);
double Bisection1D(double,double,double(*)(double,void*),double,void*);
double FalsePos(double,double,double(*)(double,void*),double,void*);
double Brent(double,double,double,double(*)(double,void*),double,void*);
double Zbrent(double,double,double(*)(double,void*),double,void*);
int ZbrentPlot(double,double,double(*)(double,void*),double,void *);
void spline(double *,double *,int,double,double,double *);
void splint(double *,double *,double *,int,double,double *,double *,int);
void outdataSpline(int, int, dataSpline,double,double *);
double Bisection1DSpline(double,double,double(*)(double,void*),double,void *);

#endif
