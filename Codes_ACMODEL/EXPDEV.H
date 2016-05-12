#ifndef _EXPDEV
#define _EXPDEV

#include "r22.h"
#include "mass.h"

void ExpDev(const char *filename,int,double,double,double,double*,MASS*,double*);
void ExpDev0(int,double,const char *filename,int,double,double,double,double*,MASS*,double*);

double Short_Tube_FlowDirec(const char *filename,HP,HP,double*);
double Short_Tube_0(const char *filename,HP,HP,double,double*);
double Short_Tube_R410A(const char *filename,HP,HP,double*);
double Short_Tube_R410A_Oil(const char *filename,HP,HP,double*);
double Short_Tube_R407C(const char *filename,HP,HP,double*);
double Short_Tube_R407C_Oil(const char *filename,HP,HP,double*);
double Short_Tube_R22(const char *filename,HP,HP,double*);
double  TxV_Nonlinear(const char *filename,HP,HP,double,double*);
double	TxV_0(const char *filename,HP,HP,double, double,double*);
#endif
