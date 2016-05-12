#ifndef _NUM3D
#define _NUM3D

int FindZero3DGolden(double X[],int(*func)(double*,double*,void*),double tol,void *P);
int FindZero3DGuess(double X[],int(*func)(double*,double*,void*),double tol,int fix,void *P);

#endif
