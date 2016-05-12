#ifndef _RefrigerantFunc
#define _RefrigerantFunc

#include "defs.h"

// Structures.
struct HP {
	double H,P;
};
struct TXP {
	double T,X,P;
};

enum {
	TSAT,
	ENTH,
	VOL,
	ENTR,
	INTENERGY,
	DENSITY
};
double PropertyTXPth(int,TXP);
enum {
	COND,
	VISC,
	SPEC,
	TENSION
};
double PropertyTXPtr(int,TXP);

// Functions.
TXP toTXP(double T,double X,double P);
HP toHP(double H, double P);
TXP HPtoTXP(HP);
HP TXPtoHP(TXP);

int Exercise_Prop();//B.S.

#endif

