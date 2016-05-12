
#include "stdafx.h"

#include <math.h>
#include <stdio.h>
#include "r22.h"
#include "ref.h"
#include "numeric.h"
#include "errorlog.h"
#include "spline1d.h"
#include "spline2d.h"


extern ErrorLog errorLog;

/* Local strucutures */
struct FlashP {
	double Gr,h0;
};
struct ChokeP {
	double Gr,h0;
};
struct FannoP {
	double Gr,h0,P;
};

/* Local functions */
double FlashFunc(double,void*);
double ChokeFunc(double,void*);
double FannoFunc(double,void*);
TXP TXPfanno(double,double,double);

/********************************************************************
Calculates pressure when subcooled liquid in fanno flow will flash
(h0=h+v^2/2).
********************************************************************/
double
Pflash(double Gr,double h0)
{
	FlashP FP;

	FP.Gr=Gr;
	FP.h0=h0;

	/* Check if Pflash>PMAXth.  If so, then return 1e20(inf). */
	double dh = FlashFunc(PMAXth,&FP);
	if(errorLog.IsError()) {
		errorLog.Add("Pflash");
		return 0;
	}
	if(dh>0) return 1e20;

	double Pflash = Zbrent(PMINth,PMAXth,FlashFunc,1e-7,&FP);
	if(errorLog.IsError()) {
		errorLog.Add("Pflash");
		return 0;
	}

	return Pflash;
}

/********************************************************************
Used by Pflash to determine flash pressure for fanno flow.
********************************************************************/
double
FlashFunc(double Pflash,void *Params)
{
	FlashP *P=(FlashP*)Params;
	TXP	TXP_prop={0,0,0};
	
	TXP_prop.P=Pflash;
	TXP_prop.X=0;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

	const double hf = PropertyTXPth(ENTH,TXP_prop);//reftplthP.h(Pflash);
	if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"hf, Pflash=%lf",Pflash);
		errorLog.Add("FlashFunc",str);
		return -1;
	}

	const double vf = PropertyTXPth(VOL,TXP_prop);//reftplthP.v(Pflash);
	if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"vf, Pflash=%lf",Pflash);
		errorLog.Add("FlashFunc",str);
		return -1;
	}
	const double arg = P->Gr*vf;
	if(arg<0.0) {
		errorLog.Add("FlashFunc","P->Gr*vf<0.0");
		return -1;
	}
	const double dh = (P->h0-hf-pow(arg,2.0)/2.0)/P->h0;
	return dh;
}

/********************************************************************
Determines the choking pressure for refrigerant in fanno flow.
********************************************************************/
double
Pchoke(double Gr,double h0)
{
	ChokeP CP;
	CP.Gr = Gr;
	CP.h0 = h0;
	double P = Brent(PMINth,(PMINth+PMAXth)/2,PMAXth,ChokeFunc,1e-8,&CP);
	if(errorLog.IsError()) {
		errorLog.Add("Pchoke");
	}

	return P;
}

/********************************************************************
Used by Pchoke to determine choking pressure for refrigerant in
fanno flow.
********************************************************************/
double
ChokeFunc(double Pc,void *Params)
{
	ChokeP *P=(ChokeP*)Params;

	TXP TXPf = TXPfanno(P->Gr,P->h0,Pc);
	if(errorLog.IsError()) {
		errorLog.Add("ChokeFunc","TXPf");
		return -1;
	}
	const double ret = -1*PropertyTXPth(ENTR,TXPf);
	if(errorLog.IsError()) {
		errorLog.Add("ChokeFunc","ret");
		return -1;
	}
	return ret;
}

/********************************************************************
Given mass flux (Gr) and stagnation enthalpy (h0), this function
calculates the complete state (TXP) at the pressure (P)
assuming the stagnation enthalpy (h0=h+V^2/2) is constant.
Conservation of stagnation enthalpy is a result of the first law
applied to any adiabatic flow with no internal generation.  This
includes fanno flow.
********************************************************************/
TXP TXPfanno(double Gr,double h0,double P)
{
	TXP	TXP_prop={0,0,0};
	TXP TXP = toTXP(0,0,0);

	FannoP FP;
	FP.Gr=Gr;
	FP.h0=h0;
	FP.P=P;

	HP HP;
	HP.P=P;

	TXP_prop.P=P;
	TXP_prop.T=TMIN;
	TXP_prop.X=0;

	const double hmin = PropertyTXPth(ENTH,TXP_prop);//refscthPT.h(P,TMIN);
	if(errorLog.IsError()) {
		errorLog.Add("TXPfanno","hmin");
		return TXP;
	}

	HP.H = Zbrent(hmin,h0,FannoFunc,1e-7,&FP);
	if(errorLog.IsError()) {
		errorLog.Add("TXPfanno","HP.H");
		return TXP;
	}

	TXP = HPtoTXP(HP);
	if(errorLog.IsError()) {
		errorLog.Add("TXPfanno","TXP");
		return TXP;
	}

	return TXP;
}

/********************************************************************
Given a guess of enthalpy (h), this function generates a residual
that conserves stagnation enthalpy when zero.  This function is
used by TXPfanno.
********************************************************************/
double
FannoFunc(double h,void *Params)
{
	FannoP *P = (FannoP*)Params;

	HP HPs;
	HPs.P = P->P;
	HPs.H = h;

	TXP TXPs = HPtoTXP(HPs);
	if(errorLog.IsError()) {
		errorLog.Add("FannoFunc");
		return 0;
	}

	double v = PropertyTXPth(VOL,TXPs);
	if(errorLog.IsError()) {
		errorLog.Add("FannoFunc");
		return 0;
	}
	
	double h0_Mid=0;// shenbo add 
	if(P->h0>=1e5) h0_Mid=P->h0;//shenbo add
	else h0_Mid=1e5;// shenbo add
/* there is a problem, when the absolute value of enthalpy is small, it is hard to converge, but the enthalpy is a relative units, 
so it does not make sense to control the exit when the enthalpy is small, so here constrain the value of enthalpy for exit, make it >=1e5, it proves to be helpful*/
	double dh = (P->h0-h-0.5*pow(P->Gr*v,2))/h0_Mid;//P->h0; shenbo add //(P->h0-h-0.5*pow(P->Gr*v,2))/P->h0; original

	return dh;
}
