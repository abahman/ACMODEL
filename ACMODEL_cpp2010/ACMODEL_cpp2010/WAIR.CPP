#include "stdafx.h"

#include <stdio.h>
#include "wair.h"
#include "numeric.h"
#include "num2d.h"
#include "errorlog.h"
#include "spline2d.h"

#include <math.h>//B.S.

extern ErrorLog errorLog;

extern WetAirPropTR wair;

#define TMIN -20
#define TMAX 50

/* Local functions */
double WHFunc(double*,double*,void*);
int WHFuncConst(double*,double*,void*);
double WPFunc(double,void*);
double HPFunc(double,void*);
double TWFunc(double,void*);
double TDFunc(double,void*);
double TZFunc(double,void*);
double THFunc(double,void*);

/* Local structures */
struct HPFuncP {
	double H,P;
};
struct WPFuncP {
	double W,P;
};
struct TWFuncP {
	double T,W;
};
struct TDFuncP {
	double T,D;
};
struct TZFuncP {
	double T,Z;
};
struct THFuncP {
	double T,H;
};

/********************************************************************
Converts format of thermodynamic state representation for moist
air from humidity ratio (W) and relative humidity (P) to the
standard format for ACMODEL: temperature (T) and relative
humidity (P).
********************************************************************/
TP
WPtoTP(double W,double P)
{
	TP TP;
	WPFuncP WP;

	WP.W=W; WP.P=P;
	TP.T=Zbrent(TMIN,TMAX,WPFunc,1e-7,&WP);
	if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"W=%lf P=%lf\n",W,P);
		errorLog.Add("WPtoTP",str);
	}
	TP.P=P;

	return TP;
}

/********************************************************************
Used by {WPtoTP} to convert thermodynamic state representations.
********************************************************************/
double
WPFunc(double T,void *Params)
{
	WPFuncP *P = (WPFuncP*)Params;

	double W = wair.HumidityRatio(T,P->P);
	if(errorLog.IsError()) {
		errorLog.Add("WPFunc");
	}

	return (W-P->W)/P->W;
}

/********************************************************************
Converts format of thermodynamic state representation for moist
air from enthalpy (H) and relative humidity (P) to the
standard format for ACMODEL: temperature (T) and relative
humidity (P).
********************************************************************/
TP
HPtoTP(double H,double P)
{
	TP TP;
	HPFuncP HP;

	HP.H = H;
	HP.P = P;
	TP.T = Zbrent(TMIN,TMAX,HPFunc,1e-7,&HP);
	if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"HP=(%lf,%lf)",H,P);
		errorLog.Add("HPtoTP",str);
	}
	TP.P = P;

	return TP;
}

/********************************************************************
Used by {\it HPtoTP} to convert thermodynamic state representations.
********************************************************************/
double
HPFunc(double T,void *Params)
{
	HPFuncP *P = (HPFuncP*)Params;
	const double H = wair.h(T,P->P);
	if(errorLog.IsError()) {
		errorLog.Add("HPFunc");
	}

	double Ref=P->H;//B.S.
	if(fabs(Ref)<1e4) Ref=1e4;//B.S.
	return (H-P->H)/Ref;//B.S.
}

/********************************************************************
Converts format of thermodynamic state representation for moist
air from temperature (T) and enthalpy (H) to the standard format
for ACMODEL: temperature (T) and relative humidity (P).
********************************************************************/
TP
THtoTP(double T,double H)
{
	TP TP;
	THFuncP TH;

	TH.T=T;
	TH.H=H;
	TP.P = Zbrent(0.05,1.0,THFunc,1e-7,&TH);
	if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"T=%lf H=%lf\n",T,H);
		errorLog.Add("THtoTP",str);
	}
	TP.T=T;

	return TP;
}

/********************************************************************
Used by {\it THtoTP} to convert thermodynamic state representations.
********************************************************************/
double
THFunc(double P,void *Params)
{
	THFuncP *Q =(THFuncP*)Params;

	double H = wair.h(Q->T,P);

	if(errorLog.IsError()) {
		errorLog.Add("THFunc");
	}

	double Ref=Q->H;//B.S.
	if(fabs(Ref)<1e4) Ref=1e4;//B.S.
	return (H-Q->H)/Ref;//B.S.
}

/********************************************************************
Converts format of thermodynamic state representation for moist
air from temperature (T) and wet bulb (Z) to the standard format
for ACMODEL: temperature (T) and relative humidity (P).
********************************************************************/
TP
TZtoTP(double T,double Z)
{
	TP TP;
	TZFuncP TZ;

	TZ.T=T;
	TZ.Z=Z;
	TP.P = Zbrent(0.05,1.0,TZFunc,1e-7,&TZ);
	if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"T=%lf Z=%lf\n",T,Z);
		errorLog.Add("TZtoTP",str);
	}
	TP.T=T;

	return TP;
}

/********************************************************************
Used by {\it TZtoTP} to convert thermodynamic state representations.
********************************************************************/
double
TZFunc(double P,void *Params)
{
	TZFuncP *Q = (TZFuncP*)Params;

	double Z = wair.WetBulb(Q->T,P);

	if(errorLog.IsError()) {
		errorLog.Add("TZFunc");
	}

	double Ref=Q->Z;//B.S.
	if(fabs(Ref)<1.0) Ref=1.0;//B.S.
	return (Z-Q->Z)/Ref;
}

/********************************************************************
Converts format of thermodynamic state representation for moist
air from temperature (T) and dew point (D) to the standard format
for ACMODEL: temperature (T) and relative humidity (P).
********************************************************************/
TP
TDtoTP(double T,double D)
{
	TP TP;
	TDFuncP TD;

	TD.T=T;
	TD.D=D;
	TP.P = Zbrent(0.05,1.0,TDFunc,1e-7,&TD);
	if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"T=%lf D=%lf\n",T,D);
		errorLog.Add("TDtoTP",str);
	}
	TP.T=T ;

	return TP ;
}

/********************************************************************
Used by {\it TDtoTP} to convert thermodynamic state representations.
********************************************************************/
double
TDFunc(double P,void *Params)
{
	TDFuncP *Q = (TDFuncP*)Params;

	double D = wair.DewPoint(Q->T,P);

	if(errorLog.IsError()) {
		errorLog.Add("TDFunc");
	}

	double Ref=Q->D;//B.S.
	if(fabs(Ref)<1.0) Ref=1.0;//B.S.

	return (D-Q->D)/Ref;
}

/********************************************************************
Converts format of thermodynamic state representation for moist
air from temperature (T) and humidity ratio (W) to the standard format
for ACMODEL: temperature (T) and relative humidity (P).
********************************************************************/
TP
TWtoTP(double T,double W)
{
	TP TP;
	TWFuncP TW;

	TW.T = T;
	TW.W = W;
	TP.P = Zbrent(0.05,1.0,TWFunc,1e-7,&TW);
	if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"T=%lf W=%lf\n",T,W);
		errorLog.Add("TWtoTP",str);
	}
	TP.T=T;

	return TP;
}

/********************************************************************
Used by {\it TWtoTP} to convert thermodynamic state representations.
********************************************************************/
double
TWFunc(double P,void *Params)
{
	TWFuncP *Q =(TWFuncP*)Params;

	double W = wair.HumidityRatio(Q->T,P);

	if(errorLog.IsError()) {
		errorLog.Add("TWFunc");
	}
	double Ref=Q->W;//B.S.
//	if(Ref<1e-3) Ref=1e-3;//B.S.
	return (W-Q->W)/Ref ;//B.S.
}

/********************************************************************
Converts format of thermodynamic state representation for moist
air from humidity ratio (W) and enthalpy (H) to the standard format
for ACMODEL: temperature (T) and relative humidity (P).
********************************************************************/
TP WHtoTP(WH WH,TP TPi)
{
	TP TPo;
	double X[2];

	X[0]=TPi.P;
	X[1]=TPi.T;
	double CE = FindZero2DConst(X,WHFunc,WHFuncConst,1e-6,&WH);

	TPo.T=X[1];
	TPo.P=X[0];
	if(errorLog.IsError() && TPo.P>0.995) {
		char str[128];
		sprintf(str,"TPo=(%lf,%lf)",TPo.T,TPo.P);
		errorLog.ClearError("WHtoTP",str);
		TPo = HPtoTP(WH.H,0.995);
	}
	if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"Can't solve to tolerance. CE=%lf",CE);
		errorLog.Add("WHtoTP",str);
	}

	return TPo;
}

/********************************************************************
Used by {WHtoTP} to convert thermodynamic state representations.
********************************************************************/
double
WHFunc(double *X,double *Y,void *Params)
{
	WH *P = (WH*)Params;

	double W = wair.HumidityRatio(X[1],X[0]);
	if(errorLog.IsError()) {
		errorLog.Add("WHFunc");
	}
	double H = wair.h(X[1],X[0]);
	if(errorLog.IsError()) {
		errorLog.Add("WHFunc");
	}

	Y[0]=(W-P->W)/P->W;
	double Ref=P->H;//B.S.
	if(fabs(Ref)<1e4) Ref=1e4;//B.S.
	Y[1]=(H-P->H)/Ref;//B.S.

	return 0;
}

/********************************************************************
Tests constraint on relative humidity when solving WHFunc.
Returns 1 if constraint is violated.
********************************************************************/
int
WHFuncConst(double *X,double* /*Y*/,void* /*Params*/)
{
	if(errorLog.IsError()) {
		errorLog.Add("WHFuncConst");
	}
	if(X[0]>1||X[0]<0.05||X[1]>TWAIRMAX||X[1]<TWAIRMIN)
		return 1;
	else
		return 0;
}
