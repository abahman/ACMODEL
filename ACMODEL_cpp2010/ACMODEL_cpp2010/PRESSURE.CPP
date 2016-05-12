#include "stdafx.h"

#include "r22.h"
#include "volume.h"
#include "pressure.h"
#include "corr.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "errorlog.h"
#include "spline1d.h"
#include "spline2d.h"

extern ErrorLog errorLog;

/* Local pressure drop calculations */
double dPdLfricSP(TXP,double,double);
double dPdLfricTP(TXP,double,double);
double PhiSqfHPSIM(double);
double PhiSqfWALLIS(double);
double dPdLfricEvapTP(TXP,double,double,double,double);
double FricFactor(TXP,double,double);
double FricFactorSP(TXP,double,double);
double FricFactorTP(TXP,double,double);

/********************************************************************
Calculates friction pressure drop for internal refrigerant flow.
********************************************************************/
double
dPdLfriction(TXP TXPi,double G,double D)
{
	TXP TXP1,TXP2;
	double y1,y2,y;
	const double X1=0.1,X2=0.9;

	if(TXPi.X>=1.0 || TXPi.X<=0) {
		y=dPdLfricSP(TXPi,G,D);
	} else if(TXPi.X<X1) {
		TXP1=toTXP(TXPi.T,0,TXPi.P);
		TXP2=toTXP(TXPi.T,X1,TXPi.P);
		y1=dPdLfricSP(TXP1,G,D);
		if(errorLog.IsError()) {
			errorLog.Add("dPdLfriction");
			return -1;
		}
		y2=dPdLfricTP(TXP2,G,D);
		y=y1+TXPi.X*(y2-y1)/X1;
	} else if(TXPi.X>X2) {
		TXP1=toTXP(TXPi.T,X2,TXPi.P);
		TXP2=toTXP(TXPi.T,1,TXPi.P);
		y1=dPdLfricTP(TXP1,G,D);
		if(errorLog.IsError()) {
			errorLog.Add("dPdLfriction");
			return -1;
		}
		y2=dPdLfricSP(TXP2,G,D);
		y=y2-(1-TXPi.X)*(y2-y1)/(1-X2);
	} else {
		y=dPdLfricTP(TXPi,G,D);
	}
	if(errorLog.IsError()) {
		errorLog.Add("dPdLfriction");
		return -1;
	}
	return y;
}

/********************************************************************
Calculates friction pressure drop for single phase internal flow.
********************************************************************/
double
dPdLfricSP(TXP TXP,double G,double D)
{
	const double mu = PropertyTXPtr(VISC,TXP);
	if(errorLog.IsError()) {
		errorLog.Add("dPdLfricSP","mu");
		return -1;
	}
	const double v = PropertyTXPth(VOL,TXP);
	if(errorLog.IsError()) {
		errorLog.Add("dPdLfricSP","v");
		return -1;
	}
	const double Re = G*D/mu;
	if(Re<0.0) {
		errorLog.Add("dPdLfricSP","Re<0.0");
		return -1;
	}
	const double f = 0.046*pow(Re,-0.2);
	if(G<0.0) {
		errorLog.Add("dPdLfricSP","G<0.0");
		return -1;
	}
	const double dPdL = -0.5e-3*f/D*v*pow(G,2);
	return dPdL;
}

/********************************************************************
Calculates friction pressure drop for two-phase internal flow.
********************************************************************/
double
dPdLfricTP(TXP TXP,double G,double D)
{
	const double dPdLf0 = dPdLfricSP(toTXP(TXP.T,0,TXP.P),G,D);
	if(errorLog.IsError()) {
		errorLog.Add("dPdLfricTP","dPdLf0");
		return -1;
	}
	const double xtt = Xtt(TXP);
	if(errorLog.IsError()) {
		errorLog.Add("dPdLfricTP","xtt");
		return -1;
	}
	const double arg = 1.0-TXP.X;
	if(arg<0.0) {
		errorLog.Add("dPdLfricTP","1.0-TXP.X<0.0");
		return -1;
	}
	const double PhiSqf0 = pow(arg,1.8)*PhiSqfWALLIS(xtt);
	//const double PhiSqf0 = pow(arg,1.8)*PhiSqfHPSIM(xtt);

	const double dPdL = PhiSqf0*dPdLf0;

	return dPdL;
}

/********************************************************************
Correlation factor for two phase frictional pressure drop.
********************************************************************/
double
PhiSqfWALLIS(double Xtt)
{
	const double n=3.5;

	double Phi = pow(1+pow(Xtt,-2/n),n);

	return Phi;
}

/********************************************************************
Calculates the pressure drop (kPa) due to a change in specific volume
usually resulting from a heat transfer process.  Assumptions include
no body forces, no viscous forces, steady flow, and constant area.
********************************************************************/
double
dPmom(double v1,double v2,double G)
{
	return 1e-3*pow(G,2)*(v2-v1);
}

/********************************************************************
Calculates the pressure drop for 180 return bend.
Assumptions: (1) smooth wall, (2) incompressible flow.
Reference:
	Idelchik, I.E.
	Handbook of Hydraulic Resistance, 2nd ed.
	Hemisphere Publishing Company, 1986.
	p.289
********************************************************************/
double
dPelbow(TXP TXP,double G,double D,double Ro)
{
	const double delta=180 ;

	double dPdL = dPdLfriction(TXP,G,D);
	if(errorLog.IsError()) {
		errorLog.Add("dPelbow");
		return -1;
	}
	double v = VolumeALL(TXP,G,D,0);
	if(errorLog.IsError()) {
		errorLog.Add("dPelbow");
		return -1;
	}
	double A1 = 0.7+0.35*delta/90;
	double B1 = 0.21/sqrt(Ro/D);
	double C1 = 1;
	double dP = 0.0175*Ro*delta*dPdL-(0.5e-3)*pow(G,2)*v*A1*B1*C1;

	return dP;
}

/********************************************************************
Calculates frictional pressure drop for internal refrigerant flow
in the evaporator.
********************************************************************/
double
dPdLfrictionEvap(TXP TXPi,double G,double D,double dx,double L)
{
	double y;
	const double X1=0.1,X2=0.9;

	if(TXPi.X>=1.0 || TXPi.X<=0) {
		y = dPdLfricSP(TXPi,G,D);
	} else if(TXPi.X<X1) {
		TXP TXP1=toTXP(TXPi.T,0,TXPi.P);
		TXP TXP2=toTXP(TXPi.T,X1,TXPi.P);
		double y1 = dPdLfricSP(TXP1,G,D);
		if(errorLog.IsError()) {
			errorLog.Add("dPdLfrictionEvap");
			return -1;
		}
		double y2 = dPdLfricEvapTP(TXP2,G,D,dx,L);
		y=y1+TXPi.X*(y2-y1)/X1;
	} else if(TXPi.X>X2) {
		TXP TXP1=toTXP(TXPi.T,X2,TXPi.P);
		TXP TXP2=toTXP(TXPi.T,1,TXPi.P);
		double y1 = dPdLfricEvapTP(TXP1,G,D,dx,L);
		if(errorLog.IsError()) {
			errorLog.Add("dPdLfrictionEvap");
			return -1;
		}
		double y2 = dPdLfricSP(TXP2,G,D);
		y = y2-(1-TXPi.X)*(y2-y1)/(1-X2);
	} else {
		y = dPdLfricEvapTP(TXPi,G,D,dx,L);
	}
	if(errorLog.IsError()) {
		errorLog.Add("dPdLfrictionEvap");
		return -1;
	}
	return y ;
}

/********************************************************************
Calculates friction pressure drop for evaporation two phase
internal flow.
********************************************************************/
double
dPdLfricEvapTP(TXP TXPm,double G,double D,double dx,double L)
{
	const double g=9.807,dx0=1e-2;
	TXP TXP_prop={0,0,0};

	TXP_prop.P=TXPm.P;
	TXP_prop.X=0;
	const double Tsat = PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(TXPm.P);
	if(errorLog.IsError()) {
		errorLog.Add("dPdLfricEvapTP","Tsat");
		return -1;
	}

	TXP_prop.P=TXPm.P;
	TXP_prop.T=Tsat;
	TXP_prop.X=0;
	const double mul = PropertyTXPtr(VISC,TXP_prop);//refsctrPT.mu(TXPm.P,Tsat);
	if(errorLog.IsError()) {
		errorLog.Add("dPdLfricEvapTP","mul");
		return -1;
	}
	const double v = PropertyTXPth(VOL,TXPm);
	if(errorLog.IsError()) {
		errorLog.Add("dPdLfricEvapTP","v");
		return -1;
	}

	TXP_prop.P=TXPm.P;
	TXP_prop.X=0;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);
	const double hl = PropertyTXPth(ENTH,TXP_prop);//reftplthP.h(TXPm.P);
	if(errorLog.IsError()) {
		errorLog.Add("dPdLfricEvapTP","hl");
		return -1;
	}

	TXP_prop.P=TXPm.P;
	TXP_prop.X=1;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

	const double hv = PropertyTXPth(ENTH,TXP_prop);//reftpvthP.h(TXPm.P);
	if(errorLog.IsError()) {
		errorLog.Add("dPdLfricEvapTP","hv");
		return -1;
	}
	const double hlv = hv-hl;
	const double Re = G*D/mul;
	const double f0 = 0.0185*pow(hlv/(Re*L*g),0.25);

	double f;
	if(dx>=dx0) {
		f = f0*pow(dx,0.25);
	} else {
		double A = f0*pow(dx0,0.25)/exp(0.25);
		double B = 1/(4*dx0);
		f = A*exp(B*dx);
	}

	double dPdL = -1e-3*f*pow(G,2)*v/D;

	return dPdL;
}

/***********************************************************************
This function provides the local and frictional losses of a return bend.
Reference:
	Idelchik, I.E.
	Handbook of Hydraulic Resistance, 2nd ed.
	Hemisphere Publishing Company, 1986.
	p.289
Inputs:
	TXP = refrigerant state
	G = refrigerant mass flux (kg/s/m^2)
	R = radius of return bend (m)
	D = inside diameter of pipe (m)
Outputs:
	(return value) return bend loss coefficient (K)
************************************************************************/
double
LossCoeffReturnBend(TXP TXP,double G,double R,double D)
{
	const double delta=180;  // angle of return bend in degrees

	// local loss coefficient
	double A1 = 0.7+0.35*delta/90;  	// delta > 100 degrees
	double B1 = 0.21/sqrt(R/D);  		// R/D > 1.0
	double C1 = 1.0;               	// circular or square x-section
	double Kloc = A1*B1*C1;

	double f = FricFactorSP(TXP,G,D);
	if(errorLog.IsError()) {
		errorLog.Add("LossCoeffReturnBend");
	}
	double Kfric = 0.0175*R/D*delta*f;

	return Kfric+Kloc;
}

/***********************************************************************
This function provides the local and frictional losses of a straight
piece of pipe.
Inputs:
	TXP = refrigerant state
	G = refrigerant mass flux (kg/s/m^2)
	L = length of pipe (m)
	D = inside diameter of pipe (m)
Outputs:
	(return value) return bend loss coefficient (K)
************************************************************************/
double
LossCoeffStraightPipe(TXP TXP,double G,double L,double D)
{
	double f = FricFactorSP(TXP,G,D);
	if(errorLog.IsError()) {
		errorLog.Add("LossCoeffStraightPipe");
	}
	double K = f*L/D;

	return K;
}


/************************************************************************
Calculates friction factor for single-phase internal refrigerant flow.
************************************************************************/
double
FricFactorSP(TXP /*TXP*/,double /*G*/,double /*D*/)
{
	// Option 1: recommended constant value
	//
	// return 0.02;

	// Option 2: smooth pipe
	//
	// char s[SL];
	// double mu = PropertyTXP(VISC,TXP,status);
	// if(*status) {sprintf(s,"FricFactorSP:\n"); bprintf(s); return -1;}
	// double Re=G*D/mu;
	// return pow(1.8*log10(Re)-1.64,-2);

	// Option 3: rough pipe
	//
	double delta=0.01;
	return pow(2*log10(3.7/delta),-2);
}


/************************************************************************
Calculates friction factor for two-phase internal refrigerant flow.
************************************************************************/
double FricFactorTP(TXP TXP,double G,double D)
{
	double fSP = FricFactorSP(toTXP(TXP.T,0,TXP.P),G,D);
	if(errorLog.IsError()) {
		errorLog.Add("FricFactorTP");
		return -1;
	}

	double xtt = Xtt(TXP);
	if(errorLog.IsError()) {
		errorLog.Add("FricFactorTP");
		return -1;
	}

	double PhiSqf0 = pow(1-TXP.X,1.8)*PhiSqfWALLIS(xtt);
	//double PhiSqf0 = pow(1-TXP.X,1.8)*PhiSqfHPSIM(xtt);

	double fTP = PhiSqf0*fSP;

	return fTP;
}


/************************************************************************
Calculates friction factor for internal refrigerant flow.
************************************************************************/
double
FricFactor(TXP TXPi,double G,double D)
{
	double y;
	const double X1=0.1,X2=0.9;

	if(TXPi.X>=1.0 || TXPi.X<=0) {
		y = FricFactorSP(TXPi,G,D);
	} else if(TXPi.X<X1) {
		TXP TXP1 = toTXP(TXPi.T,0,TXPi.P);
		TXP TXP2 = toTXP(TXPi.T,X1,TXPi.P);
		double y1 = FricFactorSP(TXP1,G,D);
		if(errorLog.IsError()) {
			errorLog.Add("FricFactor");
			return -1;
		}
		double y2 = FricFactorTP(TXP2,G,D);
		y = y1+TXPi.X*(y2-y1)/X1;
	} else if(TXPi.X>X2) {
		TXP TXP1=toTXP(TXPi.T,X2,TXPi.P);
		TXP TXP2=toTXP(TXPi.T,1,TXPi.P);
		double y1=FricFactorTP(TXP1,G,D);
		if(errorLog.IsError()) {
			errorLog.Add("FricFactor");
			return -1;
		}
		double y2=FricFactorSP(TXP2,G,D);
		y = y2-(1-TXPi.X)*(y2-y1)/(1-X2);
	} else {
		y = FricFactorTP(TXPi,G,D);
	}
	if(errorLog.IsError()) {
		errorLog.Add("FricFactor");
		return -1;
	}

	return y ;
}

/******************************
B.S. add for getting the acceleration pressure drop
  **************************/

double GET_PreAcc(double DP_ACC,
				  void *Params)

{	
	PreAcc *Preacc = (PreAcc*)Params;

	const double DP_FR=Preacc->DP_FR;
	const double G=Preacc->G;
	const double H_OUT=Preacc->H_OUT;
	const double P=Preacc->P_IN;
	const double X=Preacc->X_IN;
	TXP	TXP_prop={0,0,0};

	TXP_prop.P=P;
	TXP_prop.X=1;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

	const double DV=1.0/PropertyTXPth(VOL,TXP_prop);//1.0/reftpvthP.v(P);

	if(errorLog.IsError()) {
		errorLog.Add("PreAcc","vv");
		return 0;
	}
	
	TXP_prop.P=P;
	TXP_prop.X=0;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

	const double DL=1.0/PropertyTXPth(VOL,TXP_prop);//1.0/reftplthP.v(P);
	if(errorLog.IsError()) {
		errorLog.Add("PreAcc","vl");
		return 0;
	}
	

	const double DP_OLD=DP_ACC;
	const double P_OUT=P-DP_FR-DP_ACC;
	TXP TXP_OUT;
	HP HP_OUT;
	HP_OUT.H=H_OUT;
	HP_OUT.P=P_OUT;
	TXP_OUT=HPtoTXP(HP_OUT);

	if(errorLog.IsError()) {
		errorLog.Add("PreAcc","HPtoTXP");
		return 0;
	}
	const double X_OUT=TXP_OUT.X;
	DP_ACC=1e-3*pow(G,2e0)*(1e0/DV-1e0/DL)*(X_OUT-X);
	const double ERR_P=(DP_ACC-DP_OLD);

	return ERR_P;
	}
