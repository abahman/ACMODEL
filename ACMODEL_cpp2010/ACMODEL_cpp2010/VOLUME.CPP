#include "stdafx.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "r22.h"
#include "volume.h"
#include "numeric.h"
#include "errorlog.h"
#include "spline1d.h"
#include "spline2d.h"

extern ErrorLog errorLog;

double VolumeTP_Baroc(TXP,double,double);
double VolumeTP_Zivi(TXP,double,double);
double VolumeTP_Hugh(TXP,double,double);
double VolumeTP_ACRC(TXP,double,double);
double VolumeTP_LM(TXP,double,double);
double VolumeTP_Rigot(TXP,double,double);
double VolumeTP_Smith(TXP,double,double);
double VolumeTP_Tandon(TXP,double,double);
double VolumeTP_Thom(TXP,double,double);
double VolumeTP_Premoli(TXP,double,double);
double VolumeTP_Homo(TXP,double,double);
double VolumeTP_Rouhani(TXP,double,double);
double VolumeTP_Thome(TXP,double,double);
double VolumeTP_Taitel(TXP,double,double);
double VolumeTP_FlowPattern(TXP,double,double);

/********************************************************************
Martinelli parameter. Sqrt of the ratio of the frictional
pressure drop of liquid to vapor if they were flowing alone
in a smooth pipe.
********************************************************************/
double
Xtt(TXP TXPm)
{
	TXP TXP_prop={0,0,0};
	
	TXP_prop.P=TXPm.P;
	TXP_prop.X=0;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

	const double vl = PropertyTXPth(VOL,TXP_prop);//reftplthP.v(TXPm.P);
	if(errorLog.IsError()) {
		errorLog.Add("Xtt","vl");
		return -1;
	}

	TXP_prop.P=TXPm.P;
	TXP_prop.X=1;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

	const double vv  = PropertyTXPth(VOL,TXP_prop);//reftpvthP.v(TXPm.P);
	if(errorLog.IsError()) {
		errorLog.Add("Xtt","vv");
		return -1;
	}

	TXP_prop.P=TXPm.P;
	TXP_prop.X=0;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

	const double mul = PropertyTXPtr(VISC,TXP_prop);//refsctrPT.mu(TXPm.P,Tsat);
	if(errorLog.IsError()) {
		errorLog.Add("Xtt","mul");
		return -1;
	}
	
	TXP_prop.P=TXPm.P;
	TXP_prop.X=1;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

	const double muv = PropertyTXPtr(VISC,TXP_prop);//refshtrPT.mu(TXPm.P,Tsat);
	if(errorLog.IsError()) {
		errorLog.Add("Xtt","muv");
		return -1;
	}
	const double Xtt = pow((1-TXPm.X)/TXPm.X,0.9)*pow(vl/vv,0.5)*pow(mul/muv,0.1);

	return Xtt;
}

/********************************************************************
Solves for the void fraction (Alpha) in the Hughmark model.
It takes a guess of Alpha as an argument and returns a
residual on Alpha that must be zero.
********************************************************************/
double
VolumeTPFunc(double Alpha,void *Params)
{
	VolParams *P = (VolParams*)Params;

	// Z is the "correlating parameter"
	double Z1 = pow((P->D*P->G)/(P->ul+Alpha*(P->uv-P->ul)),0.16666667);
	double Z2 = pow((P->G*P->x)/(P->Beta*(1-P->Beta)/P->vv),2);
	double Z  = Z1 * pow(Z2/(9.81*P->D),0.125);

	// K is the flow parameter
	double K;
	if(Z<=10) {
		K = -0.1636+Z*(0.3103+Z*(-0.0352+Z*1.365e-3));
	} else {
		K = 0.7257+Z*(6.122e-3+Z*(-6.26e-5+Z*2.396e-7));
	}

	double Alpha2 = P->Beta*K;

	return 2*(Alpha-Alpha2)/(Alpha+Alpha2);
}

/********************************************************************
Computes the specific volume of refrigerant for single
and two phase.  It seperates the phase regimes and uses the
appropriate correlation to get the solution.
	TXP is the state of the refrigerant
	returns the specific volume.
	*status is non-zero if an error occurs.
********************************************************************/
double
VolumeALL(TXP TXPi,double G,double D,double q)
{
	double v;
	const double X1=0.1,X2=0.9;//0.1--original

	if(TXPi.X>0.999 || TXPi.X<0.001) {
		v = PropertyTXPth(VOL,TXPi);
	} else if(TXPi.X<X1) {
		TXP TXP1 = toTXP(TXPi.T,0,TXPi.P);
		double v1 = PropertyTXPth(VOL,TXP1);
		if(errorLog.IsError()) {
			char str[128];
			sprintf(str,"v1,TXPi.X<X1 TXP1(%lf,%lf,%lf)\n",TXP1.T,TXP1.X,TXP1.P);
			errorLog.Add("VolumeALL",str);
			return -1;
		}
		TXP TXP2 = toTXP(TXPi.T,X1,TXPi.P);
		double v2 = VolumeTP(TXP2,G,D,q);
		if(errorLog.IsError()) {
			char str[128];
			sprintf(str,"v2,TXPi.X<X1 TXP2(%lf,%lf,%lf)\n",TXP2.T,TXP2.X,TXP2.P);
			errorLog.Add("VolumeALL",str);
			return -1;
		}
		v = v1 + TXPi.X*(v2-v1)/X1;
	} else if(TXPi.X>X2) {
		TXP TXP1 = toTXP(TXPi.T,X2,TXPi.P);
		double v1 = VolumeTP(TXP1,G,D,q);
		if(errorLog.IsError()) {
			char str[128];
			sprintf(str,"v1,TXPi.X>X2 TXP1(%lf,%lf,%lf)\n",TXP1.T,TXP1.X,TXP1.P);
			errorLog.Add("VolumeALL",str);
			return -1;
		}
		TXP TXP2 = toTXP(TXPi.T,1,TXPi.P);
		double v2 = PropertyTXPth(VOL,TXP2);
		if(errorLog.IsError()) {
			char str[128];
			sprintf(str,"v2,TXPi.X>X2 TXP2(%lf,%lf,%lf)\n",TXP2.T,TXP2.X,TXP2.P);
			errorLog.Add("VolumeALL",str);
			return -1;
		}
		v = v2-(1-TXPi.X)*(v2-v1)/(1-X2);
	} else {
		v = VolumeTP(TXPi,G,D,q);
		if(errorLog.IsError()) {
			char str[128];
			sprintf(str,"v TXPi(%lf,%lf,%lf)\n",TXPi.T,TXPi.X,TXPi.P);
			errorLog.Add("VolumeALL",str);
			return -1;
		}
	}
	return v ;
}

/********************************************************************
Specific volume of two phase refrigerant.  Void fraction for two
phase flow.  This corrlation implements Hugemark model.  ASHRAE
Transactions ???, pp.309-316.
	TXP (temperature(C), quality(NU), pressure(kPa)) = thermodynamic state.
	G refrigerant mass flux (kg/m^2/s)
	D internal tube diameter.
	*status is non-zero if an error occured in this function.
	returns the specific volume in m^3/kg.
********************************************************************/
double
VolumeTP(TXP TXP,double G,double D,double q)
{
	double V=0;
	switch(1) {
			case 1: 
				V=VolumeTP_Baroc(TXP,G,D);
				break;
			case 2:
				V=VolumeTP_Zivi(TXP,G,D);
				break;
			case 3:
				V=VolumeTP_Hugh(TXP,G,D);
				break;
			case 4:
				V=VolumeTP_ACRC(TXP,G,D);
				break;
			case 5:
				V=VolumeTP_LM(TXP,G,D);
				break;
			case 6:
				V=VolumeTP_Rigot(TXP,G,D);
				break;
			case 7:
				V=VolumeTP_Smith(TXP,G,D);
				break;
			case 8:
				V=VolumeTP_Tandon(TXP,G,D);
				break;
			case 9:
				V=VolumeTP_Thom(TXP,G,D);
				break;
			case 10:
				V=VolumeTP_Premoli(TXP,G,D);
				break;
			case 11:
				V=VolumeTP_Homo(TXP,G,D);
				break;
			case 12:
			V=VolumeTP_Taitel(TXP,G,D);
				break;
			case 13:
				V=VolumeTP_Rouhani(TXP,G,D);
				break;
			case 14:
				if(fabs(q)<0.00000001) 
				{V=VolumeTP_Baroc(TXP,G,D);
				break;
				}
				V=VolumeTP_Thome(TXP,G,D,q);
				break;
			case 15:
				if(fabs(q)<0.00000001) 
				{
				V=VolumeTP_Baroc(TXP,G,D);
				break;
				}
				V=VolumeTP_FlowPattern(TXP,G,D,q);
				break;
			default: 
				V=VolumeTP_Baroc(TXP,G,D);
				break;
		};

	
	return V;
}

//B.S.------------------------------------------------
//Baroczy

double
VolumeTP_Baroc(TXP TXP1,double G,double D)
{
	VolParams P;
	double Alpha =0;
	TXP TXP_prop={0,0,0};

	// liquid refrigerant properties
	TXP_prop.P=TXP1.P;
	TXP_prop.X=0;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);
	P.vl = PropertyTXPth(VOL,TXP_prop);//reftplthP.v(TXP1.P);
	if(errorLog.IsError()) {
		errorLog.Add("VolumeTP","P.vl");
		return -1;
	}


	P.ul = PropertyTXPtr(VISC,TXP_prop);//refsctrPT.mu(TXP1.P,Tsat);
	if(errorLog.IsError()) {
		errorLog.Add("VolumeTP","P.ul");
		return -1;
	}
	//vapor refrigerant properties

	TXP_prop.P=TXP1.P;
	TXP_prop.X=1;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

	P.vv = PropertyTXPth(VOL,TXP_prop);//reftpvthP.v(TXP1.P);
	if(errorLog.IsError()) {
		errorLog.Add("VolumeTP","P.vv");
		return -1;
	}


	P.uv = PropertyTXPtr(VISC,TXP_prop);//refshtrPT.mu(TXP1.P,Tsat);
	if(errorLog.IsError()) {
		errorLog.Add("VolumeTP","P.uv");
		return -1;
	}

	// other parameters that need to be passed into the structure
	P.x=TXP1.X; P.G=G; P.D=D;

	//Baroczy model
	{int i=0, j=0;
	double PI2=0;
	double X_tt=0;

	double a[7][7]={{0.49938,-3.9952e-2,1.0212e-2,2.8821e-3,3.6251e-4,2.41e-5,6.6334e-7},
		{-2.5137e-1,-1.5633e-2,7.125e-3,1.7547e-3,1.4822e-4,4.3686e-6,0.0},
		{6.0946e-3,6.5706e-3,1.5014e-4,-6.1361e-5,-1.2632e-6,1.6624e-7,0.0},
		{1.4514e-2,2.4809e-3,-6.9250e-4,-2.2043e-4,-2.1028e-5,-7.3585e-7,0.0},
		{-6.3527e-4,-1.3860e-4,5.3116e-5,2.3740e-5,3.2202e-6,1.5895e-7,0.0},
		{-3.6424e-4,-9.3832e-5,6.0560e-6,2.4383e-6,4.7549e-8,-8.1633e-9,0.0},
		{1.6539e-5,0.0,0.0,0.0,0.0,0.0,0.0}};

	PI2=P.vl/P.vv*pow(P.ul/P.uv,0.2);
	X_tt=Xtt(TXP1);
	
	Alpha=0;
	for(int j=0;j<7;j++)
	{
	for(int i=0;i<7;i++)
	{Alpha=Alpha+a[i][j]*pow(log(X_tt),double(i))*pow(log(PI2),double(j));};

	};
	}

	if(Alpha>=1) Alpha=0.9999999999;
	if(Alpha<=0) Alpha=0.0000000001;
	return 1/(Alpha/P.vv+(1-Alpha)/P.vl);
}

//B.S.-----------------------------------------------------------
//Zivi
double
VolumeTP_Zivi(TXP TXP1,double G,double D)
{
	VolParams P;
	
	double Alpha=0;
	double PI1=0;
	double S=0;
	TXP TXP_prop={0,0,0};

	// liquid and vapor specific volumes
	TXP_prop.P=TXP1.P;
	TXP_prop.X=0;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);
	P.vl = PropertyTXPth(VOL,TXP_prop);//reftplthP.v(TXP1.P);
	if(errorLog.IsError()) {
		errorLog.Add("VolumeTP","P.vl");
		return -1;
	}

	TXP_prop.P=TXP1.P;
	TXP_prop.X=1.0;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

	P.vv = PropertyTXPth(VOL,TXP_prop);//reftpvthP.v(TXP1.P);
	if(errorLog.IsError()) {
		errorLog.Add("VolumeTP","P.vv");
		return -1;
	}
	PI1=P.vl/P.vv;
	S=pow(PI1,-0.3333333);
	Alpha=1/(1+(1-TXP1.X)/TXP1.X*P.vl/P.vv*S);
	
	if(Alpha<=0) Alpha=0.000000000001;
	if(Alpha>=1) Alpha=0.999999999999;

	return 1/(Alpha/P.vv+(1-Alpha)/P.vl);
}

//B.S.----------------------------------------------------
//Hughmark
double
VolumeTP_Hugh(TXP TXP1,double G,double D)
{
	VolParams P;
	TXP TXP_prop={0,0,0};

	// liquid and vapor specific volumes
	TXP_prop.P=TXP1.P;
	TXP_prop.X=0;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

	P.vl = PropertyTXPth(VOL,TXP_prop);//reftplthP.v(TXP1.P);
	if(errorLog.IsError()) {
		errorLog.Add("VolumeTP","P.vl");
		return -1;
	}

	TXP_prop.P=TXP1.P;
	TXP_prop.X=1.0;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

	P.vv = PropertyTXPth(VOL,TXP_prop);//reftpvthP.v(TXP1.P);
	if(errorLog.IsError()) {
		errorLog.Add("VolumeTP","P.vv");
		return -1;
	}


	// liquid and vapor viscosity
	TXP_prop.P=TXP1.P;
	TXP_prop.X=0;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

	P.ul =PropertyTXPtr(VISC,TXP_prop);//refsctrPT.mu(TXP1.P,Tsat); ;
	if(errorLog.IsError()) {errorLog.Add("VolumeTP","P.vv");
		return -1;}

	TXP_prop.P=TXP1.P;
	TXP_prop.X=1;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);
	P.uv =PropertyTXPtr(VISC,TXP_prop);//refshtrPT.mu(TXP1.P,Tsat);
	if(errorLog.IsError()) {errorLog.Add("VolumeTP","P.vv");
		return -1;}

	// other parameters that need to be passed into the structure
	P.x=TXP1.X; P.G=G; P.D=D;

	// Beta is the prediction of the homogeneous flow model.
	P.Beta = (P.x*P.vv)/(P.x*P.vv+(1-P.x)*P.vl);

	// Alpha id the void fraction.
	double Alpha =Zbrent(0,1,VolumeTPFunc,1e-7,&P);
	if(errorLog.IsError()) {errorLog.Add("VolumeTP","P.vv");
		return -1;}
	
	if(Alpha<=0) Alpha=0.001;
	if(Alpha>=1) Alpha=0.999;
	return 1/(Alpha/P.vv+(1-Alpha)/P.vl);	
}

//B.S.--------------------------------------------------
//ACRC
double
VolumeTP_ACRC(TXP TXP1,double G,double D)
{
	VolParams P;
	double Fr_ACRC=0;
	double j_g=0;
	double n=0;
	double Alpha=0;
	double X_tt=0;
	TXP TXP_prop={0,0,0};


	// liquid and vapor specific volumes
	TXP_prop.P=TXP1.P;
	TXP_prop.X=0;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

	P.vl = PropertyTXPth(VOL,TXP_prop);//reftplthP.v(TXP1.P);
	if(errorLog.IsError()) {
		errorLog.Add("VolumeTP","P.vl");
		return -1;
	}

	TXP_prop.P=TXP1.P;
	TXP_prop.X=1.0;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

	P.vv = PropertyTXPth(VOL,TXP_prop);//reftpvthP.v(TXP1.P);
	if(errorLog.IsError()) {
		errorLog.Add("VolumeTP","P.vv");
		return -1;
	}

	j_g=TXP1.X*G*P.vv;

	Fr_ACRC=pow((TXP1.X/(1-TXP1.X)),0.5)*pow(j_g*j_g/(9.8*D),0.5);

	X_tt=Xtt(TXP1);

	n=-0.321;

	Alpha=pow((1+1/Fr_ACRC+X_tt),n);

	if(Alpha>=1) Alpha=0.9999999999;
	if(Alpha<=0) Alpha=0.0000000001;

	return 1/(Alpha/P.vv+(1-Alpha)/P.vl);
}

//B.S.----------------------------------------------------
//Lockhart-Martinelli
double
VolumeTP_LM(TXP TXP1,double G,double D)
{
	VolParams P;
	double Alpha=0;
	double PI1=0;
	double S=0;
	double X_tt=0;
	TXP TXP_prop={0,0,0};

	// liquid and vapor specific volumes
	TXP_prop.P=TXP1.P;
	TXP_prop.X=0;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

	P.vl = PropertyTXPth(VOL,TXP_prop);//reftplthP.v(TXP1.P);
	if(errorLog.IsError()) {
		errorLog.Add("VolumeTP","P.vl");
		return -1;
	}
	
	TXP_prop.P=TXP1.P;
	TXP_prop.X=1.0;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

	P.vv = PropertyTXPth(VOL,TXP_prop);//reftpvthP.v(TXP1.P);
	if(errorLog.IsError()) {
		errorLog.Add("VolumeTP","P.vv");
		return -1;
	}
	PI1=P.vl/P.vv;

	X_tt=Xtt(TXP1);

	if(X_tt<=10) {Alpha=pow((1+pow(X_tt,0.8)),-0.378);}
	else {Alpha=0.823-0.157*log(X_tt);}

	if(Alpha<=0) Alpha=0.000000000001;
	if(Alpha>=1) Alpha=0.999999999999;
	return 1/(Alpha/P.vv+(1-Alpha)/P.vl);
}

//B.S.--------------------------------------------------
//Rigot
double
VolumeTP_Rigot(TXP TXP1,double G,double D)
{
	VolParams P;
	
	double Alpha=0;
	double S=0;
	TXP TXP_prop={0,0,0};

	// liquid and vapor specific volumes
	TXP_prop.P=TXP1.P;
	TXP_prop.X=0;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

	P.vl = PropertyTXPth(VOL,TXP_prop);//reftplthP.v(TXP1.P);
	if(errorLog.IsError()) {
		errorLog.Add("VolumeTP","P.vl");
		return -1;
	}

	TXP_prop.P=TXP1.P;
	TXP_prop.X=1.0;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

	P.vv =  PropertyTXPth(VOL,TXP_prop);//reftpvthP.v(TXP1.P);
	if(errorLog.IsError()) {
		errorLog.Add("VolumeTP","P.vv");
		return -1;
	}
	S=2;
	Alpha=1/(1+(1-TXP1.X)/TXP1.X*P.vl/P.vv*S);

	if(Alpha<=0) Alpha=0.000000000001;
	if(Alpha>=1) Alpha=0.999999999999;

	return 1/(Alpha/P.vv+(1-Alpha)/P.vl);

}

//B.S.-------------------------------------------------
//Smith
double 
VolumeTP_Smith(TXP TXP1,double G,double D)
{
	VolParams P;
	double Alpha=0;
	double PI1=0;
	double S=0;
	double K=0.4;
	TXP TXP_prop={0,0,0};

	// liquid and vapor specific volumes
	TXP_prop.P=TXP1.P;
	TXP_prop.X=0;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);
	P.vl = PropertyTXPth(VOL,TXP_prop);//reftplthP.v(TXP1.P);
	if(errorLog.IsError()) {
		errorLog.Add("VolumeTP","P.vl");
		return -1;
	}

	TXP_prop.P=TXP1.P;
	TXP_prop.X=1.0;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

	P.vv = PropertyTXPth(VOL,TXP_prop);//reftpvthP.v(TXP1.P);
	if(errorLog.IsError()) {
		errorLog.Add("VolumeTP","P.vv");
		return -1;
	}
	PI1=P.vl/P.vv;

	S=K+(1-K)*pow(((1/PI1+K*((1-TXP1.X)/TXP1.X))/(1+K*((1-TXP1.X)/TXP1.X))),0.5);

	Alpha=1/(1+(1-TXP1.X)/TXP1.X*P.vl/P.vv*S);
	
	if(Alpha<=0) Alpha=0.000000000001;
	if(Alpha>=1) Alpha=0.999999999999;

	return 1/(Alpha/P.vv+(1-Alpha)/P.vl);

}

//B.S.-----------------------------------------------------
//Tandon
double
VolumeTP_Tandon(TXP TXP1,double G,double D)
{
	VolParams P;
	
	double Alpha=0;
	double PI1=0;
	double S=0;
	double X_tt=0;
	double Re_L=0;
	double Tsat=0;
	double F_Xtt=0;
	TXP TXP_prop={0,0,0};

	// liquid and vapor specific volumes
	TXP_prop.P=TXP1.P;
	TXP_prop.X=0;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);
	P.vl = PropertyTXPth(VOL,TXP_prop);//reftplthP.v(TXP1.P);
	if(errorLog.IsError()) {
		errorLog.Add("VolumeTP","P.vl");
		return -1;
	}

	TXP_prop.P=TXP1.P;
	TXP_prop.X=1.0;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

	P.vv = PropertyTXPth(VOL,TXP_prop);//reftpvthP.v(TXP1.P);
	if(errorLog.IsError()) {
		errorLog.Add("VolumeTP","P.vv");
		return -1;
	}

	
	TXP_prop.P=TXP1.P;
	TXP_prop.X=0;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

	P.ul=PropertyTXPtr(VISC,TXP_prop);//refsctrPT.mu(TXP1.P,Tsat);
	if(errorLog.IsError()) {
		errorLog.Add("VolumeTP","P.vv");
		return -1;
	}
	Re_L=G*D/P.ul;

	X_tt=Xtt(TXP1);
	
	F_Xtt=0.15*(1/X_tt+2.85/pow(X_tt,0.476));

	if(Re_L<1125)
	Alpha=1-1.928*pow(Re_L,-0.315)/F_Xtt+0.9293*pow(Re_L,-0.63)/pow(F_Xtt,2);
	else
	Alpha=1-0.38*pow(Re_L,-0.088)/F_Xtt+0.0361*pow(Re_L,-0.176)/pow(F_Xtt,2);

	if(Alpha<=0) Alpha=0.000000000001;
	if(Alpha>=1) Alpha=0.999999999999;
	
	return 1/(Alpha/P.vv+(1-Alpha)/P.vl);
}

//B.S.----------------------------------------------------------
//Thome
double
VolumeTP_Thom(TXP TXP1,double G,double D)
{
	VolParams P;
	
/********************************************
Rice&Thom model
S=6.9887-481.29*PI2+14797*PI2^2-179022*PI2^3+875059*PI2^4-1.6248E+06*PI2^5+914478*PI2^6// this is for Thom&Rice
*************************************************/
	double Alpha=0;
	double PI1=0;
	double S=0;
	double PI2=0;
	double a[7]={6.9887,-481.29,14797,-179022,875059,-1.6248e06,914478};
	TXP TXP_prop={0,0,0};

	// liquid and vapor specific volumes
	TXP_prop.P=TXP1.P;
	TXP_prop.X=0;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

	P.vl = PropertyTXPth(VOL,TXP_prop);//reftplthP.v(TXP1.P);
	if(errorLog.IsError()) {
		errorLog.Add("VolumeTP","P.vl");
		return -1;
	}

	TXP_prop.P=TXP1.P;
	TXP_prop.X=1.0;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

	P.vv = PropertyTXPth(VOL,TXP_prop);//reftpvthP.v(TXP1.P);
	if(errorLog.IsError()) {
		errorLog.Add("VolumeTP","P.vv");
		return -1;
	}


	TXP_prop.P=TXP1.P;
	TXP_prop.X=0;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

	const double mul = PropertyTXPtr(VISC,TXP_prop);//refsctrPT.mu(TXP1.P,Tsat);
	if(errorLog.IsError()) {
		errorLog.Add("Xtt","mul");
		return -1;
	}

	TXP_prop.P=TXP1.P;
	TXP_prop.X=1;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);
	const double muv = PropertyTXPtr(VISC,TXP_prop);//refshtrPT.mu(TXP1.P,Tsat);
	if(errorLog.IsError()) {
		errorLog.Add("Xtt","muv");
		return -1;
	}

	PI1=P.vl/P.vv;
	PI2=PI1*pow(mul/muv,0.2);	

	for(int i=0;i<7;i++)
	{S+=a[i]*pow(PI2,i);}
	
	Alpha=1/(1+(1-TXP1.X)/TXP1.X*P.vl/P.vv*S);
	
	if(Alpha<=0) Alpha=0.000000000001;
	if(Alpha>=1) Alpha=0.999999999999;

	return 1/(Alpha/P.vv+(1-Alpha)/P.vl);

}

//B.S.----------------------------------------------
//Premolli
double
VolumeTP_Premoli(TXP TXP1,
		 double G,
		 double D)
{
	VolParams P;
	TXP TXP_prop={0,0,0};

	TXP_prop.P=TXP1.P;
	TXP_prop.X=0;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);
	const double sigma=PropertyTXPtr(TENSION,TXP_prop);//reftpltrP.Tension(TXP1.P);

	// liquid and vapor specific volumes
	TXP_prop.P=TXP1.P;
	TXP_prop.X=0;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

	P.vl = PropertyTXPth(VOL,TXP_prop);//reftplthP.v(TXP1.P);
	if(errorLog.IsError()) {
		errorLog.Add("VolumeTP","P.vl");
		return -1;
	};

	TXP_prop.P=TXP1.P;
	TXP_prop.X=1.0;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

	P.vv = PropertyTXPth(VOL,TXP_prop);//reftpvthP.v(TXP1.P);
	if(errorLog.IsError()) {
		errorLog.Add("VolumeTP","P.vv");
		return -1;
	};

	TXP_prop.P=TXP1.P;
	TXP_prop.X=0;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

	P.ul =PropertyTXPtr(VISC,TXP_prop);//refsctrPT.mu(TXP1.P,Tsat);
	if(errorLog.IsError()) {errorLog.Add("VolumeTP","P.vv");
		return -1;}

	TXP_prop.P=TXP1.P;
	TXP_prop.X=1;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

	P.uv =PropertyTXPtr(VISC,TXP_prop);//refshtrPT.mu(TXP1.P,Tsat);
	if(errorLog.IsError()) {errorLog.Add("VolumeTP","P.vv");
		return -1;}
	const double PI1=P.vl/P.vv;

	const double We=pow(G,2)*D/(sigma/P.vl);

	const double Re_L=G*D/P.ul;

	const double F_2=0.0273*We*pow(Re_L,-0.51)*pow(PI1,0.08);

	const double F_1=1.578*pow(Re_L,-0.19)*pow(PI1,-0.22);

	const double Y=(TXP1.X/(1-TXP1.X))*1/PI1;

	const double S=1+F_1*pow((Y/(1+F_2*Y)-F_2*Y),0.5);

	
	double Alpha=0;

	Alpha=1.0/(1.0+(1-TXP1.X)/TXP1.X*PI1*S);

	if(Alpha>=1) Alpha=0.9999999999;
	if(Alpha<=0) Alpha=0.0000000001;

	return 1/(Alpha/P.vv+(1-Alpha)/P.vl);
}

//B.S.------------------------------------------------
//Homogeneous
double
VolumeTP_Homo(TXP TXP1,double G,double D)
{
	VolParams P;
	TXP TXP_prop={0,0,0};

	// liquid and vapor specific volumes
	TXP_prop.P=TXP1.P;
	TXP_prop.X=0;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);
	P.vl = PropertyTXPth(VOL,TXP_prop);//reftplthP.v(TXP1.P);
	if(errorLog.IsError()) {
		errorLog.Add("VolumeTP","P.vl");
		return -1;
	}

	TXP_prop.P=TXP1.P;
	TXP_prop.X=1.0;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

	P.vv = PropertyTXPth(VOL,TXP_prop);//reftpvthP.v(TXP1.P);
	if(errorLog.IsError()) {
		errorLog.Add("VolumeTP","P.vv");
		return -1;
	}

	// homogeneous mixture specific volume
	return P.vv*TXP1.X + P.vl*(1-TXP1.X);
}

//B.S.------------------------------------------------
//Rouhani
double
VolumeTP_Rouhani(TXP TXP1,double G,double D)
{
	VolParams P;
	TXP TXP_prop={0,0,0};
	
	TXP_prop.P=TXP1.P;
	TXP_prop.X=0;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

	const double sigma=PropertyTXPtr(TENSION,TXP_prop);//reftpltrP.Tension(TXP1.P);

	// liquid and vapor specific volumes
	TXP_prop.P=TXP1.P;
	TXP_prop.X=0;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

	P.vl = PropertyTXPth(VOL,TXP_prop);//reftplthP.v(TXP1.P);
	if(errorLog.IsError()) {
		errorLog.Add("VolumeTP","P.vl");
		return -1;
	};

	TXP_prop.P=TXP1.P;
	TXP_prop.X=1.0;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

	P.vv = PropertyTXPth(VOL,TXP_prop);//reftpvthP.v(TXP1.P);
	if(errorLog.IsError()) {
		errorLog.Add("VolumeTP","P.vv");
		return -1;
	};

	TXP_prop.P=TXP1.P;
	TXP_prop.X=0.0;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

	P.ul =PropertyTXPtr(VISC,TXP_prop);//refsctrPT.mu(TXP1.P,Tsat);
	if(errorLog.IsError()) {errorLog.Add("VolumeTP","P.vv");
		return -1;}

	TXP_prop.P=TXP1.P;
	TXP_prop.X=1.0;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);
	P.uv =PropertyTXPtr(VISC,TXP_prop);//refshtrPT.mu(TXP1.P,Tsat);
	if(errorLog.IsError()) {errorLog.Add("VolumeTP","P.vv");
		return -1;}
	
	
	double Alpha=0;
	
	const double x=TXP1.X;
	const double rho_v=1/P.vv;
	const double rho_l=1/P.vl;

	Alpha=x/rho_v*1.0/((1+0.12*(1-x))*(x/rho_v+(1-x)/rho_l)+1.18*(1-x)*pow((9.8*sigma*(rho_l-rho_v)),0.25)/(G*pow(rho_l,0.5)));
	if(Alpha>=1) Alpha=0.9999999999;
	if(Alpha<=0) Alpha=0.0000000001;

	return 1/(Alpha/P.vv+(1-Alpha)/P.vl);
}

double
VolumeTP_Taitel(TXP TXP1,double G,double D)
{
	TXP TXP_prop;
	// liquid and vapor specific volumes
	TXP_prop.P=TXP1.P;
	TXP_prop.X=0;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

	const double vl = PropertyTXPth(VOL,TXP_prop);//reftplthP.v(TXP1.P);
	if(errorLog.IsError()) {
		errorLog.Add("VolumeTP","P.vl");
		return -1;
	};

	TXP_prop.P=TXP1.P;
	TXP_prop.X=1.0;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

	const double vv = PropertyTXPth(VOL,TXP_prop);//reftpvthP.v(TXP1.P);
	if(errorLog.IsError()) {
		errorLog.Add("VolumeTP","P.vv");
		return -1;
	};

	double PI1=vl/vv;
	double C_sf=1.2;
	double u_d=0.35*pow((9.8*D),0.5);
	double Alpha=1.0/(C_sf+C_sf*((1-TXP1.X)/TXP1.X)*PI1+u_d/(TXP1.X*G*vv));
	if(Alpha>=1) Alpha=0.9999999999;
	if(Alpha<=0) Alpha=0.0000000001;

	return 1/(Alpha/vv+(1-Alpha)/vl);
}

double
VolumeTP_Thome(TXP TXP1,double G,double D,double q)
{
	FlowPattern FlowPat;
	FlowPat.JudgPattern=1;
	
	TXP TXP_prop;
	// liquid and vapor specific volumes
	TXP_prop.P=TXP1.P;
	TXP_prop.X=0;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

	const double vl = PropertyTXPth(VOL,TXP_prop);//reftplthP.v(TXP1.P);
	if(errorLog.IsError()) {
		errorLog.Add("VolumeTP","P.vl");
		return -1;
	};

	TXP_prop.P=TXP1.P;
	TXP_prop.X=1.0;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

	const double vv = PropertyTXPth(VOL,TXP_prop);//reftpvthP.v(TXP1.P);
	if(errorLog.IsError()) {
		errorLog.Add("VolumeTP","P.vv");
		return -1;
	};

	if(q<0)
	Cond_FlowPattern(TXP1,G,D,q,&FlowPat);
	else
	Eva_FlowPattern(TXP1,G,D,q,&FlowPat);

	if(errorLog.IsError()) {
		errorLog.Add("VolumeTP","flow pattern");
		return -1;
	};

	double Alpha=FlowPat.epsilon;

	if(Alpha>=1) Alpha=0.9999999999;
	if(Alpha<=0) Alpha=0.0000000001;

	return 1/(Alpha/vv+(1-Alpha)/vl);
}


double
VolumeTP_FlowPattern(TXP TXP1,double G,double D,double q)
{
	FlowPattern FlowPat;
	FlowPat.JudgPattern=1;

	TXP TXP_prop;
	// liquid and vapor specific volumes
	TXP_prop.P=TXP1.P;
	TXP_prop.X=0;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

	const double vl = PropertyTXPth(VOL,TXP_prop);//reftplthP.v(TXP1.P);
	if(errorLog.IsError()) {
		errorLog.Add("VolumeTP","P.vl");
		return -1;
	};

	TXP_prop.P=TXP1.P;
	TXP_prop.X=1.0;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

	const double vv = PropertyTXPth(VOL,TXP_prop);//reftpvthP.v(TXP1.P);
	if(errorLog.IsError()) {
		errorLog.Add("VolumeTP","P.vv");
		return -1;
	};

	if(q<0)
	Cond_FlowPattern(TXP1,G,D,q,&FlowPat);
	else
	Eva_FlowPattern(TXP1,G,D,q,&FlowPat);

	if(errorLog.IsError()) {
		errorLog.Add("VolumeTP","flow pattern");
		return -1;
	};

	double V=0;

	if(FlowPat.Pattern==2)
	{
	//V=VolumeTP_ACRC(TXP1,G,D);
	V = VolumeTP_Rouhani(TXP1,G,D);
	
	}
	else
	{
	V=VolumeTP_Taitel(TXP1,G,D);}
	
	if(errorLog.IsError()) {
		errorLog.Add("VolumeTP","Flowpattern");
		return -1;
	};


	return V;
}