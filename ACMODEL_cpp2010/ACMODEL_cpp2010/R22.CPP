/********************************************************************
This file contains subroutines for evaluating the thermodynamic 
properties of refrigerants.  The particular refrigerant is 
determined by the installed lookup table.
********************************************************************/
 
#include "stdafx.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "spline.h"
#include "ref.h"
#include "numeric.h"
#include "errorlog.h"
#include "r22.h"
#include "spline1d.h"
#include "spline2d.h"

#include "outelem.h"//B.S.

extern ErrorLog errorLog;

extern RefPropTPthP reftpvthP;
extern RefPropTPthP reftplthP;

extern RefPropTPtrP reftpltrP;//B.S.
extern RefPropTPtrP reftpvtrP;//B.S.

/*#ifdef _R410A
extern RefPropSPthPT refcritthPT;//B.S.
extern RefPropSPtrPT refcrittrPT;//B.S.
extern RefPropSPthPH refcritthPH;//B.S.
#endif*/

extern RefPropSPthPT refshthPT;
extern RefPropSPthPT refscthPT;
extern RefPropSPthPH refshthPH;
extern RefPropSPthPH refscthPH;
extern RefPropSPtrPT refshtrPT;
extern RefPropSPtrPT refsctrPT;

#ifdef _RefMix
extern RefPropTPthPH refTPthPH;
extern RefPropTPthPX refTPthPX;
extern RefPropTPtrPX refTPtrPX;
#endif

//#define DEBUG_PROPERTYTXPTH

/********************************************************************
Evaluates properties given TXP data.
********************************************************************/
double
PropertyTXPth(int prop,TXP TXP)
{
	double a;
	int pos = 0;

#ifdef _R22
	if (TXP.X>=1) {
		const double Tsat = reftpvthP.Tsat(TXP.P);
		if(errorLog.IsError()) {
			char str[128];
			sprintf(str,"prop=%d T=%lf X=%lf P=%lf",prop,TXP.T,TXP.X,TXP.P);
			errorLog.Add("PropertyTXPth",str);
			return -1;
		}

		if(prop==TSAT) 
			
			{return Tsat;}

		if (TXP.T>Tsat) {
			a = refshthPT.Evaluate_y(prop-1,TXP.P,TXP.T);
			pos = 1;
		} else {
			a = reftpvthP.Evaluate(prop,TXP.P);
			pos = 2;
		}
	} else if (TXP.X<=0) {
		const double Tsat = reftplthP.Tsat(TXP.P);
		if(errorLog.IsError()) {
			char str[128];
			sprintf(str,"prop=%d T=%lf X=%lf P=%lf",prop,TXP.T,TXP.X,TXP.P);
			errorLog.Add("PropertyTXPth",str);
			return -1;
		}

		if(prop==TSAT) {return Tsat;}

		if (TXP.T<Tsat) {
			a = refscthPT.Evaluate_y(prop-1,TXP.P,TXP.T);
			pos = 3;
		} else {
			a = reftplthP.Evaluate(prop,TXP.P);
			pos = 4;
		}
	} else {
		const double al = reftplthP.Evaluate(prop,TXP.P);
		const double av = reftpvthP.Evaluate(prop,TXP.P);
		a = al + TXP.X*(av-al);
		pos = 5;
	}

	if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"(end) prop=%d T=%lf X=%lf P=%lf",prop,TXP.T,TXP.X,TXP.P);
		errorLog.Add("PropertyTXPth",str);
		return -1;
	}

	#ifdef DEBUG_PROPERTYTXPTH
	{
	static int init = 1;
	FILE* fp;
	if(init) {
		fp = fopen("r221.dbg","w");
                init = 0;
	} else {
		fp = fopen("r221.dbg","a");
	}
	if(fp) {
		fprintf(fp,"TXP=%lf,%lf,%lf a=%lf pos=%d prop=%d\n",TXP.T,TXP.X,TXP.P,a,pos,prop);
		fclose(fp);
	}
	}
	#endif

#endif

#ifdef _RefMix

	if (TXP.X>=1) {
		const double Tsat = reftpvthP.Tsat(TXP.P);
		if(errorLog.IsError()) {
			char str[128];
			sprintf(str,"prop=%d T=%lf X=%lf P=%lf",prop,TXP.T,TXP.X,TXP.P);
			errorLog.Add("PropertyTXPth",str);
			return -1;}

		if(prop==TSAT) 			
			{return Tsat;}

		if (TXP.T>Tsat) {
			a = refshthPT.Evaluate_y(prop-1,TXP.P,TXP.T);
			pos = 1;
		} else {
			a = reftpvthP.Evaluate(prop,TXP.P);
			pos = 2;
		}
	} else if (TXP.X<=0) {
		const double Tsat = reftplthP.Tsat(TXP.P);
		if(errorLog.IsError()) {
			char str[128];
			sprintf(str,"prop=%d T=%lf X=%lf P=%lf",prop,TXP.T,TXP.X,TXP.P);
			errorLog.Add("PropertyTXPth",str);
			return -1;}

		if(prop==TSAT) {return Tsat;}

		if (TXP.T<Tsat) {
			a = refscthPT.Evaluate_y(prop-1,TXP.P,TXP.T);
			pos = 3;
		} else {
			a = reftplthP.Evaluate(prop,TXP.P);
			pos = 4;
		}
	} else {
		
		if(TXP.P<4600&&TXP.P>150)//B.S.
			{
			a = refTPthPX.Evaluate_y(prop,TXP.P,TXP.X);
			}
			else
			{
			const double al = reftplthP.Evaluate(prop,TXP.P);
			const double av = reftpvthP.Evaluate(prop,TXP.P);
			a = al + TXP.X*(av-al);
			if(prop==DENSITY)
			{
			const double al = reftplthP.Evaluate(VOL,TXP.P);
			const double av = reftpvthP.Evaluate(VOL,TXP.P);
			a = al + TXP.X*(av-al);
			a=1/a;}
			}
		pos = 5;
	}

	if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"(end) prop=%d T=%lf X=%lf P=%lf",prop,TXP.T,TXP.X,TXP.P);
		errorLog.Add("PropertyTXPth",str);
		return -1;
	}

	#ifdef DEBUG_PROPERTYTXPTH
	{
	static int init = 1;
	FILE* fp;
	if(init) {
		fp = fopen("r221.dbg","w");
                init = 0;
	} else {
		fp = fopen("r221.dbg","a");
	}
	if(fp) {
		fprintf(fp,"TXP=%lf,%lf,%lf a=%lf pos=%d prop=%d\n",TXP.T,TXP.X,TXP.P,a,pos,prop);
		fclose(fp);
	}
	}
	#endif

#endif

	return a;
}

/********************************************************************
Evaluates properties given TXP data.
********************************************************************/
double
PropertyTXPtr(int prop,TXP TXP)
{
	double a;

#ifdef _R22
	if (TXP.X>=1) {
		const double Tsat = reftpvthP.Tsat(TXP.P);
		if(errorLog.IsError()) {
			char str[128];
			sprintf(str,"prop=%d T=%lf X=%lf P=%lf",prop,TXP.T,TXP.X,TXP.P);
			errorLog.Add("PropertyTXPtr",str);
			return -1;
		}
		if (TXP.T>Tsat) {//superheat vapor
			a = refshtrPT.Evaluate_y(prop,TXP.P,TXP.T);
		} else {//saturated vapor
			a = reftpvtrP.Evaluate(prop,TXP.P);//refshtrPT.Evaluate_y(prop,TXP.P,Tsat);
		}
	} else if (TXP.X<=0) {
		const double Tsat = reftplthP.Tsat(TXP.P);
		if(errorLog.IsError()) {
			char str[128];
			sprintf(str,"prop=%d T=%lf X=%lf P=%lf",prop,TXP.T,TXP.X,TXP.P);
			errorLog.Add("PropertyTXPtr",str);
			return -1;
		}
		if (TXP.T<Tsat) {//subcooled liquid
			a = refsctrPT.Evaluate_y(prop,TXP.P,TXP.T);
		} else {//saturated liquid
			a = reftpltrP.Evaluate(prop,TXP.P);//refsctrPT.Evaluate_y(prop,TXP.P,Tsat);
		}
	} else {
		errorLog.Add("PropertyTXPtr","two phase");
		return -1;
	}

	if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"(end) prop=%d T=%lf X=%lf P=%lf",prop,TXP.T,TXP.X,TXP.P);
		errorLog.Add("PropertyTXPtr",str);
		return -1;
	}
#endif

	
#ifdef _RefMix
	if (TXP.X>=1) {
		const double Tsat = reftpvthP.Tsat(TXP.P);
		if(errorLog.IsError()) {
			char str[128];
			sprintf(str,"prop=%d T=%lf X=%lf P=%lf",prop,TXP.T,TXP.X,TXP.P);
			errorLog.Add("PropertyTXPtr",str);
			return -1;}

		if (TXP.T>Tsat) {//superheat vapor
			a = refshtrPT.Evaluate_y(prop,TXP.P,TXP.T);
		
		} else {//saturated vapor
		
			a = reftpvtrP.Evaluate(prop,TXP.P);//refshtrPT.Evaluate_y(prop,TXP.P,Tsat);
		}
	} else if (TXP.X<=0) {
		const double Tsat = reftplthP.Tsat(TXP.P);
		if(errorLog.IsError()) {
			char str[128];
			sprintf(str,"prop=%d T=%lf X=%lf P=%lf",prop,TXP.T,TXP.X,TXP.P);
			errorLog.Add("PropertyTXPtr",str);
			return -1;}
		if (TXP.T<Tsat) {//subcooled liquid
			a = refsctrPT.Evaluate_y(prop,TXP.P,TXP.T);
		} else {//saturated liquid
		
			a = reftpltrP.Evaluate(prop,TXP.P);//refsctrPT.Evaluate_y(prop,TXP.P,Tsat);
		
		}
	} else {
		if(prop==TENSION)
		{
		if(TXP.P<4600)
		{
		a=refTPtrPX.Sigma(TXP.P,TXP.X);}
		else
		{
		const double al = reftpltrP.Evaluate(prop,TXP.P);
		const double av = reftpvtrP.Evaluate(prop,TXP.P);
		a = al + TXP.X*(av-al);}
		}//only the surface tension can be checked in two phase state
		else
		{
		errorLog.Add("PropertyTXPtr","two phase");
		return -1;}
	}

	if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"(end) prop=%d T=%lf X=%lf P=%lf",prop,TXP.T,TXP.X,TXP.P);
		errorLog.Add("PropertyTXPtr",str);
		return -1;
	}
#endif

	return a;
}

/********************************************************************
Convert T,X,P format to TXP format
********************************************************************/
TXP
toTXP(double T,double X,double P)
{
	static TXP txp;
	txp.T=T;
	txp.X=X;
	txp.P=P;
	return txp;
}

HP
toHP(double H,double P)
{
	static HP hp;
	hp.H = H;
	hp.P = P;
	return hp;
}

/********************************************************************
Convert HP format to TXP format.
********************************************************************/
TXP
HPtoTXP(HP HP)
{
	TXP TXP={0,0,0};

#ifdef _R22
	TXP.P = HP.P;
	double hl = reftplthP.h(HP.P);
	if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"hl(P) H=%lf P=%lf",HP.H,HP.P);
		errorLog.Add("HPtoTXP",str);
		return TXP;
	}
	double hv = reftpvthP.h(HP.P);
	if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"hv(P) H=%lf P=%lf",HP.H,HP.P);
		errorLog.Add("HPtoTXP",str);
		return TXP;
	}

	if (HP.H>hv) {
		TXP.X=1;
		TXP.T = refshthPH.T(HP.P,HP.H);
		if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"hv(P) H=%lf P=%lf",HP.H,HP.P);
		errorLog.Add("HPtoTXP",str);
		return TXP;
		}
	} else if (HP.H<hl) {
		TXP.X=0;
		TXP.T = refscthPH.T(HP.P,HP.H);
		if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"hv(P) H=%lf P=%lf",HP.H,HP.P);
		errorLog.Add("HPtoTXP",str);
		return TXP;
		}
	} else {
		TXP.T = reftplthP.Tsat(HP.P);
		TXP.X = (HP.H-hl)/(hv-hl);
	}
#endif

#ifdef _RefMix
	TXP.P = HP.P;
	double hl = reftplthP.h(HP.P);
	if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"hl(P) H=%lf P=%lf",HP.H,HP.P);
		errorLog.Add("HPtoTXP",str);
		return TXP;
	}
	double hv = reftpvthP.h(HP.P);
	if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"hv(P) H=%lf P=%lf",HP.H,HP.P);
		errorLog.Add("HPtoTXP",str);
		return TXP;
	}

	if (HP.H>hv) {
		TXP.X=1;
		TXP.T = refshthPH.T(HP.P,HP.H);
		if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"hv(P) H=%lf P=%lf",HP.H,HP.P);
		errorLog.Add("HPtoTXP",str);
		return TXP;}
	} else if (HP.H<hl) {
		TXP.X=0;
		TXP.T = refscthPH.T(HP.P,HP.H);
		if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"hv(P) H=%lf P=%lf",HP.H,HP.P);
		errorLog.Add("HPtoTXP",str);
		return TXP;}
	} else {
			
		if(TXP.P<4600)
		{
		TXP.T = refTPthPH.T(HP.P,HP.H);
		TXP.X = refTPthPH.X(HP.P,HP.H);
		}
		else
		{
		const double al = reftpltrP.Evaluate(ENTH,TXP.P);
		const double av = reftpvtrP.Evaluate(ENTH,TXP.P);
		TXP.X=(HP.H-al)/(av-al);
		const double tl=reftpltrP.Evaluate(TSAT,TXP.P);
		const double tv=reftpvtrP.Evaluate(TSAT,TXP.P);
		TXP.T=tl + TXP.X*(tv-tl);
		}
		if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"hv(P) H=%lf P=%lf",HP.H,HP.P);
		errorLog.Add("HPtoTXP",str);
		return TXP;}
	}
#endif


	return TXP;
}

/********************************************************************
Convert TXP format to HP format.
********************************************************************/
HP
TXPtoHP(TXP TXP)
{
	HP HP;

	HP.P = TXP.P;
	HP.H = PropertyTXPth(ENTH,TXP);
	if(errorLog.IsError()) {
		char str[128];
                sprintf(str,"TXP=(%lf,%lf,%lf)",TXP.T,TXP.X,TXP.P);
		errorLog.Add("TXPtoHP",str);
	}

	return HP;
}

//the following subrountine is designed for testing the property program

int Exercise_Prop()//excise property subroutine
{
	TXP TXP_prop={0,0,0};
	int i=0,j=0;

	std::ofstream outstrm; 
	outstrm.open("out_prop.xls",std::ios::app);

	for(int i=100;i<=440;i++)
	{
	TXP_prop.P=i*10;
	TXP_prop.X=0.0;
	const double Tsat=PropertyTXPth(TSAT,TXP_prop)-0.1;

	for(int j=0;j<=600;j++)
	{
	TXP_prop.X=0.0;
	TXP_prop.T=Tsat-j*3;

	if(TXP_prop.P==4830)
	{
		int shenb=0;}

	if(TXP_prop.T<-100) j=600;
	
	const double h=PropertyTXPth(ENTH,TXP_prop);
	const double v=PropertyTXPth(VOL,TXP_prop);
	const double entro=PropertyTXPth(ENTR,TXP_prop);
	const double density=PropertyTXPth(DENSITY,TXP_prop);
	outstrm<<TXP_prop.P<<"	"<<TXP_prop.X<<"	"<<TXP_prop.T<<"	"<<h<<"	"<<v<<"	"<<entro<<"	"
		<<density;
	if(errorLog.IsError()) {
		char str[128];
                sprintf(str,"TXP=(%lf,%lf,%lf)",TXP_prop.T,TXP_prop.X,TXP_prop.P);
		errorLog.Add("TXP",str);
		return 0;}

	const double cond=PropertyTXPtr(COND,TXP_prop);
	const double visc=PropertyTXPtr(VISC,TXP_prop);
	const double cp=PropertyTXPtr(SPEC,TXP_prop);
	const double tension=0;//PropertyTXPtr(TENSION,TXP_prop);
		
	outstrm<<"	"<<cond<<"	"<<visc<<"	"<<cp<<"	"<<tension<<std::endl;

	if(errorLog.IsError()) {
		char str[128];
                sprintf(str,"TXP=(%lf,%lf,%lf)",TXP_prop.T,TXP_prop.X,TXP_prop.P);
		errorLog.Add("TXP",str);
		return 0;}

	}
	}

	outstrm.close();
	return 0;
}
