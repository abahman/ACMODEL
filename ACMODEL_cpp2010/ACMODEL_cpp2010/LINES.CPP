#include "stdafx.h"

#include "r22.h"
#include "mass.h"
#include "volume.h"
#include <math.h>
#include <stdio.h>
#include "errorlog.h"
#include "corr.h"//B.S.
#include "mainfunc.h"//B.S.

#include "spline2d.h"
#include "spline1d.h"//B.S.
#include "cmine.h"//B.S.
#include "getdata.h"
#include "dataelem.h"

extern ErrorLog errorLog;
extern AirProp air;

static double pi = acos(-1.0);

/****************************************************************
Suction line model.  It is a simple model that has no pressure
drop or heat transfer.  The model only calculates the mass of
charge contained within it.
****************************************************************/
void
SuctionLine(HP HPi,double mr,HP *HPo,MASS *m)
{
	
	/*FILE *fp = fopen("./InputDoc/acmodel.dat","r");
	if(!PosComponentData(fp,SUCTIONLINE,1)) {
		errorLog.Add("SUCTIONLINE","Parameters not found");
		return;
		}
	
	double l=GetDouble(fp);//length
	double d=GetDouble(fp);//inside diameter
	double OD=GetDouble(fp);//outside diameter
	fclose(fp) ;*/ //reserved

	TXP	TXP_prop={0,0,0};

	const double L1=0.5715;
	const double D1=17.27e-3;
	
	double A = pi/4.0*D1*D1;  // cross section area
	const double V1=L1*A;  // pipe volume
	double G = mr/A;	// mass flux

	TXP TXP_1;
	TXP_1=HPtoTXP(HPi);
	if(errorLog.IsError()) {
	errorLog.Add("SuctionLine","HPtoTXP");}

	TXP_prop.P=TXP_1.P;
	TXP_prop.T=TXP_1.T;
	TXP_prop.X=TXP_1.X;
	double rho_1=1/PropertyTXPth(VOL,TXP_prop);//refshthPT.rho(TXP_1.P,TXP_1.T);
	if(errorLog.IsError()) {
	errorLog.Add("SuctionLine","rho");}

	TXP_prop.X=1;

	double MU=PropertyTXPtr(VISC,TXP_prop);//refshtrPT.mu(TXP_1.P,TXP_1.T);
	if(errorLog.IsError()) {
	errorLog.Add("SuctionLine","MU");
	}


	double RE=G*D1/MU;
	double F=0.184e0*pow(RE,(-0.2e0));
	double DP=1e-3*(pow(G,2e0)/(rho_1*2e0))*F*L1/D1;
		
	HPo->P=HPi.P+DP;
	HPo->H=HPi.H;
	
	TXP TXP_8;
	TXP_8=HPtoTXP(*HPo);
	if(errorLog.IsError()) {
	errorLog.Add("SuctionLine","HPtoTXP");}

	double v = VolumeALL(HPtoTXP(*HPo),G,D1,0);
	if(errorLog.IsError()) {
		errorLog.Add("SuctionLine","v");
	}

	m->m = V1/v;

	const double L2=4.845;
	const double D2=16.89e-3;

	A = pi/4.0*D2*D2;   // cross section area
	const double V2=L2*A;         			  // pipe volume
	G = mr/A;		// mass flux

	TXP_1=HPtoTXP(*HPo);
	if(errorLog.IsError()) {
	errorLog.Add("SuctionLine","HPtoTXP");
	}

	TXP_prop.P=TXP_1.P;
	TXP_prop.T=TXP_1.T;
	TXP_prop.X=TXP_1.X;
	rho_1=1/PropertyTXPth(VOL,TXP_prop);//refshthPT.rho(TXP_1.P,TXP_1.T);
	if(errorLog.IsError()) {
	errorLog.Add("SuctionLine","rho");
	}

	TXP_prop.X=1;

	MU=PropertyTXPtr(VISC,TXP_prop);//refshtrPT.mu(TXP_1.P,TXP_1.T);
	if(errorLog.IsError()) {
	errorLog.Add("SuctionLine","MU");
	}


	RE=G*D2/MU;
	F=0.184e0*pow(RE,(-0.2e0));
	DP=1e-3*pow(G,2e0)/(rho_1*2e0)*F*L2/D2;
	HPo->P=HPo->P+DP;
	
	const double Q_gain=114.0;//35 is the ambient temperature

	HPo->H=HPo->H-Q_gain/mr;


	TXP_8=HPtoTXP(*HPo);
	if(errorLog.IsError()) {
	errorLog.Add("SuctionLine","HPtoTXP");
	}

	v = VolumeALL(HPtoTXP(*HPo),G,D2,0);
	if(errorLog.IsError()) {
		errorLog.Add("SuctionLine","v");
	}

	const double V3=0.5924*pi/4.0*(17.27e-3)*(17.27e-3)+
		pi/4.0*(8.11e-3)*(8.11e-3)*(0.2821+0.2821+0.2402+0.2656+0.2082);//evaporator exit header inner volume
	m->m = m->m+(V2+V3)/v;
	m->V=V1+V2+V3;
}

/****************************************************************
Distribution tube model.  It is a simple model that has no pressure
drop or heat transfer.  The model only calculates the mass of
charge contained within it.
****************************************************************/
void
DistributionTubes(HP HPi,double mr,HP *HPo,MASS *m)
{
	/*FILE *fp = fopen("./InputDoc/acmodel.dat","r");
	if(!PosComponentData(fp,DISTRILINE,1)) {
		errorLog.Add("DISTRIBUTLINE","Parameters not found");
		return;
		}
	
	double l=GetDouble(fp);//length
	double d=GetDouble(fp);//inside diameter
	double OD=GetDouble(fp);//outside diameter
	fclose(fp) ;*/

	m->V=0;
	m->m=0;

	const double DP=296;//measured data from Bo Shen's R410A test.
	
	HPo->P=HPi.P;
	HPo->H=HPi.H;

}

/****************************************************************
Hot gas line model.  It is a simple model that has no pressure
drop or heat transfer.  The model only calculates the mass of
charge contained within it.
****************************************************************/
void
HotGasLine(HP HPi,double mr,HP *HPo,MASS *m,void *Params)
{
	/*FILE *fp = fopen("./InputDoc/acmodel.dat","r");
	if(!PosComponentData(fp,DISCHARGELINE,1)) {
		errorLog.Add("DISCHARGELINE","Parameters not found");
		return;
		}
	
	double l=GetDouble(fp);//length
	double d=GetDouble(fp);//inside diameter
	double OD=GetDouble(fp);//outside diameter
	fclose(fp) ;*/
	static InVars *P;

	P=(InVars*)Params;

	const double l=(1.179+0.1785);
	const double d=11.28e-3;
	const double d_o=12.7e-3;
	TXP TXP_prop={0,0,0};

	double A = pi/4.0*d*d;        // cross section area
	m->V = l*A;                      // pipe volume
	double G = mr/A;                 // mass flux

	TXP TXP_2;
	TXP_2=HPtoTXP(HPi);
	if(errorLog.IsError()) {
	errorLog.Add("HotGasLine","HPtoTXP");
	}
	
	TXP_prop.P=TXP_2.P;
	TXP_prop.T=TXP_2.T;
	TXP_prop.X=TXP_2.X;

	const double rho_2=1/PropertyTXPth(VOL,TXP_prop);//refshthPT.rho(TXP_2.P,TXP_2.T);
	if(errorLog.IsError()) {
	errorLog.Add("HotGasLine","rho");
	}

	TXP_prop.X=1;
	const double MU=PropertyTXPtr(VISC,TXP_prop);//refshtrPT.mu(TXP_2.P,TXP_2.T);
	if(errorLog.IsError()) {
	errorLog.Add("HotGasLine","MU");
	}


	double v = VolumeALL(HPtoTXP(HPi),G,d,0);
	if(errorLog.IsError()) {
		errorLog.Add("HotGasLine");
	}

	const double V_branch=pi/4.0*(8.11e-3)*(8.11e-3)*(0.08993+0.08993);
	m->V=m->V+V_branch;
	m->m=m->V/v;

	//heat transfer and pressure drop
	const double RE=G*d/MU;
	const double F=0.184e0*pow(RE,(-0.2e0));
	const double DP=1e-3*pow(G,2e0)/(2e0*rho_2)*F*l/d;

	const double rho_air=1/air.v(P->Tai);
	if(errorLog.IsError()) {
		errorLog.Add("HotGasLine","rho_air");
	}
	
	const double Ga=P->Goa*0.0004719*rho_air/(0.75*0.75); 
	const double Ta=P->Tai;
	const double ho=Coef_Hilpert(Ga,Ta,d_o);
	if(errorLog.IsError()) {
		errorLog.Add("HotGasLine","ho");}

	const double hi=ConvCoeffSP_Smooth(TXP_2,G,d);
	if(errorLog.IsError()) {
	errorLog.Add("HotGasLine","hi");}

	const double R_W=log(d_o/(d))/(2.0*pi*400.0*l);
	const double R_o=1.0/(ho*pi*l*d_o);
	const double R_i=1.0/(hi*pi*l*d);
	const double R=R_o+R_W+R_i;

	const double ma=Ga*d_o*l;

	const double Q_loss = (TXP_2.T-Ta)*CmineCrossFlow_dry(R,mr,ma,TXP_2,Ta);
	if(errorLog.IsError()) {
		errorLog.Add("HotGasLine","Q");
		return;
	}

	HPo->P=HPi.P-DP;
	HPo->H=HPi.H-Q_loss/mr;					// outlet state equals inlet state

}

/****************************************************************
Liquid line model.  It is a simple model that has no pressure
drop or heat transfer.  The model only calculates the mass of
charge contained within it.
****************************************************************/
void
LiquidLine(HP HPi,double mr,HP *HPo,MASS *m)
{
	// liquid line
	
	/*FILE *fp = fopen("./InputDoc/acmodel.dat","r");
	if(!PosComponentData(fp,LIQUIDLINE,1)) {
		errorLog.Add("liquid line","Parameters not found");
		return;
		}
	double l=GetDouble(fp);//length
	double d=GetDouble(fp);//inside diameter
	double OD=GetDouble(fp);//outside diameter
	fclose(fp) ;*/

	double l=0.4318;
	double d=8.11e-3;

	double A = pi/4.0*d*d;			/* cross section area */
	m->V = l*A;					/* pipe volume */
	double G = mr/A;						/* mass flux */
	TXP TXP_prop={0,0,0};
	
	TXP TXP_4;
	TXP_4=HPtoTXP(HPi);
	if(errorLog.IsError()) {
	errorLog.Add("LiquidLine","HPtoTXP");
	}
	
	TXP_prop.P=TXP_4.P;
	TXP_prop.T=TXP_4.T;
	TXP_prop.X=TXP_4.X;

	double rho_4=1/PropertyTXPth(VOL,TXP_prop);//refshthPT.rho(TXP_4.P,TXP_4.T);
	if(errorLog.IsError()) {
	errorLog.Add("LiquidLine","rho");
	}

	TXP_prop.X=0;
	double MU=PropertyTXPtr(VISC,TXP_prop);//refshtrPT.mu(TXP_4.P,TXP_4.T);
	if(errorLog.IsError()) {
	errorLog.Add("LiquidLine","MU");
	}


	double RE=G*d/MU;
	double F=0.184e0*pow(RE,(-0.2e0));
	double DP=1e-3*pow(G,2e0)/(2.0*rho_4)*F*l/d;

	HPo->P=HPi.P-DP;
	HPo->H=HPi.H;

	double v = VolumeALL(HPtoTXP(HPi),G,d,0);
	if(errorLog.IsError()) {
		errorLog.Add("LiquidLine");
	}
	m->m=m->V/v;

	// liquid line dryer

	/*fp = fopen("./InputDoc/acmodel.dat","r");
	if(!PosComponentData(fp,FILTER,1)) {
		errorLog.Add("filter drier","Parameters not found");
		return;
		}
	
	l=GetDouble(fp);//length
	d=GetDouble(fp);//inside diameter
	double fv=GetDouble(fp);//outside diameter
	fclose(fp) ;*/

	const double V2= 0.000185;//measured by charging oil into the drier/filter.
	G = mr/A;						/* mass flux */
	v = VolumeALL(HPtoTXP(*HPo),G,d,0);
	if(errorLog.IsError()) {
		errorLog.Add("LiquidLine");
	}
	m->m = m->m+V2/v;

	l = 4.763;					//connected liquid line length
	d = 8e-3;					//connected liquid line length

	A = pi/4.0*d*d;			// cross section area 
	G = mr/A;						//mass flux
	const double V3=l*A;

	v = VolumeALL(HPtoTXP(*HPo),G,d,0);
	if(errorLog.IsError()) {
		errorLog.Add("LiquidLine");
	}
	m->m = m->m + V3/v;
	m->V=m->V+V2+V3;

	TXP_4=HPtoTXP(*HPo);
	if(errorLog.IsError()) {
	errorLog.Add("LiquidLine","HPtoTXP");
	}
	
	TXP_prop.P=TXP_4.P;
	TXP_prop.T=TXP_4.T;
	TXP_prop.X=TXP_4.X;

	rho_4=1/PropertyTXPth(VOL,TXP_prop);//refshthPT.rho(TXP_4.P,TXP_4.T);
	if(errorLog.IsError()) {
	errorLog.Add("LiquidLine","rho");
	}

	TXP_prop.X=0;
	MU=PropertyTXPtr(VISC,TXP_prop);//refshtrPT.mu(TXP_4.P,TXP_4.T);
	if(errorLog.IsError()) {
	errorLog.Add("LiquidLine","MU");
	}

	RE=G*d/MU;
	F=0.184e0*pow(RE,(-0.2e0));
	DP=1e-3*pow(G,2e0)/(2*rho_4)*F*l/d;

	const double Q_loss=102.0;
	HPo->P=HPo->P-DP;
	HPo->H=HPo->H-Q_loss/mr;
}

