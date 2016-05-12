#include "stdafx.h"

#include "r22.h"
#include "mass.h"
#include "volume.h"
#include <math.h>
#include <stdio.h>
#include "errorlog.h"
#include "corr.h"//B.S.
#include "mainfunc.h"//B.S.
#include "pressure.h"//B.S.
#include "numeric.h"//B.S.
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
	
	FILE *fp = fopen("./InputDoc/acmodel.dat","r");
	if(!PosComponentData(fp,SUCTIONLINE,1)) {
		errorLog.Add("SUCTIONLINE","Parameters not found");
		return;
		}
	
	double l=GetDouble(fp);//length
	double d=GetDouble(fp);//inside diameter
	double OD=GetDouble(fp);//outside diameter
	double dQ = GetDouble(fp);//heat transfer rate in the suction line
	double dP = GetDouble(fp);//differential pressure drop in the suction line
	double dV = GetDouble(fp);//additional volume connected with the suction line, usually refered as the evaporator header
	fclose(fp) ;

	TXP	TXP_prop={0,0,0};

	HPo->P = HPi.P;
	HPo->H = HPi.H;

	
	const double Q_gain=dQ;

	HPo->H=HPo->H-Q_gain/mr;


	const double L1=l;
	const double D1=d;//3/4 in OD tube
	
	double A1 = pi/4.0*D1*D1;  // cross section area
	const double G = mr/A1;

	//Pressure drop
	if(dP<0)//calculate the pressure drop
	{

	TXP_prop=HPtoTXP(HPi);
	if(errorLog.IsError()) {
	errorLog.Add("SuctionLine","HPtoTXP");}

	double rho_1=1/PropertyTXPth(VOL,TXP_prop);
	if(errorLog.IsError()) {
	errorLog.Add("SuctionLine","rho_1");}

	TXP_prop.X=1;

	double MU=PropertyTXPtr(VISC,TXP_prop);
	if(errorLog.IsError()) {
	errorLog.Add("SuctionLine","MU");
	}

	double RE=G*D1/MU;
	double F=0.184e0*pow(RE,(-0.2e0));
	double DP=1e-3*(pow(G,2e0)/(rho_1*2e0))*F*L1/D1;
	HPo->P = HPo->P + DP;
	}
	else
	{
	HPo->P = HPo->P +dP;
	}

	const double V1=L1*A1;  // suction line pipe volume

	const double v = VolumeALL(HPtoTXP(*HPo),G,D1,0);
	if(errorLog.IsError()) {
		errorLog.Add("SuctionLine","v1");
	}
	

	m->m = (V1+dV)/v;
	m->V=(V1+dV);
}

void
SuctionLine_Lumped(HP HPi,double mr,HP *HPo,MASS *m)
{
	FILE *fp = fopen("./InputDoc/acmodel.dat","r");
	if(!PosComponentData(fp,SUCTIONLINE,1)) {
		errorLog.Add("SUCTIONLINE","Parameters not found");
		return;
		}
	
	double l=GetDouble(fp);//length
	double d=GetDouble(fp);//inside diameter
	double OD=GetDouble(fp);//outside diameter
	double dQ = GetDouble(fp);//heat transfer rate in the suction line
	double dP = GetDouble(fp);//differential pressure drop in the suction line
	double dV = GetDouble(fp);//additional volume connected with the suction line, usually refered as the evaporator header
	fclose(fp) ;
	
	HPo->P = HPi.P+dP;
	HPo->H = HPi.H;
	const double Q_gain=dQ;
	HPo->H=HPo->H-Q_gain/mr;
	m->m=0;
	m->V=0;
	return;

}
/****************************************************************
Distribution tube model.  It is a simple model that has no pressure
drop or heat transfer.  The model only calculates the mass of
charge contained within it.
****************************************************************/
void
DistributionTubes(HP HPi,double mr,HP *HPo,MASS *m)//only consider the pressure drop
{
	FILE *fp = fopen("./InputDoc/acmodel.dat","r");
	if(!PosComponentData(fp,DISTRILINE,1)) {
		errorLog.Add("DISTRIBUTLINE","Parameters not found");
		return;
		}
	
	double l=GetDouble(fp);//length
	double d=GetDouble(fp);//inside diameter
	double OD=GetDouble(fp);//outside diameter
	double dP = GetDouble(fp);//differential pressure drop across the distribution lines
	fclose(fp) ;

	m->V=0;
	m->m=0;

	const double DP=dP;//measured data from Bo Shen's R22 test.
	
	HPo->P=HPi.P+DP;
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
	FILE *fp = fopen("./InputDoc/acmodel.dat","r");
	if(!PosComponentData(fp,DISCHARGELINE,1)) {
		errorLog.Add("DISCHARGELINE","Parameters not found");
		return;
		}
	
	double l=GetDouble(fp);//length
	double d=GetDouble(fp);//inside diameter
	double OD=GetDouble(fp);//outside diameter
	double dQ = GetDouble(fp);//heat transfer rate in the discharge line
	double dP = GetDouble(fp);//differential pressure drop across the discharge line
	double dV = GetDouble(fp);//additional volume connected with the discharge line, usually refered to inlet header of the condenser
	fclose(fp) ;

	
	TXP TXP_prop={0,0,0};

	const double l1= l;
	const double d1= d;
	double A1 = pi/4.0*d1*d1;        // cross section area for the connection tube from compressor to four-way valve
	const double V1 = l1*A1;         // pipe volume
	double G = mr/A1;                // mass flux

	TXP_prop=HPtoTXP(HPi);
	if(errorLog.IsError()) {
	errorLog.Add("HotGasLine","HPtoTXP,1");
	}
	
	const double rho_1=1/PropertyTXPth(VOL,TXP_prop);
	if(errorLog.IsError()) {
	errorLog.Add("HotGasLine","rho_1");
	}
	
	m->m = (V1+dV) * rho_1;
	m->V = (V1+dV);

	const double Q_loss = dQ;
	
	//Pressure drop
	if(dP<0)//calculate the pressure drop
	{
	TXP_prop.X=1;

	double MU=PropertyTXPtr(VISC,TXP_prop);
	if(errorLog.IsError()) {
	errorLog.Add("SuctionLine","MU");
	}

	double RE=G*d1/MU;
	double F=0.184e0*pow(RE,(-0.2e0));
	double DP=1e-3*(pow(G,2e0)/(rho_1*2e0))*F*l1/d1;
	HPo->P = HPi.P - DP;
	}
	else
	{
	HPo->P=HPi.P-dP;
	}

	HPo->H=HPi.H-Q_loss/mr;	
}

void
HotGasLine_Lumped(HP HPi,double mr,HP *HPo,MASS *m,void *Params)
{
	FILE *fp = fopen("./InputDoc/acmodel.dat","r");
	if(!PosComponentData(fp,DISCHARGELINE,1)) {
		errorLog.Add("DISCHARGELINE","Parameters not found");
		return;
		}
	
	double l=GetDouble(fp);//length
	double d=GetDouble(fp);//inside diameter
	double OD=GetDouble(fp);//outside diameter
	double dQ = GetDouble(fp);//heat transfer rate in the discharge line
	double dP = GetDouble(fp);//differential pressure drop across the discharge line
	double dV = GetDouble(fp);//additional volume connected with the discharge line, usually refered to inlet header of the condenser
	fclose(fp) ;

	const double Q_loss1 = dQ;
	const double DP1 = dP;
	
	HPo->P=HPi.P-DP1;
	HPo->H=HPi.H-Q_loss1/mr;					

	m->m=0;
	m->V=0;
	return;
	
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
	
	FILE *fp = fopen("./InputDoc/acmodel.dat","r");
	if(!PosComponentData(fp,LIQUIDLINE,1)) {
		errorLog.Add("liquid line","Parameters not found");
		return;
		}
	double l=GetDouble(fp);//length
	double d=GetDouble(fp);//inside diameter
	double OD=GetDouble(fp);//outside diameter
	double dQ =GetDouble(fp);//heat tranfer rate in the liquild line
	double dP =GetDouble(fp);//pressure drop across the liquid line
	double dV =GetDouble(fp);//additional inner volume connected with the liquid line, usually refered to the outlet header of the condenser
	fclose(fp) ;

	TXP TXP_prop={0,0,0};

	const double A = pi/4.0*d*d;			// cross section area 
	const double G = mr/A;				//mass flux
	double V=0;
	V = l*A+dV;

	TXP_prop=HPtoTXP(HPi);
	if(errorLog.IsError()) {
	errorLog.Add("LiquidLine","HPtoTXP");
	}

	const double rho=1/PropertyTXPth(VOL,TXP_prop);
	if(errorLog.IsError()) {
	errorLog.Add("LiquidLine","rho");
	}

	m->m = V*rho;
	m->V = V;

	const double Q_loss=dQ;

	//Pressure drop
	if(dP<0)//calculate the pressure drop
	{
	
	TXP_prop.X=0;
	double MU=PropertyTXPtr(VISC,TXP_prop);
	if(errorLog.IsError()) {
	errorLog.Add("SuctionLine","MU");
	}

	double RE=G*d/MU;
	double F=0.184e0*pow(RE,(-0.2e0));
	double DP=1e-3*(pow(G,2e0)/(rho*2e0))*F*l/d;
	HPo->P = HPi.P - DP;
	}
	else
	{
	HPo->P=HPi.P-dP;
	}

	HPo->H=HPi.H-Q_loss/mr;
	
	// liquid line dryer
	fp = fopen("./InputDoc/acmodel.dat","r");
	if(!PosComponentData(fp,FILTER,1)) {
		errorLog.Add("filter drier","Parameters not found");
		return;
		}
	double fdV=GetDouble(fp);//inner volume of the drier/filter
	fclose(fp) ;

	TXP_prop=HPtoTXP(*HPo);
	if(errorLog.IsError()) {
	errorLog.Add("LiquidLine","HPtoTXP2");
	}

	const double rho_2=1/PropertyTXPth(VOL,TXP_prop);
	if(errorLog.IsError()) {
	errorLog.Add("LiquidLine","rho2");
	}
	
	m->m = m->m +fdV*rho_2;
	m->V = m->V+fdV;

	return;
}

void
LiquidLine_Lumped(HP HPi,double mr,HP *HPo,MASS *m)
{
	FILE *fp = fopen("./InputDoc/acmodel.dat","r");
	if(!PosComponentData(fp,LIQUIDLINE,1)) {
		errorLog.Add("liquid line","Parameters not found");
		return;
		}
	double l=GetDouble(fp);//length
	double d=GetDouble(fp);//inside diameter
	double OD=GetDouble(fp);//outside diameter
	double dQ =GetDouble(fp);//heat tranfer rate in the liquild line
	double dP =GetDouble(fp);//pressure drop across the liquid line
	double dV =GetDouble(fp);//additional inner volume connected with the liquid line, usually refered to the outlet header of the condenser
	fclose(fp) ;

	const double Q_loss=dQ;
	const double DP = dP;
	HPo->P=HPi.P-DP;
	HPo->H=HPi.H-Q_loss/mr;
	m->m=0;
	m->V=0;
	return;
}

/****************************************************************
Distribution capililary stube models downstream of condenser or evaporator model, 
as there are relatively small pressure drops in these capiliary tubes, there are 
regarded as connection tubes rather than expansion devices   
****************************************************************/
void
Distribution_CapTubes(HP HPi,double mr,HP *HPo,MASS *m, double L, double D, int N)
{
	const double dL = L/N;
	double DP =0;
	const double A = pi/4.0*D*D;// cross section area 
	const double G = mr/A;//mass flux
	int i=0;
	TXP TXPi={0,0,0};
	
	
	for(i=0;i<N;i++)
	{
		TXPi=HPtoTXP(HPi);
		if(errorLog.IsError()) {
		errorLog.Add("Distribution_CapTubes","HPtoTXP");
		}

		if(TXPi.X>0.05&&TXPi.X<0.99)//two-phase pressure drop
		{
			PreAcc Preacc;
			Preacc.DP_FR=FricDPTP_Smooth(TXPi,G,D, dL);
			if(errorLog.IsError()) {
			errorLog.Add("Distribution_CapTubes","FricDPTP_Smooth_0");
			}

			Preacc.G=G;
			Preacc.H_OUT=HPi.H;
			Preacc.P_IN=TXPi.P;
			Preacc.X_IN=TXPi.X;
			const double DP_ACC=Zbrent(10,-10,GET_PreAcc,1e-7,&Preacc);//two-phase acceleration pressure drop
			if(errorLog.IsError()) {
			errorLog.Add("Distribution_CapTubes", "(Preacc)");
			return;}

			HPi.P = HPi.P - FricDPTP_Smooth(TXPi,G,D, dL) - DP_ACC;
			if(errorLog.IsError()) {
			errorLog.Add("Distribution_CapTubes","FricDPTP_Smooth_1");
			}
		}
		else
		{
			HPi.P = HPi.P - FricDPSP_Smooth(TXPi,G,D, dL);
			if(errorLog.IsError()) {
			errorLog.Add("Distribution_CapTubes","FricDPSP_Smooth");
			}
		}
	}

	m->V=0;
	m->m=0;
	HPo->P=HPi.P;
	HPo->H=HPi.H;
}