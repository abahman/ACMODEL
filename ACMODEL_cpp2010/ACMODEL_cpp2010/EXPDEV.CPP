#include "stdafx.h"

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mass.h"
#include "ref.h"
#include "r22.h"
#include "expdev.h"
#include "numeric.h"
#include "fanno.h"
#include "getdata.h"
#include "errorlog.h"
#include "spline1d.h"
#include "spline2d.h"

#include "outelem.h"

extern ErrorLog errorLog;

extern RefPropTPthT reftplthT;

#define K 0.15
#define filenumbers 10

//#define DEBUG_EXPDEV
//#define DEBUG_EXPDEV0

// Bo Shen included the short tube model in the ExpDev code.  This #define needs
// to be included to use this feature.
//#define INCLUDE_SHORTTUBE_IN_EXPDEV

int expDevN=0;

// Local structures
struct GP {
	double P4,Gr2,H1,Gr,D,L;
	double mr;//B.S., mass flow rate
	double As;//B.S., short tube cross section area
	double L_c;//B.S., chamfer length
	double Ang_c;//B.S., chamfer angle
	double Dadj;
};
struct IV {
	double D,L,Dadj;
	int type;
	double L_c;//B.S., chamfer length
	double Ang_c;//B.S., chamfer angle
	double BranNum;
};
struct ExpDev0P {
	int i,nOccur;
	const char *filename;
	double *Prms,H1,P5,mr,P1;
};

struct Short_tube
{
	HP HP_up;
	HP HP_down;
	const char *filename;
	double *Prms,mr;
};

struct TxV_Para
{
	HP HP_up;
	HP HP_down;
	const char *filename;
	double Tsup;
	double *Prms,mr;
};


extern RefPropTPthP reftpvthP;//B.S. add

extern OutputDataElement gOutData;

// Local functions
double IntRho(double,double,double,double,double*,double*);
double P2Func(double,void*);
double ExpDev0Func(double,void*);
double Shorttube0Func(double,void*);
double TxV0Func(double,void*);
double VolumeALL(TXP,double,double);//B.S. 

//----------------------------------------------------------
//---------------------------------------------------------
//----------------------------------------------------------
/****************************************************
B.S., short tube expansion shenadd from O'Neal (2003)
A Mass Flowrate Correlation for Refrigerants and Refrigerant Mixtures, Journal of HVAC
***************************************************/
double Short_Tube_FlowDirec(
				const char *filename,
				HP HP_up,//upstream/inlet state
				HP HP_down,//downstream/outlet state
				double *Prms
				)

{
	static IV Q;
	TXP	TXP_prop={0,0,0};

	FILE *fp = fopen(filename,"r");
	if(!PosComponentData(fp,EXPDEV,1)) {
		errorLog.Add("Short-tube","Parameters not found");
		return 0;
		}
	Q.D=GetDouble(fp);//inside diameter
	Q.L=GetDouble(fp);//length
	Q.Dadj=GetDouble(fp);//adjusting the inside diameter
	Q.L_c=GetDouble(fp);//B.S., new parameter, chamfered length
	Q.Ang_c=GetDouble(fp);//B.S., new parameter, chamfered angle
	Q.BranNum = GetDouble(fp);//B.S., parallelled expansion device
	fclose(fp) ;
	
	const double pi = acos(-1.0);
	const double L = Q.L;//length of the short-tube
	const double L_c = Q.L_c;//chamfered length
	const double D = Q.D;//inside diameter of the short-tube
	const double Ang_c = Q.Ang_c;//chamfered angle
	const double A_s = pi/4*D*D;
	
#ifdef _R22
	const double P_c=4976;//[kPa], critical pressure of R-22
	const double T_c=369.2;//[K], critical temperature of R-22
#endif

#ifdef _R410A
	const double P_c=4963.0;// critical pressure of R410A
	const double T_c=345.65;// [K] critical temperature of R410A
#endif

	// orifice adjustment parameter
	// orifice adjustment parameter
	// orifice adjustment parameter
	double C_c = Q.Dadj;
	
	// Calculate the upstream state in TXP for HP format
	// Calculate the upstream state in TXP for HP format
	// Calculate the upstream state in TXP for HP format
	
	TXP TXP_up = HPtoTXP(HP_up);//upstream state
	if(errorLog.IsError()) {
		errorLog.Add("Short_Tube_FlowDirec", "TXP_up");
		return 0;
	}

	int stateDec = -1;

	if(TXP_up.X>0.999) 	{

		// superheated upstream state
		stateDec = 0;

		errorLog.Add("Short_Tube_FlowDirec", "Upstream state superheated.");
		return 0;

	}

	// Calculate the mass flow rate through the valve.
	// Calculate the mass flow rate through the valve.
	// Calculate the mass flow rate through the valve.


	const double P_up = HP_up.P;
	const double T_up = TXP_up.T;
	double P_sat =reftplthT.Psat(T_up);//P_sat corresponding to upstream temperature
	if(errorLog.IsError()) {
		errorLog.Add("Short_Tube_FlowDirec","P_sat");
		return 0;
	}


	TXP_prop.P=P_up;
	TXP_prop.X=0;

	const double T_sat = PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(P_up);//T_sat corresponding to the upstream pressure
	if(errorLog.IsError()) {
		errorLog.Add("Short_Tube_FlowDirec","T_sat");
		return 0;
	}


	const double T_sub=T_sat-T_up;

	TXP_prop.P=P_sat;
	TXP_prop.X=0;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

	const double rho_f=1/PropertyTXPth(VOL,TXP_prop);//1/reftplthP.v(P_sat);//upstream saturated liquid density
	if(errorLog.IsError()) {
		errorLog.Add("Short_Tube_FlowDirec","vl");
		return 0;
	}

	const double rho_g = 1/reftpvthP.v(P_sat);//upstream saturated gas density,
	if(errorLog.IsError()) {
		errorLog.Add("Short_Tube_FlowDirec","vv");
		return 0;
	}
		

	//non-dimensional groups
	const double pi_3=(P_up-P_sat)/P_c;
	const double pi_6=rho_g/rho_f;
	const double pi_9=T_sub/T_c;
	const double pi_10=L/D;

	//coeffcients
	const double a1=3.8811e-1;
	const double a2=1.1427e1;
	const double a3=-1.4194e1;
	const double a4=1.0703e0;
	const double a5=-9.1928e-2;
	const double a6=2.1425e1;
	const double a7=-5.8195e2;
	
	
	double mr =0;
	double C_tp=1;//two-phase flow rate adjustment
	//single-phase flow rate
	const double pi_1=(a1+a2*pi_3+a3*pi_9+a4*pi_6+a5*log(pi_10))/(1+a6*pi_3+a7*pi_9*pi_9);
	const double G = pi_1*pow((rho_f*P_c*1000),0.5);

	mr=G*A_s;
	
	if(TXP_up.X<0.000001) {

		// subcooled upstream state
		stateDec = 1;
		mr = mr;

	} else {

		// two-phase upstream state
		stateDec = 2;
				
		const double x_up = TXP_up.X;
		const double rho_mup=1/((1-x_up)/rho_f+x_up/rho_g);

		//non-dimensional groups
		const double tp6=rho_mup/rho_f;
		const double tp35=(P_c-P_sat)/(P_c);
		const double tp32=(P_c-P_up)/(P_c);
		const double tp27=L/D;
		const double tp34=x_up/(1-x_up)*pow((rho_f/rho_g),0.5);
		const double tp28=P_up/P_c;

		//coeffcients
		const double b1=1.1831e0;
		const double b2=-1.468e0;
		const double b3=-1.5285e-1;
		const double b4=-1.4639e1;
		const double b5=9.8401e0;
		const double b6=-1.9798e-2;
		const double b7=-1.5348e0;
		const double b8=-2.0533e0;
		const double b9=-1.7195e1;

		const double numer = (b1*tp6+b2*pow(tp6,2.0)+b3*pow(log(tp6),2.0)+b4*pow(log(tp35),2.0)+b5*pow(log(tp32),2.0)+b6*pow(log(tp27),2.0));
		const double deno = (1+b7*tp6+b8*tp34+b9*pow(tp28,3.0));
		C_tp= numer/ deno;
		if(C_tp>1) C_tp=1;//since C_tp>1 is not right.
		mr = mr*C_tp;//correct the mass flow rate by two-phase entrance
	}

	mr=mr*C_c*Prms[0];//adjust the mass flow rate via adjustment factor related with geometry
	if(mr<0) {	
		char msg[256];
		sprintf(msg,"mr<0: mr=%lf, TXP_up.X=%lf, stateDec=%d",mr,TXP_up.X,stateDec);
		errorLog.Add("Short_Tube_FlowDirec",msg);
		return 0;
	}

	mr = mr*Q.BranNum;

	return mr;
}

double Short_Tube_0(
				const char *filename,
				HP HP_up,//upstream/inlet state
				HP HP_down,//downstream/outlet state
				double mr,
				double *Prms
				)
{
	Short_tube P;
	P.HP_up = HP_up;
	P.HP_down = HP_down;
	P.Prms = Prms;
	P.filename = filename;
	P.mr = mr;

	Prms[0] = FindZero(Prms[0],Shorttube0Func,1e-7,&P);
	if(errorLog.IsError()) {
		errorLog.Add("Short_Tube_0","(1)");
		return 0;
	}

	if(errorLog.IsError())
	{
	errorLog.ClearError("Sort_tube_0","toBisection1DSpline");
	const double Xlo=0.3; const double Xhi=5.0;
	Prms[0]=Bisection1DSpline(Xlo,Xhi,Shorttube0Func,1e-7,&P);
	}

	if(errorLog.IsError()) {
		errorLog.Add("Short_Tube_0","(1)");
		return 0;
	}

	Short_Tube_FlowDirec(filename,HP_up,HP_down,Prms);
	if(errorLog.IsError()) {
		errorLog.Add("Short_Tube_0","(2)");
		return 0;
	}
	return 0;
}

double Shorttube0Func(double Adj,void *Params)
{

	Short_tube *P = (Short_tube*)Params;
	P->Prms[0] = Adj;
	const double mr=Short_Tube_FlowDirec(P->filename,P->HP_up,P->HP_down,P->Prms);

	if(errorLog.IsError()) {
		errorLog.Add("Shorttube0Func");
		return 0;
	}

	return (mr-P->mr)/fabs(mr+P->mr);
}

/****************************************************
B.S., short tube expansion shenadd from O'Neal (2003)
A Mass Flowrate Correlation for Refrigerants and Refrigerant Mixtures, Journal of HVAC
***************************************************/
double Short_Tube_R410A(
				const char *filename,
				HP HP_up,//upstream/inlet state
				HP HP_down,//downstream/outlet state
				double *Prms
				)

{
	static IV Q;
	TXP	TXP_prop={0,0,0};

	FILE *fp = fopen(filename,"r");
	if(!PosComponentData(fp,EXPDEV,1)) {
		errorLog.Add("Short-tube","Parameters not found");
		return 0;
		}
	Q.D=GetDouble(fp);//inside diameter
	Q.L=GetDouble(fp);//length
	Q.Dadj=GetDouble(fp);//adjusting the inside diameter
	Q.L_c=GetDouble(fp);//B.S., new parameter, chamfered length
	Q.Ang_c=GetDouble(fp);//B.S., new parameter, chamfered angle
	Q.BranNum = GetDouble(fp);//B.S., parallelled expansion device
	fclose(fp) ;
	
	const double pi = acos(-1.0);
	const double L = Q.L;//length of the short-tube
	const double L_c = Q.L_c;//chamfered length
	const double D = Q.D;//inside diameter of the short-tube
	const double Ang_c = Q.Ang_c;//chamfered angle
	const double A_s = pi/4*D*D;

	const double P_c=4949.65;// critical pressure of R410A
	const double T_c=345.65;// [K] critical temperature of R410A


	// orifice adjustment parameter
	// orifice adjustment parameter
	// orifice adjustment parameter
	double C_c = Q.Dadj;
	
	// Calculate the upstream state in TXP for HP format
	// Calculate the upstream state in TXP for HP format
	// Calculate the upstream state in TXP for HP format
	
	TXP TXP_up = HPtoTXP(HP_up);//upstream state
	if(errorLog.IsError()) {
		errorLog.Add("Short_Tube_FlowDirec", "TXP_up");
		return 0;
	}

	int stateDec = -1;

	if(TXP_up.X>0.999) 	{

		// superheated upstream state
		stateDec = 0;

		errorLog.Add("Short_Tube_FlowDirec", "Upstream state superheated.");
		return 0;

	}

	// Calculate the mass flow rate through the valve.
	// Calculate the mass flow rate through the valve.
	// Calculate the mass flow rate through the valve.


	const double P_up = HP_up.P;
	const double T_up = TXP_up.T;
	double P_sat =reftplthT.Psat(T_up);//P_sat corresponding to upstream temperature
	if(errorLog.IsError()) {
		errorLog.Add("Short_Tube_FlowDirec","P_sat");
		return 0;
	}

	const double P_down=HP_down.P;

	TXP_prop.P=P_up;
	TXP_prop.X=0;

	const double T_sat = PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(P_up);//T_sat corresponding to the upstream pressure
	if(errorLog.IsError()) {
		errorLog.Add("Short_Tube_FlowDirec","T_sat");
		return 0;
	}


	const double T_sub=T_sat-T_up;

	double rho_up=1/PropertyTXPth(VOL,TXP_up);
	if(errorLog.IsError()) {
		errorLog.Add("Short_Tube_FlowDirec","rho_up");
		return 0;
	}



	TXP_prop.P=P_up;//P_sat;*****
	TXP_prop.X=0;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

	const double rho_f=1/PropertyTXPth(VOL,TXP_prop);//1/reftplthP.v(P_sat);//upstream saturated liquid density
	if(errorLog.IsError()) {
		errorLog.Add("Short_Tube_FlowDirec","vl");
		return 0;
	}

	const double rho_g = 1/reftpvthP.v(P_up);//upstream saturated gas density****P_sat
	if(errorLog.IsError()) {
		errorLog.Add("Short_Tube_FlowDirec","vv");
		return 0;
	}

	if(TXP_up.X>0) rho_up=rho_f;

	//coeffcients
	double b1=0,b2=0,b3=0,b4=0,b5=0,b6=0,b7=0,b8=0,b9=0,b10=0,b11=0;
	double a1=0,a2=0,a3=0,a4=0;
	if(0)//2.2% lubricant
	{
	b1=1.050104183;
	b2=6.305986547;
	b3=0.099138818;
	b4=-0.045626106;
	b5=0.958459297;
	b6=-0.254071783;
	b7=0.137198955;
	b8=-0.276516186;
	b9=-0.014589768;
	b10=2.5121908;
	b11=0.13087558;
	a1=1.427618112;
	a2=0.530751112;
	a3=-0.365456266;
	a4=0.018669938;
	}
	else
	{
	b1=0.874805831;
	b2=3.131470913;
	b3=-0.214726407;
	b4=0.083394737;
	b5=0.901559277;
	b6=-0.020574536;
	b7=0.944446846;
	b8=-0.418400083;
	b9=-0.025322802;
	b10=2.33507746;
	b11=0.068890172;
	a1=3.693038038;
	a2=0.120175996;
	a3=0.194241638;
	a4=0.022577667;
	}

	const double PRA=P_up/P_c;
	double SUBC=T_sub/T_c;
	if(SUBC<0) SUBC=1e-20;
	const double LD=L/D;
	const double D_ref=1.524e-3;
	const double DR=D/D_ref;
	const double EVAP=(P_c-P_down)/P_c;
	
	const double P_f=P_sat*(b1+b2*pow(PRA,b3)*pow(LD,b4)*pow(SUBC,b5)+b6*pow(PRA,b7)+b8*exp(b9*DR*pow(LD,b10))+b11*EVAP);
	
	double mr =0;
	double C_tp=1;//two-phase flow rate adjustment
	//single-phase flow rate
	
	const double r_c=1.2960e10;
	const double G=pow(2*r_c*rho_up*(P_up-P_f),0.5);

	mr=G*A_s/3600;
	
	if(TXP_up.X<0.000001) {

		// subcooled upstream state
		stateDec = 1;
		C_tp=1.0;
		mr = mr;

	} else {

		// two-phase upstream state
		stateDec = 2;
		const double x_up=TXP_up.X;
		const double Y=x_up/(1-x_up)*pow(rho_f/rho_g,0.5);
	
		C_tp= 1/((1+a1*x_up)*(1+a2*pow(LD,a3)*pow(Y,a4*log(LD))));
		if(C_tp>1||C_tp<0) C_tp=1;//since C_tp>1 is not right.
		mr = mr*C_tp;//correct the mass flow rate by two-phase entrance
	}

	mr=mr*C_c*Prms[0];//adjust the mass flow rate via adjustment factor related with geometry
	if(mr<0) {	
		char msg[256];
		sprintf(msg,"mr<0: mr=%lf, TXP_up.X=%lf, stateDec=%d",mr,TXP_up.X,stateDec);
		errorLog.Add("Short_Tube_FlowDirec",msg);
		return 0;
	}

	mr = mr* Q.BranNum;

	return mr;
}

/****************************************************
B.S., short tube expansion shenadd from O'Neal (2003)
A Mass Flowrate Correlation for Refrigerants and Refrigerant Mixtures, Journal of HVAC
it is the same as the above one, just for validating the cases of lubricant presence, 
and do comparison in mainfunc_Validating()
***************************************************/
double Short_Tube_R410A_Oil(
				const char *filename,
				HP HP_up,//upstream/inlet state
				HP HP_down,//downstream/outlet state
				double *Prms
				)

{
	static IV Q;
	TXP	TXP_prop={0,0,0};

	FILE *fp = fopen(filename,"r");
	if(!PosComponentData(fp,EXPDEV,1)) {
		errorLog.Add("Short-tube","Parameters not found");
		return 0;
		}
	Q.D=GetDouble(fp);//inside diameter
	Q.L=GetDouble(fp);//length
	Q.Dadj=GetDouble(fp);//adjusting the inside diameter
	Q.L_c=GetDouble(fp);//B.S., new parameter, chamfered length
	Q.Ang_c=GetDouble(fp);//B.S., new parameter, chamfered angle
	Q.BranNum = GetDouble(fp);//B.S., parallelled expansion device
	fclose(fp) ;
	
	const double pi = acos(-1.0);
	const double L = Q.L;//length of the short-tube
	const double L_c = Q.L_c;//chamfered length
	const double D = Q.D;//inside diameter of the short-tube
	const double Ang_c = Q.Ang_c;//chamfered angle
	const double A_s = pi/4*D*D;

	const double P_c=4949.65;// critical pressure of R410A
	const double T_c=345.65;// [K] critical temperature of R410A

	// orifice adjustment parameter
	// orifice adjustment parameter
	// orifice adjustment parameter
	double C_c = Q.Dadj;
	
	// Calculate the upstream state in TXP for HP format
	// Calculate the upstream state in TXP for HP format
	// Calculate the upstream state in TXP for HP format
	
	TXP TXP_up = HPtoTXP(HP_up);//upstream state
	if(errorLog.IsError()) {
		errorLog.Add("Short_Tube_FlowDirec", "TXP_up");
		return 0;
	}

	int stateDec = -1;

	if(TXP_up.X>0.999) 	{

		// superheated upstream state
		stateDec = 0;

		errorLog.Add("Short_Tube_FlowDirec", "Upstream state superheated.");
		return 0;

	}

	// Calculate the mass flow rate through the valve.
	// Calculate the mass flow rate through the valve.
	// Calculate the mass flow rate through the valve.


	const double P_up = HP_up.P;
	const double T_up = TXP_up.T;
	double P_sat =reftplthT.Psat(T_up);//P_sat corresponding to upstream temperature
	if(errorLog.IsError()) {
		errorLog.Add("Short_Tube_FlowDirec","P_sat");
		return 0;
	}

	const double P_down=HP_down.P;

	TXP_prop.P=P_up;
	TXP_prop.X=0;

	const double T_sat = PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(P_up);//T_sat corresponding to the upstream pressure
	if(errorLog.IsError()) {
		errorLog.Add("Short_Tube_FlowDirec","T_sat");
		return 0;
	}


	const double T_sub=T_sat-T_up;

	double rho_up=1/PropertyTXPth(VOL,TXP_up);
	if(errorLog.IsError()) {
		errorLog.Add("Short_Tube_FlowDirec","rho_up");
		return 0;
	}



	TXP_prop.P=P_up;//not p_sat
	TXP_prop.X=0;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

	const double rho_f=1/PropertyTXPth(VOL,TXP_prop);//1/reftplthP.v(P_sat);//upstream saturated liquid density
	if(errorLog.IsError()) {
		errorLog.Add("Short_Tube_FlowDirec","vl");
		return 0;
	}

	const double rho_g = 1/reftpvthP.v(P_up);//upstream saturated gas density, not P_sat
	if(errorLog.IsError()) {
		errorLog.Add("Short_Tube_FlowDirec","vv");
		return 0;
	}

	if(TXP_up.X>0) rho_up=rho_f;

	//coeffcients
	double b1=0,b2=0,b3=0,b4=0,b5=0,b6=0,b7=0,b8=0,b9=0,b10=0,b11=0;
	double a1=0,a2=0,a3=0,a4=0;
	if(1)//2.2% lubricant
	{
	b1=1.050104183;
	b2=6.305986547;
	b3=0.099138818;
	b4=-0.045626106;
	b5=0.958459297;
	b6=-0.254071783;
	b7=0.137198955;
	b8=-0.276516186;
	b9=-0.014589768;
	b10=2.5121908;
	b11=0.13087558;
	a1=1.427618112;
	a2=0.530751112;
	a3=-0.365456266;
	a4=0.018669938;
	}
	else
	{
	b1=0.874805831;
	b2=3.131470913;
	b3=-0.214726407;
	b4=0.083394737;
	b5=0.901559277;
	b6=-0.020574536;
	b7=0.944446846;
	b8=-0.418400083;
	b9=-0.025322802;
	b10=2.33507746;
	b11=0.068890172;
	a1=3.693038038;
	a2=0.120175996;
	a3=0.194241638;
	a4=0.022577667;
	}

	const double PRA=P_up/P_c;
	double SUBC=T_sub/T_c;
	if(SUBC<0) SUBC=1e-20;
	const double LD=L/D;
	const double D_ref=1.524e-3;
	const double DR=D/D_ref;
	const double EVAP=(P_c-P_down)/P_c;
	
	const double P_f=P_sat*(b1+b2*pow(PRA,b3)*pow(LD,b4)*pow(SUBC,b5)+b6*pow(PRA,b7)+b8*exp(b9*DR*pow(LD,b10))+b11*EVAP);
	
	double mr =0;
	double C_tp=1;//two-phase flow rate adjustment
	//single-phase flow rate
	
	const double r_c=1.2960e10;
	const double G=pow(2*r_c*rho_up*(P_up-P_f),0.5);

	mr=G*A_s/3600;
	
	if(TXP_up.X<0.000001) {

		// subcooled upstream state
		stateDec = 1;
		C_tp=1.0;
		mr = mr;

	} else {

		// two-phase upstream state
		stateDec = 2;
		const double x_up=TXP_up.X;
		const double Y=x_up/(1-x_up)*pow(rho_f/rho_g,0.5);
	
		C_tp= 1/((1+a1*x_up)*(1+a2*pow(LD,a3)*pow(Y,a4*log(LD))));
		if(C_tp>1||C_tp<0) C_tp=1;//since C_tp>1 is not right.
		mr = mr*C_tp;//correct the mass flow rate by two-phase entrance
	}

	mr=mr*C_c*Prms[0];//adjust the mass flow rate via adjustment factor related with geometry
	if(mr<0) {	
		char msg[256];
		sprintf(msg,"mr<0: mr=%lf, TXP_up.X=%lf, stateDec=%d",mr,TXP_up.X,stateDec);
		errorLog.Add("Short_Tube_FlowDirec",msg);
		return 0;
	}

	mr = mr*Q.BranNum;

	return mr;
}


/****************************************************
B.S., nonlinear TxV model from Haorong Li paper
***************************************************/
double TxV_Nonlinear(
				const char *filename,
				HP HP_up,//upstream/inlet state
				HP HP_down,//downstream/outlet state
				double Tsup,
				double *Prms
				)
{
	static IV Q;//borrowed
	TXP TXP_prop={0,0,0};

	FILE *fp = fopen(filename,"r");
	if(!PosComponentData(fp,EXPDEV,1)) {
		errorLog.Add("TxV","Parameters not found");
		return 0;
		}
	
	Q.D = GetDouble(fp);
	const double Tsh_static = Q.D;//4;
	Q.L = GetDouble(fp);
	const double Tsh_max = Q.L;//6;
	Q.L_c = GetDouble(fp);
	const double C= Q.L_c;//1.2656;
	Q.Dadj = GetDouble(fp);
	const double Adj= Q.Dadj;//0.7630;
	double A=(Tsup-Tsh_static)/Tsh_max;
	if(A>1.0) A=1.0;
	fclose(fp) ;

	TXP TXP_up = HPtoTXP(HP_up);//upstream state
	if(errorLog.IsError()) {
		errorLog.Add("TxV", "TXP_up");
		return 0;
	}

	int stateDec = -1;

	if(TXP_up.X>0.999) 	
	{

	// superheated upstream state
	stateDec = 0;
	errorLog.Add("TxV", "Upstream state superheated.");
	return 0;
	}
	else if(TXP_up.X>0.0)
	{
	// saturated upstream state
	stateDec = 0;
	}
	else
	{
	// subcooled upstream state
	stateDec = -1;
	}
	
//	const double rho_up=1/PropertyTXPth(VOL,TXP_up);;//upstream saturated liquid density
	TXP_prop.P=TXP_up.P;
	TXP_prop.X=0;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);
	const double rho_up=1/PropertyTXPth(VOL,TXP_prop);//1/reftplthP.v(TXP_up.P);

	if(errorLog.IsError()) {
		errorLog.Add("TxV", "rho_f");
		return 0;
	}

	const double P_up=HP_up.P;
	const double P_down=HP_down.P;

	double mr = C*(2*A-A*A)*pow(rho_up*(P_up-P_down)*1e-2,0.5)*1e-3;
	
	mr=mr*Adj*Prms[0];//adjust the mass flow rate via adjustment factor related with geometry
	if(mr<0) {	
		char msg[256];
		sprintf(msg,"mr<0: mr=%lf, TXP_up.X=%lf, stateDec=%d",mr,TXP_up.X,stateDec);
		errorLog.Add("TxV",msg);
		return 0;
	}

	return mr;
}

double TxV_0(	const char *filename,
				HP HP_up,//upstream/inlet state
				HP HP_down,//downstream/outlet state
				double Tsup,
				double mr,
				double *Prms
				)
{
	TxV_Para P;
	P.HP_up = HP_up;
	P.HP_down = HP_down;
	P.Prms = Prms;
	P.filename = filename;
	P.Tsup=Tsup;
	P.mr = mr;

	Prms[0] = FindZero(Prms[0],TxV0Func,1e-7,&P);
	if(errorLog.IsError()) {
		errorLog.Add("TxV_0","(1)");
		return 0;
	}

	TxV_Nonlinear(filename,HP_up,HP_down,Tsup,Prms);
	if(errorLog.IsError()) {
		errorLog.Add("TxV_0","(2)");
		return 0;
	}
	return 0;
}

double TxV0Func(double Adj,void *Params)
{

	TxV_Para *P = (TxV_Para*)Params;
	P->Prms[0] = Adj;
	const double mr=TxV_Nonlinear(P->filename,P->HP_up,P->HP_down,P->Tsup,P->Prms);

	if(errorLog.IsError()) {
		errorLog.Add("TxV0Func");
		return 0;
	}
	return (mr-P->mr)/fabs(mr+P->mr);
}

//--------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
//----------------------------------------------------------------------------------

//the following subroutines are not used currently

/********************************************************************
Fixed orifice expansion device model.  The model predicts inlet
pressure (P1) from stagnation enthalpy (H1), outlet pressure (P5),
and mass flow rate (mr).  The file (fn) contains physical dimensions
of the expansion device.
States:   (1) inlet state
		  (2) state after enterence losses
		  (3) flash point (liquid / two phase boundary)
		  (4) outlet pressure
		  (5) exit back pressure
********************************************************************/
void ExpDev(const char *filename,int nOccur,double H1,double P5,double mr,double *P1,MASS *m,double *Prms)
{
	int i,j;
	double v,P2,P4;
	double Pvisc,mu,Re,f,Irho,L1,L2,rhoa,rhob,a;
	static int Occur[filenumbers];
	static IV Q[filenumbers];
	static GP GP;
	TXP TXP_prop={0,0,0};

	j=-1;
	for(int i=0;i<expDevN;i++) if(nOccur==Occur[i]) {j=i; break;}
	if (j<0) {
		FILE *fp = fopen(filename,"r");
		if(!PosComponentData(fp,EXPDEV,nOccur)) {
			errorLog.Add("ExpDev","Parameters not found");
			return;
		}
		j=expDevN;
		Q[j].D=GetDouble(fp);
		Q[j].L=GetDouble(fp);
		Q[j].Dadj=GetDouble(fp);
		fclose(fp) ;
		Occur[expDevN]=nOccur;
		expDevN++;

		#ifdef DEBUG_EXPDEV
		{
		char fn[128];
		sprintf(fn,"ExpDev%dDebug.txt",nOccur);
		FILE* fp = fopen(fn,"w");
		if(fp) {
			fprintf(fp,"Starting...\n");
			fclose(fp);
		}
		}
		#endif
	}

	#ifdef DEBUG_EXPDEV
	{
	char fn[128];
	sprintf(fn,"ExpDev%dDebug.txt",nOccur);
	FILE* fp = fopen(fn,"a");
	if(fp) {
		fprintf(fp,"H1=%.1lf, P5=%.1lf, mr=%.6lf\n",H1,P5,mr);
		fclose(fp);
	}
	}
	#endif

	GP.D = Q[j].D * Q[j].Dadj * Prms[0];
	GP.L = Q[j].L;

	// Inside cross-sectional area of device.
	static double pi = acos(-1.0);
	const double Ax = pi/4*pow(GP.D,2);

	// refrigerant mass flux.
	const double Gr = mr/Ax;
	const double Gr2 = Gr*Gr;

	// exit pressure.
	const double Pc = Pchoke(Gr,H1);
	if(errorLog.IsError()) {
		errorLog.Add("ExpDev","Pc");
		return;
	}
	if(P5<Pc) P4=Pc; else P4=P5;

	// flash pressure.
	const double P3 = Pflash(Gr,H1);
	if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"P3 (P3=%lf Gr=%lf H1=%lf)",P3,Gr,H1);
		errorLog.Add("ExpDev",str);
		return;
	}

	if(P3>PMAXth) {
		/* two phase enters and exits the orifice. */
		GP.Gr = Gr;
		GP.Gr2 = Gr2;
		GP.H1 = H1;
		GP.P4 = P4;
		P2 = Zbrent(P4,PMAXth,P2Func,1e-7,&GP);
		if(errorLog.IsError()) {
			char str[128];
			sprintf(str,"(P2) P4=%lf P2=%lf",P4,P2);
			errorLog.Add("ExpDev",str);
			return;
		}
		v = PropertyTXPth(VOL,TXPfanno(Gr,H1,P2));
		if(errorLog.IsError()) {
			errorLog.Add("ExpDev","v(1)");
			return;
		}
		*P1 = P2+0.5e-3*(1+K)*Gr2*v;
	} else if(P3>P4) {
		// two phase exits.
		// length of two phase section assuming it starts at P3.
		// Transport property data has limited range.  Approximating
		// viscosity at highest pressure for which transport
		// properties are known if it is greater than the maximum.
		Pvisc = (P3+P4)/2;
		if(Pvisc>PMAXtr) Pvisc=PMAXtr;

		TXP_prop.P=Pvisc;
		TXP_prop.X=0;

		const double Tsat_visc = PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(Pvisc);
		
		TXP_prop.P=Pvisc;
		TXP_prop.T=Tsat_visc;
		TXP_prop.X=0;

		mu = PropertyTXPtr(VISC,TXP_prop);//refsctrPT.mu(Pvisc,Tsat_visc);
		if(errorLog.IsError()) {
			errorLog.Add("ExpDev","mu(1)");
			return;
		}
		Re = Gr*GP.D/mu;
		f = 0.775/pow(Re,0.5)*exp(1/2.4);
		Irho = IntRho(P3,P4,Gr,H1,&rhoa,&rhob);
		if(errorLog.IsError()) {
			char str[128];
			sprintf(str,"Irho (P3=%lf P4=%lf Gr=%lf H1=%lf)",P3,P4,Gr,H1);
			errorLog.Add("ExpDev",str);
			return;
		}
		a = rhoa/rhob;
		if(a<0) {
			errorLog.Add("ExpDev","log domain error");
			return;
		}
		L2 = GP.D/(2*f)*(Irho/Gr2-log(a));
		if(L2>GP.L) {
			// two phase flow through entire tube.
			GP.Gr = Gr;
			GP.Gr2 = Gr2;
			GP.H1 = H1;
			GP.P4 = P4;
			P2 = Zbrent(P4,PMAXth,P2Func,1e-7,&GP);

			// There is a small numerical problem when P2 is extremely
			// close to P4.  The convergence criteria is relaxed a bit
			// to avoid errors in this case.
			if(errorLog.IsError() && (P2-P4)/P4<1e-3) {
				errorLog.ClearError("ExpDev","P2Func(1)");
				if(P2Func(P2,&GP)<1e-3) {
					errorLog.ClearError("ExpDev","P2Func(2)");
				} else {
					errorLog.Add("ExpDev","P2Func(3)");
				}
			}

			if(errorLog.IsError()) {
				char str[128];
				sprintf(str,"(2) H1,P5,mr=(%lf,%lf,%lf) P4=%lf P2=%lf",H1,P5,mr,P4,P2);
				errorLog.Add("ExpDev",str);
				return;
			}
			v = PropertyTXPth(VOL,TXPfanno(Gr,H1,P2));
			if(errorLog.IsError()) {
				errorLog.Add("ExpDev","v(2)");
				return;
			}
			*P1=P2+0.5e-3*(1+K)*Gr2*v;
		} else {
			// two phase flow begins in tube at Pf.
			L1 = GP.L-L2;
			// Transport property data has limited range.  Approximating
			// viscosity at highest pressure for which transport
			// properties are known if it is greater than the maximum.
			Pvisc = P3;
			if(Pvisc>PMAXtr) Pvisc=PMAXtr;

			TXP_prop.P=Pvisc;
			TXP_prop.X=0;

			const double Tsat_visc = PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(Pvisc);

			TXP_prop.P=Pvisc;
			TXP_prop.T=Tsat_visc;
			TXP_prop.X=0;

			mu = PropertyTXPtr(VISC,TXP_prop);//refsctrPT.mu(Pvisc,Tsat_visc);
			if(errorLog.IsError()) {
				errorLog.Add("ExpDev","mu(2)");
				return;
			}

			TXP_prop.P=P3;
			TXP_prop.X=0;
			TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

			v = PropertyTXPth(VOL,TXP_prop);//reftplthP.v(P3);
			if(errorLog.IsError()) {
				errorLog.Add("ExpDev","v(3)");
				return;
			}
			Re = Gr*GP.D/mu;
			if(Re>2000) {
				// f = 0.046*pow(Re,-0.2) ;
				f = pow(1.8*log10(Re)-1.64,-2)/4;
			} else {
				f = 16/Re;
			}
			P2 = P3+2e-3*Gr2*f*v*L1/GP.D;
			*P1 = P2+0.5e-3*(1+K)*Gr2*v;
		}
	} else {
		// liquid exits.
		// Transport property data has limited range.  Approximating
		// viscosity at highest pressure for which transport
		// properties are known if it is greater than the maximum.
		Pvisc = P4;
		if(Pvisc>PMAXtr) Pvisc=PMAXtr;
		TXP_prop.P=Pvisc;
		TXP_prop.X=0;
		const double Tsat_visc = PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(Pvisc);

		TXP_prop.P=Pvisc;
		TXP_prop.T=Tsat_visc;
		TXP_prop.X=0;

		mu = PropertyTXPtr(VISC,TXP_prop);//refsctrPT.mu(Pvisc,Tsat_visc);
		if(errorLog.IsError()) {
			errorLog.Add("ExpDev","mu(4)");
			return;
		}

		TXP_prop.P= P4;
		TXP_prop.X=0;
		TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

		v = PropertyTXPth(VOL,TXP_prop);//reftplthP.v(P4);
		if(errorLog.IsError()) {
			errorLog.Add("ExpDev","v(5)");
			return;
		}
		Re = Gr*GP.D/mu ;
		if(Re>2000) {
			//f=0.046*pow(Re,-0.2) ;
			f=pow(1.8*log10(Re)-1.64,-2)/4;
		} else {
			f = 16/Re;
		}
		P2 = P4+2e-3*Gr2*f*v*GP.L/GP.D;
		*P1 = P2+0.5e-3*(1+K)*Gr2*v;
	}

	/* There is no charge in the expansion device. */
	m->m=0 ;
	m->V=0 ;

	#ifdef DEBUG_EXPDEV
	{
	char fn[128];
	sprintf(fn,"ExpDev%dDebug.txt",nOccur);
	FILE* fp = fopen(fn,"a");
	if(fp) {
		fprintf(fp,"P1=%.1lf\n\n",*P1);
		fclose(fp);
	}
	}
	#endif
}


/********************************************************************
Used with a numerical method in ExpDev() to get the inlet state
pressure (P2) ginve the outlet state pressure (P4).
********************************************************************/
double P2Func(double P2,void *Params)
{
	GP* P = (GP*)Params;
	TXP TXP_prop={0,0,0};

	// Transport property data has limited range.  Approximating
	// viscosity at highest pressure for which transport properties
	// are known if it is greater than the maximum.
	double Pvisc = (P2+P->P4)/2;
	if(Pvisc>PMAXtr) Pvisc=PMAXtr;

	TXP_prop.P=Pvisc;
	TXP_prop.X=0;
	const double Tsat_visc = PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(Pvisc);

	TXP_prop.P=Pvisc;
	TXP_prop.T= Tsat_visc;
	TXP_prop.X=0;

	const double mu = PropertyTXPtr(VISC,TXP_prop);//refsctrPT.mu(Pvisc,Tsat_visc);

	if(errorLog.IsError()) {
		errorLog.Add("P2Func");
		return 0;
	}

	const double Re = P->Gr*P->D/mu;
	const double f = 0.775/pow(Re,0.5)*exp(1/2.4);

	double rho2,rho4;
	const double Irho = IntRho(P2,P->P4,P->Gr,P->H1,&rho2,&rho4);
	if(errorLog.IsError()) {
		errorLog.Add("P2Func");
		return 0;
	}

	const double a = rho2/rho4;
	if(a<0) {
		errorLog.Add("P2Func","log domain error");
		return 0;
	}

	const double E = 2*f*P->L/P->D+log(a)-Irho/P->Gr2;

	return E ;
}

/********************************************************************
Integrates the refrigerant density (Rho) as a function of pressure
from P1 to P2.
********************************************************************/
double IntRho(double P1,double P2,double Gr,double H1,double *rho1,double *rho2)
{
	const int n=100;

	if(P1<P2) {
		errorLog.Add("IntRho","P1<P2");
		return 1.0;
	}

	// pressure step size
	double dP = (P1-P2)/(n-1);

	// inlet density
	TXP TXP = TXPfanno(Gr,H1,P1);
	if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"(rho1) Gr=%lf H1=%lf P1=%lf",Gr,H1,P1);
		errorLog.Add("IntRho",str);
		return -1;
	}

	*rho1 = PropertyTXPth(DENSITY,TXP);
	if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"(rho1) T=%lf X=%lf P=%lf",TXP.T,TXP.X,TXP.P);
		errorLog.Add("IntRho",str);
		return -1;
	}

	/* outlet density */
	TXP = TXPfanno(Gr,H1,P2);
	if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"(rho2) Gr=%lf H1=%lf P2=%lf",Gr,H1,P2);
		errorLog.Add("IntRho",str);
		return -1;
	}
	*rho2 = PropertyTXPth(DENSITY,TXP);
	if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"(rho2) T=%lf X=%lf P=%lf",TXP.T,TXP.X,TXP.P);
		errorLog.Add("IntRho",str);
		return -1;
	}

	/* integral density from P1 to P2 */
	double Irho = (*rho1+*rho2)/2;
	for(int i=2;i<n;i++) {
		double Px = P1-(i-1)*dP;
		TXP=TXPfanno(Gr,H1,Px);
		if(errorLog.IsError()) {
			char str[128];
			sprintf(str,"(Irho) Gr=%lf H1=%lf Px=%lf",Gr,H1,Px);
			errorLog.Add("IntRho",str);
			return -1;
		}
		Irho += PropertyTXPth(DENSITY,TXP);
		if(errorLog.IsError()) {
			char str[128];
			sprintf(str,"(Irho) T=%lf X=%lf P=%lf",TXP.T,TXP.X,TXP.P);
			errorLog.Add("IntRho",str);
			return -1;
		}
	}
	Irho*=1e3*dP;

	return Irho;
}

/********************************************************************
Solves the expansion device model iteratively.  It takes the inlet
pressure as an argument and returns the diameter of the orifice.
********************************************************************/
void ExpDev0(int i,double P1,const char *filename,int nOccur,double H1,double P5,double mr,double *Pb,MASS *m,double *Prms)
{
	ExpDev0P P;
	P.i = i;
	P.Prms = Prms;
	P.filename = filename;
	P.nOccur = nOccur;
	P.H1 = H1;
	P.P5 = P5;
	P.mr = mr;
	P.P1 = P1;

	Prms[i] = FindZero(Prms[i],ExpDev0Func,1e-7,&P);
	if(errorLog.IsError()) {
		errorLog.Add("ExpDev0","(1)");
		return;
	}

	ExpDev(filename,nOccur,H1,P5,mr,Pb,m,Prms);
	if(errorLog.IsError()) {
		errorLog.Add("ExpDev0","(2)");
		return;
	}
}

/********************************************************************
Function used by ExpDev0() to iterate on.
********************************************************************/
double ExpDev0Func(double Adj,void *Params)
{
	ExpDev0P *P = (ExpDev0P*)Params;

	P->Prms[P->i] = Adj;
	double P1;
	MASS m;
	ExpDev(P->filename,P->nOccur,P->H1,P->P5,P->mr,&P1,&m,P->Prms);
	if(errorLog.IsError()) {
		errorLog.Add("ExpDev0Func");
		return 0;
	}

	#ifdef DEBUG_EXPDEV0
	printf("Adj=%.4lf, P1=%.1lf, P1(goal)=%.1lf (H1=%.1lf, P5=%.1lf, mr=%.4lf).\n",Adj,P1,P->P1,P->H1,P->P5,P->mr);
	#endif

	return (P->P1-P1)/P->P1;
}
