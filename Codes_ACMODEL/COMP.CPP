#include "stdafx.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ref.h"
#include "r22.h"
#include "mass.h"
#include "numeric.h"
#include "comp.h"
#include "getdata.h"
#include "errorlog.h"
#include "spline1d.h"
#include "spline2d.h"

extern ErrorLog errorLog;

extern RefPropTPthT reftpvthT;
extern RefPropTPthT reftplthT;
extern LubriProp	lubriPropComp;

int compInit=1;

// Local structures
struct CompP101 {
	HP HP1;
	double v1,P2,npAdj;
	double H2,n;
	// parameters
	double rpm,dispVolume,np;
	double PowerP[10],VolumeEffP[10],Adj[3];
	double InnerVolume,OilCharge,OilTemp;//parameters for free inner volume of the compressor, oil charge amount, and the oil temperature
	int OilType;//oil type that can be related to different oil solubility equation
};

struct CompP102 {
	HP HP1;
	double v1,P2,npAdj;
	double H2,n;
	// parameters
	double np;
	double powerP[11],massFlowP[11],Adj[3];
	double InnerVolume,OilCharge,OilTemp;//parameters for free inner volume of the compressor, oil charge amount, and the oil temperature
	int OilType;//oil type that can be related to different oil solubility equation
};

//B.S.------------------------
//add for ARI standard compressor equations in the form of polynomial with ten terms

struct CompP103 {
	HP HP1;
	double v1,P2,npAdj;
	double H2,n;
	// parameters
	double np;
	//B.S. the following two parameters have no use for ARI compressor equations, but place them here in case.
	double rpm,dispVolume;
	double powerP[10],massFlowP[10],Adj[3];
	double InnerVolume,OilCharge,OilTemp;//parameters for free inner volume of the compressor, oil charge amount, and the oil temperature
	int OilType;//oil type that can be related to different oil solubility equation
};
//--------------------------B.S.

struct Comp0P {
	double H0,P1,P2,Toa,*Prms;
};

// Local functions
double CompFunc101(double,void*);
void Compressor101(CompP101& Q,double H1,double P1,double P2,double Toa,double *mr,HP *HPo,HP *HPi,double *Ei,MASS *m,double *Prms);
double CompFunc102(double,void*);
void Compressor102(CompP102& Q,double H1,double P1,double P2,double Toa,double *mr,HP *HPo,HP *HPi,double *Ei,MASS *m,double *Prms);

// B.S. the following is added for ARI standard
void Compressor103(CompP103& Q,double H1,double P1,double P2,double Toa,double *mr,HP *HPo,HP *HPi,double *Ei,MASS *m,double *Prms);
double Dissolved_Oil(double,double,double,int);//accounting refrigerant dissolved in the oil

/********************************************************************
Hermetically sealed reciprocating compressor model.  Process
modeled as polytropic compression.

States:  (1)-inlet to compressor,
			(2)-outlet of compressor.

Inputs:	H1 (enthalpy state 1, J/kg),
	P1 (pressure state 1, kPa),
	P2 (pressure state 2, kPa),
	Toa (ambient temperature, C), and
	Prms (pointer to array of 3 multiplicative fitting factors for
	volumetric, motor electrical, and polytropic efficiencies
	respectively).  These should all be 1 except when used
	for fitting experimental data.

Outputs:mr (refrigerant mass flow rate, kg/s),
	HPo (refrigerant state at compressor outlet, J/kg, kPa),
	HPi (refrigerant state at compressor inlet, J/kg, kPa),
	Ei (compressor motor power consumption, W),
	m (mass of charge in compressor and total volume, kg, m^3), and
	status (set non-zero to indicate error),
********************************************************************/
void Compressor(const char *filename,
				double H1,
				double P1,
				double P2,
				double Toa,
				double *mr,
				HP *HPo,
				HP *HPi,
				double *Ei,
				MASS *m,
				double *Prms)
{
	static int type;
	static CompP101 Q101;
	static CompP102 Q102;
	// B.S. the following is added for ARI standard
	static CompP103 Q103;
	static int endSign;

	// Load model parameters from file "comp.dat"
	if(compInit) {
		FILE *fp = fopen(filename,"r");
		if(!(type=PosComponentData(fp,COMPRESSOR,1))) {
			errorLog.Add("Compressor","Parameters not found");
			return;
		}
		switch(type) {
			int i;
			case 101: // physical model
				Q101.rpm=GetDouble(fp);
				Q101.dispVolume=GetDouble(fp);
				Q101.np=GetDouble(fp);
				for(i=0;i<5;i++) Q101.PowerP[i]=GetDouble(fp);
				for(i=0;i<4;i++) Q101.VolumeEffP[i]=GetDouble(fp);
				for(i=0;i<3;i++) Q101.Adj[i]=GetDouble(fp);
				Q101.InnerVolume=GetDouble(fp);
				Q101.OilCharge=GetDouble(fp);
				Q101.OilTemp=GetDouble(fp);
				Q101.OilType = (int)GetDouble(fp);
				break;
			case 102: // compressor map model
				Q102.np = GetDouble(fp);
				for(i=0;i<11;i++) Q102.powerP[i] = GetDouble(fp);
				for(i=0;i<11;i++) Q102.massFlowP[i] = GetDouble(fp);
				for(i=0;i<3;i++)  Q102.Adj[i] = GetDouble(fp);
				Q102.InnerVolume=GetDouble(fp);
				Q102.OilCharge=GetDouble(fp);
				Q102.OilTemp=GetDouble(fp);
				Q102.OilType = (int)GetDouble(fp);
				break;
//B.S.------------------------
//add for ARI standard compressor equations in the form of polynomial with ten terms
				case 103:
				Q103.rpm = GetDouble(fp);//B.S. no use, put it here in case
				Q103.dispVolume=GetDouble(fp);//B.S. no use, put it here in case
				Q103.np=GetDouble(fp);
				for(i=0;i<10;i++) Q103.powerP[i] = GetDouble(fp);
				for(i=0;i<10;i++) Q103.massFlowP[i] = GetDouble(fp);
				for(i=0;i<3;i++)  Q103.Adj[i] = GetDouble(fp);
				Q103.InnerVolume=GetDouble(fp);//compressor free inner volume
				Q103.OilCharge=GetDouble(fp);//oil chargre amount
				//lubricant temperature
				Q103.OilTemp=GetDouble(fp);//oil temperature
				//---------------------------
				Q103.OilType = (int)GetDouble(fp);// oil type for getting the correct solubility equation
				break;
//--------------------------------------B.S.
				default:
				break;

		};

		endSign = (int)GetDouble(fp);//-256
		fclose(fp);

		if(endSign!=-256)
		{
		errorLog.Add("Compressor","Parameter input wrong");
		return;
		}

		compInit=0;
	}

	switch(type) {
		case 101:
			Compressor101(Q101,H1,P1,P2,Toa,mr,HPo,HPi,Ei,m,Prms);
			break;
		case 102:
			Compressor102(Q102,H1,P1,P2,Toa,mr,HPo,HPi,Ei,m,Prms);
			break;
//B.S.------------------------
//add for ARI standard compressor equations in the form of polynomial with ten terms
		case 103:
			Compressor103(Q103,H1,P1,P2,Toa,mr,HPo,HPi,Ei,m,Prms);
//--------------------------------------B.S.
			break;
	};

}

void Compressor101(CompP101& Q,double H1,double P1,double P2,double Toa,double *mr,HP *HPo,HP *HPi,double *Ei,MASS *m,double *Prms)
{
	HP HPmin;
	TXP TXP2;
	double neAdj,nvAdj;
	TXP TXP_prop={0,0,0};

	nvAdj = Q.Adj[0]*Prms[0];
	neAdj = Q.Adj[1]*Prms[1];
	Q.npAdj = Q.Adj[2]*Prms[2];

	// Set up parameters needed to solve for v2
	Q.HP1.H=H1;
	Q.HP1.P=P1;
	Q.P2=P2;
	TXP TXP1 = HPtoTXP(Q.HP1);
	if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"(TXP1) H=%lf P=%lf",Q.HP1.H,Q.HP1.P);
		errorLog.Add("Compressor101",str);
		return;
	}

	Q.v1 = PropertyTXPth(VOL,TXP1);
	if(errorLog.IsError()) {
		errorLog.Add("Compressor101","Q.v1");
		return;
	}

	// Bound v2 between v2min and v2max
	TXP_prop.P=P2;
	TXP_prop.T=TMAX;
	TXP_prop.X=1;
	const double v2max = PropertyTXPth(VOL,TXP_prop);//refshthPT.v(P2,TMAX);
	if(errorLog.IsError()) {
		errorLog.Add("Compressor101","v2max");
		return;
	}

	HPmin.H=H1;
	HPmin.P=P2;

	double v2min = PropertyTXPth(VOL,HPtoTXP(HPmin));
	if(errorLog.IsError()) {
		errorLog.Add("Compressor101","v2min");
		return;
	}

	// Iteratively solve for v2
	Zbrent(v2min,v2max,CompFunc101,1e-7,&Q);
	if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"v2: H1=%lf P1=%lf P2=%lf",H1,P1,P2);
		errorLog.Add("Compressor101",str);
		return;
	}

	// output state
	HPo->H=Q.H2;
	HPo->P=P2;

	TXP2 = HPtoTXP(*HPo);
	if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"(TXP2) H=%lf P=%lf",HPo->H,HPo->P);
		errorLog.Add("Compressor101",str);
		return;
	}

	/* mass flow rate */
	double dT = Toa-TXP1.T;
	double rP = P2/P1;
	double T = dT<=0 ? 1e-6 : dT;
	double *D = Q.VolumeEffP;
	const double arg = 1.8*T;
	if(arg<0.0) {
		errorLog.Add("Compressor101","1.8*T<0.0");
		return;
	}
	double nv = nvAdj * (1/(D[0]+rP*D[1])-D[2]*pow(arg,D[3]));
	double w = Q.rpm/60;
	double Vs = Q.dispVolume*0.0254*0.0254*0.0254;
	*mr = nv * w * Vs / Q.v1 ;

	/* input power */
	double W = *mr * (Q.H2-H1);
	D = Q.PowerP;
	if(rP<0.0) {
		errorLog.Add("Compressor101","rP<0.0");
		return;
	}
	*Ei = neAdj * (D[0]+W*(D[1]+W*D[2])+D[3]*pow(rP,D[4]));

	//input state
	HPi->H=H1;
	HPi->P=P1;

	/* charge mass */
	const double h=0.3302,d=0.2183,fv=0.7;
	m->V = fv*h*acos((-1)/4.0*d*d);
	m->m = m->V/Q.v1;
	const double m_dis_oil = Dissolved_Oil(P1,Q.OilCharge,Q.OilTemp,Q.OilType);//calculate the refrigerant dissolved in the compressor oil

	m->m = m->m+m_dis_oil;
}

// ARI Compressor map model
// INPUTS
// H1 = inlet enthalpy (J/kg)
// P1 = inlet pressure (Pa)
// P2 = outlet pressure (Pa)
// Toa = ambient air temperature around compressor (C)
// Prms = normalized parameter adjustments
//        Prms[0] = volumetric efficiency
//        Prms[1] = electrical efficiency
//        Prms[2] = polytropic efficiency
// Q = compressor parameter data structure
//        np = polytropic efficiency
//        powerP[11] = ARI coefficients for power
//        massFlowP[11] = ARI coefficients for mass flow rate
//        Adj[3] = normalized parameter adjustments (see above)
// OUTPUTS
// mr = refrigerant mass flow rate (kg/s)
// HPo = outlet state (h in J/kg, P in Pa)
// HPi = inlet state (h in J/kg, P in Pa)
// Ei = energy consumption (kW)
// m = mass of stored refrigerant (kg)
void Compressor102(CompP102& Q,double H1,double P1,double P2,double Toa,double *mr,HP *HPo,HP *HPi,double *Ei,MASS *m,double *Prms)
{
	double nvAdj = Q.Adj[0]*Prms[0];
	double neAdj = Q.Adj[1]*Prms[1];
	Q.npAdj = Q.Adj[2]*Prms[2];

	// Set up parameters needed to solve for v2
	Q.HP1.H=H1;
	Q.HP1.P=P1;
	Q.P2=P2;

	TXP TXP_prop={0,0,0};
	TXP TXP1 = HPtoTXP(Q.HP1);
	if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"(TXP1) H=%lf P=%lf",Q.HP1.H,Q.HP1.P);
		errorLog.Add("Compressor102",str);
		return;
	}

	Q.v1 = PropertyTXPth(VOL,TXP1);
	if(errorLog.IsError()) {
		errorLog.Add("Compressor102");
		return;
	}

	/* Bound v2 between v2min and v2max */
	TXP_prop.P=P2;
	TXP_prop.T=TMAX;
	TXP_prop.X=1;
	double v2max = PropertyTXPth(VOL,TXP_prop);//refshthPT.v(P2,TMAX);
	if(errorLog.IsError()) {
		errorLog.Add("Compressor102","v2max");
		return;
	}

	HP HPmin;
	HPmin.H=H1;
	HPmin.P=P2;

	double v2min = PropertyTXPth(VOL,HPtoTXP(HPmin));
	if(errorLog.IsError()) {
		errorLog.Add("Compressor102","v2min");
		return;
	}

	// Iteratively solve for v2
	Zbrent(v2min,v2max,CompFunc102,1e-7,&Q);
	if(errorLog.IsError()) {
		errorLog.Add("Compressor102", "Zbrent");
		return;
	}

	// output state
	HPo->H=Q.H2;
	HPo->P=P2;
	//TXP TXP2 = HPtoTXP(*HPo);
	if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"(TXP2) H=%lf P=%lf",HPo->H,HPo->P);
		errorLog.Add("Compressor102",str);
		return;
	}

	// mass flow rate
	double* map;
	double Psuc = P1;
	double Pdisc = P2;

	map = Q.massFlowP;
	*mr =
		map[0]+
		Psuc*(map[1]+Psuc*(map[2]+Psuc*map[3]))+
		Pdisc*(map[4]+Pdisc*(map[5]+Pdisc*map[6]))+
		Psuc*Pdisc*(map[7]+Pdisc*map[8]+Psuc*map[9]+Psuc*Pdisc*map[10]);
	*mr *= nvAdj;

	map = Q.powerP;
	*Ei =
		map[0]+
		Psuc*(map[1]+Psuc*(map[2]+Psuc*map[3]))+
		Pdisc*(map[4]+Pdisc*(map[5]+Pdisc*map[6]))+
		Psuc*Pdisc*(map[7]+Pdisc*map[8]+Psuc*map[9]+Psuc*Pdisc*map[10]);
	*Ei *= neAdj;

	//input state
	HPi->H=H1;
	HPi->P=P1;

	// charge mass
	const double h=0.3302,d=0.2183,fv=0.7;
	m->V = fv*h*acos((-1)/4.0*d*d);
	m->m = m->V/Q.v1;
}

//B.S.------------------------------------------
// ARI Compressor map model with ten terms
// INPUTS
// H1 = inlet enthalpy (J/kg)
// P1 = inlet pressure (Pa)
// P2 = outlet pressure (Pa)
// Toa = ambient air temperature around compressor (C)
// Prms = normalized parameter adjustments
//        Prms[0] = volumetric efficiency
//        Prms[1] = electrical efficiency
//        Prms[2] = polytropic efficiency
// Q = compressor parameter data structure
//        np = polytropic efficiency
//        powerP[11] = ARI coefficients for power
//        massFlowP[11] = ARI coefficients for mass flow rate
//        Adj[3] = normalized parameter adjustments (see above)
// OUTPUTS
// mr = refrigerant mass flow rate (kg/s)
// HPo = outlet state (h in J/kg, P in Pa)
// HPi = inlet state (h in J/kg, P in Pa)
// Ei = energy consumption (kW)
// m = mass of stored refrigerant (kg)

void Compressor103(CompP103& Q,//struct contained the compressor parameters
				   double H1,//suction enthalpy, initial guess
				   double P1,//suction pressure, initial guess
				   double P2,//discharge pressure, intial guess
				   double Toa,//air ambient temperature
				   double *mr,//refrigerant mass flow rate, for output
				   HP *HPo,//discharge state of the compressor
				   HP *HPi,//suction state of the compressor
				   double *Ei,//compressor power consumption
				   MASS *m,//for compressor inner volume and charge amount
				   double *Prms)//adjustment parameters
{
	double nvAdj = Q.Adj[0]*Prms[0];
	double neAdj = Q.Adj[1]*Prms[1];
	double Vol20=0,Tevap=0,Tcond=0;
	Q.npAdj = Q.Adj[2]*Prms[2];

	// Set up parameters needed to solve for v2
	Q.HP1.H=H1;
	Q.HP1.P=P1;
	Q.P2=P2;

	
	TXP TXP1 = HPtoTXP(Q.HP1);
	if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"(TXP1) H=%lf P=%lf",Q.HP1.H,Q.HP1.P);
		errorLog.Add("Compressor103",str);
		return;
	}

	TXP TXP_prop={0,0,0};

	TXP_prop.X=1;//dew point
	TXP_prop.P=P1;
	Tevap=PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(P1);
	TXP_prop.X=1;//dew point
	TXP_prop.P=P2;
	Tcond=PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(P2);

	
	TXP_prop.P=TXP1.P;
	TXP_prop.T=TXP1.T;
	TXP_prop.X=TXP1.X;


	Q.v1 = PropertyTXPth(VOL,TXP_prop);
	if(errorLog.IsError()) {
		errorLog.Add("Compressor103", "v1");
		return;	}

	/* Calc vol @ comp inlet w/ 20F superheat (ARI) */
	TXP VolARI;
	VolARI.T=Tevap+11.11;//11.11 is defined as the standard superheat degree for getting the compressor map
	VolARI.X=1;
	VolARI.P=TXP1.P;
	Vol20=PropertyTXPth(VOL,VolARI);
	if(errorLog.IsError()) {
	errorLog.Add("Compressor103", "Vol20");
	return;
	}


	// mass flow rate
	double* D;
	double Pow_Comp=0;
	double Psuc = P1;
	double Pdisc = P2;

	D=Q.massFlowP;

	double TeF,TcF,mf;

	TeF=Tevap*1.8+32;
	TcF=Tcond*1.8+32;

	mf=D[0]+D[1]*TeF+D[3]*pow(TeF,2)+D[6]*pow(TeF,3)+
          (D[2]+D[4]*TeF+D[7]*pow(TeF,2))*TcF+(D[5]+D[8]*TeF)*pow(TcF,2)
			 +D[9]*pow(TcF,3) ;
	
	
	*mr=mf*1.26e-4*(1+0.75*(Vol20/Q.v1-1.0));//Rice suction density correction
	
	*mr *= nvAdj; 

	D=Q.powerP;

	Pow_Comp=D[0]+D[1]*TeF+D[3]*pow(TeF,2)+D[6]*pow(TeF,3)+
          (D[2]+D[4]*TeF+D[7]*pow(TeF,2))*TcF+(D[5]+D[8]*TeF)*pow(TcF,2)
			 +D[9]*pow(TcF,3);
	
	*Ei=neAdj*Pow_Comp; // B.S. Without adjusting the compressor power predictions, since on some occasions, the suction density correction does not do good.

	//new--------------- output state
	const double HeatLossR = 1.0;
	HPo->H=H1+neAdj*Pow_Comp*HeatLossR/(*mr);
	HPo->P=P2;
	//new---------------------


	//input state
	HPi->H=H1;
	HPi->P=P1;

	// charge mass, existing in compressor inner volume and dissolved in the compressor oil
	m->V = Q.InnerVolume;
	const double T_oil = Prms[3];
	const double m_dis_oil = 0;//Dissolved_Oil(P1,Q.OilCharge,Q.OilTemp,Q.OilType);//calculate the refrigerant dissolved in the compressor oil
	const double m_gas_com=m->V/Q.v1; //calculate the refrigerant gas amount in the compressor free inner volume
	m->m = m_gas_com+m_dis_oil;//refrigerant mass in the compressor is the addition of the refrigerant gas and the refrigerant dissovled in the oil
}


//---------------------------------------------------B.S.
/**********************************************************
Martz, W. L., Burton, C.M., and Jacobi, A.M., 1996 "Local composition modelling
of the thermodynamic properties of refrigerant and oil mixtures," 
International Journal of Refrigerantion, Vol. 19, pp. 25-33.
or 
Get solubility by interpolating manufacturer's data
*************************************************************/
double 
Dissolved_Oil(double P1,
			  double OilCharge,
			  double OilTemp,
			  int OilType)
{
	const double T_bottom=OilTemp;
	OilTemp =T_bottom;//oil temperature is assumed to be the compressor bottom temperature

	double epsilon=0;
	const double T=T_bottom;
	const double P=P1;
	switch(OilType)
	{
	case 1:
		epsilon=P1/(773.0+18.5*OilTemp+0.525*pow(OilTemp,2.0));     //Martz et al.
		//epsilon=lubriPropComp.Solubility(OilTemp,P1);
		
		break;
	case 2:
		epsilon=lubriPropComp.Solubility(OilTemp,P1);
	
		break;
	default:
		epsilon=P1/(773.0+18.5*OilTemp+0.525*pow(OilTemp,2.0));     //Martz et al.
		break;
	}
	
	if(epsilon<0||epsilon>1.0) 
	{errorLog.Add("Dissolved_Oil", "solubility");
	return -1;
	}

	if(errorLog.IsError()) {
		errorLog.Add("Dissolved_Oil", "solubility");
		return -1;}

	const double m_dis_oil=(OilCharge)*(epsilon/(1-epsilon)); 
	return m_dis_oil;
}

/********************************************************************
Generates a residual used by a numerical method to solve for the
compressor outlet state (v2).  The residual is zero when the guess
(v2a) is correct.

Inputs: v2a (guess value for specific volume of state 2, m^3/kg),
	Params->HP1 (refrigerant state at 1, J/kg, kPa),
	Params->v1 (specific volume at 1, m^3/kg),
	Params->P2 (pressure at 2, kPa),
	Params->npAdj (multiplicative fitting factor for polytropic
		       efficiency).  Should normally be 1, except
		       when fitting experimental data.

Outputs: Params->H2 (enthalpy at 2, J/kg),
	 params->n (polytropic exponent, dimensionless),
	 status (set non-zero to indicate error),
	 CompFunc (normalized specific volume residual, m^3/kg).
********************************************************************/
double CompFunc101(double v2a, void *Params)
{
	CompP101 *P = (CompP101*)Params;

	double np = P->np*P->npAdj;

	P->n = log(P->HP1.P/P->P2)/log(v2a/P->v1);

	double rn = (P->n-1)/P->n;
	double rP = P->P2/P->HP1.P;

	HP HP2;
	if(rP<0.0) {
		errorLog.Add("CompFunc101","rP<0.0");
		return -1;
	}
	HP2.H = P->H2 = P->HP1.H + 1000*P->HP1.P*P->v1*(pow(rP,rn)-1)/(rn*np);
	HP2.P = P->P2;

	TXP TXP2 = HPtoTXP(HP2);
	if(errorLog.IsError()) {
		errorLog.Add("CompFunc101","TXP2");
		return -1;
	}

	double v2b = PropertyTXPth(VOL,TXP2);
	if(errorLog.IsError()) {
		errorLog.Add("CompFunc101","v2b");
		return -1;
	}

	double dv = (v2a-v2b)/(v2a+v2b);

	return dv;
}

double
CompFunc102(double v2a,void *Params)
{
	CompP102* P = (CompP102*)Params;

	double np = P->np*P->npAdj;

	P->n = log(P->HP1.P/P->P2)/log(v2a/P->v1);

	double rn = (P->n-1)/P->n;
	double rP = P->P2/P->HP1.P;

	HP HP2;
	if(rP<0.0) {
		errorLog.Add("CompFunc102","rP<0.0");
		return -1;
	}
	HP2.H = P->H2=P->HP1.H+1000*P->HP1.P*P->v1*(pow(rP,rn)-1)/(rn*np);
	HP2.P = P->P2;

	TXP TXP2 = HPtoTXP(HP2);
	if(errorLog.IsError()) {
		errorLog.Add("CompFunc102","TXP2");
		return -1;
	}

	const double v2b = PropertyTXPth(VOL,TXP2);
	if(errorLog.IsError()) {
		errorLog.Add("CompFunc102","v2b");
		return -1;
	}
	const double dv = (v2a-v2b)/((v2a+v2b)/2);
	return dv ;
}

/********************************************************************
Compares the calculated and actual values of compressor power and
flow rate.  The actual data comes from a file in the DATA directory
called comptest.dat which was made from comp1b.ees
********************************************************************/
void CompressorMap(const char *filename)
{
	double Te,Tc,P2,mr,Ei,mrM,EiM;
	double CompPrms[3]={1,1,1};
	double d[30];
	double dE,dM,Emax=-1e20,Mmax=-1e20;
	int i;
	TXP TXP1;
	HP HP1,HPso,HPsi;
	MASS m;

	FILE *fp = fopen("comptest.dat","r");
	for(i=0;i<19;i++) fscanf(fp,"%lf",&d[i]);
	do {
		errorLog.ClearError("CompressorMap");

		Te=(d[0]-32)/1.8;
		Tc=(d[1]-32)/1.8;
		EiM=d[2];
		mrM=d[3]*1.26e-4;

		TXP1.T = Te+20/1.8;
		TXP1.X = 1;
		TXP1.P = reftpvthT.Psat(Te);
		HP1 = TXPtoHP(TXP1);
		P2 = reftplthT.Psat(Tc);

		Compressor(filename,HP1.H,HP1.P,P2,35,&mr,&HPso,&HPsi,&Ei,&m,CompPrms);

		dE=100*(Ei-EiM)/EiM;
		if(fabs(dE)>Emax) Emax=fabs(dE);
		dM=100*(mr-mrM)/mrM;
		if(fabs(dM)>Mmax) Mmax=fabs(dM);
		printf("%.0lf %.0lf %lf %lf %lf %lf %lf %lf \n",Te,Tc,EiM,Ei,dE, mrM,mr,dM);

		for(i=0;i<19;i++) fscanf(fp,"%lf",&d[i]);

	} while(!feof(fp));

	fclose(fp);
	printf("%lf %lf\n",Emax,Mmax);
}
