#include "stdafx.h"

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "air.h"
#include "r22.h"
#include "wair.h"
#include "mass.h"
#include "comp.h"
#include "lines.h"
#include "evap.h"
#include "cond.h"
#include "expdev.h"
#include "num2d.h"
#include "num3d.h"
#include "numeric.h"
#include "mainfunc.h"
#include "outelem.h"
#include "mtime.h"
#include "subcool.h"
#include "compstat.h"
#include "ref.h"
#include "errorlog.h"
#include "defs.h"
#include "spline1d.h"
#include "spline2d.h"//shenbo add for calculating superheat state
#include "defs.h"


//#define DEBUG_MAINFUNC

// Bo Shen did not want the added complication of the restriction model.  Include this
// #define to incorporate the restriction model, otherwise it is ommitted.
//#define INCLUDE_RESTRICTION	

// Indicates which component (e.g. compressor, condenser, etc.) is
// curently being solved.
int cyclePosition=0;

// Cycle position labels
#define NCYCLEPOS 20
const char* cyclePositionLabels[NCYCLEPOS] = {
	"",                   // 0
	"Compressor(1)",      // 1
	"Compressor(2)",
	"Oil Separator",
	"Condenser",          // 4
	"Receiver",
	"Subcooler",          // 6
	"Accumulator",
	"Evaporator",         // 8
	"Compressor",         // 9
	"Expansion Device",   // 10
	"Restriction",        // 11
	"",
	"",
	"",
	"",
	"",
	"",
	"",
	""
};

// Tracks the amount of time (in milliseconds) that is spent
// computing each component.
unsigned long int cycleTime[NCYCLEPOS]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

extern OutputDataElement gOutData;
extern ErrorLog errorLog;
extern WetAirPropTR wair;

/*********************************************************************
Vapor compression cycle model.  This function combines the components
to make the cycle.  It that three inputs and returns three residuals
which are all zero for the correct operating state.

Inputs:

	X[0] = h1 (compressor inlet enthalpy)
	X[1] = P1 (compressor inlet pressure)
	X[2] = P2 (compressor outlet pressure)

	P->Tai    (condenser air inlet temperature)
	P->TPi.T  (evaporator air inlet temperature)
	P->TPi.P  (evaporator air inlet relative humidity)
	P->Goa    (condenser air mass flux)
	P->Gra    (evaporator air mass flux)
	P->Charge (refrigerant charge)
	P->Dadj   (liquid line restriction)
	P->nvAdj  (compressor valve leak)

Outputs:

	Y[0] = dP (normalized pressure residual)
	Y[1] = dh (normalized enthalpy residual)
	Y[2] = dm (normalized charge residual)

	return = CE (convergence error)

States

	1 = compressor inlet, suction line outlet
	2 = compressor outlet
	***
	3 = hot gas line outlet
	4 = condenser outlet, liquid line inlet
	5 = liquid line outlet
	***
	10 = restriction inlet, liquid line outlet
	9 = expansion device inlet, restriction outlet
	6 = distribution tubes inlet, expansion device outlet
	7 = evaporator inlet
	8 = suction line inlet

*********************************************************************/

/*********************************************************************
Main function for system type #1.  This is a vapor compression
cycle with a short orifice expansion device.

Returns 0 is successful and non-zero error code if a problem occurred.
*********************************************************************/
int MainFunc1(double *X,double *Y,void *Params)
{
	static double mr,Ei,P1,P2,Tao,Aflow,CE=-1;
	static double Ones[5]={1,1,1,1,1};
	static double CompPrms[4]={1,1,1,50};//CompPrms[3]=50, for lubricant temperature
	static double RstPrms[1]={1};
	static MASS m[9];
	static HP HPs[9];
	static TXP TXPs[9];
	static double Tsat[9];
	static TP TPo;
	static InVars *P;
	static CGP Cond_struc;//B.S. condenser construct
	static ETdim Evap_struc;//B.S. evaporator construct
	static double mr_exp=0;//B.S. mass flow rate through the expansion device
	static HP HP_Sub={0,0};//enthalpy ahead of expansion device
	static TXP TXP_sub;//state ahead of expansion device
	static double Hevpex,Pevpex;//evaporator exit enthalpy and pressure
	TXP TXP_prop={0,0,0};

	P=(InVars*)Params;

	//
	// Keep track of how much time is spent crunching each component.
	// This provides the ability to determine which component should
	// be made more efficient.
	static int init=1;
	static unsigned long int saveTime;
	if(init) {
		saveTime = MilliTime();
		init=0;

		#ifdef DEBUG_MAINFUNC
		{
		FILE* fp = fopen("MainFuncDebug.txt","w");
		if(fp) {
			fprintf(fp,"Starting...\n");
			fclose(fp);
		}
		}
		#endif

	}
	unsigned long int currentTime = MilliTime();
	cycleTime[cyclePosition] += currentTime - saveTime;
	saveTime = currentTime;

	#ifdef DEBUG_MAINFUNC
	{
	FILE* fp = fopen("MainFuncDebug.txt","a");
	if(fp) {
		fprintf(fp,"Guess: Hin=%.1lf,Pin=%.1lf,Pout=%.1lf\n",X[0],X[1],X[2]);
		fprintf(fp,"Tai=%.1lf,TPi.T=%.1lf,TPi.P=%.1lf\n",P->Tai,P->TPi.T,P->TPi.P);
		fprintf(fp,"Gca=%.1lf,Gea=%.1lf\n",P->Goa,P->Gra);
		fprintf(fp,"Charge=%.4lf,nvAdj=%.4lf,Dadj=%.4lf\n",P->Charge,P->nvAdj,P->Dadj);
		fprintf(fp,"filename=[%s]\n",P->filename);
		fclose(fp);
	}
	}
	#endif

//	gOutData.ClearAll();

	gOutData.Charge_target = P->Charge;
	gOutData.Gca = P->Goa;
	gOutData.Gea = P->Gra;
	gOutData.nvAdj = P->nvAdj;
	gOutData.Dadj = P->Dadj;

	CompPrms[0]=P->nvAdj;
	CompPrms[3]=P->CompT;
	cyclePosition = 9;
	Compressor(P->filename,X[0],X[1],X[2],P->Tai,&mr,&HPs[2],&HPs[1],&Ei,&m[0],CompPrms);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc1_Fwd","Compressor");
		return 1;
	}

	// compressor power and flow rate
	gOutData.Ecmp = Ei;
	gOutData.mr[1] = mr;
	gOutData.m[0] = m[0].m;

	// state 1 - compressor inlet
	TXPs[1] = HPtoTXP(HPs[1]);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc1_Fwd","Tsat[1]");
		return 2;
	}

	TXP_prop.P=HPs[1].P;
	TXP_prop.X=1.0;
	Tsat[1] = PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(HPs[1].P);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc1_Fwd","TXPs[1]");
		return 3;
	}
	gOutData.H[1] = HPs[1].H;
	gOutData.T[1] = TXPs[1].T;
	gOutData.X[1] = TXPs[1].X;
	gOutData.P[1] = TXPs[1].P;
	gOutData.Tsat[1] = Tsat[1];

	// state 2 - compressor outlet
	TXP_prop.P=HPs[2].P;
	TXP_prop.X=1.0;
	Tsat[2] = PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(HPs[2].P);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc1_Fwd","Tsat[2]");
		return 4;
	}
	TXPs[2] = HPtoTXP(HPs[2]);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc1_Fwd","TXPs[2]");
		return 5;
	}
	gOutData.H[2] = HPs[2].H;
	gOutData.T[2] = TXPs[2].T;
	gOutData.X[2] = TXPs[2].X;
	gOutData.P[2] = TXPs[2].P;
	gOutData.Tsat[2] = Tsat[2];

	HotGasLine(HPs[2],mr,&HPs[3],&m[1],P);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc1_Fwd");
		return 6;
	}
	
	//state 3-condenser inlet
	gOutData.H[3] = HPs[3].H;
	gOutData.P[3] = HPs[3].P;
	gOutData.m[1] = m[1].m;

	SuctionLine(HPs[1],mr,&HPs[8],&m[2]);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc1_Fwd");
		return 7;
	}
	//state 8-evaporator exit
	Hevpex=HPs[8].H;
	Pevpex = HPs[8].P;
	TXP_prop.P=HPs[8].P;
	TXP_prop.X=1;
	Tsat[8] = PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(HPs[8].P);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc1_Fwd","Tsat[8]");
		return 7;
	}
	TXPs[8] = HPtoTXP(HPs[8]);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc1_Fwd","TXPs[8]");
		return 7;
	}

	gOutData.T[8] = TXPs[8].T;
	gOutData.X[8] = TXPs[8].X;
	gOutData.P[8] = TXPs[8].P;
	gOutData.Tsat[8] = Tsat[8];

	gOutData.H[8] = HPs[8].H;
	gOutData.P[8] = HPs[8].P;

	cyclePosition = 4;
	while(Condenser(P->filename,mr,HPs[3],P->Tai,P->Goa,&HPs[4],&Tao,&m[3],&Cond_struc,Ones)) {	}
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc1_Fwd","Condenser");
		return 8;
	}
	P->Cond = Cond_struc;//B.S. condenser struct

	// condenser air side states
	gOutData.cfmca = Cond_struc.cfma; // cond air cfm
	gOutData.mc = P->Goa;
	gOutData.Tci = P->Tai;
	gOutData.Tco = Tao;
	gOutData.m[3] = m[3].m;

	// state 4 - condenser outlet
	TXP_prop.P=HPs[4].P;
	TXP_prop.X=0;
	Tsat[4] = PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(HPs[4].P);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc1_Fwd","Tsat[4]");
		return 9;
	}
	TXPs[4] = HPtoTXP(HPs[4]);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc1_Fwd","TXPs[4]");
		return 10;
	}
	gOutData.H[4] = HPs[4].H;
	gOutData.T[4] = TXPs[4].T;
	gOutData.X[4] = TXPs[4].X;
	gOutData.P[4] = TXPs[4].P;
	gOutData.Tsat[4] = Tsat[4];
	
	
	LiquidLine(HPs[4],mr,&HPs[5],&m[5]);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc1_Fwd");
		return 11;
	}

	//state 5-liquid line exit
	TXPs[5]= HPtoTXP(HPs[5]);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc1_Fwd","TXPs[5]");
		return 12;
	}

	TXP_prop.P=HPs[5].P;
	TXP_prop.X=0;
	Tsat[5] = PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(HPs[5].P);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc1_Fwd","Tsat[5]");
		return 13;
	}

	gOutData.H[5] = HPs[5].H;
	gOutData.T[5] = TXPs[5].T;
	gOutData.X[5] = TXPs[5].X;
	gOutData.P[5] = TXPs[5].P;
	gOutData.Tsat[5] = Tsat[5];
	gOutData.m[5] = m[5].m;

	gOutData.Twb_ei = wair.WetBulb(P->TPi.T,P->TPi.P);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc1_Fwd","Twb_ei");
		return 14;
	}


	// This evaporator model is just 1/4 of the entire coil.  There are 4 independent sections stacked
	// on top of eachother, each with their own expansion device.

	cyclePosition = 8;

	HPs[7].H=HPs[5].H;//evaporator inlet enthalp equal to liquid line outlet
	Evap_struc.REV=P->REV_CAL;
	while(Evaporator(P->filename,mr,&HPs[8],P->Gra,P->TPi,&HPs[7],&TPo,&m[4],&Aflow,&Evap_struc,Ones)) { }
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc1_Fwd","Evaporator");
		return 15;
	}

	P->Evap = Evap_struc;//B.S. evaporator struct

	// evaporator air side states
	gOutData.cfmea = Evap_struc.cfma; // evap air cfm
	gOutData.me = P->Gra;
	gOutData.Tei = P->TPi.T;
	gOutData.Teo = TPo.T;
	gOutData.HumPei = P->TPi.P;
	gOutData.HumPeo = TPo.P;
	gOutData.m[4] = m[4].m;

	// state 7 - evaporator inlet
	TXP_prop.P=HPs[7].P;
	TXP_prop.X=1;
	Tsat[7] = PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(HPs[7].P);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc1_Fwd","Tsat[7]");
		return 16;
	}
	TXPs[7] = HPtoTXP(HPs[7]);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc1_Fwd","TXPs[7]");
		return 17;
	}
	gOutData.H[7] = HPs[7].H;
	gOutData.T[7] = TXPs[7].T;
	gOutData.X[7] = TXPs[7].X;
	gOutData.P[7] = TXPs[7].P;
	gOutData.Tsat[7] = Tsat[7];

	DistributionTubes(HPs[7],mr,&HPs[6],&m[6]);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc1_Fwd","Distribution Tubes");
		return 18;
	}
	
	//state 6-after expansion device
	gOutData.H[6] = HPs[6].H;
	gOutData.P[6] = HPs[6].P;
	gOutData.m[6] = m[6].m;

	cyclePosition = 10;

	if(P->PISTON)//short-tube
	{

		if(P->SPEC_SUB == 2)//specify the subcooling state ahead of expansion device
		{
			if(P->Tsub<0.0)//two-phase upstream 
			{
			TXP_sub.P=TXPs[5].P; TXP_sub.X=-P->Tsub; TXP_sub.T=PropertyTXPth(TSAT,TXP_sub);
			}
			else
			{
			TXP_sub.P=TXPs[5].P; TXP_sub.X=0; TXP_sub.T=Tsat[5]-P->Tsub;
			}				
			HP_Sub = TXPtoHP(TXP_sub);
			if(errorLog.IsError()) {
			errorLog.Add("MainFunc1_Fwd","TXP_sub to HP_Sub");
			return 18;}
		}
		else//get the calculated subcooling state ahead of expansion device
		{
			HP_Sub = HPs[5];
		}

#ifdef _R22
		mr_exp = Short_Tube_R22(P->filename,HP_Sub,HPs[6],Ones);
#endif

#ifdef _R410A
		HPs[6].P=HPs[1].P+375.09;//B.S.
		mr_exp = Short_Tube_R410A(P->filename,HP_Sub,HPs[6],Ones);
#endif

#ifdef _R407C
		mr_exp = Short_Tube_FlowDirec(P->filename,HP_Sub,HPs[6],Ones);
#endif
		
		gOutData.mr_exp = mr_exp;
		gOutData.P2 = HPs[5].P;
		m[7].m=0;
		m[7].V=0;
	}
	else//apparent TxV
	{
	mr_exp = TxV_Nonlinear(P->filename,HPs[5],HPs[6],TXPs[8].T-Tsat[8],Ones);
	gOutData.mr_exp = mr_exp;
	gOutData.P2 = HPs[5].P;
	m[7].m=0;
	m[7].V=0;
	}

	if(errorLog.IsError()) {
		errorLog.Add("MainFunc1_Fwd","Expansion Device");
		return 19;
	}
	m[7].m=0;
	
#ifdef INCLUDE_RESTRICTION	
	// B.S. assume there is no such a flow restriction device in the liquid line
	RstPrms[0]=P->Dadj;
	cyclePosition = 11;
	ExpDev(P->filename,2,HPs[5].H,P2,mr,&P1,&m[8],RstPrms);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc1_Fwd","Restriction");
		return 20;
	}

	gOutData.P1 = P1;
#else
	gOutData.P1 = P2;
#endif

	m[8].m=0;
	cyclePosition = 0;
	
	gOutData.m[7] = m[7].m;
	gOutData.m[8] = m[8].m;
	//B.S.------------------------------
	//calculate cooling and condensing capacity
	gOutData.Qc=mr*(HPs[3].H-HPs[4].H);
	gOutData.Qe=mr*(HPs[7].H-HPs[1].H);
	gOutData.Tsub4=-1*(gOutData.T[4]-gOutData.Tsat[4]); 
	gOutData.Tsub5=-1*(gOutData.T[5]-gOutData.Tsat[5]); 
	gOutData.Tsup1=gOutData.T[1]-gOutData.Tsat[1];
	//----------------------------------B.S.


	// calculate normalized residuals
	Y[0] =2*(mr-mr_exp)/(mr+mr_exp);//calcuate the mass flow rate residual
	if(P->REV_CAL)//backforward evaporator scheme
	{
	Y[1] =  2*(HPs[5].H-HPs[6].H)/(HPs[5].H+HPs[6].H);
	}
	else{//forward evaporator scheme
		if(P->EvaSign==0) Y[1] = 2*(HPs[8].H-Hevpex)/(HPs[8].H+Hevpex);//Newton iteration
		else Y[1] = (HPs[8].P-Pevpex);//iteration by updating lumped models
	}

	const double Chg  =  m[0].m+m[1].m+m[2].m+m[3].m+m[5].m+m[8].m+(m[4].m+m[6].m+m[7].m);
	P->ChargMass = Chg;//calculated charge inventory, P->Charge, specified charge inventory

	if(P->SPEC_SUB)//specify subcooling, P->SPEC_SUB=1 or 2
	{
		if(P->Tsub<0.0)//two-phase state 
		{
			if(P->SPEC_SUB == 2)//specify upstream state of expansion device
			{
			TXP_sub.P=TXPs[5].P; TXP_sub.X=-P->Tsub; TXP_sub.T=PropertyTXPth(TSAT,TXP_sub);
			}
			else{//specify condenser exit state
			TXP_sub.P=TXPs[4].P; TXP_sub.X=-P->Tsub; TXP_sub.T=PropertyTXPth(TSAT,TXP_sub);
			}
			HP_Sub.H=PropertyTXPth(ENTH,TXP_sub);
			if(errorLog.IsError()) {
			errorLog.Add("MainFunc1_Fwd","H_sub");
			return 21;}
		}
		else//subcooled state
		{
			if(P->SPEC_SUB == 2)
			{
			TXP_sub.P=TXPs[5].P; TXP_sub.X=0; TXP_sub.T=Tsat[5]-P->Tsub;
			}
			else{
			TXP_sub.P=TXPs[4].P; TXP_sub.X=0; TXP_sub.T=Tsat[4]-P->Tsub;
			}
			HP_Sub.H=PropertyTXPth(ENTH,TXP_sub);
			if(errorLog.IsError()) {
			errorLog.Add("MainFunc1_Fwd","H_sub");
			return 21;}
		}

		if(P->SPEC_SUB == 2){//specify expansion device upstream state
			if(P->EvaSign==1) {//using simplied models for calculations
				if(P->Tsub>0.0) Y[2]=(Tsat[5]-TXPs[5].T)-P->Tsub;
				else Y[2]=(TXPs[5].X+P->Tsub)*10;
				}
			else{//Newton method for iteration
				Y[2]=2*(HPs[5].H-HP_Sub.H)/(HPs[5].H+HP_Sub.H);
				}
			}
		else{//specify condenser exit state
			if(P->EvaSign==1) {
				if(P->Tsub>0.0) Y[2]=(Tsat[4]-TXPs[4].T)-P->Tsub;
				else Y[2]=(TXPs[4].X+P->Tsub)*10;
				}
			else{
				Y[2]=2*(HPs[4].H-HP_Sub.H)/(HPs[4].H+HP_Sub.H);
				}			
			}
	}
	else//specify charge
	{
	Y[2] =- (Chg-P->Charge)/(P->Charge);
	}

	// calculate convergence error
	CE = sqrt((Y[0]*Y[0]+Y[1]*Y[1]+Y[2]*Y[2])/3);

	gOutData.CE = CE;

	//B.S.-------------------------------------
	//output data
	gOutData.devMr=Y[0];
	gOutData.devH=Y[1];
	gOutData.devCharge=Y[2];
	gOutData.CE=CE;
	gOutData.inputH1=X[0];
	gOutData.inputP1=X[1];
	gOutData.inputP2=X[2];
	gOutData.mr_TP = mr;
	gOutData.Charge = Chg;
	gOutData.LiqLCond=Cond_struc.LiqL;
	gOutData.TPLCond=Cond_struc.TPL;
	gOutData.VapLCond=Cond_struc.VapL;
	gOutData.LiqLEvap=Evap_struc.LiqL;
	gOutData.TPLEvap=Evap_struc.TPL;
	gOutData.VapLEvap=Evap_struc.VapL;
	//-------------------------------B.S.

	printf("MainFunc1_Fwd: CE=%1.8lf (%1.8lf,%1.8lf,%1.8lf)\n",CE,Y[0],Y[1],Y[2]);//B.S.


	#ifdef DEBUG_MAINFUNC
	{
	FILE* fp = fopen("MainFuncDebug.txt","a");
	if(fp) {
		fprintf(fp,"ExpDev pressure drop: HPs[6].P=%.1lf, P2=%.1lf, P1=%.1lf (Hsat=%.1lf mr/4=%.4lf)\n",HPs[6].P,P1,P2,HPs[5].H,mr/4);
		fprintf(fp,"CE=%.4lf, Y[0]=%.4lf, Y[1]=%.4lf, Y[2]=%.4lf\n\n",CE,Y[0],Y[1],Y[2]);
		fclose(fp);
	}
	}
	#endif

	return 0;
}

/*********************************************************************
Main function for system type #2.  This is a vapor compression
cycle with a simple TxV that maintains constant superheat.
*********************************************************************/
int MainFunc2(double *X,double *Y,void *Params)
{
	static double mr,Ei,Tao,Aflow,Chg,CE=-1;
	static double Ones[5]={1,1,1,1,1},CompPrms[4]={1,1,1,50};//CompPrms[3]=50, for lubricant temperature
	static double p1,p2,h1,Tsat1;
	static MASS m[9];
	static HP HPs[9];
	static TXP TXPs[9];
	static double Tsat[9];
	static TP TPo;
	static InVars *P;
	static TXP TXP1; 
	static HP HP1;
	static CGP Cond_struc;//B.S.keep the condenser parameters
	static ETdim Evap_struc;//B.S. keep the evaporator parameters 
	static double H_Sub=0;//enthalpy ahead of expansion device
	static TXP TXP_sub;//state ahead of expansion device
	static double H_evpex, P_evpex;//for get the tolerance of evaporator exit enthalpy and pressure
	TXP TXP_prop={0,0,0};

	// Keep track of how much time is spent crunching each component.
	// This provides the ability to determine which component should
	// be made more efficient.
	static int init=1;
	static unsigned long int saveTime;
	if(init) {
		saveTime = MilliTime();
		init=0;
	}
	unsigned long int currentTime = MilliTime();
	cycleTime[cyclePosition] += currentTime - saveTime;
	saveTime = currentTime;
	//
	//

	//gOutData.ClearAll();
	P=(InVars*)Params;//move here from behind, shenbo

	p1=X[0];
	p2=X[1];

	TXP_prop.P=p1;
	TXP_prop.X=1.0;
	Tsat1 = PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(p1);
	if(P->Tsup<0)
	{
	TXP_prop.T=Tsat1;
	h1=PropertyTXPth(ENTH,TXP_prop)+PropertyTXPtr(SPEC,TXP_prop)*P->Tsup;
	HP1=toHP(h1,p1);
	TXP1=HPtoTXP(HP1);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc2_Fwd","h1");
		return 1;
	}
	}
	else{
	TXP1=toTXP(Tsat1+P->Tsup,1,p1);
	h1= PropertyTXPth(ENTH,TXP1);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc2_Fwd","h1");
		return 1;
	}
	}

//	P=(InVars*)Params;

	cyclePosition = 9;
	CompPrms[0]=P->nvAdj;
	CompPrms[3]=P->CompT;//B.S.new
	Compressor(P->filename,h1,p1,p2,P->Tai,&mr,&HPs[2],&HPs[1],&Ei,&m[0],CompPrms);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc2_Fwd");
		return 2;
	}

	// compressor power and flow rate
	gOutData.Ecmp = Ei;
	gOutData.mr[1] = mr;
	gOutData.m[0] = m[0].m;

	// state 1 - compressor inlet
	TXP_prop.P=HPs[1].P;
	TXP_prop.X=1;
	Tsat[1] = PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(HPs[1].P);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc2_Fwd","Tsat[1]");
		return 3;
	}
	TXPs[1] = HPtoTXP(HPs[1]);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc2_Fwd","TXPs[1]");
		return 4;
	}
	gOutData.H[1] = HPs[1].H;
	gOutData.T[1] = TXPs[1].T;
	gOutData.X[1] = TXPs[1].X;
	gOutData.P[1] = TXPs[1].P;
	gOutData.Tsat[1] = Tsat[1];

	// state 2 - compressor outlet
	TXP_prop.P=HPs[2].P;
	TXP_prop.X=1.0;
	Tsat[2] = PropertyTXPth(TSAT,TXP_prop);//PropertyTXPth();//reftplthP.Tsat(HPs[2].P);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc2_Fwd","Tsat[2]");
		return 5;
	}
	TXPs[2] = HPtoTXP(HPs[2]);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc2_Fwd","TXPs[2]");
		return 6;
	}
	gOutData.H[2] = HPs[2].H;
	gOutData.T[2] = TXPs[2].T;
	gOutData.X[2] = TXPs[2].X;
	gOutData.P[2] = TXPs[2].P;
	gOutData.Tsat[2] = Tsat[2];

	HotGasLine(HPs[2],mr,&HPs[3],&m[1],P);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc2_Fwd");
		return 7;
	}
	
	//state 3-condenser inlet
	gOutData.H[3] = HPs[3].H;
	gOutData.P[3] = HPs[3].P;
	gOutData.m[1] = m[1].m;
	
	SuctionLine(HPs[1],mr,&HPs[8],&m[2]);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc2_Fwd");
		return 8;
	}
	
	//state 8-evaporator exit
	H_evpex=HPs[8].H;
	P_evpex = HPs[8].P;
	TXP_prop.P=HPs[8].P;
	TXP_prop.X=1;
	Tsat[8] = PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(HPs[8].P);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc2_Fwd","Tsat[8]");
		return 8;
	}
	TXPs[8] = HPtoTXP(HPs[8]);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc2","TXPs[8]");
		return 8;
	}

	gOutData.T[8] = TXPs[8].T;
	gOutData.X[8] = TXPs[8].X;
	gOutData.P[8] = TXPs[8].P;
	gOutData.Tsat[8] = Tsat[8];

	gOutData.H[8] = HPs[8].H;
	gOutData.P[8] = HPs[8].P;
	gOutData.m[2] = m[2].m;

	cyclePosition = 4;

	Condenser(P->filename,mr,HPs[3],P->Tai,P->Goa,&HPs[4],&Tao,&m[3],&Cond_struc,Ones);
		if(errorLog.IsError()) {
			errorLog.Add("MainFunc2_Fwd");
			return 9;
		}
	P->Cond = Cond_struc;//B.S. condenser struct
	

	// condenser air side states
	gOutData.mc = P->Goa;
	gOutData.Tci = P->Tai;
	gOutData.Tco = Tao;

	// state 4 - condenser outlet state
	TXP_prop.P=HPs[4].P;
	TXP_prop.X=0.0;
	Tsat[4] = PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(HPs[4].P);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc2_Fwd","Tsat[4]");
		return 11;
	}

	TXPs[4] = HPtoTXP(HPs[4]);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc2_Fwd","TXPs[4]");
		return 12;
	}
	gOutData.H[4] = HPs[4].H;
	gOutData.T[4] = TXPs[4].T;
	gOutData.X[4] = TXPs[4].X;
	gOutData.P[4] = TXPs[4].P;
	gOutData.Tsat[4] = Tsat[4];
	gOutData.m[3] = m[3].m;

	LiquidLine(HPs[4],mr,&HPs[5],&m[5]);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc2_Fwd");
		return 12;
	}
	
	//state 5-liquid line exit
	TXPs[5]= HPtoTXP(HPs[5]);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc2_Fwd","TXPs[5]");
		return 13;
	}

	TXP_prop.P=HPs[5].P;
	TXP_prop.X=0;
	Tsat[5] = PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(HPs[5].P);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc2_Fwd","Tsat[5]");
		return 14;
	}

	gOutData.H[5] = HPs[5].H;
	gOutData.T[5] = TXPs[5].T;
	gOutData.X[5] = TXPs[5].X;
	gOutData.P[5] = TXPs[5].P;
	gOutData.Tsat[5] = Tsat[5];
	gOutData.m[5] = m[5].m;

	gOutData.Twb_ei = wair.WetBulb(P->TPi.T,P->TPi.P);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc2_Fwd","Twb_ei");
		return 15;
	}

	cyclePosition = 8;

	HPs[7].H=HPs[5].H;//evaporator inlet enthalp equal to liquid line outlet
	Evap_struc.REV=P->REV_CAL;
	Evaporator(P->filename,mr,&HPs[8],P->Gra,P->TPi,&HPs[7],&TPo,&m[4],&Aflow,&Evap_struc,Ones);
		if(errorLog.IsError()) {
			errorLog.Add("MainFunc2_Fwd");
			return 16;
		}
	P->Evap = Evap_struc;//B.S. evaporator struct
	
	// evaporator air side states
	gOutData.me = P->Gra;
	gOutData.Tei = P->TPi.T;
	gOutData.Teo = TPo.T;
	gOutData.HumPei = P->TPi.P;
	gOutData.HumPeo = TPo.P;
	gOutData.m[4] = m[4].m;

	// state 7 - evaporator inlet
	TXP_prop.P=HPs[7].P;
	TXP_prop.X=1;
	Tsat[7] = PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(HPs[7].P);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc2_Fwd","Tsat[7]");
		return 17;
	}
	TXPs[7] = HPtoTXP(HPs[7]);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc2_Fwd","TXPs[7]");
		return 18;
	}

	gOutData.H[7] = HPs[7].H;
	gOutData.T[7] = TXPs[7].T;
	gOutData.X[7] = TXPs[7].X;
	gOutData.P[7] = TXPs[7].P;
	gOutData.Tsat[7] = Tsat[7];


	DistributionTubes(HPs[7],mr,&HPs[6],&m[6]);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc2_Fwd");
		return 19;
	}
	
	//state 6-after expansion device
	gOutData.H[6] = HPs[6].H;
	gOutData.P[6] = HPs[6].P;
	gOutData.m[6] = m[6].m;

	cyclePosition=0;

	//B.S-----------------------
	gOutData.Qc=mr*(HPs[3].H-HPs[4].H);
	gOutData.Qe=mr*(HPs[7].H-HPs[8].H);
	gOutData.mr_TP=mr;
	gOutData.Tsub4=-1*(gOutData.T[4]-gOutData.Tsat[4]); 
	gOutData.Tsub5=-1*(gOutData.T[5]-gOutData.Tsat[5]); 
	gOutData.Tsup1=gOutData.T[1]-gOutData.Tsat[1];
	gOutData.LiqLCond=Cond_struc.LiqL;
	gOutData.TPLCond=Cond_struc.TPL;
	gOutData.VapLCond=Cond_struc.VapL;
	gOutData.LiqLEvap=Evap_struc.LiqL;
	gOutData.TPLEvap=Evap_struc.TPL;
	gOutData.VapLEvap=Evap_struc.VapL;
	//---------------------B.S.

	// calculate normalized residuals
	if(P->REV_CAL)//backward evaporator scheme
	{
	Y[0] =  2*(HPs[5].H-HPs[6].H)/(HPs[5].H+HPs[6].H);
	}
	else{//forward evaporator scheme
		if(P->EvaSign==1) {//cycle of simplified models
		Y[0] = (HPs[8].P-P_evpex);
		}
		else{
		Y[0] = 2*(HPs[8].H-H_evpex)/(HPs[8].H+H_evpex);//evaporator exit state
		}
	}

	Chg  =  m[0].m+m[1].m+m[2].m+m[3].m+m[5].m+(m[4].m+m[6].m);
	P->ChargMass = Chg;//calculated charge inventory, P->Charge, specified charge inventory

	if(P->SPEC_SUB)//specify subcooling, P->SPEC_SUB=1 or 2
	{
		if(P->Tsub<0.0)//two-phase 
		{
			if(P->SPEC_SUB ==2)//specify the state upstream of expansion device
			{
			TXP_sub.P=TXPs[5].P; TXP_sub.X=-P->Tsub; TXP_sub.T=PropertyTXPth(TSAT,TXP_sub);
			}
			else//specify the state upstream of condenser exit
			{
			TXP_sub.P=TXPs[4].P; TXP_sub.X=-P->Tsub; TXP_sub.T=PropertyTXPth(TSAT,TXP_sub);
			}
			H_Sub=PropertyTXPth(ENTH,TXP_sub);
			if(errorLog.IsError()) {
			errorLog.Add("MainFunc2_Fwd","H_sub");
			return 21;}
		}
		else//subcooled state
		{	
			if(P->SPEC_SUB ==2)//specify the state upstream of expansion device
			{
			TXP_sub.P=TXPs[5].P; TXP_sub.X=0; TXP_sub.T=Tsat[5]-P->Tsub;
			}
			else//specify the state at condenser exit
			{
			TXP_sub.P=TXPs[4].P; TXP_sub.X=0; TXP_sub.T=Tsat[4]-P->Tsub;
			}
			H_Sub=PropertyTXPth(ENTH,TXP_sub);
			if(errorLog.IsError()) {
			errorLog.Add("MainFunc2_Fwd","H_sub");
			return 20;
			}
		}

		if(P->SPEC_SUB ==2){//specify expansion device upstream state
			if(P->EvaSign==1){//cycle of simplified models
				if(P->Tsub>0.0) Y[1]=(Tsat[5]-TXPs[5].T)-P->Tsub;
				else Y[1]=(TXPs[5].X+P->Tsub)*10;
			}
			else{//Newton iteration method
				Y[1]=2*(HPs[5].H-H_Sub)/(HPs[5].H+H_Sub);
			}
		}
		else{//specify condenser exit
			if(P->EvaSign==1){//cycle of simplified models
				if(P->Tsub>0.0) Y[1]=(Tsat[4]-TXPs[4].T)-P->Tsub;
				else Y[1]=(TXPs[4].X+P->Tsub)*10;
			}
			else{
				Y[1]=2*(HPs[4].H-H_Sub)/(HPs[4].H+H_Sub);
			}
		}
	}
	else
	{
	Y[1] =- (Chg-P->Charge)/(P->Charge*1000);
	}


	// calculate convergence error
	CE=sqrt(pow(Y[0],2)+pow(Y[1],2));

	//B.S--------------------------
	P->Qual = TXPs[5].X;//output the condenser exit quality
	P->HP5 = HPs[5];//state upstream the orifice
	P->HP6 = HPs[6];//state downstream the orifice
	P->HP6.P=X[0]+375.09;//B.S.
	P->mr = mr;//mass flow rate calculated with the compressor model
	P->Tsup8 = TXPs[8].T-Tsat[8];//superheat at the evaporator exit

	if(!P->SPEC_SUB)
	{
	P->Tsub =(Tsat[4]-TXPs[4].T);//output condenser exit subcooling degree
	}
	//-------------------------B.S.


	//B.S. ---------------------------------------
	gOutData.devMr=0;
	gOutData.devH=Y[0];
	gOutData.devCharge=Y[1];
	gOutData.CE=CE;
	gOutData.inputH1=0;
	gOutData.inputP1=X[0];
	gOutData.inputP2=X[1];	
	gOutData.Charge=Chg;
	//--------------------------------------------B.S.

	printf("MainFunc2_Fwd: CE=%1.8lf (%1.8lf,%1.8lf)\n",CE,Y[0],Y[1]);//B.S.

	//printf("MainFunc2_Fwd: CE=%1.8lf\n",CE);

	return 0;
}


// Main function for system type #3.  This is a vapor compression
// cycle chiller at NIST with two stages of compression and a
// subcooler.
//
// X[0] = P8  (low stage inlet pressure)
// X[1] = P11 (low stage outlet and high stage inlet pressures)
// X[2] = P1  (high stage outlet pressure)
//
// STATES
// 1  - compressor high pressure stage outlet
// 2  - condenser inlet
// 3  - condenser outlet
// 4  - desuperheater/subcooler high pressure inlet
// 5  - main TxV inlet
//    - desuperheater/subcooler high pressure outlet
// 6  - evaporator inlet
// 7  - evaporator outlet
// 8  - compressor low stage inlet
// 9  - desuperheater/subcooler low pressure outlet
// 10 - compressor low stage outlet
// 11 - compressor high stage inlet
// 12 - desuperheater/subcooler low pressure inlet
int MainFunc3(double *X,double *Y,void *Params)
{
	static double mr1,mr2,mr3;
	static double Ones[5]={1,1,1,1,1};
	static MASS m[20];
	static HP HPs[20];
	static TXP TXPs[20];
	static double Tsat[20];
	TXP TXP_prop={0,0,0};

	InVars* P = (InVars*)Params;
	static CGP Cond_struc;
	static ETdim Evap_struc;
	//
	// Keep track of how much time is spent crunching each component.
	// This provides the ability to determine which component should
	// be made more efficient.
	static int init=1;
	static unsigned long int saveTime;
	if(init) {

		#ifdef DEBUG_MAINFUNC
		{
		FILE *fp = fopen("mainfunc.dbg","w");
		fprintf(fp,"Opening...\n");
		fclose(fp);
		}
		{
		FILE *fp = fopen("mainfunc.tst","w");
		fprintf(fp,"Opening...\n");
		fclose(fp);
		}
		#endif

		saveTime = MilliTime();
		init=0;
	}

	//gOutData.ClearAll();

	{
	// atan() function maps -inf to +inf to +pi/2 to -pi/2.  Therefore,
	// it is used to allow the numerical method to propose valid guesses
	// for -inf to +inf.  The contraint is Pevap < Pint < Pcond.
	double pi = acos(-1.0);
	HPs[8].P = PMINth + (PMAXth-PMINth)*(atan(X[0])+pi/2)/pi;
	HPs[9].P = HPs[8].P + (PMAXth-HPs[8].P)*(atan(X[1])+pi/2)/pi;
	HPs[1].P = HPs[9].P + (PMAXth-HPs[9].P)*(atan(X[2])+pi/2)/pi;
	}

	#ifdef DEBUG_MAINFUNC
	{
	FILE *fp = fopen("mainfunc.tst","a");
	fprintf(fp,"P8=%lf P9=%lf P1=%lf\n",HPs[8].P,HPs[9].P,HPs[1].P);
	fclose(fp);
	}
	#endif

	HPs[10].P = HPs[9].P;
	HPs[11].P = HPs[9].P;


	// low stage inlet enthalpy (TxV 1)
	{
		TXP_prop.P=HPs[8].P;
		TXP_prop.X=1;
		Tsat[8] = PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(HPs[8].P);
		if(errorLog.IsError()) {
			errorLog.Add("MainFunc3","Tsat[8]");
			return 0;
		}
		TXPs[8] = toTXP(Tsat[8]+TSH_LP,1,HPs[8].P);
		HPs[8] = TXPtoHP(TXPs[8]);
		if(errorLog.IsError()) {
			char str[128];
			sprintf(str,"HPs[8], (TXPs[8]=%lf,%lf,%lf)",TXPs[8].T,TXPs[8].X,TXPs[8].P);
			errorLog.Add("MainFunc3",str);
			return 0;
		}
		gOutData.H[8] = HPs[8].H;
		gOutData.T[8] = TXPs[8].T;
		gOutData.X[8] = TXPs[8].X;
		gOutData.P[8] = TXPs[8].P;
		gOutData.Tsat[8] = Tsat[8];
		#ifdef DEBUG_MAINFUNC
			{
			FILE *fp = fopen("mainfunc.dbg","a");
			fprintf(fp,"\n");
			fprintf(fp,"T8sat  = %8.3lf\n",Tsat[8]);
			fprintf(fp,"T8     = %8.3lf\n",Tsat[8]+11.1);
			fprintf(fp,"HTXP8  = %8.3lf,%8.3lf,%8.3lf,%8.3lf\n\n",0.001*HPs[8].H,TXPs[8].T,TXPs[8].X,TXPs[8].P);
			fclose(fp);
			}
		#endif
	}

	// low stage compressor
	{
		{
		unsigned long int currentTime = MilliTime();
		cycleTime[cyclePosition] += currentTime - saveTime;
		saveTime = currentTime;
                }

		cyclePosition=1;

		double Ei;
		CompressorCCLowStage(P->filename,HPs[8].H,HPs[8].P,HPs[10].P,HPs[1].P,HPs[10].H,mr1,Ei,m[1],Ones);
		if(errorLog.IsError()) {
			errorLog.Add("MainFunc3","CompressorCCLowStage");
			return 0;
		}

		TXP_prop.P=HPs[10].P;
		TXP_prop.X=0;
		Tsat[10] = PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(HPs[10].P);
		if(errorLog.IsError()) {
			errorLog.Add("MainFunc3","Tsat[8]");
			return 0;
		}
		TXPs[10] = HPtoTXP(HPs[10]);
		if(errorLog.IsError()) {
			errorLog.Add("MainFunc3","TXPs[10]");
			return 0;
		}
		gOutData.Ecmp = Ei;
		gOutData.mr[1] = mr1;
		gOutData.H[10] = HPs[10].H;
		gOutData.T[10] = TXPs[10].T;
		gOutData.X[10] = TXPs[10].X;
		gOutData.P[10] = TXPs[10].P;
		gOutData.Tsat[10] = Tsat[10];
		#ifdef DEBUG_MAINFUNC
			{
			FILE *fp = fopen("mainfunc.dbg","a");
			fprintf(fp,"Low stage compressor.\n");
			fprintf(fp,"mr     = %8.3lf\n",mr1);
			fprintf(fp,"Ei     = %8.3lf\n",Ei);
			fprintf(fp,"HP8    = %8.3lf,%8.3lf\n",HPs[8].H/1000,HPs[8].P);
			fprintf(fp,"HTXP10 = %8.3lf,%8.3lf,%8.3lf,%8.3lf\n\n",0.001*HPs[10].H,TXPs[10].T,TXPs[10].X,TXPs[10].P);
			fclose(fp);
			}
		#endif
	}

	// high stage inlet enthalpy (TxV 2)
	{

		TXP_prop.P=HPs[11].P;
		TXP_prop.X=0;
		Tsat[11] = PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(HPs[11].P);
		if(errorLog.IsError()) {
			errorLog.Add("MainFunc3","5");
			return 0;
		}
		TXPs[11] = toTXP(Tsat[11]+TSH_HP,1,HPs[11].P);
		HPs[11] = TXPtoHP(TXPs[11]);
		if(errorLog.IsError()) {
			errorLog.Add("MainFunc3","6");
			return 0;
		}
		gOutData.H[11] = HPs[11].H;
		gOutData.T[11] = TXPs[11].T;
		gOutData.X[11] = TXPs[11].X;
		gOutData.P[11] = TXPs[11].P;
		gOutData.Tsat[11] = Tsat[11];
	}
	//
	// high stage compressor
	{
		{
		unsigned long int currentTime = MilliTime();
		cycleTime[cyclePosition] += currentTime - saveTime;
		saveTime = currentTime;
		}

		cyclePosition=2;

		CompressorCCHighStage(P->filename,HPs[11].H,HPs[11].P,HPs[1].P,HPs[1].H);
		if(errorLog.IsError()) {
			errorLog.Add("MainFunc3","CompressorCCHighStage");
			return 0;
		}

		TXP_prop.P=HPs[1].P;
		TXP_prop.X=1;
		Tsat[1] = PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(HPs[1].P);
		if(errorLog.IsError()) {
			errorLog.Add("MainFunc3","Tsat[1]");
			return 0;
		}
		TXPs[1] = HPtoTXP(HPs[1]);
		if(errorLog.IsError()) {
			errorLog.Add("MainFunc3","TXPs[1]");
			return 0;
		}
		gOutData.H[1] = HPs[1].H;
		gOutData.T[1] = TXPs[1].T;
		gOutData.X[1] = TXPs[1].X;
		gOutData.P[1] = TXPs[1].P;
		gOutData.Tsat[1] = Tsat[1];
		#ifdef DEBUG_MAINFUNC
			{
			FILE *fp = fopen("mainfunc.dbg","a");
			fprintf(fp,"High stage compressor.\n");
			fprintf(fp,"HP11   = %8.3lf,%8.3lf\n",HPs[11].H/1000,HPs[11].P);
			fprintf(fp,"HTXP1  = %8.3lf,%8.3lf,%8.3lf,%8.3lf\n\n",0.001*HPs[1].H,TXPs[1].T,TXPs[1].X,TXPs[1].P);
			fclose(fp);
			}
		#endif
	}
	//
	// volume flow rate of low stage is twice the high stage
	{
		double v8  = PropertyTXPth(VOL,TXPs[8]);
		if(errorLog.IsError()) {
			char str[128];
			sprintf(str,"v8: P(8,9,1)=(%lf,%lf,%lf) TXP[8]=(%lf,%lf,%lf)",HPs[8].P,HPs[9].P,HPs[1].P,TXPs[8].T,TXPs[8].X,TXPs[8].P);
			errorLog.Add("MainFunc3",str);
			return 0;
		}
		double v11 = PropertyTXPth(VOL,TXPs[11]);
		if(errorLog.IsError()) {
			errorLog.Add("MainFunc3","v11");
			return 0;
		}
		mr2 = 0.5*mr1*v8/v11;
		gOutData.mr[2] = mr2;
	}
	//
	// Refrigerant entering intermediate stage from subcooler.
	// Use mass and energy balance.
	{
		mr3 = mr2 - mr1;
		gOutData.mr[3] = mr3;
		HPs[9].H = (mr2*HPs[11].H-mr1*HPs[10].H)/mr3;

		TXP_prop.P=HPs[9].P;
		TXP_prop.X=0;
		Tsat[9] = PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(HPs[9].P);
		if(errorLog.IsError()) {
			errorLog.Add("MainFunc3","Tsat[9]");
			return 0;
		}
		TXPs[9] = HPtoTXP(HPs[9]);
		if(errorLog.IsError()) {
			char str[128];
			sprintf(str,"TXPs[9] P8=%lf P10=%lf P1=%lf, H11=%lf H10=%lf m1=%lf m2=%lf m3=%lf\n",HPs[8].P,HPs[10].P,HPs[1].P,HPs[11].H,HPs[10].H,mr1,mr2,mr3);
			errorLog.Add("Mainfunc3",str);
			return 0;
		}
		gOutData.H[9] = HPs[9].H;
		gOutData.T[9] = TXPs[9].T;
		gOutData.X[9] = TXPs[9].X;
		gOutData.P[9] = TXPs[9].P;
		gOutData.Tsat[9] = Tsat[9];
		#ifdef DEBUG_MAINFUNC
			{
			FILE *fp = fopen("mainfunc.dbg","a");
			fprintf(fp,"Intermediate stage.\n");
			fprintf(fp,"mr2 = %lf\n",mr2);
			fprintf(fp,"mr3 = %lf\n",mr3);
			fprintf(fp,"HTXP9  = %8.3lf,%8.3lf,%8.3lf,%8.3lf\n\n",0.001*HPs[9].H,TXPs[9].T,TXPs[9].X,TXPs[9].P);
			fclose(fp);
			}
		#endif
	}
	//
	// oil separator
	{
		{
		unsigned long int currentTime = MilliTime();
		cycleTime[cyclePosition] += currentTime - saveTime;
		saveTime = currentTime;
		}

		cyclePosition=3;

		HPs[2] = HPs[1];
		m[3].m = m[3].V = 0.0;
	}
	//
	// condenser
	{
		{
		unsigned long int currentTime = MilliTime();
		cycleTime[cyclePosition] += currentTime - saveTime;
		saveTime = currentTime;
		}

		cyclePosition=4;

		double Tao;
		Condenser(P->filename,mr2,HPs[2],P->Tai,P->Goa,&HPs[3],&Tao,&m[4],&Cond_struc,Ones);
		if(errorLog.IsError()) {
			errorLog.Add("MainFunc3","Condenser");
			return 1;
		}

		gOutData.mc = P->Goa;
		gOutData.Tci = P->Tai;
		gOutData.Tco = Tao;

		TXP_prop.P=HPs[2].P;
		TXP_prop.X=0;
		Tsat[2] = PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(HPs[2].P);
		if(errorLog.IsError()) {
			errorLog.Add("MainFunc3","Tsat[2]");
			return 0;
		}
		TXPs[2] = HPtoTXP(HPs[2]);
		if(errorLog.IsError()) {
			errorLog.Add("MainFunc3","TXPs[2]");
			return 0;
		}
		gOutData.H[2] = HPs[2].H;
		gOutData.T[2] = TXPs[2].T;
		gOutData.X[2] = TXPs[2].X;
		gOutData.P[2] = TXPs[2].P;
		gOutData.Tsat[2] = Tsat[2];

		TXP_prop.P=HPs[3].P;
		TXP_prop.X=0;
		Tsat[3] = PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(HPs[3].P);
		if(errorLog.IsError()) {
			errorLog.Add("MainFunc3","Tsat[3]");
			return 0;
		}
		TXPs[3] = HPtoTXP(HPs[3]);
		if(errorLog.IsError()) {
			errorLog.Add("MainFunc3","TXPs[3]");
			return 0;
		}
		gOutData.H[3] = HPs[3].H;
		gOutData.T[3] = TXPs[3].T;
		gOutData.X[3] = TXPs[3].X;
		gOutData.P[3] = TXPs[3].P;
		gOutData.Tsat[3] = Tsat[3];
		#ifdef EXTRA_REPORTING
			printf("Condenser.\n");
			printf("HTXP2  = %8.3lf,%8.3lf,%8.3lf,%8.3lf\n",0.001*HPs[2].H,TXPs[2].T,TXPs[2].X,TXPs[2].P);
			printf("HTXP3  = %8.3lf,%8.3lf,%8.3lf,%8.3lf\n\n",0.001*HPs[3].H,TXPs[3].T,TXPs[3].X,TXPs[3].P);
		#endif

		gOutData.Qc = mr2*(HPs[2].H-HPs[3].H);
	}
	//
	// receiver
	{
		{
		unsigned long int currentTime = MilliTime();
		cycleTime[cyclePosition] += currentTime - saveTime;
		saveTime = currentTime;
		}

		cyclePosition=5;

		HPs[4] = HPs[3];

		TXP_prop.P=HPs[4].P;
		TXP_prop.X=0;

		Tsat[4] = PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(HPs[4].P);
		if(errorLog.IsError()) {
			errorLog.Add("MainFunc3","Tsat[4]");
			return 0;
		}
		TXPs[4] = HPtoTXP(HPs[4]);
		if(errorLog.IsError()) {
			errorLog.Add("MainFunc3","TXPs[4]");
			return 0;
		}
		gOutData.H[4] = HPs[4].H;
		gOutData.T[4] = TXPs[4].T;
		gOutData.X[4] = TXPs[4].X;
		gOutData.P[4] = TXPs[4].P;
		gOutData.Tsat[4] = Tsat[4];

		m[5].m = m[5].V = 0.0;
	}
	//
	// subcooler
	{
		{
		unsigned long int currentTime = MilliTime();
		cycleTime[cyclePosition] += currentTime - saveTime;
		saveTime = currentTime;
		}

		cyclePosition=6;

		Subcooler(P->filename,mr1,HPs[4],mr3,HPs[9],HPs[5],HPs[12],m[6]);
		if(errorLog.IsError()) {
			errorLog.Add("MainFunc3","Subcooler");
			return 0;
		}

		TXP_prop.P=HPs[5].P;
		TXP_prop.X=0;

		Tsat[5] = PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(HPs[5].P);
		if(errorLog.IsError()) {
			errorLog.Add("MainFunc3","Tsat[5]");
			return 0;
		}
		TXPs[5] = HPtoTXP(HPs[5]);
		if(errorLog.IsError()) {
			errorLog.Add("MainFunc3","TXPs[5]");
			return 0;
		}
		gOutData.H[5] = HPs[5].H;
		gOutData.T[5] = TXPs[5].T;
		gOutData.X[5] = TXPs[5].X;
		gOutData.P[5] = TXPs[5].P;
		gOutData.Tsat[5] = Tsat[5];

		TXP_prop.P=HPs[12].P;
		TXP_prop.X=0;
		Tsat[12] = PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(HPs[12].P);
		if(errorLog.IsError()) {
			errorLog.Add("MainFunc3","Tsat[12]");
			return 0;
		}
		TXPs[12] = HPtoTXP(HPs[12]);
		if(errorLog.IsError()) {
			errorLog.Add("MainFunc3","TXPs[12]");
			return 0;
		}
		gOutData.H[12] = HPs[12].H;
		gOutData.T[12] = TXPs[12].T;
		gOutData.X[12] = TXPs[12].X;
		gOutData.P[12] = TXPs[12].P;
		gOutData.Tsat[12] = Tsat[12];
		#ifdef EXTRA_REPORTING
			printf("Subcooler.\n");
			printf("HTXP5  = %8.3lf,%8.3lf,%8.3lf,%8.3lf\n",0.001*HPs[5].H,TXPs[5].T,TXPs[5].X,TXPs[5].P);
			printf("HTXP12 = %8.3lf,%8.3lf,%8.3lf,%8.3lf\n\n",0.001*HPs[12].H,TXPs[12].T,TXPs[12].X,TXPs[12].P);
		#endif
	}

	// accumulator
	{
		{
		unsigned long int currentTime = MilliTime();
		cycleTime[cyclePosition] += currentTime - saveTime;
		saveTime = currentTime;
		}

		cyclePosition=7;

		HPs[7] = HPs[8];

		TXP_prop.P=HPs[7].P;
		TXP_prop.P=0;
		Tsat[7] = PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(HPs[7].P);
		if(errorLog.IsError()) {
			errorLog.Add("MainFunc3","Tsat[7]");
			return 0;
		}
		TXPs[7] = HPtoTXP(HPs[7]);
		if(errorLog.IsError()) {
			errorLog.Add("MainFunc3","TXPs[7]");
			return 0;
		}
		gOutData.H[7] = HPs[7].H;
		gOutData.T[7] = TXPs[7].T;
		gOutData.X[7] = TXPs[7].X;
		gOutData.P[7] = TXPs[7].P;
		gOutData.Tsat[7] = Tsat[7];

		m[7].m = m[7].V = 0.0;
	}
	//
	// evaporator
	{
		{
		unsigned long int currentTime = MilliTime();
		cycleTime[cyclePosition] += currentTime - saveTime;
		saveTime = currentTime;
		}

		cyclePosition=8;

		double Tewo;
		Evaporator(P->filename,mr1,HPs[7],P->Gra,P->TPi.T,HPs[6],Tewo,m[8],Ones);
		if(errorLog.IsError()) {
			errorLog.Add("MainFunc3","Evaporator");
			return 2;
		}

		gOutData.me = P->Gra;
		gOutData.Tei = P->TPi.T;
		gOutData.Teo = Tewo;

		TXP_prop.P=HPs[6].P;
		TXP_prop.X=0;
		Tsat[6] = PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(HPs[6].P);
		if(errorLog.IsError()) {
			errorLog.Add("MainFunc3","Tsat[6]");
			return 0;
		}
		TXPs[6] = HPtoTXP(HPs[6]);
		if(errorLog.IsError()) {
			errorLog.Add("MainFunc3","TXPs[6]");
			return 0;
		}
		gOutData.H[6] = HPs[6].H;
		gOutData.T[6] = TXPs[6].T;
		gOutData.X[6] = TXPs[6].X;
		gOutData.P[6] = TXPs[6].P;
		gOutData.Tsat[6] = Tsat[6];
		#ifdef EXTRA_REPORTING
			printf("Evaporator.\n");
			printf("HTXP6  = %8.3lf,%8.3lf,%8.3lf,%8.3lf\n",0.001*HPs[6].H,TXPs[6].T,TXPs[6].X,TXPs[6].P);
			printf("Tewi/o = %8.3lf,%8.3lf\n",Tewi,Tewo);
		#endif
		gOutData.Qe = mr1*(HPs[7].H-HPs[6].H);
	}

	// residuals

	// dh=0 across TxV1
	Y[0]=(HPs[5].H-HPs[6].H)/(HPs[5].H+HPs[6].H);

	// dh=0 across TxV2
	Y[1]=(HPs[4].H-HPs[12].H)/(HPs[4].H+HPs[12].H);

	// charge is conserved
	{
	if(fabs(P->Charge)<1e-20) {
		errorLog.Add("MainFunc3","Specified charge level is invalid");
		return 3;
	}
	double Charge = 0;
	for(int i=1;i<=8;i++) Charge += m[i].m;
	gOutData.Charge = Charge;
	Y[2]=-(Charge-P->Charge)/P->Charge;
	}

	return 0;
}

//validate the component models with measured data
int MainFunc_Validating(void *Params)
{
	double mr,Ei,Tao,Aflow;
	double Ones[5]={1,1,1,1,1};
	double CompPrms[4]={1,1,1,50};//CompPrms[3]=50, for lubricant temperature
	MASS m[9];
	HP HPs[9];
	TP TPo;
	InVars *P;
	CGP Cond_struc;//B.S. condenser construct
	ETdim Evap_struc;//B.S. evaporator construct
	double X[3]={0,0,0};
	TXP TXP_prop={0,0,0};
	static int FirstTime=0;
	int i=0;

	std::ofstream outstrm; 
	outstrm.open("./outputDoc/components_out.xls",std::ios::app);

	if(FirstTime==0)
	{
		outstrm	<<"Tci"
		<<"	"<<"Tei"<<"	"<<"HumPi"<<"	"<<"Gca"<<"	"<<"Gea"<<"	"<<"Charge"<<"	"<<"Tsup1"<<"	"<<"Tsub5"<<"	"<<"mr"
		<<"	"<<"dev_mrcom"<<"	"<<"dev_Ei"<<"	"<<"h2_m"<<"	"<<"h2_c"<<"	"<<"dev_h2"
		<<"	"<<"DPhotgas_m"<<"	"<<"DPhotgas_c"
		<<"	"<<"DQhotgas_m"<<"	"<<"DQhotgas_c"
		<<"	"<<"DPsuc_m"<<"	"<<"DPsuc_c"
		<<"	"<<"DQsuc_m"<<"	"<<"DQsuc_c"
		<<"	"<<"DPcond_m"<<"	"<<"DPcond_c"<<"	"<<"Dev_DPcond"
		<<"	"<<"DHcond_m"<<"	"<<"DHcond_c"<<"	"<<"Dev_DHcond"<<"	"<<"Tsub4_c"<<"	"<<"Tco_c"
		<<"	"<<"DPliq_m"<<"	"<<"DPliq_c"
		<<"	"<<"DQliq_m"<<"	"<<"DQliq_c"<<"	"<<"Tsub5_c"
		<<"	"<<"DPevap_m"<<"	"<<"DPevap_c"<<"	"<<"Dev_DPevap"
		<<"	"<<"DHevap_m"<<"	"<<"DHevap_c"<<"	"<<"Dev_DHevap"<<"	"<<"Tsup8_c"<<"	"<<"Teo_c"<<"	"<<"HumPo_c";
		for(int i=0;i<4;i++){outstrm<<"	"<<"DPE_"<<i;}
		for(int i=0;i<4;i++) {outstrm<<"	"<<"Tsup8_"<<i;}

		outstrm<<"	"<<"DP5_m"<<"	"<<"DP5_c"
		<<"	"<<"dev_mrR410Apure"<<"	"<<"dev_mrMulRef"<<"	"<<"dev_mrR410Aoil"<<"	"<<"Chg"<<"	"<<"dev_charge"<<"	"<<"LiqL"
		<<"	"<<"Ecmp"<<"	"<<"Qc"<<"	"<<"Qe"<<"	"<<"dev_hotgas"<<"	"<<"dev_liqline"<<"	"<<"dev_suc"<<std::endl;
		FirstTime=1;
	}


	P=(InVars*)Params;

	outstrm	<<P->Tai
	<<"	"<<P->TPi.T<<"	"<<P->TPi.P<<"	"<<P->Goa<<"	"<<P->Gra<<"	"<<P->Charge<<"	"<<P->Tsup<<"	"<<P->Tsub<<"	"<<P->mr<<"	";


	X[1]=P->P[1];
	X[2]=P->P[2];

	if(P->X[0]<1) //P->X[0] means the quality at the compressor suction
	{
	TXP_prop.X=P->X[0];
	TXP_prop.P=P->P[1];
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);
	}
	else 
	{
	TXP_prop.X=1;
	TXP_prop.P=P->P[1];
	TXP_prop.T=P->T[1];
	}

	X[0] = PropertyTXPth(ENTH,TXP_prop);
	const double h1=X[0];

	if(errorLog.IsError()) {
		errorLog.Add("MainFunc1_Components","Compressor_prop");
		return 1;
	}

	cyclePosition = 9;
	CompPrms[3]=P->CompT;
	Compressor(P->filename,X[0],X[1],X[2],P->Tai,&mr,&HPs[2],&HPs[1],&Ei,&m[0],CompPrms);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc1_Components","Compressor");
		return 1;
	}

	TXP_prop.X=1;
	TXP_prop.P=P->P[2];
	TXP_prop.T=P->T[2];
	const double h2= PropertyTXPth(ENTH,TXP_prop);
	
	const double dev_mrcom=(mr-P->mr)/P->mr*100;
	const double dev_Ei=(Ei-P->mr*(h2-X[0]))/(P->mr*(h2-X[0]))*100;
	const double dev_h2 = (h2-HPs[2].H)/h2*100;
	outstrm<<dev_mrcom<<"	"<<dev_Ei<<"	"<<h2<<"	"<<HPs[2].H<<"	"<<dev_h2<<"	";

	HPs[2].P=P->P[2];
	HPs[2].H=h2;

	HotGasLine(HPs[2],P->mr,&HPs[3],&m[1],P);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc1_Components","Hotgasline");
		return 6;
	}
	
	TXP_prop.X=1;
	TXP_prop.P=P->P[3];
	TXP_prop.T=P->T[3];
	const double h3=PropertyTXPth(ENTH,TXP_prop);

	const double dev_P_hotgasline=((HPs[2].P-HPs[3].P)-(P->P[2]-P->P[3]))/(P->P[2]-P->P[3])*100;
	const double dev_h_hotgasline=((HPs[2].H-HPs[3].H)-(h2-h3))/(h2-h3)*100;
	
	outstrm<<P->P[2]-P->P[3]<<"	"<<(HPs[2].P-HPs[3].P)<<"	";
	outstrm<<(h2-h3)*P->mr<<"	"<<(HPs[2].H-HPs[3].H)*P->mr<<"	";

	HPs[1].H=h1;
	HPs[1].P=P->P[1];

	SuctionLine(HPs[1],P->mr,&HPs[8],&m[2]);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc1_Components","suctionline");
		return 7;
	}

	if(P->X[3]<1) //P->X[3] quality at the exit of the evaporator
	{
	TXP_prop.X=P->X[3];
	TXP_prop.P=P->P[8];
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);
	}
	else 
	{
	TXP_prop.X=1;
	TXP_prop.P=P->P[8];
	TXP_prop.T=P->T[8];
	}

	const double h8 = PropertyTXPth(ENTH,TXP_prop);
	
	const double dev_P_suctionline=((HPs[8].P-HPs[1].P)-(P->P[8]-P->P[1]))/(P->P[8]-P->P[1])*100;
	const double dev_h_suctionline=((HPs[1].H-HPs[8].H)-(h1-h8))/(h1-h8)*100;
	outstrm<<P->P[8]-P->P[1]<<"	"<<(HPs[8].P-HPs[1].P)<<"	";
	outstrm<<(h1-h8)*P->mr<<"	"<<(HPs[1].H-HPs[8].H)*P->mr<<"	";

	cyclePosition = 4;

	Ones[0] = P->Cadj;

	while(Condenser(P->filename,P->mr,HPs[3],P->Tai,P->Goa,&HPs[4],&Tao,&m[3],&Cond_struc,Ones)) {	}
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc1_Components","Condenser");
		return 8;
	}

	Ones[0]=1;

	if(P->X[1]>0) //P->X[1], quality at the condenser exit
	{
	TXP_prop.X=P->X[1];
	TXP_prop.P=P->P[4];
	TXP_prop.T=P->T[4];
	}
	else 
	{
	TXP_prop.X=1;
	TXP_prop.P=P->P[4];
	TXP_prop.T=P->T[4];
	}
		

	TXP TXP4;
	TXP4 = HPtoTXP(HPs[4]);
	
	TXP_prop=TXP4;
	TXP_prop.X=0;

	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

	const double Tsub4=TXP_prop.T-TXP4.T;
		
	LiquidLine(HPs[4],P->mr,&HPs[5],&m[5]);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc1_Components","liquidline");
		return 11;
	}

	if(P->X[2]>0) //P->X[2], quality at the exit of the liquid line
	{
	TXP_prop.X=P->X[2];
	TXP_prop.P=P->P[4];
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);
	}
	else 
	{
	TXP_prop.X=0;
	TXP_prop.P=P->P[4];
	TXP_prop.T=P->T[4];
	}
	const double h4 = PropertyTXPth(ENTH,TXP_prop);
	const double h5 = h4;

	const double dev_P_cond=((HPs[3].P-HPs[4].P)-(P->P[3]-P->P[4]))/(P->P[3]-P->P[4])*100;
	const double dev_h_cond=((HPs[2].H-HPs[4].H)-(h2-h4))/(h2-h4)*100;
	
	outstrm<<P->P[3]-P->P[4]<<"	"<<(HPs[3].P-HPs[4].P)<<"	"<<dev_P_cond<<"	";
	outstrm<<(h2-h4)<<"	"<<(HPs[2].H-HPs[4].H)<<"	"<<dev_h_cond<<"	"<<Tsub4<<"	"<<Tao<<"	";


	TXP TXP5;
	TXP5 = HPtoTXP(HPs[5]);
	
	TXP_prop=TXP5;
	TXP_prop.X=0;

	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

	const double Tsub5=TXP_prop.T-TXP5.T;

	const double dev_P_liquidline=((HPs[4].P-HPs[5].P)-(P->P[4]-P->P[5]))/(P->P[4]-P->P[5])*100;
	//const double dev_h_liquidline=((HPs[4].H-HPs[5].H)-(h4-h5))/(h4-h5)*100;
	outstrm<<P->P[4]-P->P[5]<<"	"<<(HPs[4].P-HPs[5].P)<<"	";
	outstrm<<(h4-h5)*P->mr<<"	"<<(HPs[4].H-HPs[5].H)*P->mr<<"	"<<Tsub5<<"	";

	cyclePosition = 8;
	
	HPs[7].H=h5;//evaporator inlet enthalp equal to liquid line outlet
	HPs[8].P=HPs[8].P;//P->P[8];
	HPs[8].H=HPs[8].H;//h8;

	Evap_struc.REV=P->REV_CAL;
	Ones[0] = P->Eadj;
	while(Evaporator(P->filename,P->mr,&HPs[8],P->Gra,P->TPi,&HPs[7],&TPo,&m[4],&Aflow,&Evap_struc,Ones)) { }
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc1_Components","Evaporator");
		return 15;
	}
	
	Ones[0]=1;


	TXP TXP8;
	TXP8 = HPtoTXP(HPs[8]);
	
	TXP_prop=TXP8;
	TXP_prop.X=1;

	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

	double Tsup8=TXP8.T-TXP_prop.T;

	const double dev_P_Evap=((Evap_struc.DPr[0]/1.135)-(P->P[7]-P->P[8]))/(P->P[7]-P->P[8])*100;
	const double dev_h_Evap=((HPs[8].H-HPs[7].H)-(h1-h5))/(h1-h5)*100;
	outstrm<<P->P[7]-P->P[8]<<"	"<<(HPs[7].P-HPs[8].P)<<"	"<<dev_P_Evap<<"	";
	outstrm<<(h1-h5)<<"	"<<(HPs[8].H-HPs[7].H)<<"	"<<dev_h_Evap<<"	"<<Tsup8<<"	"<<TPo.T<<"	"<<TPo.P<<"	";

	for(int i=0;i<4;i++){outstrm<<Evap_struc.DPr[i]<<"	";}
	double sup8[4];
	for(int i=0;i<4;i++)
	{
	HP HP_out;
	HP_out.H=Evap_struc.Hout8[i];
	HP_out.P=P->P[8];
	TXP8 = HPtoTXP(HP_out);
	TXP_prop=TXP8;
	TXP_prop.X=1;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);
	Tsup8=TXP8.T-TXP_prop.T;
	sup8[i]=Tsup8;
	outstrm<<Tsup8<<"	";
	}

	
	HPs[7].P=P->P[7];
	HPs[7].H=h5;
	DistributionTubes(HPs[7],P->mr,&HPs[6],&m[6]);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc1_Components","Distribution Tubes");
		return 18;
	}
	
	const double dev_P_distributor=((HPs[6].P-HPs[7].P)-(P->P[6]-P->P[7]))/(P->P[6]-P->P[7])*100;
	outstrm<<P->P[6]-P->P[7]<<"	"<<(HPs[6].P-HPs[7].P)<<"	";
	
	cyclePosition = 10;

	HPs[5].H=h5;
	HPs[5].P=P->P[5];
	HPs[6].H=h5;
	HPs[6].P=P->P[1]+375.09;

	mr = Short_Tube_R410A(P->filename,HPs[5],HPs[6],Ones);
	const double dev_mrR410Apure=(mr-P->mr)/P->mr*100;
	Ones[0] = 0.9975;
	mr = Short_Tube_R410A_Oil(P->filename,HPs[5],HPs[6],Ones);
	const double dev_mrR410Aoil=(mr-P->mr)/P->mr*100;
	Ones[0]=1.026;
	mr = Short_Tube_FlowDirec(P->filename,HPs[5],HPs[6],Ones);
	const double dev_mrMulRef=(mr-P->mr)/P->mr*100;


	const double Chg  =  m[0].m+m[1].m+m[2].m+m[3].m+m[4].m+m[5].m+m[6].m;
	const double dev_charge=(P->Charge-Chg);
	outstrm<<dev_mrR410Apure<<"	"<<dev_mrMulRef<<"	"<<dev_mrR410Aoil
		<<"	"<<Chg<<"	"<<dev_charge<<"	"<<Cond_struc.LiqL<<"	";

	const double Qc=P->mr*(h2-h4);
	const double Qe=P->mr*(h1-h5);

	outstrm<<Ei<<"	"<<Qc<<"	"<<Qe<<"	"<<P->mr*(h2-h3)/Qc*100<<"	"<<P->mr*(h4-h5)/Qc*100<<"	"<<P->mr*(h1-h8)/Qe*100<<std::endl;//actual value

	outstrm.close();
	printf("MainFunc_Components: dev=%1.8lf \n",dev_mrR410Apure/*dev_h_conddev_h_conddev_mrexp*/);//B.S.
	printf("MainFunc_Components: dev=%1.8lf \n",dev_h_cond/*dev_h_conddev_h_conddev_mrexp*/);//B.S.

	return 0;
}

//Tune the system model to get condenser and evaporator heat transer adjustment and system charge
int MainFunc_Tuning(void *Params)
{
	double mr,Ei,Tao,Aflow;
	double Prms[10] = {1,1,1,1,1,1,1,1,1,1};
	double CompPrms[4]={1,1,1,50};//CompPrms[3]=50, for a default lubricant temperature
	double Tsub4=0,Tsub5=0;//subcooling 
	MASS m[9];
	HP HPs[9];
	TP TPo;
	InVars *P;
	CGP Cond_struc;//B.S. condenser construct
	ETdim Evap_struc;//B.S. evaporator construct
	double X[3]={0,0,0};
	TXP TXP_prop={0,0,0};
	static int FirstTime=0;

	std::ofstream outstrm; 
	outstrm.open("./outputDoc/Tune_out.xls",std::ios::app);
	
	//print title
	if(FirstTime==0)
	{
		outstrm<<"OutdT"<<"	"<<"IndT"<<"	"<<"IndRH"<<"	"<<"Goa"<<"	"<<"Gia"
			<<"	"<<"Charge"<<"	"<<"Tsup"<<"	"<<"Tsub"<<"	"<<"mr"<<"	";
		outstrm<<"Coeff_com"<<"	"<<"dev_Ei"<<"	"<<"dev_h2"<<"	";
		outstrm<<"Coeff_cond"<<"	";
		outstrm<<"Coeff_evap"<<"	";
		outstrm<<"mr_exp"<<"	"<<"Coeff_Exp"<<"	"<<"Chg"<<"	"<<"dev_charge"<<"	";
		outstrm<<"m_com"<<"	"<<"m_hotgas"<<"	"<<"m_suc"<<"	"<<"m_con"	
			<<"	"<<"m_eva"<<"	"<<"m_liq"<<"	"<<"m_distri"<<"	"<<"LiqL";
		outstrm<<std::endl;
		FirstTime=1;
	}


	P=(InVars*)Params;

	outstrm	<<P->Tai
	<<"	"<<P->TPi.T<<"	"<<P->TPi.P<<"	"<<P->Goa<<"	"<<P->Gra<<"	"
	<<P->Charge<<"	"<<P->Tsup<<"	"<<P->Tsub<<"	"<<P->mr<<"	";


	X[1]=P->P[1];
	X[2]=P->P[2];

	if(P->X[0]<1) //P->X[0] means the quality at the compressor suction
	{
	TXP_prop.X=P->X[0];
	TXP_prop.P=P->P[1];
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);
	}
	else 
	{
	TXP_prop.X=1;
	TXP_prop.P=P->P[1];
	TXP_prop.T=P->T[1];
	}

	X[0] = PropertyTXPth(ENTH,TXP_prop);
	const double h1=X[0];

	if(errorLog.IsError()) {
		errorLog.Add("MainFunc1_Tuning","Compressor_prop");
		return 1;
	}

	CompPrms[3]=P->CompT;//Lubricant temperature inputted from a data file
	Compressor(P->filename,X[0],X[1],X[2],P->Tai,&mr,&HPs[2],&HPs[1],&Ei,&m[0],CompPrms);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc1_Tuning","Compressor");
		return 1;
	}

	TXP_prop.X=1;
	TXP_prop.P=P->P[2];
	TXP_prop.T=P->T[2];
	const double h2= PropertyTXPth(ENTH,TXP_prop);
	
	const double Coeff_com=mr/P->mr;
	const double dev_Ei=(Ei-P->mr*(h2-X[0]))/(P->mr*(h2-X[0]))*100;
	const double dev_h2 = (h2-HPs[2].H)/h2*100;
	outstrm<<Coeff_com<<"	"<<dev_Ei<<"	"<<dev_h2<<"	";

	HPs[2].P=P->P[2];
	HPs[2].H=h2;

	HotGasLine(HPs[2],P->mr,&HPs[3],&m[1],P);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc1_Tuning","Hotgasline");
		return 6;
	}
	

	HPs[1].H=h1;
	HPs[1].P=P->P[1];

	SuctionLine(HPs[1],P->mr,&HPs[8],&m[2]);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc1_Tuning","suctionline");
		return 7;
	}


	const double h8 = HPs[8].H;

	if(P->X[2]>0) //P->X[1], quality at the condenser exit
	{
	Tsub4 = -1*P->X[2];
	TXP_prop.X= P->X[2];
	TXP_prop.P=P->P[4];
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);
	}
	else 
	{
	TXP_prop.X=0;
	TXP_prop.P=P->P[4];//TXP_prop.P=P->P[5];
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);
	Tsub4=TXP_prop.T-P->T[4];//Tsub4=TXP_prop.T-P->T[5];
	}
	
	const double h4 = PropertyTXPth(ENTH,TXP_prop);

	//B.S.
	Condenser0(Tsub4,P->filename,P->mr,HPs[3],P->Tai,P->Goa,&HPs[4],&Tao,&m[3],Prms);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc1_Components","Condenser");
		return 8;
	}
	outstrm<<Prms[0]<<"	";
	
	while(Condenser(P->filename,P->mr,HPs[3],P->Tai,P->Goa,&HPs[4],&Tao,&m[3],&Cond_struc,Prms)) {	}
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc1_Tuning","Condenser");
		return 8;
	}

	Prms[0]=1;

	LiquidLine(HPs[4],P->mr,&HPs[5],&m[5]);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc1_Tuning","liquidline");
		return 11;
	}

	if(P->X[2]>0) //P->X[2], quality at the exit of the liquid line
	{
	TXP_prop.X=P->X[2];
	TXP_prop.P=P->P[5];
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);
	}
	else 
	{
	TXP_prop.X=0;
	TXP_prop.P=P->P[5];
	TXP_prop.T=P->T[5];
	}

	const double h5 = PropertyTXPth(ENTH,TXP_prop);

	TXP TXP5;
	TXP5 = HPtoTXP(HPs[5]);
	
	TXP_prop=TXP5;
	TXP_prop.X=0;
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

	Tsub5=TXP_prop.T-TXP5.T;

	
	HPs[7].H=h4;//h5;//evaporator inlet enthalp equal to liquid line outlet
	HPs[8].P=HPs[8].P;//P->P[8];
	HPs[8].H=HPs[8].H;//h8;

	double hi=0;
	if(P->REV_CAL)
	{
	hi=h4;//h5
	}
	else{
	hi = HPs[8].H;	
	}

	Evaporator0(0,hi,P->filename,P->mr,HPs[8],P->Gra,P->TPi,&HPs[7],&TPo,&m[4],&Aflow,Prms,P->REV_CAL);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc1_Components","Evaporator");
		return 15;
	}
	
	Evap_struc.REV=P->REV_CAL;
	while(Evaporator(P->filename,P->mr,&HPs[8],P->Gra,P->TPi,&HPs[7],&TPo,&m[4],&Aflow,&Evap_struc,Prms)) { }
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc1_Tuning","Evaporator");
		return 15;
	}

	outstrm<<Prms[0]<<"	";

	Prms[0]=1;

	TXP TXP8;
	TXP8 = HPtoTXP(HPs[8]);
	
	TXP_prop.P=TXP8.P;
	TXP_prop.X=1;

	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);

	double Tsup8=TXP8.T-TXP_prop.T;


	HPs[7].P=P->P[7];
	HPs[7].H=h5;
	DistributionTubes(HPs[7],P->mr,&HPs[6],&m[6]);
	if(errorLog.IsError()) {
		errorLog.Add("MainFunc1_Tuning","Distribution Tubes");
		return 18;
	}
	

	HPs[5].H=h5;
	HPs[5].P=P->P[5];
	HPs[6].H=h5;
	HPs[6].P=P->P[1]+375.09;//P->P[6];//B.S.,

	mr = Short_Tube_R410A(P->filename,HPs[5],HPs[6],Prms);
//	mr = Short_Tube_FlowDirec(P->filename,HPs[5],HPs[6],Ones);
	const double Coeff_Exp=mr/P->mr;
	const double Chg  =  m[0].m+m[1].m+m[2].m+m[3].m+m[4].m+m[5].m+m[6].m;
	const double dev_charge=(P->Charge-Chg);
	outstrm<<mr<<"	"<<Coeff_Exp<<"	"<<Chg<<"	"<<dev_charge<<"	"
	<<m[0].m<<"	"<<m[1].m<<"	"<<m[2].m<<"	"<<m[3].m<<"	"
	<<m[4].m<<"	"<<m[5].m<<"	"<<m[6].m<<"	"<<Cond_struc.LiqL<<std::endl;

	outstrm.close();

	return 0;
}
/*********************************************************************
Print time spent crunching each component
*********************************************************************/
void PrintCycleTimes()
{
	FILE *fp = fopen("times","w");
	if(fp==NULL) return;
	int i;

	time_t t;
	time(&t);
	fprintf(fp,"%s",ctime(&t));

	// Add up all the time spent crunching
	double totalSecs = 0.0;
	for(int i=1;i<NCYCLEPOS;i++) totalSecs += 0.001*(double)cycleTime[i];

	// Print the time crunching each component and the percentage of
	// the total.
	for(int i=1;i<NCYCLEPOS;i++) {
		if(cycleTime[i]==0) continue;
		double secs = 0.001*(double)(cycleTime[i]);
		fprintf(fp,"%2d: %20s %8.3lf seconds, ",i,cyclePositionLabels[i],secs);
		if(fabs(totalSecs)>1e-10)
			fprintf(fp,"%5.1lf%%\n",100.0*secs/totalSecs);
		else
			fprintf(fp,"\n");
	}

	fprintf(fp,"Total seconds = %8.3lf\n",totalSecs);
	fclose(fp);
}

//---------------------------------------------------------------
//---------------------------------------------------------------
//-----------------------------------------------------------------
