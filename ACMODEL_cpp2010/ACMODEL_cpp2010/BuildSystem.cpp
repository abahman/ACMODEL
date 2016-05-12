//these subfuctions are for setting up and running the simplified models. But they are not in use now.

#include "stdafx.h"
#include <math.h>
#include "compu.h"
#include "errorlog.h"
#include "spline1d.h"
#include "expdev.h"//B.S.


extern RefPropTPthP reftplthP;
extern RefPropTPthP reftpvthP;
extern RefPropTPthT reftpvthT;//new

extern ErrorLog errorLog;


struct WholEquip
{
	InVars I;
	double p1;
	double h1;
	double p2;
	double Tsub;
	double Qual;
	int twophase;
	double Tsup;
	double LiqL;
	double p5;
	double h5;
	double m;
	int Lumped;
	double mr;
	TP TPi;
};

double Com_Con(double,void*);
double EvaSide_Cal(double,void*);
double WholEquip_Cal(double,void*);
double InitialGuess(double*, InVars*);


//B.S-----------------------------------------------
double Com_Con(double p2,void*Params)//conduct the compressor and condenser simplified model analysis
{
	WholEquip* ComCond=(WholEquip*)Params;//Equipment parameters

	static double CompPrms[3]={1,1,1};
	InVars P =ComCond->I;//Equipment
	CGP Cond_struc;//condenser struct
	HP HP5;
	HP HP4;
	HP HP3;
	HP HP2;
	HP HP1;
	double Ei =0;
	MASS m ;
	double mr;
	double Tsub;

	CompPrms[0]=P.nvAdj;
	
	Compressor(P.filename,ComCond->h1,ComCond->p1,p2,P.Tai,&mr,&HP2,&HP1,&Ei,&m,CompPrms);//compressor model
	if(errorLog.IsError()) {
		errorLog.Add("ComAndCon","compressor");
		return 0;
	}

	ComCond->m = m.m;

	HotGasLine(HP2,mr,&HP3,&m,&P);//hotgas line model
	if(errorLog.IsError()) {
		errorLog.Add("ComAndCon","HotGasLine");
		return 0;
	}
	
	ComCond->m = ComCond->m+m.m;

	Cond_struc = P.Cond;//condenser struct
	if(ComCond->Lumped==2) Condenser_Simple(mr, HP3,P.Tai,&HP4,&m,&Cond_struc);//moving boundary model analysis
	else Condenser_Lumped(mr, HP3,P.Tai,&HP4,&m,&Cond_struc);//lumped model analysis
	if(errorLog.IsError()) {
		errorLog.Add("ComAndCon","Condenser_Simple");
		return 0;
	}
	ComCond->m = ComCond->m+m.m;
	//state 5-liquid line exit
	TXP TXP4;
	TXP TXP5;
	double Tsat5;

	TXP4= HPtoTXP(HP4);
	if(errorLog.IsError()) {
		errorLog.ClearError("ComAndCon","TXP4");
		HP5.H = HP4.H;
		HP5.P = HP4.P;
		Tsub = 50;
		goto end_Con;}

	LiquidLine(HP4,mr,&HP5,&m);//liquid line model
	if(errorLog.IsError()) {
		errorLog.ClearError("ComAndCon","Liquidline");
		HP5.H = HP4.H;
		HP5.P = HP4.P;
		Tsub = 50;//set a maximum subcooling degree
		goto end_Con;
	}
	ComCond->m = ComCond->m+m.m;

	//state 5-liquid line exit
	TXP5= HPtoTXP(HP5);
	if(errorLog.IsError()) {
		errorLog.Add("ComAndCon","TXP5");
		return 0;
	}

	Tsat5 = reftplthP.Tsat(HP5.P);
	if(errorLog.IsError()) {
		errorLog.Add("ComAndCon","Tsat5");
		return 0;
	}
	Tsub = Tsat5-TXP5.T;

end_Con:;
	ComCond->mr = mr;//mass flow rate
	ComCond->p2 = p2;//condensing pressure
	ComCond->h5 = HP5.H;//condenser exit enthalpy
	ComCond->p5 = HP5.P;//condenser exit pressure
	ComCond->LiqL = Cond_struc.LiqL;//condenser liquid length

	double residual = 0;
	if(ComCond->twophase==0)//subcooled exit
	{
	residual= (Tsub - ComCond->Tsub)/1e3;
	}
	else//two-phaes exit
	{
	residual = (TXP5.X-ComCond->Qual)/1;
	}
	return residual;
}

double EvaSide_Cal(double T_in,void*Params)//conduct evaporator side simplified model analysis
{
	WholEquip* Eva=(WholEquip*)Params;
	ETdim Evap_struc;
	const double mr = Eva->mr;//mass flow rate
	MASS m;
	HP HP8, HP7,HP1;
	InVars P =Eva->I;
	Eva->TPi.T = T_in;//temperature guess
	Eva->TPi.P = P.TPi.P;//keep the original humidity ratio
	HP1.H = Eva->h1;//suction enthalpy
	HP1.P = Eva->p1;//suction pressure
	
	SuctionLine(HP1,mr,&HP8,&m);//suction line
	if(errorLog.IsError()) {
		errorLog.Add("EvaSide_Cal","SuctionLine");
		return 0;
	}
	Eva->m = m.m;

	Evap_struc=P.Evap;//evaporator struct
	if(Eva->Lumped==2)//conduct the moving boundary model analysis 
	Evaporator_Simple(mr/4, HP8,Eva->TPi,&HP7,&m,&Evap_struc);
	else 
	Evaporator_Lumped(mr/4, HP8,Eva->TPi,&HP7,&m,&Evap_struc);//conduct the lumped model analysis

	if(errorLog.IsError()) {
		errorLog.ClearError("EvaSide_Cal","Evaporator_Simple");
		HP7.H=-1e10;
	}
	
	Eva->m = Eva->m+m.m*4;//important for charge calculation

	return (HP7.H - Eva->h5)/1e5;//fabs(HP7.H + Eva->h5);
}

double WholEquip_Cal(double p1,void*Params)//conduct simplified system model analysis
{
	WholEquip* Equip=(WholEquip*)Params;
	InVars P;
	P= Equip->I;//can't use static variable here
	
	Equip->p1 = p1;
	
	const double P_max = 2900;//maximum condensing pressure
	const double P_min =reftpvthT.Psat(P.Tai);//minimum evaporating pressure
	if(errorLog.IsError()) {
	errorLog.Add("WholEquip_Cal","P_min");
	return 0;
	}

	const double Tsat1 = reftplthP.Tsat(p1);//suction pressure
	if(errorLog.IsError()) {
	errorLog.Add("WholEquip_Cal","Tsat1");
	return 0;
	}
	const double Tsup = Equip->I.Tsup;//superheat degree

	TXP TXP1;
	TXP1=toTXP(Tsat1+Tsup,1,p1);//suction state at the compressor entrance
	if(errorLog.IsError()) {
	errorLog.Add("WholEquip_Cal","TXP1");
	return 0;
	}

	Equip->h1= PropertyTXPth(ENTH,TXP1);//suction enthalpy
	if(errorLog.IsError()) {
	errorLog.Add("WholEquip_Cal","h1");
	return 0;
	}
 
	Equip->m =0;
	Zbrent(P_max,P_min,Com_Con,1e-7,Equip);//conduct the condenser side analysis

	const double m_cond = Equip->m;

	if(errorLog.IsError()) {
	errorLog.Add("WholEquip_Cal","Comcond_Zbrent");
	return 0;
	}
	
	
	const double T_max =80;
	const double T_min = Tsat1;
	
	Equip->m =0;
	Zbrent(T_max,T_min,EvaSide_Cal,1e-7,Equip);//conduct the evaporator side analysis
	const double m_eva = Equip->m;

	if(errorLog.IsError()) {
	errorLog.Add("WholEquip_Cal","Eva_Zbrent");
	return 0;
	}

	Equip->m = m_cond + m_eva;
	double residual = (Equip->TPi.T-P.TPi.T)/fabs(Equip->TPi.T+P.TPi.T);
	return residual;
}

double InitialGuess(double *X, InVars* P)
{
	InVars I;
	InVars I1;//Equipment state backup
	double X1[3];//initial guess backup
//	double Y[3];//iteration deviation
	double Ini = 1;//signal variable
	int Num = 0;
	I = *P;

	const double ConvTol = 1e-6;

	errorLog.ClearError ("beginning");

	I.Tsup = P->Tsup;//this is a real superheat degree for runing the model

	
	//keep the results from previous run
	X1[0] = X[0];//store the initial guess
	X1[1] = X[1];//store the initial guess
	I1 = I;//store the intial equipment state
	
	I.Lumped = 1;//use the lumped model to get the initial guess
	I.EvaSign =0;//use the detailed evaporator model
	I.CondSign = 0;//use the detailed condenser model

	while(BuildSystem(X,&I)){	};//guess the initial guess with the lumped model
	
	if(reftpvthP.Tsat(X[0])+I.Tsup>I.TPi.T) 
	{
	errorLog.Add("evaporating pressure guess is too large");
	}

	if(errorLog.IsError())
	{
	printf("---------------Lumped model guess error---------------\n");
	errorLog.ClearError("Lumped model guess");
	X[0] = X1[0];//restore the initial guess
	X[1] = X1[1];//restore the initial guess
	I.Evap = I1.Evap;//restore the evaporator state from the initial equipment state, but the condenser state can be used
	}
	
	if(I1.CE<I.CE) {I=I1; X[0] = X1[0]; X[1] = X1[1];}//restore
	

	return 0;
}


//B.S.---------------------------------------------new
int BuildSystem(double*X1, void *Params)//get the initial guesses with lumped model and moving boundary model
{
	InVars *P;//for storing the equipment parameters
	P=(InVars*)Params;
	int Num=0;//iteration number
	InVars P1;
	
	const double ConvTol =1e-6;
	double Y1[3];


	X1[0]=663;//give an intial heat exchanger state
	X1[1] = 2760;//provide an intial heat exchanger state
	P1=*P;
	P->Tsup = 17.1;
	while(MainFunc2(X1,Y1,P)) { };

	if(errorLog.IsError()) {
		errorLog.Add("BuildSystem","state initialized");
		return 0;
	}
	P->Tsup = P1.Tsup;
	
	WholEquip Equip;//class for the simplified model
	Equip.I = *P;
	Equip.Lumped = P->Lumped;//conduct the lumped model analysis "Lumped = 1" or moving boundary model analysis "Lumped = 2"

	const double P_max = 900;//maximum evaporating pressure
	const double P_min = 300;//minimum evaporating pressure
	
	double Tsubmax = 20;//maximum subcooling degree
	double Tsubmin = 0.1;//0.0000001;//minimum subcooling degree
	double residual = 1;//residual for determing the iteration direction
	
	Equip.twophase = 0;//assuming the subcooling exit

	while(fabs(Tsubmax-Tsubmin)> 1e-7){
	Equip.I = *P;
	Equip.Tsub = (Tsubmax+Tsubmin)/2.0;//subcooling initial guess
	
	P->subset = Equip.Tsub;


	Zbrent(P_max,P_min,WholEquip_Cal,1e-7,&Equip);//conduct the simplified model analysis
	Equip.Lumped=1;
	if(errorLog.IsError()) {
	errorLog.ClearError("BuildSystem","Zbrent");
	Equip.m = 1000;
	Equip.Lumped=0;
	}


	if(reftpvthT.Psat(P->TPi.T-P->Tsup)>Equip.p1) 
	{X1[0] = Equip.p1;}//initial guess
	else 
	{
	double T_diff = 12+(P->Tsup - 0.1)/(25)*(0.1-12);
	X1[0] = reftpvthT.Psat(P->TPi.T-P->Tsup-1);//initial guess for the evaporating pressure
	}
	X1[1] = Equip.p2;//initial guess

	if(Equip.Lumped==1)//conduct the lumped model analysis 
	{
	int signal =1;
	int count=0;
	while(signal)
	{
	while(MainFunc2(X1,Y1,P)) {	};//calculate the system charge
	if(errorLog.IsError()||P->Evap.LiqL>0) {
		errorLog.ClearError("BuildSystem","Calculated Charge");
		X1[0] = X1[0] + 0.5;
	}
	else if(P->Evap.TPL<2.0)
	{X1[0] = X1[0] - 0.5;}
	else signal =0;

	count = count +1;
	if(count>100)  
	{
	signal = 0;
	X1[0] = X1[0] - 0.25;
	}
	}
	Equip.Tsub = P->Tsub;//subcooling degree from the detailed model
	Equip.m = P->ChargMass;//system charge calculated from the detailed model according to the subcooling degree	
	}
	
  	residual= P->Charge- Equip.m;//charge residual
	if(residual>0) Tsubmin = Equip.Tsub;//modify the subcooling guess
	else Tsubmax = Equip.Tsub;//modify the subcooling guess

	if(Tsubmax<0.1) {Tsubmax=0; Tsubmin =0;}

	Num = Num+1;
	if(Num>=40)// limit the maximum iteration number
	{Tsubmax = (Tsubmax+Tsubmin)/2; Tsubmin = Tsubmax;}

	}

	double Qualmax = 0.3;//maximum quality
	double Qualmin = 0.00001;//minimum quality
	Num = 0;//renew the iteration number

	if(Tsubmax<=0.101){//conduct the analysis for the two-phase exit of the condenser
	
	residual =1;
	Equip.twophase = 1;//assuming the two-phase exit of the condenser
	
	while(fabs(Qualmax-Qualmin)> 1e-6){
	Equip.I = *P;
	Equip.Qual = (Qualmax+Qualmin)/2.0;//get the guess for the two-phase exit

	Zbrent(P_max,P_min,WholEquip_Cal,1e-7,&Equip);
	if(errorLog.IsError()) {
		errorLog.Add("BuildSystem","Zbrent");
		return 0;
	}

	if(reftpvthT.Psat(P->TPi.T-P->Tsup)>Equip.p1) 
	{X1[0] = Equip.p1;}//initial guess
	else 
	{
	double T_diff = 12+(P->Tsup - 0.1)/(25)*(0.1-12);
	X1[0] = reftpvthT.Psat(P->TPi.T-P->Tsup-1);//initial guess for the evaporating pressure
	}
	X1[1] = Equip.p2;//initial guess

	if(Equip.Lumped==1)//conduct the analysis for lumped model 
	{
	double signal = 1;
	int count =0;
	while(signal)
	{
	while(MainFunc2(X1,Y1,P)) {	};//calculate the system charge
	if(errorLog.IsError()||P->Evap.LiqL>0) {
		errorLog.ClearError("BuildSystem","Calculated Charge");
		X1[0] = X1[0] + 0.5;
	}
	else if(P->Evap.TPL<2.0)
	{X1[0] = X1[0] - 0.5;}
	else signal =0;

	count = count +1;
	if(count>100)  
	{
	signal = 0;
	X1[0] = X1[0] - 0.25;
	}
	}
	Equip.Qual = P->Qual;
	Equip.m = P->ChargMass;
	}

  	residual= P->Charge- Equip.m;
	if(residual<0) Qualmin = Equip.Qual;
	else Qualmax = Equip.Qual;

	Num = Num+1;
	if(Num>=40)//limit the maximum iteration number 
	{
	Qualmax = (Qualmax+Qualmin)/2; Qualmin = Qualmax;}

	}
	}
	
	return 0;
}
//*************************************************************************************
//*************************************************************************************
//*************************************************************************************