//these subfuctions are for setting up and running the simplified models. But they are not in use now.

#include "stdafx.h"
#include <math.h>
#include "compu.h"
#include "errorlog.h"
#include "spline1d.h"
#include "expdev.h"

extern RefPropTPthP reftplthP;
extern RefPropTPthP reftpvthP;
extern RefPropTPthT reftpvthT;
extern ErrorLog errorLog;


struct WholEquip//equipment parameters for lumped model analysis
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
	double P7;
	double h7;
	double m;
	int Lumped;
	double mr;
	TP TPi;
};

double Com_Con(double,void*);//condenser side simplied model analysis
double EvaSide_Cal(double,void*);//evaporator side simplied model analysis
double WholEquip_Cal(double,void*);//simplied system model analysis
//*************************************************************************************
//*************************************************************************************
//*************************************************************************************
//conduct the compressor and condenser simplified model analysis
//the key for the lumped model analysis is to make it pure mathematical, never have the refrigerant property calculation involved in the calculation
double Com_Con(double p2,void*Params)
{
	WholEquip* ComCond=(WholEquip*)Params;//transfer the input parameters
	static double CompPrms[4]={1,1,1,50};//adjustment coefficients for compressor mass flow rate, powe consumption, and another, CompPrms[3] = 50 C, for compressor lubricant temperature
	InVars P =ComCond->I;//Equipment parameters
	CGP Cond_struc;//condenser struct parameters
	HP HP1={0,0};//compressor suction state
	HP HP2={0,0};//compressor dicharge state
	HP HP3={0,0};//condenser inlet state
	HP HP4={0,0};//condenser outlet state
	HP HP5={0,0};//inlet state before the expansion	
	
	double Ei =0;//compressor power consumption
	MASS m ;//struct parameter for keeping charge and volume inside one component
	double mr;//compressor mass flow rate
	TXP TXP_prop = {0,0,0};//TXP struct for obtaining any properties
	
	CompPrms[0]=P.nvAdj;//adjustment parameter for compressor mass flow rate	
	Compressor(P.filename,ComCond->h1,ComCond->p1,p2,P.Tai,&mr,&HP2,&HP1,&Ei,&m,CompPrms);//compressor model
	if(errorLog.IsError()) {
		errorLog.Add("ComAndCon","compressor");
		return 0;}

	ComCond->m = m.m;//charge inside the compressor

	HotGasLine_Lumped(HP2,mr,&HP3,&m,&P);//lumped hotgas line model, without refrigerant property involved
	if(errorLog.IsError()) {
		errorLog.Add("ComAndCon","HotGasLine_Lumped");
		return 0;}
	
	ComCond->m = ComCond->m+m.m;//add charge in the discharge line

	Cond_struc = P.Cond;//condenser struct which was generated from the detailed model analysis
	
	Condenser_Lumped(mr, HP3,P.Tai,&HP4,&m,&Cond_struc);//lumped model analysis
	if(errorLog.IsError()) {
		errorLog.Add("ComAndCon","Condenser_Lumped");
		return 0;}
	ComCond->m = ComCond->m+m.m;//add charge in the condenser

	TXP_prop.P= HP4.P;//condenser exit pressure from this run
	const double Tsat4 = reftplthP.Tsat(HP4.P);//condenser exit saturation temperature from this run
	if(ComCond->twophase==0)//subcooled exit
	{TXP_prop.T = Tsat4 - ComCond->Tsub; TXP_prop.X=0;}
	else//two phase exit
	{TXP_prop.T = Tsat4; TXP_prop.X = ComCond->Qual;} 

	HP HP_set4 = TXPtoHP(TXP_prop);//obtain the required condenser exit enthalpy by the condensing pressure of this run and the specify condenser exit subcooling degree

	LiquidLine_Lumped(HP4,mr,&HP5,&m);//lumped liquid line model
	if(errorLog.IsError()) {
		errorLog.Add("ComAndCon","LiquidLine");
		return 0;}

	ComCond->m = ComCond->m+m.m;//add charge in the liquid line

	TXP_prop.P= HP5.P;//liquidline exit pressure from this run
	const double Tsat5 = reftplthP.Tsat(HP5.P);//liquidline exit saturation tempertaure in this run
	//here is an alternative for specifying the liquid line exit state
	if(ComCond->twophase==0)//subcooled exit state
	{TXP_prop.T = Tsat5 - ComCond->Tsub; TXP_prop.X=0;}
	else//two-phase exit state
	{TXP_prop.T = Tsat5; TXP_prop.X = ComCond->Qual;} 

	HP HP_set5 = TXPtoHP(TXP_prop);//can be used for specifying the liquid line exit state

	//output parameters
	ComCond->mr = mr;//mass flow rate
	ComCond->p2 = p2;//condensing pressure
	ComCond->h5 = HP5.H;//condenser exit enthalpy for evaporator inlet
	ComCond->p5 = HP5.P;//condenser exit pressure for evaporator inlet
	ComCond->LiqL = Cond_struc.LiqL;//condenser liquid length

	double residual = 1.0;
	if(P.SPEC_SUB == 2)//specify the state before expansion
	{
		residual = (HP5.H-HP_set5.H)/1e6;
	}
	else{//specify the state at the condenser exit
		residual = (HP4.H-HP_set4.H)/1e6;//residual between the specified condenser outlet state and the calculated outlet state
	}
	 return residual;
}

//*************************************************************************************
//*************************************************************************************
//*************************************************************************************
//conduct evaporator side simplified model analysis,
//the key for the lumped model analysis is to make it pure mathematical, never have the refrigerant property calculation involved in the calculation
double EvaSide_Cal(double T_in,void*Params)
{
	WholEquip* Eva=(WholEquip*)Params;//transfering the input parameters
	ETdim Evap_struc;//evaporator state parameters
	MASS m;//struct parameter for keeping charge and volume inside one component
	HP HP8={0,0}, HP7={0,0},HP1={0,0};//evaporator exit, inlet, and compressor suction state
	InVars P =Eva->I;//equipment parameters
	Eva->TPi.T = T_in;//change to the indoor temperature guess
	Eva->TPi.P = P.TPi.P;//keep the original humidity ratio
	HP1.H = Eva->h1;//suction enthalpy
	HP1.P = Eva->p1;//suction pressure
	const double mr = Eva->mr;//mass flow rate obtained from the condenser side analysis
	
	SuctionLine_Lumped(HP1,mr,&HP8,&m);//lumped suction line
	if(errorLog.IsError()) {
		errorLog.Add("EvaSide_Cal","SuctionLine_Lumped");
		return 0;}
	Eva->m = m.m;//charge inside the suction line

	Evap_struc=P.Evap;//evaporator struct generated from the detailed model analysis
	Evaporator_Lumped(mr, HP8,Eva->TPi,&HP7,&m,&Evap_struc);//conduct the lumped model analysis

	if(errorLog.IsError()) {
		errorLog.Add("EvaSide_Cal","Evaporator_Lumped");
		return 0;}

	//output parameters	
	Eva->m = Eva->m+m.m;//add charge inside the evaporator model
	Eva->h7 = HP7.H;//evaporator inlet enthalpy
	Eva->P7 = HP7.P;//evaportor inlet pressure

	return (HP7.H - Eva->h5)/1e5;//residual between the evaporator inlet enthalpies by this run and the value from the condenser side lumped model analysis
}

//*************************************************************************************
//*************************************************************************************
//*************************************************************************************
double WholEquip_Cal(double p1,void*Params)//conduct lumped system model analysis
{
	WholEquip* Equip=(WholEquip*)Params;//transfer the input parameters
	TXP TXP_prop = {0,0,0};//TXP format for getting property
	TXP TXP1 = {0,0,0};//TXP for the compressor suction state
	HP HP1={0,0};//HP for the compressor suction state
	InVars P;//for keeping the parameters inside Equip->I, for it can be changed in the iteration later
	P= Equip->I;//Attention: can't use static variable here
	
	Equip->p1 = p1;//compressor suction initial guess
	
	const double P_max = 2300;//maximum condensing pressure, can be changed
	const double P_min =reftpvthT.Psat(P.Tai);//minimum condensing pressure
	if(errorLog.IsError()) {
	errorLog.Add("WholEquip_Cal","P_min");
	return 0;}

	const double Tsat1 = reftpvthP.Tsat(p1);//compressor suction saturated pressure
	if(errorLog.IsError()) {
	errorLog.Add("WholEquip_Cal","Tsat1");
	return 0;}

	TXP_prop.T = Tsat1; TXP_prop.X=1; TXP_prop.P = p1;//prepare for getting the compressor suction saturation state
	const double Tsup = Equip->I.Tsup;//compressor suction superheat degree

	if(Tsup<0)//minus suction superheat degree introduce the two-phase suction enthalpy like the below equation
	{
	Equip->h1=PropertyTXPth(ENTH,TXP_prop)+PropertyTXPtr(SPEC,TXP_prop)*Tsup;//two-phase suction enthalpy
	HP1=toHP(Equip->h1,p1);
	TXP1=HPtoTXP(HP1);
	if(errorLog.IsError()) {
		errorLog.Add("WholEquip_Cal","h1");
		return 0;}
	}
	else{
	TXP1=toTXP(Tsat1+Tsup,1,p1);
	if(errorLog.IsError()) {
		errorLog.Add("WholEquip_Cal","TXP1");
		return 1;}
	Equip->h1= PropertyTXPth(ENTH,TXP1);//superheated compressor suction enthalpy
	if(errorLog.IsError()) {
	errorLog.Add("WholEquip_Cal","h1");
	return 0;}
	}
	
 
	Equip->m =0;
	//search the condensing pressure by the given suction pressure and superheat degree, then give the compressor mass flow rate, evaporator inlet state 
	Zbrent(P_max,P_min,Com_Con,1e-7,Equip);//conduct the condenser side analysis, to search the specified condenser or liquid line exit state (by subcooling degree or quality)

	const double m_cond = Equip->m;//condenser side charge, not really used for lumped model analysis

	if(errorLog.IsError()) {
	errorLog.Add("WholEquip_Cal","Comcond_Zbrent");
	return 0;
	}
	
	
	const double T_max =50;//maximum indoor temperature. can be adjusted for some occasions
	const double T_min = Tsat1;//minimum indoor temperature
	
	Equip->m =0;
	//search the required indoor temperature by the given suction pressure, superheat degree, and the compressor mass flow rate, evaporator inet state obtained from the condenser side lumped model analysis
	Zbrent(T_max,T_min,EvaSide_Cal,1e-7,Equip);//conduct the evaporator side lumped analysis
	const double m_eva = Equip->m;//evaporator side charge, not really used for lumped model analysis

	if(errorLog.IsError()) {
	errorLog.Add("WholEquip_Cal","Eva_Zbrent");
	return 0;
	}

	Equip->m = m_cond + m_eva;//the whole system charge, not really used for lumped model analysis
	double residual = (Equip->TPi.T-P.TPi.T);//residual between the actual indoor temperature and the required indoor temperature by the given suction pressure, superheat degree and the condenser exit state
	return residual;
}

//*************************************************************************************
//*************************************************************************************
//*************************************************************************************
//interation loop based on simplified model analysis, and updating the simplified models from the detailed model analysis
//for 2-D system having a specified expansion device
int ENGINE::BuildSystem_SubC_2D(double*X1, void *Params)
{
	InVars *P;//equipment parameters
	P=(InVars*)Params;//transfer the equipment parameters
	InVars P1;//for backup
	
	double Y1[3];//residul outputs from the detailed system model analysis
	if(X1[0]<100)  X1[0]=440;//if no proper initial guess, give one
	if(X1[1]<100)  X1[1] = 1500;//if no proper initial guess, give one
	P1=*P;//backup the given equipment parameters

	MainFunc2(X1,Y1,P);//get the initial states of HXs
	if(errorLog.IsError()) {
		errorLog.Add("BuildSystem","state initialized");
		return 0;}

	P->Tsup = P1.Tsup;//recover the given compressor entrance state
	P->Tsub = P1.Tsub;//recover the given condenser exit state
	
	WholEquip Equip;//class for the simplified model
	Equip.Lumped = P->Lumped;//conduct the lumped model analysis "Lumped = 1" or moving boundary model analysis "Lumped = 2", but not exactly used now

	const double P_max = reftpvthT.Psat(P->TPi.T);//maximum evaporating pressure
	const double P_min = 250;//minimum evaporating pressure

	if(P->Tsub>0)//having subcooled condenser or liquid line outlet state
	{Equip.twophase = 0;//judgement for subcooling exit
	Equip.Tsub = P->Tsub;//subcooling degree
	}
	else//having two-phase condenseror liquid line outlet state 
	{Equip.twophase = 1;//two-phase exit
	Equip.Qual = -P->Tsub;//condenser exit quality
	}


	double p1=0,p2=0;//suction and discharge pressure
	int Num = 0;

	do{
		Equip.I = *P;//update the parameters of HXs for simplied model analysis
		p1=X1[0];p2=X1[1];//keep the initial guesses from the last run
		Zbrent(P_max,P_min,WholEquip_Cal,1e-7,&Equip);//conduct the simplified model analysis
		if(errorLog.IsError()) {
		errorLog.Add("BuildSystem_SubC_2D","Zbrent");
		return 0;}
		//initial guess from the simplified model analysis
		X1[0] = Equip.p1;//compressor suction pressure
		X1[1] = Equip.p2;//compressor discharge pressure
		
		//get the new initial guesses for detailed model analysis, using damping factor
		const double R=0.5;//damping factor
		X1[0] = p1+R*(X1[0]-p1);//compressor suction pressure
		X1[1] = p2+R*(X1[1]-p2);//compressor discharge pressure
	
		int ErrLog=0;//error log
		int CountN=0;//number of the error occurances

		do
		{
			ErrLog=0;//no error at the beginning
			MainFunc2(X1,Y1,P);//calculate the system charge, and outputs from detailed models
			if(errorLog.IsError()) {
			errorLog.ClearError("BuildSystem_2D","MainFunc2");
			//have the initial guesses go back a bit, using the damping factor
			X1[0] = p1+R*(X1[0]-p1);
			X1[1] = p2+R*(X1[1]-p2);
			ErrLog =1; CountN=CountN+1;//set the error log and keep the number of error occurances
			}

			if(CountN>5)//if too many error occurances, jump out 
			{errorLog.Add("BuildSystem_2D","MainFunc2"); CountN=0; return 0;}

		}while(ErrLog==1);//keeping doing this if there is an error

		Num = Num+1;
		if(Num>50)//if the iteration times exceed 50, jump out and return
		{errorLog.Add("BuildSystem_2D","Lumped Analysis"); Num=0; return 0;}
		
	}while(fabs(Y1[0])>5e-1||fabs(Y1[1])>1e-1);
	//tolerance for suction pressure less than 0.5 kPa, tolerance for subcooling degree is 0.1 K

	return 0;
}

//*************************************************************************************
//*************************************************************************************
//*************************************************************************************
//interation loop based on simplified model analysis, and updating the simplified models from the detailed model analysis
//for 3-D system having a specified expansion device
int ENGINE::BuildSystem_SubC_3D(double*X1, void *Params)
{
	InVars *P;//for inputting the equipment parameters
	P=(InVars*)Params;//transfering the input parameters
	InVars P1;//for backup the input parameters, especially keeping superheat and subcooling
	double Y1[3];
	
	P1=*P;//backup the input parameters

	if(X1[1]<100)  X1[1]=437;//if no correct given values, using this one
	if(X1[2]<100)  X1[2] = 1542;//if no correct given values, using this one
	
	if(P->Tsup<0) P->Tsup = 5.0;//if no correct given values, using this one
	X1[0]=X1[1];X1[1]=X1[2];//use MainFunc2() to generate the initial heat exchanger effectiveness

	MainFunc2(X1,Y1,P);//generate initial heat exchanger effectiveness
	if(errorLog.IsError()) {
		errorLog.Add("BuildSystem_SubC_3D","Initialize HXs");
		return 0;}

	X1[2]=X1[1];X1[1]=X1[0];//prepare for using MainFunc3()
	P->Tsup = P1.Tsup;//restore the input superheat degree
	P->Tsub = P1.Tsub;//restore the input subcooling degree or condenser quality

	WholEquip Equip;//class for keeping the input parameters to run the simplified model

	const double P_max = reftpvthT.Psat(P->TPi.T);//maximum compressor suction pressure
	const double P_min = 400;//minimum compressor suction pressure
	
	if(P->Tsub>0)//subcooled condenser exit
	{Equip.twophase = 0;//subcooled exit
	Equip.Tsub = P->Tsub;//specify the subcooling degree for lumped model analysis
	}
	else //two-phase condenser exit, P->Tsub is a minus value
	{Equip.twophase = 1;//two-phase exit
	Equip.Qual = -P->Tsub;//P->Tsub gives the minor condenser exit quality
	}

	int Num =0;
	double p1=0,p2=0;
	double dev3=0;
	

	double Tsup_max = 20, Tsup_min = -10;//maximum or minimum superheat degree for the bisetion root search

	do{

		Num=0;

		do{//iteration for searching the roots with specifying both the subcooling degree and superheat degree
			P->Tsup = (Tsup_max+Tsup_min)/2.0;//specify the superheat degree
			P->Tsub = P1.Tsub;//specify the subcooling degree
			Equip.I = *P;//input the system parameters for the lumped model analysis
			p1=X1[1];p2=X1[2];//keep the initial guesses from the last run
			
			Zbrent(P_max,P_min,WholEquip_Cal,1e-7,&Equip);//conduct the lumped model analysis
			if(errorLog.IsError()) {
			errorLog.Add("BuildSystem_SubC_3D","Zbrent");
			return 0;}
			
			//get the new initial guesses from the lumped model analysis
			X1[0] = Equip.h1;//initial guess of suction enthalpy
			X1[1] = Equip.p1;//suction pressure
			X1[2] = Equip.p2;//discharge pressure
	
			const double R=0.5;//damping factor
			X1[1] = p1+R*(X1[1]-p1);//slow the searching step using damping factor
			X1[2] = p2+R*(X1[2]-p2);//try to improve the convergence
			 
	
			int ErrLog=0;//error log for running the detailed model
			int CountN=0;//number for error occurances

			do{
				ErrLog=0;//no error at the beginning
				MainFunc1(X1,Y1,P);//calculate the system charge, and HX effectivenesses from detailed models
				if(errorLog.IsError()) {
				errorLog.ClearError("BuildSystem","MainFunc1");
				//if error, clear the error, make the initial guess coming back a bit and try again
				X1[1] = p1+R*(X1[1]-p1);
				X1[2] = p2+R*(X1[2]-p2);
				ErrLog =1; CountN=CountN+1;//set the error sign and keep the times for error occurances
				}

				if(CountN>5)//if the error occurances exceed five, jump out and return
				{errorLog.Add("BuildSystem","MainFunc1"); CountN=0; return 0;}

				}while(ErrLog==1);//if there is an error, repeat the calculation with corrected initial guesses
			
			Num = Num +1;
			if(Num>50)//if the iteration times exceed 50, jump out and return
			{errorLog.Add("BuildSystem","Lumped Analysis"); Num=0; return 0;}

		}while(fabs(Y1[1])>5||fabs(Y1[2])>2e-1);//tolerance for suction pressure less than 5 kPa, tolerance for subcooling degree is 0.2 K

				
		dev3= Y1[0];//mass flow rate deviation between the compressor and the expansion device correlation
			
		if(dev3<0) Tsup_max = P->Tsup;//if the mass flow rate of the expansion device correlation larger than the compressor mass flow rate, the specified superheat degree is too large
		else Tsup_min = P->Tsup;

	}while(fabs(dev3)>0.5e-3&&fabs(Tsup_max-Tsup_min)>0.01);//do the iteration, while the mass flow deviation is larger than 0.5 g/s and the superheat degree deviation larger than 0.1 K.
	
	return 0;
}
//*************************************************************************************
//*************************************************************************************
//*************************************************************************************
//interation loop based on simplified model analysis, and updating the simplified models from the detailed model analysis
//forsystem simulation by specifying the charge inventory
int ENGINE::BuildSystem_Charg(double*X1, void *Params)
{
	InVars *P;//for storing the equipment parameters
	P=(InVars*)Params;//transfer  the equipment parameters
	double Y1[3] = {0,0,0};
	const double Charg = P->Charge;//tranfer the specified system charge

	P->SPEC_SUB =1;//running in the mode of specifying the condenser or liquid line exit state
	double Tsubmax = 22;//max subcooling degree, minus value means quality, can be changed
	double Tsubmin = 0;//min subcooling degree, minus means quality, can be changed
	double dev1=0,dev2=0;//deviation for subcooling degree and charge inventory
	
	do
	{
		P->Tsub = (Tsubmax+Tsubmin)/2.0;//specify the subcooling degree in bi-section method
			
		switch(systemType) {
		case 1:
			BuildSystem_SubC_3D(X1, P);//solve the system by specifying subcooling degree and superheat degree
			if(errorLog.IsError()) {
			errorLog.Add("BuildSystem_Charg","BuildSystem_2D_SubC");
			return 0;}
		
			MainFunc3(X1,Y1,P);//get the system charge
			if(errorLog.IsError()) {
			errorLog.Add("BuildSystem_Charg","MainFunc2");
			return 0;}			
		break;
		
		case 2:
			BuildSystem_SubC_2D(X1, P);//solve the system by specifying subcooling degree and superheat degree
			if(errorLog.IsError()) {
			errorLog.Add("BuildSystem_Charg","BuildSystem_2D_SubC");
			return 0;}
		
			MainFunc2(X1,Y1,P);//get the system charge
			if(errorLog.IsError()) {
			errorLog.Add("BuildSystem_Charg","MainFunc2");
			return 0;}
		break;

		default:
			errorLog.Add("BuildSystem_Charg","Mode not found");
		break;
		}
		
		if(P->ChargMass>Charg) Tsubmax = P->Tsub;
		else Tsubmin = P->Tsub;

		dev2 =(P->ChargMass-Charg);//deviation for charge inventory
		dev1 = Tsubmax-Tsubmin;//deviation for subcooling degree
		if(P->Tsub<0) dev1 = dev1*100;//if two-phase exit, need a more strict evaluation

	}while(fabs(dev1)>0.1&&fabs(dev2)>0.02);
	//residule for charge mass is 0.02 kg, for subcooling degree is 0.l K

	P->SPEC_SUB =0;//recover the judgement digit for specifying system charge
	return 0;
}

//*************************************************************************************
//*************************************************************************************
//*************************************************************************************
//*****************************