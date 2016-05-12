#include "stdafx.h"

#include <stdio.h>
#include <stdio.h>
#include <math.h>
#include "mass.h"
#include "wair.h"
#include "air.h"
#include "r22.h"
#include "evap.h"
#include "pressure.h"
#include "volume.h"
#include "corr.h"
#include "numeric.h"
#include "num2d.h"
#include "getdata.h"
#include "cmine.h"
#include "ref.h"
#include "errorlog.h"
#include "spline2d.h"
#include "spline1d.h"
#include "outelem.h"//B.S.

extern OutputDataElement gOutData;//B.S.


extern ErrorLog errorLog;

extern WetAirPropTR wair;
extern AirProp air;

//#define EXTRA_REPORTING

#define EVAPORATOR_PRESSURE_DROP

int evapInit=1;


struct Evap0P {
	int i;
	int REV;
	double Ga,mr,Hi,*Prms;
	HP HPo;
	HP HPi;//B.S.
	TP TPi;
	const char *filename;
};
struct Evap1P {
	double mr,Hin,Pout,Ga,*Prms;
	TP TPi;
	const char *filename;
};
struct ELRP {
	HP HPo;
	double q,K,G,D,vo;
	double hRefAdj,Api,Ro;
	double ma,mr;
	TP TPi;
};

//B.S.-------------------------------------
//add for evaporative heat transfer iteration, since most of the evporating heat transfer calculation contains the nucleate boiling component, which needs to know the tube wall temperature at first
struct EVA_Get_Q
{//shenbo add
	double ma;//air mass flow rate
	double mr;//refrigerant mass flow rate
	TXP TXPo;//refrigerant state in the segment
	TP TPi;//air inlet state of the segment
	double Gr;//mass flow flux
	double q;//heat flux
	double T_w;//inside tube wall temperature
	ETdim P;//evaporator struct
	double W;//outlet air humidity 
	double Cmin;//Cmine
	double hi;//inside refrigerant heat transfer coefficent
	HP HPo;//inlet state of the refrigerant
	double T_S_O;//effective tube surface temperature
	};
//------------------------------------------------B.S.

//B.S.------------------------------------------
struct Phase_eva
{
	ETdim P;
	HP HPo;
	double mr;
	TP TPi;
	HP HPi;
	double L;
	double m;
};

double GasCal_Eva(double, void *);//B.S., new
double TPCal_Eva(double, void *);//B.S., new
double LiqCal_Eva(double, void *);//B.S., new, empty now

//-------------------------------------------B.S.


// Local functions
#ifdef EVAPORATOR_PRESSURE_DROP
double Get_Q_EVA(double,void*);//B.S., new
double Get_Q_Single(double, void*);//B.S., new
double Get_Q_Single_For(HP, void*);//B.S., new
#else
static int EvapTubeLnoDP(double,HP*,double,TP,WH*,MASS*,ETdim*);
#endif

//not used
/****#ifdef EVAPORATOR_PRESSURE_DROP
static int EvapTubeL(double,HP*,double,TP,WH*,MASS*,ETdim*);
static double CompLossResid(double*,double*,void*);
static int HPLimitsConst(double*,double*,void*);
#endif*********/

static double Evap0Func(double,void*);
static double Evap1Func(double,void*);

/********************************************************************
Model of a direct expansion evaporator  coil. The model is comprised 
of small pieces of single finned tubes put together to form rows.  
Rows of finned tubes are put together using elbows and manifolds to 
make an entire coil. This function defines the coil's circuit pattern.
Model parameters are stored in the file "evap.dat".

INPUTS
mr - refrigerant mass flow rate (kg/s).
HPo - outlet refrigerant thermodynamic state (J/kg,kPa). 
Ga - air mass flux (kg/s*m^2).
TPi - inlet air thermodynamic state (C,NU).
Prms - array of three multiplicative correlation correction factors.
	Prms[0] - corrects air side convection coefficient (dry coil).
	Prms[1] - corrects refrigerant side convection coefficient.
	Prms[2] - corrects refrigerant side pressure drop correlation.
	Prms[3] - corrects air side convection coefficient (wet coil).
	Prms[4] - reduction in air flow rate (0-1) to simulate fouling.
	Set these parameters to 1 (Prms={1,1,1}) for no corrections.

OUTPUTS
HPi - inlet refrigerant thermodynamic state (J/kg,kPa).
TPo - outlet air thermodynamic state (C,NU).
Sm - mass of charge in condenser (kg) and volume of charge (m^3).
Aflow = air flow area (ma=Ga*Aflow) (m^2).
status - error flag.  Set to 1 to indicate error has occurred.
********************************************************************/
int Evaporator(const char *filename,
			   double mr,//refrigerant mass flow rate
			   HP *HPo,//refrigerant outlet state
			   double GaI,//air mass flux
			   TP TPi,//inlet air state
			   HP *HPi,//inlet refrigerant state of the segment
			   TP *TPo,//outlet air state
			   MASS *Sm,//mass and volume in the evaporator
			   double *Aflow,//air flow area
			   ETdim *Evap_struc,
			   double *Prms,
			   int NSeg)//adjustment parameters
{
	static ETdim D;
	static double hAirAdj,hRefAdj,PRefAdj,WAirAdj;
	static double FoulFac,GaNom;
	static double pi = acos(-1.0);
	static int endSign;

	if(evapInit || NSeg>0) {
		FILE *fp = fopen(filename,"r");
		if(!(D.type=PosComponentData(fp,EVAPORATOR,1))) {
			errorLog.Add("Evaporator","Parameters not found");
			return 0;
		}
		const double PI=4e0*atan(1e0);

		D.Di=GetDouble(fp);             // tube inside diameter
		D.L=GetDouble(fp);              // tube length 
		//printf("evap tube length = %lf\n",D.L);
		D.xp=GetDouble(fp);             // tube thickness
		D.Df=GetDouble(fp);
		//get by equal fin surface area
		D.z=GetDouble(fp);//B.S., it is used as the fin pitch
		D.th=GetDouble(fp);
		D.vsp=GetDouble(fp);//B.S, actually it is P_t
		D.P_l = GetDouble(fp);//B.S. spacing between tubs in the longitudual direction (m)
		D.NSeg=(int)GetDouble(fp);
		if(NSeg>0) D.NSeg = NSeg;
		D.Brad=GetDouble(fp);
		GaNom=GetDouble(fp);

		D.NBranchs = (int)GetDouble(fp); // number of equivalent branches
		D.NBraTube = (int)GetDouble(fp); // number of tube in each branch
		D.Nrows=(int)GetDouble(fp);		// number of rows high
		D.Ndeep=(int)GetDouble(fp);		// number of rows deep
		D.Frontal_A = GetDouble(fp);	// frontal area, B.S.

		hAirAdj=GetDouble(fp);
		hRefAdj=GetDouble(fp);
		PRefAdj=GetDouble(fp);
		WAirAdj=GetDouble(fp);
		FoulFac=GetDouble(fp);

		//B.S. ------------------------
		//new input parameter for condenser
		D.microfin=(int)GetDouble(fp); //microfin type, 0=smooth tube, 1=helical, 2=cross-grooved, 3=herringbone
		D.finN = GetDouble(fp); //fin number in a micro-fin tube
		D.gama = GetDouble(fp);//fin apex angle in a micro-fin tube
		D.beta = GetDouble(fp);	//fin helix angle in a micro-fin tube
		D.finH = GetDouble(fp); //fin height in a micro-fin tube
		D.w_b= GetDouble(fp); //base width of a single fin
		D.w_e= GetDouble(fp); //top width of a single fin
		D.w_z= GetDouble(fp); //base distance between two neighboring fins
		// airside fin type
		D.airFin = (int)GetDouble(fp);//1-plain, 2-corrugated, 3-slit, 4-louvered, 5-convex louvered, 6-smooth wavy, 7-spine
		D.sub1= GetDouble(fp);//fin substructures, for the possible use later, such as the waffle height of the wavy fin, louver number of the louverd fin
		D.sub2= GetDouble(fp);//fin substructures
		D.sub3= GetDouble(fp);//fin substructures
		D.sub4= GetDouble(fp);//fin substructures
		D.sub5= GetDouble(fp);//fin substructures
		D.K_T = GetDouble(fp);//400
		//this is the conductance factor of copper, for calculating tube resistance
		D.K_F= GetDouble(fp);//237;
		// conductance factor of fin
		//------------------------------B.S.
	
		endSign=(int)GetDouble(fp);//-256
		fclose(fp);
		//close input data file

		if(endSign!=-256)
		{
		errorLog.Add("Evaporator","parameter input wrong");
		return -1;
		}

		evapInit = 0;
		D.Ls=D.L/(double)D.NSeg;
		
		//B.S.-----------------------------------
		// tube inside geometry
		if(D.microfin)
		{
		D.D_b = D.Di + 2*D.finH;//tube diameter at the base of the fin
		D.Do=D.D_b+2*D.xp; //Pipe outside diameter.
		D.D_m = (D.Di + D.D_b)/2; //mean diameter of the micro-fin tube
		D.P_H = D.finN*(D.w_z+2*pow((pow(D.finH,2.0)+pow((0.5*(D.w_b-D.w_e)),2.0)),0.5)+D.w_e);// the hydraulical circumference
		D.Acs = PI*pow((D.D_b/2.0),2.0)-D.finN*0.5*D.finH*(D.w_b+D.w_e);//cross area of the micro-fin tube, this is the actual cross-section area
		D.Dh_i = 4*D.Acs/D.P_H;//inside hydraulical diameter 
		D.Ax = PI/4.0*D.Di*D.Di;// Inside pipe cross sectional area, based on fin tips
		D.Api=PI*D.Di*D.Ls ;// Inside pipe surface area (one tube segment), based on fin tips
		}
		else 
		{
		D.finH=0; //without fins
		D.D_b = D.Di;//tube diameter at the base of the fin
		D.Do=D.D_b+2*D.xp; //Pipe outside diameter.
		D.D_m = D.Di ; //averge diameter of the micro-fin tube
		D.Ax = PI/4.0*D.Di*D.Di;// Inside pipe cross sectional area
		D.Api=PI*D.Di*D.Ls ;// Inside pipe surface area (one tube segment)
		D.P_H = PI*D.D_b;// the hydraulical circumference
		D.Acs = D.Ax;//cross area of the micro-fin tube, this is the actual cross-section area
		D.Dh_i = D.Di;//inside hydraulical diameter
		}
		//-----------------------------------------B.S.
		D.Vs=D.Acs*D.Ls;
		D.Blen=pi*D.Brad;
		D.BVs=D.Blen*D.Ax;

		/*External Tube Parameters*/
		//B.S.-----------------------------
		//get the fin reference length, according to Schimidt equation
		const double M=D.vsp/2e0;
		const double R_O=D.Do/2;
		const double L=0.5e0*pow((pow(M,2e0)+pow(D.P_l,2e0)),0.5);
		const double BETA=L/M;
		const double PHI=M/R_O;
		const double R_EFF=R_O*1.27e0*PHI*pow((BETA-0.3e0),0.5e0);
		D.Df = 2.0*R_EFF;//confirm the outside reference fin diameter
		D.L_F=(R_EFF-R_O)*(1e0+0.35e0*log(R_EFF/R_O));//get the fin reference length for calculating the fin efficiency, according to Schimidt equation
		//--------------------------------------B.S.

		D.y=(D.Df-D.Do)/2;//Distance from outside of pipe to edge of fin. 
		D.N=(int)(D.Ls/(D.z));//B.S., change the original, (int)(D.Ls/(D.z+D.th));	
		D.Apo=pi*D.Do*(D.Ls-D.N*D.th);//B.S., correct pi*D.Do*(D.Ls);

		const double x = D.Df + D.th/2.0;
		D.Af = D.N*pi/2.0*(x*x-D.Do*D.Do);//fin has double faces
		D.Aflow=(D.vsp-D.Do)*(D.Ls-D.N*D.th);//minimum airflow area per segment
		// printf("Minimum free flow area = %lf (m^2)\n",D.Aflow);

		D.Dh=4*D.Aflow/(D.Af+D.Apo)*D.vsp;     /* hydraulic diameter */
	
	}

	D.hAirAdj = Prms[0]*hAirAdj;
	D.hRefAdj = Prms[1]*hRefAdj;
	D.PRefAdj = Prms[2]*PRefAdj;
	D.WAirAdj = Prms[3]*WAirAdj;

	const double Ga = GaI<0 ? FoulFac*GaNom*Prms[4] : GaI;
	
	// Calculate air cfm after Ga is established.
	const double A_total = D.Aflow*D.NSeg*(double)D.Nrows; //total area air flows throguh for entire coil (m^2)
	const double v_air = air.v(TPi.T); // specific volume of air
	const double vfra = Ga*A_total*v_air; // nominal volumetric flow rate air
	D.cfma = vfra*(60.0/((12.0*0.0254)*(12.0*0.0254)*(12.0*0.0254))); // nominal cfm of air flow

	//B.S., prepare to set up the simple model
	D.V_TOT =0; D.V_TP=0; D.V_Liq =0; D.V_Vap = 0;
	D.m_TOT = 0; D.m_TP =0; D.m_Liq =0; D.m_Vap =0;
	D.mr = 0; D.ma_TOT =0;
	D.UA_TOT=0; D.UA_Liq =0;D.UA_Vap = 0; D.UA_TP =0;
	D.UAw_TOT=0; D.UAw_Liq =0;D.UAw_Vap = 0; D.UAw_TP =0;
	D.LiqL = 0; D.TPL = 0; D.VapL = 0;
	D.L_dry = 0; D.L_wet =0;
	D.count1=0;D.count2=0;
	D.HP_TP1.H=0; D.HP_TP1.P=0;
	D.HP_TP2.H=0; D.HP_TP2.P=0;
	D.Qtp_dry =0; D.Qtp_wet = 0;
	//-----------------------------------------B.S.
	
	D.GetP =0;//don't calculate the airside pressure drop in the following procedure
	
	D.REV= Evap_struc->REV;//iteration loop direction
	EvapCircuit(D.type,mr,HPo,Ga,TPi,HPi,TPo,Sm,Aflow,&D);
	if(errorLog.IsError()) {
		errorLog.Add("Evaporator","EvapCircuit");
		return 0;
	}

	const double DP_circuit = Circuit_DP(Ga,TPi,*TPo,&D);//B.S., calculate the airside pressure drop across the circuit
	gOutData.Eva_AirDP = DP_circuit;//for output

	if(errorLog.IsError()) {//the error in calculating the airside pressure drop should not influence the main iteration loop
		errorLog.ClearError("Circuit_DP","Evap");
	}


	//B.S.---------------------------------------
	D.HP_out = *HPo;//suciton site of the evaporator
	D.HP_in = *HPi;//entrance of the evaporator
	D.mr = mr; D.Ga = Ga; D.TPi=TPi; D.ma_TOT = Ga*D.Aflow*D.NSeg*(D.Nrows);
	Build_Eva_Lumped(&D);//correlate lumped evaporator model
	*Evap_struc = D;//output the evaporator construct
	//---------------------------------------B.S.

	return 0;
}

/********************************************************************
Solves the evaporator model iteratively.  It takes input enthalpy
as an arguement and solves for the appropriate value of Prms[i].
********************************************************************/
void Evaporator0(int i,double Hi,const char *filename,double mr,HP HPo,double Ga,TP TPi, HP *HPi,TP *TPo,MASS *Sm,double *Aflow,double *Prms,int REV)
{
	Evap0P P;
	P.filename = filename;
	P.i=i;
	P.Ga=Ga;
	P.mr=mr;
	P.HPo=HPo;
	P.HPi = *HPi;
	P.TPi=TPi;
	P.Prms=Prms;
	P.REV=REV;

	P.Hi=Hi;
	
	Prms[i] = FindZero(Prms[i],Evap0Func,1e-4,&P);//1e-7

	if(errorLog.IsError())
	{
	errorLog.ClearError("Evaporator0","toBisection1DSpline");
	const double Xlo=0.3; const double Xhi=5.0;
	Prms[i]=Bisection1DSpline(Xlo,Xhi,Evap0Func,1e-7,&P);
	}

	if(errorLog.IsError()) {
		errorLog.Add("Evaporator0");
		return;
	}

	ETdim Evap_struc;
	Evap_struc.REV=REV;
	Evaporator(filename,mr,&HPo,Ga,TPi,HPi,TPo,Sm,Aflow,&Evap_struc,Prms);
	if(errorLog.IsError()) {
		errorLog.Add("Evaporator0");
		return;
	}
}

/********************************************************************
Function used by Evaporator0() during its iteration.
********************************************************************/
static double Evap0Func(double Adj,void *Params)
{
	Evap0P *P = (Evap0P*)Params;
	P->Prms[P->i] = Adj;

	HP HPi;

	HPi=P->HPi;//B.S.
	TP TPo;
	MASS Sm;
	double Aflow;
	ETdim Evap_struc;
	Evap_struc.REV=P->REV;
	Evaporator(P->filename,P->mr,&P->HPo,P->Ga,P->TPi,&HPi,&TPo,&Sm,&Aflow,&Evap_struc,P->Prms);
	if(errorLog.IsError()) {
		errorLog.Add("Evap0Func");
		return 0;
	}
	
	if(P->REV)
	{
	return (HPi.H-P->Hi)/P->Hi;
	}
	else
	{
	return (P->HPo.H-P->Hi)/P->Hi;
	}
}


/********************************************************************
Models an evaporator tube return bend.
********************************************************************/
void EvapTubeBend(double Gr,HP *HPo,MASS *m,ETdim *P)
{
	TXP TXPo = HPtoTXP(*HPo);
	if(errorLog.IsError()) {
		errorLog.Add("EvapTubeBend");
		return;
	}

	double vo = VolumeALL(TXPo,Gr,P->Di,0);
	if(errorLog.IsError()) {
		errorLog.Add("EvapTubeBend");
		return;
	}

	double dP = dPelbow(TXPo,Gr,P->Di,P->Brad);
	if(errorLog.IsError()) {
		errorLog.Add("EvapTubeBend");
		return;
	}

	HP HPm;
	HPm.H=HPo->H;
	HPm.P=HPo->P-dP;
	TXP TXPm=HPtoTXP(HPm);
	if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"dP=%lf HPm.H=%lf HPm.P=%lf\n",dP,HPm.H,HPm.P);
		errorLog.Add("EvapTubeBend");
		return;
	}

	double vm = VolumeALL(TXPm,Gr,P->Di,0);
	if(errorLog.IsError()) {
		errorLog.Add("EvapTubeBend");
		return;
	}
	dP = dPmom(vo,vm,Gr);
	if(errorLog.IsError()) {
		errorLog.Add("EvapTubeBend");
		return;
	}

	HPo->H=HPm.H;
	HPo->P=HPm.P-dP;
	if(errorLog.IsError()) {
		errorLog.Add("EvapTubeBend");
		return;
	}

	m->V=P->BVs; m->m=P->BVs/vm;

	
	//B.S.----------------------------------------------------
	if(TXPm.X==1)//B.S., the following caculate the inner volume and mass in the tube bend	
	{
	P->m_Vap = P->m_Vap+m->V/vm;//B.S. record the vapor mass
	P->V_Vap = P->V_Vap + m->V;//B.S. record the vapor volume
	}
	else if(TXPm.X<1&&TXPm.X>0.00)
	{
	P->m_TP = P->m_TP+m->V/vm;//B.S. record the two-phase mass
	P->V_TP = P->V_TP + m->V;//B.S. record the two-phase volume
	}
	else 
	{
	P->m_Liq = P->m_Liq+m->V/vm;//B.S. record the liquid mass
	P->V_Liq = P->V_Liq + m->V;//B.S. record the liquid volume
	}
	//---------------------------------------------B.S.
}

void EvapTubeBend_Fwd(double Gr,HP *HPo,MASS *m,ETdim *P)
{
	TXP TXPo = HPtoTXP(*HPo);
	if(errorLog.IsError()) {
		errorLog.Add("EvapTubeBend");
		return;}

	double vo = VolumeALL(TXPo,Gr,P->Di,0);
	if(errorLog.IsError()) {
		errorLog.Add("EvapTubeBend");
		return;}

	double dP = dPelbow(TXPo,Gr,P->Di,P->Brad);
	if(errorLog.IsError()) {
		errorLog.Add("EvapTubeBend");
		return;}

	HP HPm;
	HPm.H=HPo->H;
	HPm.P=HPo->P+dP;
	TXP TXPm=HPtoTXP(HPm);
	if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"dP=%lf HPm.H=%lf HPm.P=%lf\n",dP,HPm.H,HPm.P);
		errorLog.Add("EvapTubeBend");
		return;}

	double vm = VolumeALL(TXPm,Gr,P->Di,0);
	if(errorLog.IsError()) {
		errorLog.Add("EvapTubeBend");
		return;}

	dP = dPmom(vo,vm,Gr);
	if(errorLog.IsError()) {
		errorLog.Add("EvapTubeBend");
		return;}

	HPo->H=HPm.H;
	HPo->P=HPm.P-dP;
	if(errorLog.IsError()) {
		errorLog.Add("EvapTubeBend");
		return;}

	m->V=P->BVs; m->m=P->BVs/vm;

	
	//B.S.----------------------------------------------------
	if(TXPm.X==1)//B.S., the following caculate the inner volume and mass in the tube bend	
	{
	P->m_Vap = P->m_Vap+m->V/vm;//B.S. record the vapor mass
	P->V_Vap = P->V_Vap + m->V;//B.S. record the vapor volume
	}
	else if(TXPm.X<1&&TXPm.X>0.00)
	{
	P->m_TP = P->m_TP+m->V/vm;//B.S. record the two-phase mass
	P->V_TP = P->V_TP + m->V;//B.S. record the two-phase volume
	}
	else 
	{
	P->m_Liq = P->m_Liq+m->V/vm;//B.S. record the liquid mass
	P->V_Liq = P->V_Liq + m->V;//B.S. record the liquid volume
	}
	//---------------------------------------------B.S.

}

#ifdef EVAPORATOR_PRESSURE_DROP
/*********************************************************************
Evaporator tube segment model (reversed in flow direction) which considers pressure drops
Inputs:
	Gr = refrigerant mass flux (kg/s/m^2)
	HPi = refrigerant inlet state (h,P)
	Ga = air mass flux (kg/s/m^2)
	Tai = air inlet temperature (C)
Outputs:
	HPi = refrigerant outlet state (h,P)
	hao = air outlet enthalpy (J/kg)
	m = mass of charge in return bend (kg)
*********************************************************************/
int EvapTubeL_Rev(double Gr,//refrigerant mass flux
					HP *HPo,//refrigerant outlet and inlet state
					double Ga,//air mass flux
					TP TPi,//air inlet state
					WH *WHo,//air outlet state
					MASS *m,//charge and inner volume in the evaporator
					ETdim *P)//evaporator struct
{
	double v,hai,ma,mr;
	TXP TXPo;
	EVA_Get_Q EVA_Get_Q;//B.S., this struct is for storing some parameters for iteration
	double y=0;
	double hi=0;
	double Ri=0,R=0;
	double q=0;

	/* Calculate air side resistance */
	P->wet=0;
	const double ho = P->hAirAdj*ConvCoeffAir_EVA(TPi,Ga,P);//B.S., the airside dry heat transfer coefficient can be got for different fin types
	if(errorLog.IsError()) {
		errorLog.Add("EvapTubeL", "get ho");
		return 0;
	}

	const double phi = FinEffect_Schmidt(ho,233,P->th,P->y,P->Do);//B.S. calculate the fin efficiency with the Schmidt equation
	if(errorLog.IsError()) {
		errorLog.Add("EvapTubeL", "FinEffect_Schmidt");
		return 0;
	}

	P->Ro = 1/(ho*(P->Apo+phi*P->Af));
	P->ho=ho;
	P->Ga=Ga;

	
	// convert inlet state into TXP format
	TXPo = HPtoTXP(*HPo);
	if(errorLog.IsError()) {
		errorLog.Add("EvapTubeL_new", "(TXPo,1)");
		return 0;
	}

	//B.S.--------------------------------------
	TXP TXP_bak;
	TXP_bak = TXPo;//backup the outlet state of this segment
	//---------------------------------------B.S.

	// Mass flow rates of air and refrigerant
	ma=Ga*P->Aflow;//B.S. Ga is the maximum airflow flux, P->Aflow is the minimun air flow cross-sectional area
	mr=Gr*P->Ax;

	// Calculate heat transfered per unit mass of refrigerant.

	const double hi_max=10e4;//largest possible refrigerant side heat tranfer coeffcient
	const double hi_min=100;//minimum possible refrigerant side heat transfer coeffcient

	//B.S.-----------------------------
	//store the parameters for iteration
	EVA_Get_Q.ma=ma;
	EVA_Get_Q.mr=mr;
	EVA_Get_Q.TPi=TPi;
	EVA_Get_Q.Gr=Gr;
	EVA_Get_Q.P=*P;
	EVA_Get_Q.q=0;
	EVA_Get_Q.W=wair.HumidityRatio(TPi.T,TPi.P);
	if(errorLog.IsError()) {
		errorLog.Add("EvapTubeL_new", "(W)");
		return 0;
	}
	//---------------------------------------B.S.

	const double X1=0.05,X2=0.95;//	const double X1=0.01,X2=0.99;
	//for calculating the heat transfer and pressure drop in single-phase or approximated single-phase region

	if(TXPo.X>X2||TXPo.X<X1) {

		if(TXPo.X>=1 || TXPo.X<=0) {//Pure single phase region

			y = ConvCoeffSP(TXPo,Gr,P);//get the single-phase heat transfer coefficient
			if(errorLog.IsError()) {
				errorLog.Add("EvapTubeLNew","(hi_SP)");
				return 0;
			}
		} else if(TXPo.X<X1) {
	
			TXP TXP1 = toTXP(TXPo.T,0.0,TXPo.P);
			TXP TXP2 = toTXP(TXPo.T,X1,TXPo.P);

			double y1 = ConvCoeffSP(TXP1,Gr,P);
			if(errorLog.IsError()) {
				errorLog.Add("EvapTubeLNew","(hi_SP)");
				return 0;
			}
		
			EVA_Get_Q.TXPo=TXP2;//local refrigerant state
			hi = Zbrent(hi_max,hi_min,Get_Q_EVA,1e-5,&EVA_Get_Q);//B.S., get the heat transfer coefficient at the two-phase region
			if(errorLog.IsError()) {
				errorLog.Add("EvapTubeLNew","Zbrent1");
				return 0;
			}

			double y2 = hi/P->hRefAdj;//B.S., calculate the two-phase heat transfer coefficient at this region
			
			if(y1>y2) y1=y2;//B.S., make the single phase flow less than two phase
			y = y1+TXPo.X*(y2-y1)/X1;//B.S., get the heat transfer coefficient in this region with intropolation
		} else {

			TXP TXP1 = toTXP(TXPo.T,X2,TXPo.P);
			TXP TXP2 = toTXP(TXPo.T,1,TXPo.P);
	
			EVA_Get_Q.TXPo = TXP1;
			hi = Zbrent(hi_max,hi_min,Get_Q_EVA,1e-5,&EVA_Get_Q);//two-phase refrigerant heat transfer coefficient
			//B.S., get the heat transfer coefficient at the two-phase region
			if(errorLog.IsError()) {
				char msg[256];
				sprintf(msg,"Zbrent2: TXPo.T=%.2lf, T_max=%.2lf, T_min=%.2lf",TXPo.T,hi_max,hi_min);
				errorLog.Add("EvapTubeLNew",msg);
				return 0;
			}

			double y1 = hi/P->hRefAdj;//ConvCoeffEvapTP_microfin(TXP1,Gr,P);//B.S. get two-phase heat transfer coefficient at this region	
		
			double y2 = ConvCoeffSP(TXP2,Gr,P);//B.S., single-phase heat transfer coefficient
			if(errorLog.IsError()) {
				errorLog.Add("EvapTubeLNew","(hi_SP)");
				return 0;
			}

			if(y2>y1) y2=y1; //B.S., for making the single phase flow less than two phase
			y = y2-(1-TXPo.X)*(y2-y1)/(1-X2);//B.S., get the heat transfer coefficient in this region with intropolation
		} 

		hi = P->hRefAdj*y;//B.S., heat transfer coefficient at this region
		Ri = 1/(hi*P->Api);//B.S., inside thermal resistance
		const double R_W=log(P->Do/(P->Do-2.0*P->xp))/(2.0*3.1415*P->K_T*P->Ls);
		R = P->Ro+R_W+Ri;//overall thermal resistance
	
		//B.S., prepare the parameters for iteration
		EVA_Get_Q.TXPo = TXPo;
		EVA_Get_Q.Cmin = CmineCrossFlow_dry(R,mr,ma,TXPo,TPi.T);
		EVA_Get_Q.HPo = *HPo;
		
		if(errorLog.IsError()) {
				errorLog.Add("EvapTubeL_new", "T_sat");
				return 0;
			}

		if(0/*TXPo.T<T_sat+1*/) {
			// TXPo.T is actually the outlet refrigerant temperature
			// however TXPo.T = TXPi.T in the two phase region and
			// TXPo.T ~ to TXPi.T in the superheated region if steps are small

			q = (TPi.T-TXPo.T)/mr*CmineCrossFlow_dry(R,mr,ma,TXPo,TPi.T);

			if(errorLog.IsError()) {
				errorLog.Add("EvapTubeL_new", "(q_singlephase)");
				return 0;
			}	

		} else {//when it is close to the saturated state, the function zbrent can not converge

			const double H_max=HPo->H;
			const double H_min=HPo->H-4e4;
			EVA_Get_Q.hi=hi;//input the refrigerant side heat transfer coefficient
			Zbrent(H_max,H_min,Get_Q_Single,1e-7,&EVA_Get_Q);//B.S., to calculate the single-phase heat transfer and pressure drop at this region
			q=EVA_Get_Q.q;
		}

		if(errorLog.IsError()) {
			errorLog.ClearError("EvapTubeL_new", "(q_singlephase)");
			q = (TPi.T-TXPo.T)/mr*CmineCrossFlow_dry(R,mr,ma,TXPo,TPi.T);
		}/**///don't delete this

		//when the superheat is setted too high, there might exist a problem
	/*	if(q<=0) {//B.S., to remove the wrong result
			q=0;
		}*///don't delete this

		if(errorLog.IsError()) {
			errorLog.Add("EvapTubeL_new", "(q_singlephase)");
			return 0;
		}	

	
		//B.S.---------------------------------------------
	
		if(TXPo.X>=0.9999999) {//B.S., the following caculate the whole thermal resistances separately, prepared to adjust them separately.	
			R = (P->Ro+Ri);
			P->UA_Vap=P->UA_Vap+1/R;//B.S., vapor heat transfer conductance
			P->VapL = P->VapL+ P->Ls;//B.S., vapor length
		} else if(TXPo.X<=0.00000001) {
			R = (P->Ro + Ri);
			P->LiqL = P->LiqL+ P->Ls;//B.S., keep the liquid length in the evaporator
			P->UA_Liq=P->UA_Liq+1/R;//B.S. liquid heat transfer conductance	
		} else {
			R = P->Ro+Ri;
			P->UA_TP=P->UA_TP+1/R;//B.S., dry thermal resistance of the two-phase heat transfer
			P->TPL = P->TPL+ P->Ls;//B.S., length of two-phase flow
		}
	//-----------------------------------------B.S.

	} else {//calculating the heat transfer in the two-phase region


		EVA_Get_Q.TXPo=TXPo;
		hi = Zbrent(hi_max,hi_min,Get_Q_EVA,1e-5,&EVA_Get_Q);//iterate the heat transfer coefficient

		if(errorLog.IsError()) {
			// plot error vs. T to see whgat function looks like that caused error
		//	ZbrentPlot(hi_max,hi_min,Get_Q_EVA,1e-4,&EVA_Get_Q);//Haorong change from -7 to-2
			char msg[2048];
			sprintf(msg,"(Zbrent3) T=%.4lf, T_min=%.4lf, T_max=%.4lf",hi,hi_min,hi_max);
			errorLog.Add("EvapTubeLnew", msg);
			return 0;
		}; 

		q=EVA_Get_Q.q;

		//B.S.-----------------------------------
		P->Qtp_wet = P->Qtp_wet+EVA_Get_Q.P.Qtp_wet;//wet heat transfer amount
		P->Qtp_dry =P->Qtp_dry+EVA_Get_Q.P.Qtp_dry;//dry heat transfer amount
		P->UAw_TP = EVA_Get_Q.P.UAw_TP+P->UAw_TP;//heat conductance of two-phase wet heat transfer
		P->UA_TP = EVA_Get_Q.P.UA_TP+P->UA_TP;//heat conductance of two-phase dry heat transfer
		P->L_wet = EVA_Get_Q.P.L_wet + P->L_wet;//wet heat transfer length in the two-phase region
		P->L_dry = EVA_Get_Q.P.L_dry + P->L_dry;//wet heat transfer length in the two-phase region
		P->TPL = P->TPL+ P->Ls;//length of the two-phase heat transfer
		//------------------------------------------B.S.
	}

	HPo->H-=q;

	// calculate pressure drop
	double P_out;

	Gr = mr/P->Acs;//this mass flux is for calculating the pressure drop
	const double DP_FR = FricDP(TXPo, Gr, q, P);
	
	if(errorLog.IsError()) {
		errorLog.Add("EvapTubeL_new", "DP_FR");
		return 0;
	}
	

	if(TXPo.X>0.05&&TXPo.X<0.99)
	{
	//B.S., prepared to calculate the two-phase acceleration pressure drop
	PreAcc Preacc;
	Preacc.DP_FR=-1*DP_FR;
	Preacc.G=Gr;
	Preacc.H_OUT=HPo->H;
	Preacc.P_IN=TXPo.P;
	Preacc.X_IN=TXPo.X;
	
	const double DP_ACC=Zbrent(10,-10,GET_PreAcc,1e-7,&Preacc);//B.S., calculate the acceleration pressure drop 
	if(errorLog.IsError()) {
		errorLog.Add("EvapTubeL_new", "(Preacc)");
		return 0;
	}
	P_out=TXPo.P+DP_FR-DP_ACC;
	}
	else {
	P_out=TXPo.P+DP_FR;
	}

	HPo->P=P_out;

	// Determine mass of charge.  It is based on the inlet state.
	// The specific volume is recalculated so that a different model
	// can be used from the one used to calculate the pressure drop.

	TXPo = HPtoTXP(*HPo);
	if(errorLog.IsError()) {
		errorLog.Add("EvapTubeLnew", "(TXPo,2)");
		return 0;
	}

	v = VolumeALL(TXPo,Gr,P->Di,mr*q/P->Api);		// seperate flow model
	// v = PropertyTXP(VOL,TXPo);		// homogeneous flow model
	if(errorLog.IsError()) {
		errorLog.Add("EvapTubeLnew", "(v)");
		return 0;
	}

	m->V=P->Ls*P->Acs;//B.S., use the real cross-sectional area to calculate the inner volume
	m->m=m->V/v;

	// Calculate output air state
	hai = wair.h(TPi.T,TPi.P);
	if(errorLog.IsError()) {
		errorLog.Add("EvapTubeLnew", "(hai)");
		return 0;
	}

	const double W_I=wair.HumidityRatio(TPi.T,TPi.P);
	if(errorLog.IsError()) {
		errorLog.Add("EvapTubeLnew", "(W_I)");
		return 0;
	}

	const double HF_water=(EVA_Get_Q.T_S_O*4.1877+0.0594)*1e3;

	WHo->H=hai-q*mr/ma-HF_water*(W_I-EVA_Get_Q.W);

	WHo->W = EVA_Get_Q.W;//B.S. the outlet air humidity is calculated with the fuction Get_Q_EVA below.
	if(errorLog.IsError()) {
		errorLog.Add("EvapTubenew");
		return 0;
	}

	
	//B.S.----------------------------------------
	if(TXPo.X>=0.999999999999)//B.S., the following caculate the whole thermal resistances separately, prepared to adjust them separately.	
	{
	P->m_Vap = P->m_Vap+m->V/v;//B.S.vapor mass
	P->V_Vap = P->V_Vap + m->V;//B.S.vapor volume
	}
	else if(TXPo.X<=0.0000001)
	{
	P->m_Liq = P->m_Liq+m->V/v;//B.S., liquid mass
	P->V_Liq = P->V_Liq + m->V;//B.S., liquid volume
	
	if(TXP_bak.X>0.0000001)
	{
	P->HP_TP1.H=P->HP_TP1.H+mr*HPo->H;//B.S., corresponding to the first segment that involves with two-phase heat transfer, counted from the entrance of the evaporator
	P->HP_TP1.P=P->HP_TP1.P+mr*HPo->P;//B.S., corresponding to the first segment that involves with two-phase heat transfer, counted from the entrance of the evaporator
	P->count1 = P->count1+mr;//B.S., this addition is for counting how many evaporator circuits that involved with heat transfer calculation
	}

	}
	else 
	{
	P->m_TP = P->m_TP+m->V/v;//B.S., mass of two-phase flow
	P->V_TP = P->V_TP + m->V;//B.S., volume of two-phase flow

	if(TXP_bak.X>=0.999999999999)
	{
	P->HP_TP2.H=P->HP_TP2.H+mr*HPo->H;//B.S. corresponding to the last segment that involves with two-phase heat transfer, counted from the entrance of the evaporator
	P->HP_TP2.P=P->HP_TP2.P+mr*HPo->P;//B.S. corresponding to the last segment that involves with two-phase heat transfer, counted from the entrance of the evaporator
	P->count2 = P->count2+mr;//B.S., this addition is for counting how many evaporator circuits that involved with heat transfer calculation
	}
	}
	//---------------------------------------------B.S.

	return 0;
}

/*********************************************************************
Evaporator tube segment model (forward in flow direction) which considers pressure drops
Inputs:
	Gr = refrigerant mass flux (kg/s/m^2)
	HPi = refrigerant inlet state (h,P)
	Ga = air mass flux (kg/s/m^2)
	Tai = air inlet temperature (C)
Outputs:
	HPi = refrigerant outlet state (h,P)
	hao = air outlet enthalpy (J/kg)
	m = mass of charge in return bend (kg)
*********************************************************************/
int EvapTubeL_Fwd(double Gr,//refrigerant mass flux
					HP *HPo,//refrigerant outlet and inlet state
					double Ga,//air mass flux
					TP TPi,//air inlet state
					WH *WHo,//air outlet state
					MASS *m,//charge and inner volume in the evaporator
					ETdim *P)//evaporator struct
{
	double v,hai,ma,mr;
	TXP TXPo;
	EVA_Get_Q EVA_Get_Q;//B.S., this struct is for storing some parameters for iteration
	double y=0;
	double hi=0;
	double Ri=0,R=0;
	double q=0;

	/* Calculate air side resistance */
	P->wet=0;//without considering wet heat tranfer adjustment

	const double ho = P->hAirAdj*ConvCoeffAir_EVA(TPi,Ga,P);//B.S., the airside dry heat transfer coefficient can be got for different fin types
	if(errorLog.IsError()) {
		errorLog.Add("EvapTubeL", "get ho");
		return 0;}

	const double phi = FinEffect_Schmidt(ho,233,P->th,P->y,P->Do);//B.S. calculate the fin efficiency with the Schmidt equation
	if(errorLog.IsError()) {
		errorLog.Add("EvapTubeL", "FinEffect_Schmidt");
		return 0;}

	P->Ro = 1/(ho*(P->Apo+phi*P->Af));
	P->ho=ho;
	P->Ga=Ga;

	
	// convert inlet state into TXP format
	TXPo = HPtoTXP(*HPo);
	if(errorLog.IsError()) {
		errorLog.Add("EvapTubeL_new", "(TXPo,1)");
		return 0;}

	TXP TXP_bak;
	TXP_bak = TXPo;//backup the inlet state of this segment

	// Mass flow rates of air and refrigerant
	ma=Ga*P->Aflow;//B.S. Ga is the maximum airflow flux, P->Aflow is the minimun air flow cross-sectional area
	mr=Gr*P->Ax;

	// Calculate heat transfered per unit mass of refrigerant.

	const double hi_max=10e4;//largest possible refrigerant side heat tranfer coeffcient
	const double hi_min=100;//minimum possible refrigerant side heat transfer coeffcient

	//store the parameters for iteration
	EVA_Get_Q.ma=ma;
	EVA_Get_Q.mr=mr;
	EVA_Get_Q.TPi=TPi;
	EVA_Get_Q.Gr=Gr;
	EVA_Get_Q.P=*P;
	EVA_Get_Q.q=0;
	EVA_Get_Q.W=wair.HumidityRatio(TPi.T,TPi.P);
	if(errorLog.IsError()) {
		errorLog.Add("EvapTubeL_new", "(W)");
		return 0;}

	const double X1=0.05,X2=0.95;//	const double X1=0.01,X2=0.99;
	//for calculating the heat transfer and pressure drop in single-phase or approximated single-phase region

	if(TXPo.X>X2||TXPo.X<X1) {

		if(TXPo.X>=1 || TXPo.X<=0) {//Pure single phase region

			y = ConvCoeffSP(TXPo,Gr,P);
			//get the single-phase heat transfer coefficient
			if(errorLog.IsError()) {
				errorLog.Add("EvapTubeLNew","(hi_SP)");
				return 0;}
		} 
		else if(TXPo.X<X1) {
	
			TXP TXP1 = toTXP(TXPo.T,0.0,TXPo.P);
			TXP TXP2 = toTXP(TXPo.T,X1,TXPo.P);

			double y1 = ConvCoeffSP(TXP1,Gr,P);
			if(errorLog.IsError()) {
				errorLog.Add("EvapTubeLNew","(hi_SP)");
				return 0;}
		
			EVA_Get_Q.TXPo=TXP2;//local refrigerant state
			hi = Zbrent(hi_max,hi_min,Get_Q_EVA,1e-5,&EVA_Get_Q);//B.S., get the heat transfer coefficient at the two-phase region
			if(errorLog.IsError()) {
				errorLog.Add("EvapTubeLNew","Zbrent1");
				return 0;}
		
			double y2 = hi/P->hRefAdj;//ConvCoeffEvapTP_microfin(TXP2,Gr,P)
		
			if(y1>y2) y1=y2;//B.S., make the single phase flow less than two phase
			y = y1+TXPo.X*(y2-y1)/X1;//B.S., get the heat transfer coefficient in this region with intropolation
		} 
		else {

			TXP TXP1 = toTXP(TXPo.T,X2,TXPo.P);
			TXP TXP2 = toTXP(TXPo.T,1,TXPo.P);
	
			EVA_Get_Q.TXPo = TXP1;
			hi = Zbrent(hi_max,hi_min,Get_Q_EVA,1e-5,&EVA_Get_Q);//two-phase refrigerant heat transfer coefficient
			//B.S., get the heat transfer coefficient at the two-phase region
			if(errorLog.IsError()) {
				char msg[256];
				sprintf(msg,"Zbrent2: TXPo.T=%.2lf, T_max=%.2lf, T_min=%.2lf",TXPo.T,hi_max,hi_min);
				errorLog.Add("EvapTubeLNew",msg);
				return 0;}

			double y1 = hi/P->hRefAdj;//B.S. get two-phase heat transfer coefficient at this region	
			
			double y2 = ConvCoeffSP(TXP2,Gr,P);//B.S., single-phase heat transfer coefficient

			if(y2>y1) y2=y1; //B.S., for making the single phase flow less than two phase
			y = y2-(1-TXPo.X)*(y2-y1)/(1-X2);//B.S., get the heat transfer coefficient in this region with intropolation
		} 

		hi = P->hRefAdj*y;//B.S., heat transfer coefficient at this region
		Ri = 1/(hi*P->Api);//B.S., inside thermal resistance
		const double R_W=log(P->Do/(P->Do-2.0*P->xp))/(2.0*3.1415*P->K_T*P->Ls);
		R = P->Ro+R_W+Ri;//overall thermal resistance
	
		//B.S., prepare the parameters for iteration
		EVA_Get_Q.TXPo = TXPo;
		EVA_Get_Q.Cmin = CmineCrossFlow_dry(R,mr,ma,TXPo,TPi.T);
		EVA_Get_Q.HPo = *HPo;
	
		if(0/*TXPo.T<T_sat+1*/) {
			// TXPo.T is actually the outlet refrigerant temperature
			// however TXPo.T = TXPi.T in the two phase region and
			// TXPo.T ~ to TXPi.T in the superheated region if steps are small

			q = (TPi.T-TXPo.T)/mr*CmineCrossFlow_dry(R,mr,ma,TXPo,TPi.T);

			if(errorLog.IsError()) {
				errorLog.Add("EvapTubeL_new", "(q_singlephase)");
				return 0;
			}	

		} 
		else {//when it is close to the saturated state, the function zbrent can not converge
			EVA_Get_Q.hi=hi;//input the refrigerant side heat transfer coefficient
			Get_Q_Single_For(*HPo, &EVA_Get_Q);
			if(errorLog.IsError()) {
			errorLog.Add("EvapTubeL_new", "(q_singlephase_For)");
			return 0;}	
			q=EVA_Get_Q.q;
		}
	
	
		if(TXPo.X>=0.9999999) {//B.S., the following caculate the whole thermal resistances separately, prepared to adjust them separately.	
			R = (P->Ro+R_W+Ri);
			P->UA_Vap=P->UA_Vap+1/R;//B.S., vapor heat transfer conductance
			P->VapL = P->VapL+ P->Ls;//B.S., vapor length
		} else if(TXPo.X<=0.00000001) {
			R = (P->Ro + R_W+ Ri);
			P->LiqL = P->LiqL+ P->Ls;//B.S., keep the liquid length in the evaporator
			P->UA_Liq=P->UA_Liq+1/R;//B.S. liquid heat transfer conductance	
		} else {
			R = P->Ro+ R_W+ Ri;
			P->UA_TP=P->UA_TP+1/R;//B.S., dry thermal resistance of the two-phase heat transfer
			P->TPL = P->TPL+ P->Ls;//B.S., length of two-phase flow
		}


	} 
	else {//calculating the heat transfer in the two-phase region


		EVA_Get_Q.TXPo=TXPo;
		hi = Zbrent(hi_max,hi_min,Get_Q_EVA,1e-5,&EVA_Get_Q);//iterate the heat transfer coefficient

		if(errorLog.IsError()) {//try again
		errorLog.ClearError("EvapTubeLNew","zrent2");
		EVA_Get_Q.TXPo.X = 	EVA_Get_Q.TXPo.X+0.05;
		const double hi = Zbrent(hi_max,hi_min,Get_Q_EVA,1e-5,&EVA_Get_Q);
		EVA_Get_Q.TXPo.X = 	EVA_Get_Q.TXPo.X-0.05;
		Get_Q_EVA(hi,&EVA_Get_Q);
		}

		if(errorLog.IsError()) {
			// plot error vs. T to see whgat function looks like that caused error
		//	ZbrentPlot(hi_max,hi_min,Get_Q_EVA,1e-4,&EVA_Get_Q);//Haorong change from -7 to-2
			char msg[2048];
			sprintf(msg,"(Zbrent3) h=%.4lf, h_min=%.4lf, h_max=%.4lf",hi,hi_min,hi_max);
			errorLog.Add("EvapTubeLnew", msg);
			return 0;
		}; 

		q=EVA_Get_Q.q;

		P->Qtp_wet = P->Qtp_wet+EVA_Get_Q.P.Qtp_wet;//wet heat transfer amount
		P->Qtp_dry =P->Qtp_dry+EVA_Get_Q.P.Qtp_dry;//dry heat transfer amount
		P->UAw_TP = EVA_Get_Q.P.UAw_TP+P->UAw_TP;//heat conductance of two-phase wet heat transfer
		P->UA_TP = EVA_Get_Q.P.UA_TP+P->UA_TP;//heat conductance of two-phase dry heat transfer
		P->L_wet = EVA_Get_Q.P.L_wet + P->L_wet;//wet heat transfer length in the two-phase region
		P->L_dry = EVA_Get_Q.P.L_dry + P->L_dry;//wet heat transfer length in the two-phase region
		P->TPL = P->TPL+ P->Ls;//length of the two-phase heat transfer
	}

	HPo->H+=q;

	// calculate pressure drop
	double P_out;

	Gr = mr/P->Acs;//this mass flux is for calculating the pressure drop
	const double DP_FR = FricDP(TXPo, Gr, q, P);
	
	if(errorLog.IsError()) {
		errorLog.Add("EvapTubeL_new", "DP_FR");
		return 0;}
	

	if(TXPo.X>0.05&&TXPo.X<0.99)
	{
	//B.S., prepared to calculate the two-phase acceleration pressure drop
	PreAcc Preacc;
	Preacc.DP_FR=DP_FR;
	Preacc.G=Gr;
	Preacc.H_OUT=HPo->H;
	Preacc.P_IN=TXPo.P;
	Preacc.X_IN=TXPo.X;
	
	const double DP_ACC=Zbrent(10,-10,GET_PreAcc,1e-7,&Preacc);//B.S., calculate the acceleration pressure drop 
	if(errorLog.IsError()) {
		errorLog.Add("EvapTubeL_new", "(Preacc)");
		return 0;}
	P_out=TXPo.P-(DP_FR+DP_ACC)*P->PRefAdj;
	}
	else {
	P_out=TXPo.P-DP_FR*P->PRefAdj;
	}

	HPo->P=P_out;

	// Determine mass of charge.  It is based on the inlet state.
	// The specific volume is recalculated so that a different model
	// can be used from the one used to calculate the pressure drop.

	v = VolumeALL(TXPo,Gr,P->Di,mr*q/P->Api);		// seperate flow model
	// v = PropertyTXP(VOL,TXPo);		// homogeneous flow model
	if(errorLog.IsError()) {
		errorLog.Add("EvapTubeLnew", "(v)");
		return 0;}

	m->V=P->Ls*P->Acs;//B.S., use the real cross-sectional area to calculate the inner volume
	m->m=m->V/v;

	TXPo = HPtoTXP(*HPo);
	if(errorLog.IsError()) {
		errorLog.Add("EvapTubeLnew", "(TXPo,2)");
		return 0;}

	// Calculate output air state
	hai = wair.h(TPi.T,TPi.P);
	if(errorLog.IsError()) {
		errorLog.Add("EvapTubeLnew", "(hai)");
		return 0;}

	const double W_I=wair.HumidityRatio(TPi.T,TPi.P);
	if(errorLog.IsError()) {
		errorLog.Add("EvapTubeLnew", "(W_I)");
		return 0;}

	const double HF_water=(EVA_Get_Q.T_S_O*4.1877+0.0594)*1e3;

	WHo->H=hai-q*mr/ma-HF_water*(W_I-EVA_Get_Q.W);

	WHo->W = EVA_Get_Q.W;//B.S. the outlet air humidity is calculated with the fuction Get_Q_EVA below.
	if(errorLog.IsError()) {
		errorLog.Add("EvapTubenew");
		return 0;}

	
	//B.S.----------------------------------------
	if(TXP_bak.X>=0.999999999999)//B.S., the following caculate the whole thermal resistances separately, prepared to adjust them separately.	
	{
	P->m_Vap = P->m_Vap+m->V/v;//B.S.vapor mass
	P->V_Vap = P->V_Vap + m->V;//B.S.vapor volume
	}
	else if(TXP_bak.X<=0.0000001)
	{
	P->m_Liq = P->m_Liq+m->V/v;//B.S., liquid mass
	P->V_Liq = P->V_Liq + m->V;//B.S., liquid volume
	
	if(TXPo.X>0.0000001)
	{
	P->HP_TP1.H=P->HP_TP1.H+mr*HPo->H;//B.S., corresponding to the first segment that involves with two-phase heat transfer, counted from the entrance of the evaporator
	P->HP_TP1.P=P->HP_TP1.P+mr*HPo->P;//B.S., corresponding to the first segment that involves with two-phase heat transfer, counted from the entrance of the evaporator
	P->count1 = P->count1+mr;//B.S., this addition is for counting how many evaporator circuits that involved with heat transfer calculation
	}

	}
	else 
	{
	P->m_TP = P->m_TP+m->V/v;//B.S., mass of two-phase flow
	P->V_TP = P->V_TP + m->V;//B.S., volume of two-phase flow

	if(TXPo.X>=0.999999999999)
	{
	P->HP_TP2.H=P->HP_TP2.H+mr*HPo->H;//B.S. corresponding to the last segment that involves with two-phase heat transfer, counted from the entrance of the evaporator
	P->HP_TP2.P=P->HP_TP2.P+mr*HPo->P;//B.S. corresponding to the last segment that involves with two-phase heat transfer, counted from the entrance of the evaporator
	P->count2 = P->count2+mr;//B.S., this addition is for counting how many evaporator circuits that involved with heat transfer calculation
	}
	}
	//---------------------------------------------B.S.

	return 0;
}
#endif

/********************************************************************
Iteration to get evaporative heat transfer rate at two-phase region
wet coil calculation from 
Braun, J. E., Klein, S.A., and Michell, J.W., 1989,"effectiveness models for cooling towers
and cooling coils", ASHRAE Transactions, Vol. 95-2, pp. 164-174
***********************************************************************/
double Get_Q_EVA(double hi_0,//this function is for getting the refrigerant side heat transfer coefficent 
				 void *Params)//this struct contains the necessary parameters for iteration

{
	EVA_Get_Q* EVA_Q=(EVA_Get_Q*)Params;
	static double pi = acos(-1.0);

	double hi=0,Ri=0;
	double T_w=0;
	double Q=0;
	
	//B.S.--------------------------
	EVA_Q->P.Qtp_wet = 0;//new
	EVA_Q->P.Qtp_dry = 0;//new
	EVA_Q->P.UAw_TP = 0;
	EVA_Q->P.UA_TP = 0;
	EVA_Q->P.L_wet = 0;
	EVA_Q->P.L_dry = 0;
	//---------------------------B.S.
	
	hi=hi_0;

	const double W_I=wair.HumidityRatio(EVA_Q->TPi.T,EVA_Q->TPi.P);//B.S., inlet air humidity
	if(errorLog.IsError()) {
		errorLog.Add("Get_Q_EVA","wair.HumidityRatio1");
		return 0;
	}

	
	
	ETdim Eva_dim;
	Eva_dim = EVA_Q->P;//B.S. get the evaporator struct parameters
	

	//dry condition
	Ri = 1/(hi*EVA_Q->P.Api);//B.S., inside thermal resistance
	const double R_W=log(EVA_Q->P.Do/(EVA_Q->P.Do-2.0*EVA_Q->P.xp))/(2.0*pi*EVA_Q->P.K_T*EVA_Q->P.Ls);
	const double R = EVA_Q->P.Ro+R_W+Ri;//B.S., external thermal resistance under dry condition

	const double CP_M=wair.Cp(EVA_Q->TPi.T,EVA_Q->TPi.P);
	if(errorLog.IsError()) {
		errorLog.Add("Get_Q_EVA","CP_M");
		return 0;
	}
	

//begin with wet condition
		const double H_A_I=wair.h(EVA_Q->TPi.T,EVA_Q->TPi.P);//inlet air enthalpy
		if(errorLog.IsError()) {
		errorLog.Add("Get_Q_EVA","H_A_I");
		return 0;
		}
		
		const double H_SAT=wair.h(EVA_Q->TXPo.T,0.999);
		if(errorLog.IsError()) {
		errorLog.Add("Get_Q_EVA","H_SAT");
		return 0;
		}

		
		const double enthal1=wair.h(EVA_Q->TXPo.T,0.999);
		const double enthal2=wair.h(EVA_Q->TXPo.T-0.1,0.999);

		if(errorLog.IsError()) {
		errorLog.Add("Get_Q_EVA","get C_S");
		return 0;
		}

		const double C_S=(enthal1-enthal2)/0.1e0;

		
		const double M_DOT_A=EVA_Q->ma;
		EVA_Q->P.wet=1;
		const double H_W= EVA_Q->P.hAirAdj*ConvCoeffAir_EVA(EVA_Q->TPi,EVA_Q->P.Ga,&EVA_Q->P);
		EVA_Q->P.wet=0;
		//calculate wet fin efficiency
		const double M_F_W=pow((2e0*H_W*C_S/(EVA_Q->P.K_F*EVA_Q->P.th*CP_M)),0.5e0);//for calculating the wet fin efficiency, including both the heat and mass transfer
		const double ETA_F_W=tanh(M_F_W*EVA_Q->P.L_F)/(M_F_W*EVA_Q->P.L_F);//wet fin efficiency
		
		const double A_F=EVA_Q->P.Af;//fin surface area
		const double A_T=EVA_Q->P.Af+EVA_Q->P.Apo;//the whole external area
		const double ETA_O_W=1e0-A_F/A_T*(1e0-ETA_F_W);//overall external fin efficiency
		const double R_SURF_W=1/(ETA_O_W*H_W*A_T);//airside thermal resistance under wet condition
		const double NTU_O_W=1/(R_SURF_W*M_DOT_A*CP_M);//for calculating the effective saturated air enthalpy at T_S_O

		const double NTU_O=1/(EVA_Q->P.Ro*M_DOT_A*CP_M);//for calculating the outlet air temperature, which is involved with sensible heat transfer

		const double UA_W=1e0/(C_S*Ri+C_S*R_W+CP_M*R_SURF_W);//UA under wet condition, considering both the heat and mass transfer
		const double NTU_W=1*UA_W/(M_DOT_A);//NTU under wet condition, for calculating heat transfer
		const double EPSILON_W=1e0-exp(-NTU_W);//epsilon under wet condition

		Q=EPSILON_W*M_DOT_A*(H_A_I-H_SAT);//heat transfer under wet condition
		EVA_Q->q =Q/EVA_Q->mr;

		//get outlet humidity and effective surface temperature
		const double H_A_O=H_A_I-Q/(M_DOT_A);//outlet air enthalpy, actually the energy balance shoudl include the water flowing away, like -(83.84e3+(T_S-20)*4.183e3)*(W_I-W_O), since this part is small and need iteration, so it is ignored
		double H_S_S_O=H_A_I-(H_A_I-H_A_O)/(1-exp(-NTU_O_W));//saturated air enthalpy at T_S_S_O, effective surface temperature
		

		TP TP_S_O;
//		const double H_zero=wair.h(2,0.99);//fixed the minimum possible H_S_S_O
//		if(H_S_S_O<=H_zero) H_S_S_O=H_zero;
		
		TP_S_O=HPtoTP(H_S_S_O,0.999);
		if(errorLog.IsError()) {
		errorLog.Add("Get_Q_EVA","TP_S_O");
		return 0;
		}
		const double T_S_O=TP_S_O.T;//effective surface temperature
		const double T_O=T_S_O+(EVA_Q->TPi.T-T_S_O)*exp(-NTU_O);//outlet air temperature
		EVA_Q->T_S_O= T_S_O;

		TP TP_dew;
		TP_dew=WPtoTP(W_I,0.999);//get the dew point corresponding to the inlet humidity ratio
		if(errorLog.IsError()) {
		errorLog.Add("Get_Q_EVA","TP_dew");
		return 0;
		}
		const double T_dew=TP_dew.T;//dew temperature corresponding to the inlet air enthalpy
	
		if(T_S_O>=T_dew) //dry condition
		{

	// TXPo.T is actually the outlet refrigerant temperature
	// however TXPo.T = TXPi.T in the two phase region and
	// TXPo.T ~ to TXPi.T in the superheated region if steps are small
		const double q_dry=(EVA_Q->TPi.T-EVA_Q->TXPo.T)/EVA_Q->mr*CmineCrossFlow_wet(R,EVA_Q->mr,EVA_Q->ma,EVA_Q->TXPo,EVA_Q->TPi.T,CP_M);
		if(T_S_O>=T_dew+0.2)
		{
		EVA_Q->q = (EVA_Q->TPi.T-EVA_Q->TXPo.T)/EVA_Q->mr*CmineCrossFlow_wet(R,EVA_Q->mr,EVA_Q->ma,EVA_Q->TXPo,EVA_Q->TPi.T,CP_M);
		}
		else
		{
		EVA_Q->q=(q_dry-Q/EVA_Q->mr)/(0.2)*(T_S_O-T_dew) + Q/EVA_Q->mr;
		}//B.S.

		if(errorLog.IsError()) {
		errorLog.Add("Get_Q_EVA","q_dry");
		return 0;
		}

		Q=EVA_Q->q*EVA_Q->mr;//heat transfer amount under dry condition
		
		EVA_Q->W=W_I;//without dehumidifying
		T_w=Q*Ri+EVA_Q->TXPo.T;
		
		EVA_Q->T_w=T_w;
		//B.S.---------------------------
		EVA_Q->P.Qtp_dry = Q;//dry heat transfer
		EVA_Q->P.UA_TP = 1/R;//dry heat transfer conductance of this segment
		EVA_Q->P.L_dry = EVA_Q->P.Ls;//dry heat transfer tube length of this segment
		//--------------------------------B.S.

		Eva_dim.T_w=T_w;
		Eva_dim.q_flux=Q/EVA_Q->P.Api;
		hi =EVA_Q->P.hRefAdj*ConvCoeffEvapTP_microfin(EVA_Q->TXPo,EVA_Q->Gr,&Eva_dim);//B.S., based on the tube wall temperature to get the refrigerant side heat transfer coefficient
		
		if(errorLog.IsError()) {
		errorLog.Add("Get_Q_EVA","hi");
		return 0;
		}

		const double hi_dry = (hi-hi_0)/(hi+hi_0);
		return hi_dry;
		}

		TP TP_O;
		TP_O.T=T_O;
		const double HH_min = wair.h(T_O,0.05);
		const double HH_max = wair.h(T_O,0.999);

		if(errorLog.IsError()) {
		errorLog.Add("Get_Q_EVA","HH_min&HH_max");
		return 0;
		}
		
		if(H_A_O>=HH_max)
		{
		TP_O.T=T_O;
		TP_O.P=0.999;
		}
		else if(H_A_O<=HH_min)
		{
		TP_O.T=T_O;
		TP_O.P=0.05;
		}
		else
		{
		TP_O=THtoTP(T_O,H_A_O);//outlet air state
		if(errorLog.IsError()) {
		errorLog.ClearError("Get_Q_EVA","TP_O");
		TP_O.T=T_O;
		TP_O.P=0.999;
		}
		}

		double W_O=wair.HumidityRatio(T_O,TP_O.P);//outlet air humidity ratio
		if(errorLog.IsError()) {
		errorLog.ClearError("Get_Q_EVA","W_O");
		W_O=W_I;
		}

		EVA_Q->W=W_O;
		T_w=Q*Ri+EVA_Q->TXPo.T;
		EVA_Q->T_w=T_w;

		//B.S.---------------------------
		EVA_Q->P.Qtp_wet = EVA_Q->mr*EVA_Q->q;//wet heat transfer
		EVA_Q->P.UAw_TP = UA_W;//wet heat transfer conductance of this segment
		EVA_Q->P.L_wet = EVA_Q->P.Ls;//wet heat transfer tube length of this segment
		//--------------------------------B.S.

		Eva_dim.T_w=T_w;
		Eva_dim.q_flux=Q/EVA_Q->P.Api;
		hi =EVA_Q->P.hRefAdj*ConvCoeffEvapTP_microfin(EVA_Q->TXPo,EVA_Q->Gr,&Eva_dim);//B.S., based on the tube wall temperature to get the refrigerant side heat transfer coefficient

		if(errorLog.IsError()) {
		errorLog.Add("Get_Q_EVA","hi");
		return 0;
		}

	double dhi=(hi-hi_0)/(hi+hi_0);
	
	return dhi;
}

//B.S.------------------------------------------------------------------------
//for iteration to get evaporative heat transfer at single-phase region (backward in flow direction), use the inlet refrigerant state as the reference state

double Get_Q_Single(double H_in, //inlet refrigerant temperature in the segment
					void *Params)//this struct stores the parameters for iteration
{
	EVA_Get_Q* EVA_Q=(EVA_Get_Q*)Params;
	static double pi = acos(-1.0);
	TP TPi;
	HP HPo;
	TXP TXPo;
	const double mr=EVA_Q->mr;
	double q=0,Q=0;
	double Ri=0;
	const double Gr=EVA_Q->Gr;

	ETdim Eva_dim;
	Eva_dim = EVA_Q->P;//B.S. get the evaporator struct parameters

	HPo=EVA_Q->HPo;//refrigerant outlet state
	const double H_out=HPo.H;//backup the outlet enthalpy
	TPi = EVA_Q->TPi;//air inlet state
	TXPo = EVA_Q->TXPo;//refrigerant outlet state
	
	q=HPo.H-H_in;
	
	const double DP_SP = FricDP(TXPo, Gr, q, &Eva_dim);//single-phase pressure drop
	const double P_in=TXPo.P+DP_SP;
	HPo.P = P_in;
	HPo.H=H_in;

	TXPo = HPtoTXP(HPo);
	if(errorLog.IsError()) {
		errorLog.Add("Get_Q_Single", "(TXPo)");
		return 0;
	}

	const double W_I=wair.HumidityRatio(EVA_Q->TPi.T,EVA_Q->TPi.P);//B.S., inlet air humidity
	if(errorLog.IsError()) {
		errorLog.Add("Get_Q_Single","wair.HumidityRatio1");
		return 0;
	}

	
	//dry condition
	Ri = 1/(EVA_Q->hi*EVA_Q->P.Api);//B.S., inside thermal resistance
	const double R_W=log(EVA_Q->P.Do/(EVA_Q->P.Do-2.0*EVA_Q->P.xp))/(2.0*pi*EVA_Q->P.K_T*EVA_Q->P.Ls);
	const double R = EVA_Q->P.Ro+ R_W +Ri;//B.S., external thermal resistance under dry condition

	const double CP_M=wair.Cp(EVA_Q->TPi.T,EVA_Q->TPi.P);
	if(errorLog.IsError()) {
		errorLog.Add("Get_Q_Single","CP_M");
		return 0;
	}

//begin with wet condition
		const double H_A_I=wair.h(EVA_Q->TPi.T,EVA_Q->TPi.P);//inlet air enthalpy
		if(errorLog.IsError()) {
		errorLog.Add("Get_Q_Single","H_A_I");
		return 0;
		}

		const double H_SAT=wair.h(TXPo.T,0.999);
		if(errorLog.IsError()) {
		errorLog.Add("Get_Q_Single","H_SAT");
		return 0;
		}

		
		const double enthal1=wair.h(TXPo.T,0.999);
		const double enthal2=wair.h(TXPo.T-0.1,0.999);

		if(errorLog.IsError()) {
		errorLog.Add("Get_Q_Single","get C_S");
		return 0;
		}

		const double C_S=(enthal1-enthal2)/0.1e0;

		
		const double M_DOT_A=EVA_Q->ma;
		EVA_Q->P.wet=1;
		const double H_W=EVA_Q->P.hAirAdj*ConvCoeffAir_EVA(EVA_Q->TPi,EVA_Q->P.Ga,&EVA_Q->P);
		EVA_Q->P.wet=0;
		//calculate wet fin efficiency
		const double M_F_W=pow((2e0*H_W*C_S/(EVA_Q->P.K_F*EVA_Q->P.th*CP_M)),0.5e0);//for calculating the wet fin efficiency, including both the heat and mass transfer
		const double ETA_F_W=tanh(M_F_W*EVA_Q->P.L_F)/(M_F_W*EVA_Q->P.L_F);//wet fin efficiency
		
		const double A_F=EVA_Q->P.Af;//fin surface area
		const double A_T=EVA_Q->P.Af+EVA_Q->P.Apo;//the whole external area
		const double ETA_O_W=1e0-A_F/A_T*(1e0-ETA_F_W);//overall external fin efficiency
		const double R_SURF_W=1/(ETA_O_W*H_W*A_T);//airside thermal resistance under wet condition
		const double NTU_O_W=1/(R_SURF_W*M_DOT_A*CP_M);//for calculating the effective saturated air enthalpy at T_S_O

		const double NTU_O=1/(EVA_Q->P.Ro*M_DOT_A*CP_M);//for calculating the outlet air temperature, which is involved with sensible heat transfer

		const double UA_W=1e0/(C_S*Ri+C_S*R_W+CP_M*R_SURF_W);//UA under wet condition, considering both the heat and mass transfer
		const double NTU_W=1*UA_W/(M_DOT_A);//NTU under wet condition, for calculating heat transfer
		const double EPSILON_W=1e0-exp(-NTU_W);//epsilon under wet condition

		Q=EPSILON_W*M_DOT_A*(H_A_I-H_SAT);//heat transfer under wet condition
		EVA_Q->q =Q/EVA_Q->mr;
		q=Q/EVA_Q->mr;

		//get outlet humidity and effective surface temperature
		const double H_A_O=H_A_I-Q/(M_DOT_A);//outlet air enthalpy, actually the energy balance shoudl include the water flowing away, like -(83.84e3+(T_S-20)*4.183e3)*(W_I-W_O), since this part is small and need iteration, so it is ignored
		double H_S_S_O=H_A_I-(H_A_I-H_A_O)/(1-exp(-NTU_O_W));//saturated air enthalpy at T_S_S_O, effective surface temperature
		

		TP TP_S_O;
//		const double H_zero=wair.h(2,0.99);//fixed the minimum possible H_S_S_O
//		if(H_S_S_O<=H_zero) H_S_S_O=H_zero;
		
		TP_S_O=HPtoTP(H_S_S_O,0.999);
		if(errorLog.IsError()) {
		errorLog.Add("Get_Q_Single","TP_S_O");
		return 0;
		}
		const double T_S_O=TP_S_O.T;//effective surface temperature
		const double T_O=T_S_O+(EVA_Q->TPi.T-T_S_O)*exp(-NTU_O);//outlet air temperature
		EVA_Q->T_S_O= T_S_O;

		TP TP_dew;
		TP_dew=WPtoTP(W_I,0.999);//get the dew point corresponding to the inlet humidity ratio
		if(errorLog.IsError()) {
		errorLog.Add("Get_Q_Single","TP_dew");
		return 0;
		}
		const double T_dew=TP_dew.T;//dew temperature corresponding to the inlet air enthalpy
				
		if(T_S_O>T_dew) //dry condition
		{

		EVA_Q->q = (EVA_Q->TPi.T-TXPo.T)/EVA_Q->mr*CmineCrossFlow_wet(R,EVA_Q->mr,EVA_Q->ma,EVA_Q->TXPo,EVA_Q->TPi.T,CP_M);
		
		if(errorLog.IsError()) {
		errorLog.Add("Get_Q_Single","q_dry");
		return 0;
		}
		q= EVA_Q->q;
		Q=EVA_Q->q*EVA_Q->mr;//heat transfer amount under dry condition
		
		EVA_Q->W=W_I;//without dehumidifying
		}

		TP TP_O;
		TP_O.T=T_O;
		const double HH_min = wair.h(T_O,0.05);
		const double HH_max = wair.h(T_O,0.999);

		if(errorLog.IsError()) {
		errorLog.Add("Get_Q_Single","HH_min&HH_max");
		return 0;
		}
		
		if(H_A_O>=HH_max)
		{
		TP_O.T=T_O;
		TP_O.P=0.999;
		}
		else if(H_A_O<=HH_min)
		{
		TP_O.T=T_O;
		TP_O.P=0.05;
		}
		else
		{
		TP_O=THtoTP(T_O,H_A_O);//outlet air state
		if(errorLog.IsError()) {
		errorLog.ClearError("Get_Q_Single","TP_O");
		TP_O.T=T_O;
		TP_O.P=0.999;
		}
		}

		double W_O=wair.HumidityRatio(T_O,TP_O.P);//outlet air humidity ratio
		if(errorLog.IsError()) {
		errorLog.ClearError("Get_Q_Single","W_O");
		W_O=W_I;
		}

		EVA_Q->W=W_O;
	
		HPo.H=H_out - q;
		const double dH=(HPo.H-H_in)/(HPo.H+H_in);
		return dH;
}
//--------------------------------------------------B.S.

//for iteration to get evaporative heat transfer at single-phase region (forward in flow direction), use the inlet refrigerant state as the reference state

double Get_Q_Single_For(HP HPi, //inlet refrigerant temperature in the segment
					void *Params)//this struct stores the parameters for iteration
{
	EVA_Get_Q* EVA_Q=(EVA_Get_Q*)Params;
	static double pi = acos(-1.0);
	TP TPi;
	TXP TXPi;
	const double mr=EVA_Q->mr;
	double Q=0;
	double Ri=0;
	const double Gr=EVA_Q->Gr;

	ETdim Eva_dim;
	Eva_dim = EVA_Q->P;//B.S. get the evaporator struct parameters

	TPi = EVA_Q->TPi;//air inlet state
	

	TXPi = HPtoTXP(HPi);
	if(errorLog.IsError()) {
		errorLog.Add("Get_Q_Single", "(TXPo)");
		return 0;
	}

	const double W_I=wair.HumidityRatio(EVA_Q->TPi.T,EVA_Q->TPi.P);//B.S., inlet air humidity
	if(errorLog.IsError()) {
		errorLog.Add("Get_Q_Single","wair.HumidityRatio1");
		return 0;
	}

	
	//dry condition
	Ri = 1/(EVA_Q->hi*EVA_Q->P.Api);//B.S., inside thermal resistance
	const double R_W=log(EVA_Q->P.Do/(EVA_Q->P.Do-2.0*EVA_Q->P.xp))/(2.0*pi*EVA_Q->P.K_T*EVA_Q->P.Ls);
	const double R = EVA_Q->P.Ro+ R_W +Ri;//B.S., external thermal resistance under dry condition

	const double CP_M=wair.Cp(EVA_Q->TPi.T,EVA_Q->TPi.P);
	if(errorLog.IsError()) {
		errorLog.Add("Get_Q_Single","CP_M");
		return 0;
	}

//begin with wet condition
		const double H_A_I=wair.h(EVA_Q->TPi.T,EVA_Q->TPi.P);//inlet air enthalpy
		if(errorLog.IsError()) {
		errorLog.Add("Get_Q_Single","H_A_I");
		return 0;
		}

		const double H_SAT=wair.h(TXPi.T,0.999);
		if(errorLog.IsError()) {
		errorLog.Add("Get_Q_Single","H_SAT");
		return 0;
		}

		
		const double enthal1=wair.h(TXPi.T,0.999);
		const double enthal2=wair.h(TXPi.T-0.1,0.999);

		if(errorLog.IsError()) {
		errorLog.Add("Get_Q_Single","get C_S");
		return 0;
		}

		const double C_S=(enthal1-enthal2)/0.1e0;

		
		const double M_DOT_A=EVA_Q->ma;
		EVA_Q->P.wet=1;
		const double H_W=EVA_Q->P.hAirAdj*ConvCoeffAir_EVA(EVA_Q->TPi,EVA_Q->P.Ga,&EVA_Q->P);
		EVA_Q->P.wet=0;
		//calculate wet fin efficiency
		const double M_F_W=pow((2e0*H_W*C_S/(EVA_Q->P.K_F*EVA_Q->P.th*CP_M)),0.5e0);//for calculating the wet fin efficiency, including both the heat and mass transfer
		const double ETA_F_W=tanh(M_F_W*EVA_Q->P.L_F)/(M_F_W*EVA_Q->P.L_F);//wet fin efficiency
		
		const double A_F=EVA_Q->P.Af;//fin surface area
		const double A_T=EVA_Q->P.Af+EVA_Q->P.Apo;//the whole external area
		const double ETA_O_W=1e0-A_F/A_T*(1e0-ETA_F_W);//overall external fin efficiency
		const double R_SURF_W=1/(ETA_O_W*H_W*A_T);//airside thermal resistance under wet condition
		const double NTU_O_W=1/(R_SURF_W*M_DOT_A*CP_M);//for calculating the effective saturated air enthalpy at T_S_O

		const double NTU_O=1/(EVA_Q->P.Ro*M_DOT_A*CP_M);//for calculating the outlet air temperature, which is involved with sensible heat transfer

		const double UA_W=1e0/(C_S*Ri+C_S*R_W+CP_M*R_SURF_W);//UA under wet condition, considering both the heat and mass transfer
		const double NTU_W=1*UA_W/(M_DOT_A);//NTU under wet condition, for calculating heat transfer
	
		//const double EPSILON_W=1e0-exp(-NTU_W);//epsilon under wet condition, at two-phase region
		//pay attention to heat and mass transfer at single-phase section
		TXP TXP_prop;
		TXP_prop = TXPi;
		TXP_prop.X=1;
		const double CP_r = PropertyTXPtr(SPEC,TXP_prop);
		const double Cratio = M_DOT_A/(EVA_Q->mr*CP_r/C_S);
		const double l = -1*Cratio*(1-exp(-1*NTU_W));
		const double EPSILON_W = (1/Cratio)*(1-exp(l));
		//pay attention to heat and mass transfer at single-phase section

		Q=EPSILON_W*M_DOT_A*(H_A_I-H_SAT);//heat transfer under wet condition
		EVA_Q->q =Q/EVA_Q->mr;
	
		//get outlet humidity and effective surface temperature
		const double H_A_O=H_A_I-Q/(M_DOT_A);//outlet air enthalpy, actually the energy balance shoudl include the water flowing away, like -(83.84e3+(T_S-20)*4.183e3)*(W_I-W_O), since this part is small and need iteration, so it is ignored
		double H_S_S_O=H_A_I-(H_A_I-H_A_O)/(1-exp(-NTU_O_W));//saturated air enthalpy at T_S_S_O, effective surface temperature
		

		TP TP_S_O;
//		const double H_zero=wair.h(2,0.99);//fixed the minimum possible H_S_S_O
//		if(H_S_S_O<=H_zero) H_S_S_O=H_zero;
		
		TP_S_O=HPtoTP(H_S_S_O,0.999);
		if(errorLog.IsError()) {
		errorLog.Add("Get_Q_Single","TP_S_O");
		return 0;
		}
		const double T_S_O=TP_S_O.T;//effective surface temperature
		const double T_O=T_S_O+(EVA_Q->TPi.T-T_S_O)*exp(-NTU_O);//outlet air temperature
		EVA_Q->T_S_O= T_S_O;

		TP TP_dew;
		TP_dew=WPtoTP(W_I,0.999);//get the dew point corresponding to the inlet humidity ratio
		if(errorLog.IsError()) {
		errorLog.Add("Get_Q_Single","TP_dew");
		return 0;
		}
		const double T_dew=TP_dew.T;//dew temperature corresponding to the inlet air enthalpy
				
		if(T_S_O>T_dew) //dry condition
		{

		EVA_Q->q = (EVA_Q->TPi.T-TXPi.T)/EVA_Q->mr*CmineCrossFlow_wet(R,EVA_Q->mr,EVA_Q->ma,EVA_Q->TXPo,EVA_Q->TPi.T,CP_M);
		
		if(errorLog.IsError()) {
		errorLog.Add("Get_Q_Single","q_dry");
		return 0;
		}
		Q=EVA_Q->q*EVA_Q->mr;//heat transfer amount under dry condition
		
		EVA_Q->W=W_I;//without dehumidifying
		return 0;
		}

		TP TP_O;
		TP_O.T=T_O;
		const double HH_min = wair.h(T_O,0.05);
		const double HH_max = wair.h(T_O,0.999);

		if(errorLog.IsError()) {
		errorLog.Add("Get_Q_Single","HH_min&HH_max");
		return 0;
		}
		
		if(H_A_O>=HH_max)
		{
		TP_O.T=T_O;
		TP_O.P=0.999;
		}
		else if(H_A_O<=HH_min)
		{
		TP_O.T=T_O;
		TP_O.P=0.05;
		}
		else
		{
		TP_O=THtoTP(T_O,H_A_O);//outlet air state
		if(errorLog.IsError()) {
		errorLog.ClearError("Get_Q_Single","TP_O");
		TP_O.T=T_O;
		TP_O.P=0.999;
		}
		}

		double W_O=wair.HumidityRatio(T_O,TP_O.P);//outlet air humidity ratio
		if(errorLog.IsError()) {
		errorLog.ClearError("Get_Q_Single","W_O");
		W_O=W_I;
		}

		EVA_Q->W=W_O;
	
		return 0;
}
//--------------------------------------------------B.S.

//------------------ Water chiller evaporator --------------------
//------------------ Water chiller evaporator --------------------
//------------------ Water chiller evaporator --------------------
//
// Shell-and-tube evaporator with refrigerant flowing through
// the tubes.
// Given:
//    1. refrigerant flow rate and outlet state
//    2. water flow rate and inlet state
// Calculate:
//    1. refrigerant inlet state
//    2. water outlet state
//    3. mass of refrigerant in tubes
int Evaporator(const char* filename,double mr,HP HPo,double mw,double Twi,HP& HPi,double& Two,MASS& m,double* /*Prms*/)
{
	// Load model parameters from parameter file
	static double UA = 1100;
	static double Vtotal = 0.15;
	static int nPasses = 30;

	if(evapInit) {
		FILE *fp = fopen(filename,"r") ;
		int type = PosComponentData(fp,EVAPORATOR,1);
		if(type!=303) {
			errorLog.Add("Evaporator(shell-and-tube","Parameters not found");
			return 0;
		}
		UA = GetDouble(fp);
		Vtotal = GetDouble(fp);
		nPasses = (int)GetDouble(fp);
		fclose(fp);
		evapInit=0;
	}

	// no refrigerant pressure drop
	HPi.P = HPo.P;

	// From Incropera & DeWitt p.658-660
	//
	// q = Cmin*e*(Thi-Tci)
	// q = Ch*(Thi-Tho)
	// q = Cc*(Tco-Tci)
	//
	// Therefore (derived):
	//    q = (Tho-Tci)*Ch*Cmin*e/(Ch-Cmin*e)
	//    q = (Thi-Tco)*Cc*Cmin*e/(Cc-Cmin*e)
	// For Cc very large (two phase refrigerant):
	//    q = (Thi-Tco)*Cmin*e
	//
	// Know: Tco (refrigerant) and Thi (water)

	// heat transfer resistance
	double R = 1/UA;
	double Cw = mw*4200.0;

	// refrigerant mass
	m.m = m.V = 0;

	#ifdef EXTRA_REPORTING
	{
		FILE* fp = fopen("evap.dbg","w");
		if(fp) {
			fprintf(fp,"UA=%.3lf HPo=%.3lf,%.3lf\n",UA,0.001*HPo.H,HPo.P);
			fclose(fp);
		}
	}
	#endif

	double q;
	for(int i=0;i<nPasses;i++) {

		TXP TXPo = HPtoTXP(HPo);
		if(errorLog.IsError()) {
			char str[256];
			sprintf(str,"HPo(%.3lf,%.3lf)->TXPo(pass=%d)",0.001*HPo.H,HPo.P,i);
			errorLog.Add("Evaporator(shell-and-tube)",str);
			return 0;
		}

		double Cmine = CmineShellAndTube(nPasses*R,mr,mw,TXPo,Twi);
		if(errorLog.IsError()) {
			errorLog.Add("Evaporator(shell-and-tube)","2");
			return 0;
		}

		if(IsTwoPhase(TXPo.X)) {
			q = (Twi-TXPo.T)*Cmine;
		} else {
			double Cr = mr*PropertyTXPtr(SPEC,TXPo);
			if(errorLog.IsError()) {
				errorLog.Add("Evaporator(shell-and-tube)","3");
				return 0;
			}
			q = (Twi-TXPo.T)*Cr*Cmine/(Cr-Cmine);
		}

		Twi -= q/Cw;
		HPo.H -= q/mr;

		double Vpass = Vtotal/(double)nPasses;
		double v = PropertyTXPth(VOL,TXPo);   // homogeneous flow model
		if(errorLog.IsError()) {
			errorLog.Add("Evaporator(shell-and-tube)","4");
			return 0;
		}
		m.V += Vpass;
		m.m += Vpass/v;

		#ifdef EXTRA_REPORTING
		{
			FILE* fp = fopen("evap.dbg","a");
			if(fp) {
				fprintf(fp,"q = %.3lf, Tw = %.3lf, HTXP = %.3lf,%.3lf,%.3lf,%.3lf\n",q,Twi,0.001*HPo.H,TXPo.T,TXPo.X,TXPo.P);
				fclose(fp);
			}
		}
		#endif
	}
	Two = Twi;
	HPi.H = HPo.H;

	return 0;
}

//------------------------------------------------------------------------
//-----------------------------------------------------------------------
//--------------------------------------------------------------------------
//------------------------------------------------------------------------
//-----------------------------------------------------------------------
//--------------------------------------------------------------------------
//the followings are some old programs, not used currently

/********************************************************************
Solves the evaporator model iteratively.  It takes known values of
the outlet pressure and the inlet enthalpy and solves for the inlet
pressure and the outlet enthalpy.
********************************************************************/
void Evaporator1(double Hout,const char *filename,double mr,double Hin,double Pout,double Ga,TP TPi, HP *HPi,HP *HPo,TP *TPo,MASS *m,double *Aflow,double *Prms)
{
	Evap1P P;
	P.filename=filename;
	P.mr=mr;
	P.Hin=Hin;
	P.Pout=Pout;
	P.Ga=Ga;
	P.TPi=TPi;
	P.Prms=Prms;

	double Ho = FindZero(Hout,Evap1Func,1e-7,&P);
	if(errorLog.IsError()) {
		errorLog.Add("Evaporator1");
		return;
	}

	HPo->H = Ho;
	HPo->P = Pout;
	
	ETdim Evap_struc;
	Evaporator(filename,mr,HPo,Ga,TPi,HPi,TPo,m,Aflow,&Evap_struc,Prms);
	if(errorLog.IsError()) {
		errorLog.Add("Evaporator1");
		return;
	}
}

/********************************************************************
Function used by Evaporator1() during its iteration.
********************************************************************/
static double Evap1Func(double Hout,void *Params)
{
	Evap1P* P=(Evap1P*)Params;

	HP HPo;
	HPo.H = Hout;
	HPo.P = P->Pout;

	HP HPi;
	TP TPo;
	MASS m;
	double Aflow;

	ETdim Evap_struc;
	Evaporator(P->filename,P->mr,&HPo,P->Ga,P->TPi,&HPi,&TPo,&m,&Aflow,&Evap_struc,P->Prms);
	if(errorLog.IsError()) {
		errorLog.Add("Evap1Func");
		return 0;
	}

	double dH = (HPi.H-P->Hin)/P->Hin;

	return dH;
}


#ifndef EVAPORATOR_PRESSURE_DROP
/*********************************************************************
Evaporator tube segment model which neglects pressure drops
Inputs:
	Gr = refrigerant mass flux (kg/s/m^2)
	HPi = refrigerant inlet state (h,P)
	Ga = air mass flux (kg/s/m^2)
	Tai = air inlet temperature (C)
Outputs:
	HPi = refrigerant outlet state (h,P)
	hao = air outlet enthalpy (J/kg)
	m = mass of charge in return bend (kg)
*********************************************************************/
static int EvapTubeLnoDP(double Gr,HP *HPo,double Ga,TP TPi,WH *WHo,MASS *m,ETdim *P)
{
	double v,hai,ma,mr;
	TXP TXPo;
	double hi,R,Ri,q;

	// convert inlet state into TXP format
	TXPo = HPtoTXP(*HPo);
	if(errorLog.IsError()) {
		errorLog.Add("EvapTubeLnoDP (TXPo,1)");
		return 0;
	}

	// Mass flow rates of air and refrigerant
	ma=Ga*P->Aflow;
	mr=Gr*P->Ax;

	// Calculate heat transfered per unit mass of refrigerant.
	hi = P->hRefAdj*ConvCoeffInside(TXPo,Gr,P->Di);
	if(errorLog.IsError()) {
		errorLog.Add("EvapTubeLnoDP (hi)");
		return 0;
	}

	Ri = 1/(hi*P->Api);
	R = P->Ro+Ri;

	// TXPo.T is actually the outlet refrigerant temperature
	// however TXPo.T = TXPi.T in the two phase region and
	// TXPo.T ~ to TXPi.T in the superheated region if steps are small

	q = (TPi.T-TXPo.T)/mr*CmineCrossFlow(R,mr,ma,TXPo,TPi.T);
	if(errorLog.IsError()) {
		errorLog.Add("EvapTubeLnoDP (q)");
		return 0;
	}

	HPo->H-=q;

	// Determine mass of charge.  It is based on the inlet state.
	// The specific volume is recalculated so that a different model
	// can be used from the one used to calculate the pressure drop.

	TXPo = HPtoTXP(*HPo);
	if(errorLog.IsError()) {
		errorLog.Add("EvapTubeLnoDP (TXPo,2)");
		return 0;
	}

	v = VolumeALL(TXPo,Gr,P->Di);		// seperate flow model
	// v = PropertyTXP(VOL,TXPo);		// homogeneous flow model
	if(errorLog.IsError()) {
		errorLog.Add("EvapTubeLnoDP (v)");
		return 0;
	}

	m->V=P->Ls*P->Ax;
	m->m=m->V/v;

	// Calculate output air state
	hai = wair.h(TPi.T,TPi.P);
	if(errorLog.IsError()) {
		errorLog.Add("EvapTubeLnoDP (hai)");
		return 0;
	}

	WHo->H=hai-q*mr/ma;

	WHo->W = wair.HumidityRatio(TPi.T,TPi.P);
	if(errorLog.IsError()) {
		errorLog.Add("EvapTubeLnoDP (WHo)");
		return 0;
	}

	return 0;
}

/********************************************************************
Used to solve for the change of state considering the energy and
momentum conservation equations for compressible flow.
********************************************************************/
static double CompLossResid(double *X,double *Y,void *Params)
{
	HP* HPi = (HP*)X;
	ELRP* P = (ELRP*)Params;

	// inlet specific volume and entropy
	TXP TXPi = HPtoTXP(*HPi);
	if(errorLog.IsError()) {
		errorLog.Add("CompLossResid");
		return 0;
	}

	/* Calculate heat transfered per unit mass of refrigerant. */
	double hi = P->hRefAdj*ConvCoeffInside(TXPi,P->G,P->D);
	if(errorLog.IsError()) {
		errorLog.Add("CompLossResid");
		return 0;
	}
	double Ri = 1/(hi*P->Api);
	double R = P->Ro+Ri;
	P->q = (P->TPi.T-TXPi.T)/P->mr*CmineCrossFlow_dry(R,P->mr,P->ma,TXPi,P->TPi.T);
	if(errorLog.IsError()) {
		errorLog.Add("CompLossResid");
		return 0;
	}

	/* specific volume */
	double vi=VolumeALL(TXPi,P->G,P->D);	         // seperate flow model
//	double vi = PropertyTXP(VOL,TXPi);			// homogeneous flow model
	if(errorLog.IsError()) {
		errorLog.Add("CompLossResid");
		return 0;
	}

	// inlet and outlet stagnation enthalpy
	double h0o = P->HPo.H + 0.5*pow(P->vo*P->G,2);	   	// with kinetic energy terms
	double h0i = HPi->H   + 0.5*pow(vi*P->G,2);

	// momentum conservation
	double t1 = pow(P->G,2)*(vi+P->vo)*P->K/4;
	double t2 = 1000*(P->HPo.P - HPi->P);
	double t3 = pow(P->G,2)*(P->vo-vi);

	// stagnation enthalpy residual
	Y[0]=(h0o-h0i-P->q)/h0o;
	Y[1]=(t1+t2+t3)/t1;

	return 0;
}

/********************************************************************
Provides constraints used in solving for the outlet state of the
manifold.  The numerical solver updates guesses of enthalpy and
pressure.  This functions states a new guess to ensure that is
within the property tables.
********************************************************************/
static int HPLimitsConst(double *X,double* /*Y*/,void* /*Params*/)
{
	HP* HPo = (HP*)X;

	if(HPo->P<PMINth || HPo->P>PMAXth) return 1;

	double h = PropertyTXP(ENTH,toTXP(TMIN,0,HPo->P));
	if(errorLog.IsError()) {
		errorLog.Add("HPLimitsConst");
		return 0;
	}

	if(HPo->H<h) return 1;

	h = PropertyTXP(ENTH,toTXP(TMAX,1,HPo->P));
	if(errorLog.IsError()) {
		errorLog.Add("HPLimitsConst");
		return 0;
	}

	if(HPo->H>h) return 1;

	return 0;
}

/*********************************************************************
Evaporator tube segment model.
Inputs:
	Gr = refrigerant mass flux (kg/s/m^2)
	HPi = refrigerant inlet state (h,P)
	Ga = air mass flux (kg/s/m^2)
	Tai = air inlet temperature (C)
Outputs:
	HPi = refrigerant outlet state (h,P)
	hao = air outlet enthalpy (J/kg)
	m = mass of charge in return bend (kg)
*********************************************************************/
static int EvapTubeL(double Gr,HP *HPo,double Ga,TP TPi,WH *WHo,MASS *m,ETdim *P)
{
	double v,hai,ma,mr;
	TXP TXPo;
	ELRP Q;

	/* convert inlet state into TXP format */
	TXPo = HPtoTXP(*HPo);
	if(errorLog.IsError()) {
		errorLog.Add("EvapTubeL");
		return 0;
	}

	// Q.vo = VolumeALL(TXPo,Gr,P->Di);	/* seperate flow model */
	Q.vo = PropertyTXP(VOL,TXPo);			/* homogeneous flow model */
	if(errorLog.IsError()) {
		errorLog.Add("EvapTubeL");
		return 0;
	}

	/* Mass flow rates of air and refrigerant. */
	ma = Ga*P->Aflow;
	mr = Gr*P->Ax;

	Q.K = LossCoeffStraightPipe(TXPo,Gr,P->Ls,P->Di);
	if(errorLog.IsError()) {
		errorLog.Add("EvapTubeL");
		return 0;
	}

	Q.D=P->Di;
	Q.HPo=*HPo;
	Q.G=Gr;
	Q.hRefAdj = P->hRefAdj;
	Q.Api = P->Api;
	Q.Ro = P->Ro;
	Q.TPi = TPi;
	Q.ma = ma;
	Q.mr = mr;

	/* Solve for outlet state. */
	FindZero2DConst((double*)HPo,CompLossResid,HPLimitsConst,1e-6,&Q);
	if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"G=%lf HPo=(%lf,%lf)\n",Gr,HPo->H,HPo->P);
		errorLog.Add("EvapTubeL",str);
		return 0;
	}

	// Determine mass of charge.  It is based on the inlet state.
	//	The specific volume is recalculated so that a different model
	// can be used from the one used to calculate the pressure drop.

	TXPo = HPtoTXP(*HPo);
	if(errorLog.IsError()) {
		errorLog.Add("EvapTubeL");
		return 0;
	}

	v = VolumeALL(TXPo,Gr,P->Di);		/* seperate flow model */
	// v = PropertyTXP(VOL,TXPo);		/* homogeneous flow model */
	if(errorLog.IsError()) {
		errorLog.Add("EvapTubeL");
		return 0;
	}

	m->V=P->Ls*P->Ax;
	m->m=m->V/v;

	/* Calculate output air state. */
	hai = wair.h(TPi.T,TPi.P);
	WHo->H = hai-Q.q*mr/ma ;
	WHo->W = wair.HumidityRatio(TPi.T,TPi.P);

	return 0;
}
#endif

/********************************************************************
Models an entire evaporator tube.
Gr = refrigerant mass flux (kg/s/m^2)
Ga = air mass flux (kg/s/m^2)
********************************************************************/
int EvapTube(double Gr,//refrigerant mass flux
			 HP *HPo,//refrigerant outlet and inlet state
			 double Ga,//air mass flux
			 TP TPi,// air inlet state
			 WH *WHo,//air airlet state
			 MASS *m,// mass and inner volume of the evaporator tube
			 ETdim *P)//evaporator struct
{
	WH WHl;
	int i;
	double ho,phi;
	MASS mi;

	/* Calculate air side resistance */
	ho = P->hAirAdj*ConvCoeffAir_EVA(TPi,Ga,P);//B.S., the airside dry heat transfer coefficient can be got for different fin types
	if(errorLog.IsError()) {
		errorLog.Add("EvapTube", "get ho");
		return 0;
	}

	phi = FinEffect_Schmidt(ho,233,P->th,P->y,P->Do);//B.S. calculate the fin efficiency with the Schmidt equation
	
	P->Ro = 1/(ho*(P->Apo+phi*P->Af));
	P->ho=ho;

	m->m=0;
	m->V=0;
	WHo->W=0;
	WHo->H=0;

	for(i=0;i<P->NSeg;i++) {

		#ifdef EVAPORATOR_PRESSURE_DROP
			EvapTubeL_Rev(Gr,HPo,Ga,TPi,&WHl,&mi,P);//B.S., new code accounts for pressure drop and wet condition
		#else
			EvapTubeLnoDP(Gr,HPo,Ga,TPi,&WHl,&mi,P);
		#endif

		//printf("evap seg outlet seg=%d/%d, H=%.1lf, P=%.1lf\n",i+1,P->NSeg,HPo->H,HPo->P);

		if(errorLog.IsError()) {
			char msg[256];
			sprintf(msg,"EvapTube (segment %d/%d)",i+1,P->NSeg);
			errorLog.Add(msg);
			return 0;
		}

		m->m+=mi.m;
		m->V+=mi.V;
		WHo->W+=WHl.W;
		WHo->H+=WHl.H;
	}

	WHo->W/=(double)P->NSeg; WHo->H/=(double)P->NSeg;

	return 0;
}


//B.S.--------------------------------------------------
//*******************************************************************************
//*******************************************************************************
//*******************************************************************************
//the follow subroutines are for correlating the moving boundary and lumped heat exchanger models
void Build_Eva_Lumped(ETdim* P)
{
	const double cp_a= wair.Cp(P->TPi.T,P->TPi.P);//air specific heat
	if(errorLog.IsError()) {
		errorLog.Add("Build_Eva", "cp_a");
		return;
	}
	
	TXP TXP_prop = {0,0,0};
	TXP TXP8;//outlet state of the evaporator
	TXP8= HPtoTXP(P->HP_out);
	if(errorLog.IsError()) {
		errorLog.Add("Build_Eva", "TXP4");
		return;
	}

	TXP_prop.P=TXP8.P;
	TXP_prop.X=1;

	const double Tsat8 = PropertyTXPth(TSAT,TXP_prop);//reftpvthP.Tsat(TXP4.P);//calculate the vapor specific heat with the outlet state of the evaporator
	if(errorLog.IsError()) {
		errorLog.Add("Build_Eva", "Tsat8");
		return;
	}
	
	
	//lumped model-----------------------------
	const double Q_TOT = P->mr*(P->HP_out.H-P->HP_in.H);//overall heat transfer amount

	P->DP_TOT = (P->HP_in.P-P->HP_out.P);//overall pressure drop
	const double epsilonT_dry = Q_TOT/(P->ma_TOT*cp_a)/(P->TPi.T-Tsat8);//dry heat transfer effectiveness
	const double NTUT_dry = -log(1-epsilonT_dry);//dry heat transfer NTU
	P->UA_TOT = NTUT_dry*P->ma_TOT*cp_a;//dry heat transfer overall conductance

	/*const double H_SAT_E=wair.h(Tsat4,0.999);//saturated air enthalpy at the outlet refrigerant temperature
	if(errorLog.IsError()) {
	errorLog.Add("Build_Eva","H_SAT_E");
	return;
	}
	const double epsilonT_wet = Q_TOT/P->ma_TOT/(H_A_I-H_SAT_E);//wet heat transfer effectiveness
	const double NTUT_wet = -log(1-epsilonT_wet);
	P->UAw_TOT = NTUT_wet*P->ma_TOT;//wet heat transfer overall conductance*/

	//lumped model ended
	return;
}

void Build_Eva_MovBound(ETdim* P)
{
	static ETdim P1;//for keeping the constructs calculated from last run

	TXP TXP_prop={0,0,0};

	//average density, pressure drop, cross-sectional area of each phase-------
	if(P->count1>0)
	{
	P->HP_TP1.H=P->HP_TP1.H/P->count1;//get the average enthalpy of the first two-phase segment, counted from the entrance of the evaporator
	P->HP_TP1.P=P->HP_TP1.P/P->count1;//pressure of the first two-phase segment in the evaporator, counted from the entrance of the evaporator
	}
	else
	{
	P->HP_TP1 = P->HP_in;//otherwise. the first segment of two-phase heat transfe is just at the entrance of the evaporator
	}

	if(P->count2>0)
	{
	P->HP_TP2.H=P->HP_TP2.H/P->count2;//get the average enthalpy of the last two-phase segment, counted from the entrance of the evaporator
	P->HP_TP2.P=P->HP_TP2.P/P->count2;//pressure of the last two-phase segment in the evaporator, counted from the entrance of the evaporator
	}
	else
	{
	P->HP_TP2 = P->HP_out;//otherwise. the last segment of two-phase heat transfe is just at the exit of the evaporator
	}
	
	P->V_TOT = P->V_Liq+P->V_Vap+P->V_TP;//total inner volume
	P->L_TOT = P->VapL+P->TPL+P->LiqL;//total tube length
	P->A_TOT=P->V_TOT/P->L_TOT;//total heat transfer outside surface

	
	//-----------------------------the followings are average density, pressure drop, cross-sectional area of each phase


	//prepare heat transfer calculation---------------
	P->Ga_meanL = P->ma_TOT/P->L_TOT;//average air flow rate for per tube length

	const double ma_v = P->Ga_meanL * P->VapL;//air flow rate across the vapor-phase region
	const double ma_tp = P->Ga_meanL * P->TPL;//air flow rate across the two-phase region
	const double ma_l = P->Ga_meanL * P->LiqL;//air flow rate across the liquid-phase region
	
	const double cp_a= wair.Cp(P->TPi.T,P->TPi.P);//air specific heat
	if(errorLog.IsError()) {
		errorLog.Add("Build_Eva", "cp_a");
		return;
	}	
	
	TXP TXP1;//inlet state of the evaporator
	TXP1= HPtoTXP(P->HP_in);
	if(errorLog.IsError()) {
		errorLog.Add("Build_Eva", "TXP1");
		return;
	}
	TXP_prop.P=TXP1.P;
	TXP_prop.X=1;


	const double Tsat1 = PropertyTXPth(TSAT,TXP_prop);//reftpvthP.Tsat(TXP1.P);//saturation temperature corresponding to the inlet state of the evaporator
	if(errorLog.IsError()) {
		errorLog.Add("Build_Eva", "Tsat1");
		return;
	}
	
	TXP_prop.P=TXP1.P;
	TXP_prop.T=Tsat1;
	TXP_prop.X=0;

	const double cp_rl = PropertyTXPtr(SPEC,TXP_prop);//refsctrPT.Cp(TXP1.P,Tsat1);//refrigerant liquid specific heat corresponding to the inlet state of the evaporator
	if(errorLog.IsError()) {
		errorLog.Add("Build_Eva", "cp_rl");
		return;
	}


	TXP TXP2;//the first two-phase state point counted from the entrance of the evaporator

	TXP2= HPtoTXP(P->HP_TP1);
	if(errorLog.IsError()) {
		errorLog.Add("Build_Cond", "TXP3");
		return;
	}


	TXP TXP4;//outlet state of the evaporator
	TXP4= HPtoTXP(P->HP_out);
	if(errorLog.IsError()) {
		errorLog.Add("Build_Eva", "TXP4");
		return;
	}

	TXP TXP3;//the last two-phase state point counted from the entrance ofthe evaporator

	TXP3= HPtoTXP(P->HP_TP2);
	if(errorLog.IsError()) {
		errorLog.Add("Build_Eva", "TXP3");
		return;
	}

	TXP_prop.P=TXP3.P;
	TXP_prop.X=1;

	const double Tsat3 = PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(TXP3.P);//saturation temperature corresponding to the last two-phase segment of the evaporator
	if(errorLog.IsError()) {
		errorLog.Add("Eva_Cond", "Tsat3");
		return;
	}

	TXP_prop.P=P->HP_TP2.P;
	TXP_prop.T=Tsat3;
	TXP_prop.X=1;
	const double Hsat2= PropertyTXPth(ENTH,TXP_prop);//refshthPT.h(P->HP_TP2.P,Tsat3);//saturation enthalpy corresponding to the last two-phase segment of the evaporator
	if(errorLog.IsError()) {
		errorLog.Add("Eva_Cond","Hsat2");
		return;
	}
	if(P->TPL>0) P->H2_residual = Hsat2 - P->HP_TP2.H;//prepared for maintaining the consistency in the moving boundary model

	TXP_prop.P=TXP4.P;
	TXP_prop.X=1;

	const double Tsat4 = PropertyTXPth(TSAT,TXP_prop);//reftpvthP.Tsat(TXP4.P);//calculate the vapor specific heat with the outlet state of the evaporator
	if(errorLog.IsError()) {
		errorLog.Add("Build_Eva", "Tsat3");
		return;
	}
	
	TXP_prop.P=TXP4.P;
	TXP_prop.T=Tsat4;
	TXP_prop.X=1;

	const double cp_rv = PropertyTXPtr(SPEC,TXP_prop);//refshtrPT.Cp(TXP4.P,Tsat4);//calculate the vapor specific heat with the outlet state of the evaporator
	if(errorLog.IsError()) {
		errorLog.Add("Build_Eva", "cp_rv");
		return;
	}

	//---------------------------Prepare heat transfer calculation
	
	
	//-----------------------------------------heat transfer of the vapor phase
	double cmin_v=0;
	double Cr_v = 0;
	const double Q_v = P->mr*(P->HP_out.H-P->HP_TP2.H);//actual gas phase heat transfer
	if(ma_v>0&&Q_v>0)
	{
	
	Cr_v =(P->mr*cp_rv)/(ma_v*cp_a);//calculate Cr_v, cmin_v
	if(ma_v*cp_a>P->mr*cp_rv) {cmin_v = P->mr*cp_rv; Cr_v = Cr_v;}
	else {cmin_v = ma_v*cp_a;Cr_v = 1/Cr_v;}

	const double Qv_max = cmin_v*(P->TPi.T -/*Tsat4*/TXP4.T);//based on the exit gas temperature, this is for avoiding the additional iteration in the later calculation
	const double e_v = Q_v/Qv_max;//actual effectiveness
	const double NTU_v = P->UA_Vap/cmin_v;
	double e_sup =0;	
	if(ma_v*cp_a>P->mr*cp_rv) {//cmax unmixed
	const double l = 1-exp(-Cr_v*NTU_v);
	e_sup = 1 - exp(-l/Cr_v);//theoretical effectiveness
	}
	else {//cmin unmixed
	const double l = -Cr_v*(1-exp(-NTU_v));
	e_sup = (1-exp(l))/Cr_v;//theoretical effectiveness
	}

	P->r_v = e_v/e_sup;//parameter for adjusting the effectiveness
	P->U_Vap = P->UA_Vap/(P->VapL);//get the average heat transfer coefficent
	P->rho_Vap = P->m_Vap/P->V_Vap;//averge gas-phase density
	P->A_Vap=P->V_Vap/P->VapL;//average gas-phase heat transfer area
	P->DP_Vap = (P->HP_TP2.P-P->HP_out.P)/P->V_Vap;//average gas-phase pressure drop gradient
	}
	else//use the parameter stored at last time
	{
	P->r_v = P1.r_v;
	P->U_Vap = P1.U_Vap;
	P->rho_Vap = P1.rho_Vap;
	P->A_Vap=P1.A_Vap;
	P->DP_Vap = P1.DP_Vap;
	}
	//--------------------------------------------------------heat transfer of vapor phase



	//heat transfer of liquid phase--------------------
	double cmin_l=0;
	double Cr_l = 0;
		
	if(ma_l>0)
	{
	const double Q_l = P->mr*(P->HP_TP1.H-P->HP_in.H);//liquid phase heat transfer amount

	Cr_l = (P->mr*cp_rl)/(ma_l*cp_a);//calculate Cr_l and cmin_l with the theoritical method
	if(ma_l*cp_a>P->mr*cp_rl) {cmin_l = P->mr*cp_rl; Cr_l=Cr_l;}
	else {
		cmin_l = ma_l*cp_a; Cr_l =1/Cr_l;}

	const double Ql_max = cmin_l*(P->TPi.T-TXP2.T);//calculate the effectiveness with the theoritical method
	const double e_l = Q_l/Ql_max;//actual liquid-phase heat transfer effectiveness


	const double NTU_l=P->UA_Liq/cmin_l;
	
	double e_sub =0;//theoritical two-phase heat transfer effectiveness
	
	if(ma_l*cp_a>P->mr*cp_rl) {//cmax unmixed 
	const double l = 1-exp(-Cr_l*NTU_l);
	e_sub = 1 - exp(-l/Cr_l);
	}
	else {//cmin unmixed
	const double l = -Cr_l*(1-exp(-NTU_l));
	e_sub = (1-exp(l))/Cr_l;
	}

	P->r_l = e_l/e_sub;//parameter for adjusting the theoritical heat transfer effectiveness
	P->U_Liq = P->UA_Liq/(P->LiqL);//get the averge heat transfer conductance
	P->rho_Liq = P->m_Liq/P->V_Liq;//averge liquid phase density
	P->A_Liq=P->V_Liq/P->LiqL;//averge liquid-phase heat transfer surface area
	P->DP_Liq = (P->HP_in.P-P->HP_TP1.P)/P->V_Liq;//average liquid-phase pressure drop gradient
	}
	else//use the parameters stored last run
	{
	P->r_l = P->r_v;//use the parameters of gas for the liquid-phase
	P->U_Liq = P->U_Vap;
	}
	//--------------------------------------heat transfer of liquid-phase
	

	//heat transfer of two-phase------------------------
	const double H_A_I=wair.h(P->TPi.T,P->TPi.P);//inlet air enthalpy
	if(errorLog.IsError()) {
		errorLog.Add("Build_eva","H_A_I");
		return ;
		}

	const double H_SAT=wair.h(TXP3.T,0.999);//air saturation enthalpy at the refrigerant temperature of the last two-phase segment, counted from the entrance of the evaporator
		if(errorLog.IsError()) {
		errorLog.Add("Build_Eva","H_SAT");
		return;
		}

	if(P->TPL>0)
	{
	const double Q_tp = P->mr*(P->HP_TP2.H-P->HP_TP1.H);//two-phase heat transfer amount
	P->Uw_TP = P->UAw_TP/P->L_wet;//average wet two-phase heat transfer conductance
	P->U_TP = P->UA_TP/(P->TPL-P->L_wet);//averge dry two-phase heat transfer conductance
	
	//dry heat transfer
	double cmintp_dry = ma_tp*cp_a;	 

	double UA_TP = P->U_TP*P->TPL;
	const double NTU_tp = UA_TP/cmintp_dry;
	const double e_dry = 1-exp(-NTU_tp);
	const double Qtp_dry = e_dry*cmintp_dry*(P->TPi.T-TXP3.T);

	//wet heat transfer
	double UAw_TP = P->Uw_TP*P->TPL;
	const double NTUw_tp = UAw_TP/ma_tp;
	const double e_wet = 1-exp(-NTUw_tp);
	

	const double Qtp_wet = e_wet*ma_tp*(H_A_I- H_SAT);
	
	if(Qtp_wet>Qtp_dry) 
	{ P->r_tp= Q_tp/Qtp_wet;}
	else {
		P->r_tp= Q_tp/Qtp_dry;}
	
	P->rho_TP = P->m_TP/P->V_TP;//average two-phase density
	P->A_TP=P->V_TP/P->TPL;//averge two-phase heat transfer surface area
	P->DP_TP = (P->HP_TP1.P-P->HP_TP2.P)/P->V_TP;//average two-phase pressure drop gradient
	}
	else//get the parameters from last run
	{
	P->rho_TP = P1.rho_TP;
	P->A_TP=P1.A_TP;	
	P->DP_TP = P1.DP_TP;
	P->r_tp= P1.r_tp;
	P->Uw_TP = P1.Uw_TP;//wet
	P->U_TP = P1.U_TP;//dry
	}
	//-----------------------------heat transfer of two-phase
	

	P1=*P;//store the information from this run

	return;	
}

//******************************************************************************
//******************************************************************************
//******************************************************************************
//the following subrountines are for moving boundary models and lumped models analysis. they are not used now.

double Evaporator_Simple(
			double mr,//refrigerant mass flow rate
			   HP HPo,//refrigerant outlet state
			   TP TPi,//inlet air state
			   HP *HPi,//inlet refrigerant state of the segment
			   MASS *Sm,//mass and volume in the evaporator
			   ETdim *P)//adjustment parameters
{
	//gas phase
	Phase_eva Gas;
	Gas.P = *P;
	Gas.mr = mr;
	Gas.HPo = HPo;
	Gas.TPi = TPi;

	double residual=0;

	residual = GasCal_Eva(P->L_TOT, &Gas);
	if(residual >0)
	{
	*HPi = Gas.HPi;
	Sm->m = Gas.m;
	return 0;
	}

	Zbrent(P->L_TOT,0.0000001,GasCal_Eva,1e-7,&Gas);

	if(errorLog.IsError()) {
		errorLog.Add("Evaporator_Simple", "(Zbrent1)");
		return 0;
	};
 
	Phase_eva TP;
	TP.P = *P;
	TP.mr = mr;
	TP.HPo = Gas.HPi;
	TP.TPi = TPi;
	
	
	const double L_left = P->L_TOT - Gas.L;
	TPCal_Eva(L_left,&TP);

	if(errorLog.IsError()) {
		errorLog.Add("Evaporator_Simple", "TPCal");
		return 0;
	};

	*HPi = TP.HPi;
	Sm->m = Gas.m+TP.m;

	return 0;
}

double GasCal_Eva(double L, void *Params)
{
	Phase_eva* D=(Phase_eva*)Params;
	TXP TXP_prop={0,0,0};

	double cmin_v=0;
	double Cr_v = 0;
	const double ma_v = L*D->P.Ga_meanL;
	
	const double cp_a= wair.Cp(D->TPi.T,D->TPi.P);
	if(errorLog.IsError()) {
		errorLog.Add("GasCal_Eva", "cp_a");
		return 0;
	}

	TXP TXP4;
	TXP4= HPtoTXP(D->HPo);
	if(errorLog.IsError()) {
		errorLog.Add("GasCal_Eva", "TXP4");
		return 0;
	}

	if(TXP4.X<1)
	{
	D->L =0;
	D->HPi=D->HPo;
	D->m =0;
	return 0;
	}

	TXP_prop.P=D->HPo.P;
	TXP_prop.X=1;
	
	const double Tsat4 = PropertyTXPth(TSAT,TXP_prop);//reftpvthP.Tsat(D->HPo.P);
	if(errorLog.IsError()) {
		errorLog.Add("GasCal_Eva", "Tsat4");
		return 0;
	}
	
	TXP_prop.P=TXP4.P;
	TXP_prop.T=Tsat4;
	TXP_prop.X=1;

	const double cp_rv = PropertyTXPtr(SPEC,TXP_prop);//refshtrPT.Cp(TXP4.P,Tsat4);
	if(errorLog.IsError()) {
		errorLog.Add("GasCal_Eva", "cp_rv");
		return 0;
	}

	Cr_v = (D->mr*cp_rv)/(ma_v*cp_a);

	if(ma_v*cp_a>D->mr*cp_rv) {cmin_v = D->mr*cp_rv; Cr_v = Cr_v;}
	else {cmin_v = ma_v*cp_a;Cr_v = 1/Cr_v;}

	double Qv_max = cmin_v*(D->TPi.T- TXP4.T/*Tsat4*/);//attention for avoiding additional iteration
	
	if(Qv_max<=0) 
	{
	Qv_max = 0;}

	const double NTU_v = D->P.U_Vap*L/cmin_v;
	
	double e_sup =0;
	if(ma_v*cp_a>D->mr*cp_rv) {
	const double l = 1-exp(-Cr_v*NTU_v);
	e_sup = 1 - exp(-l/Cr_v);
	}
	else {
	const double l = -Cr_v*(1-exp(-NTU_v));
	e_sup = (1-exp(l))/Cr_v;
	}

	e_sup =e_sup*D->P.r_v;
	
	const double Q_gas = Qv_max*e_sup;

	D->HPi.H = D->HPo.H - Q_gas/D->mr;

	const double V_Vap = L*D->P.A_Vap;
	const double DP_Vap = V_Vap*D->P.DP_Vap;

	D->HPi.P = DP_Vap+D->HPo.P;
	D->m = D->P.rho_Vap*V_Vap;
	D->L = L;

	TXP_prop.P=D->HPi.P;
	TXP_prop.X=1;

	const double Tsat3 = PropertyTXPth(TSAT,TXP_prop);//reftpvthP.Tsat(D->HPi.P);
	if(errorLog.IsError()) {
		errorLog.Add("GasCal_Eva", "Tsat3");
		return 0;
	}

	TXP_prop.P=D->HPi.P;
	TXP_prop.T=Tsat3;
	TXP_prop.X=1;

	const double Hsat= PropertyTXPth(ENTH,TXP_prop) - D->P.H2_residual;;//refshthPT.h(D->HPi.P,Tsat3) - D->P.H2_residual;
	if(errorLog.IsError()) {
		errorLog.Add("GasCal_Eva","Hsat");
		return 0;
	}

	return (D->HPi.H- Hsat)/1e6;//fabs(Hsat + D->HPi.H);

}

double TPCal_Eva(double L, void *Params)
{
	Phase_eva* D=(Phase_eva*)Params;
	
	const double ma_tp = L*D->P.Ga_meanL;

	const double cp_a= wair.Cp(D->TPi.T,D->TPi.P);
	if(errorLog.IsError()) {
		errorLog.Add("TPCal_Eva", "cp_a");
		return 0;
	}

	TXP TXP3;
	TXP3= HPtoTXP(D->HPo);
	if(errorLog.IsError()) {
		errorLog.Add("TPCal_Eva", "TXP3");
		return 0;
	}
	
	//dry heat transfer
	double cmintp_dry = ma_tp*cp_a;	

	const double UA_TP = D->P.U_TP*L;
	const double NTU_tp = UA_TP/cmintp_dry;
	const double e_dry = 1-exp(-NTU_tp);

	const double Qtp_dry = e_dry*cmintp_dry*(D->P.TPi.T-TXP3.T);

	//wet heat transfer
	const double UAw_TP = D->P.Uw_TP*L;
	const double NTUw_tp = UAw_TP/ma_tp;
	const double e_wet = 1-exp(-NTUw_tp);
	
	const double H_A_I=wair.h(D->TPi.T,D->TPi.P);//inlet air enthalpy
	if(errorLog.IsError()) {
		errorLog.Add("TPCal_eva","H_A_I");
		return 0;
		}

	const double H_SAT=wair.h(TXP3.T,0.999);
		if(errorLog.IsError()) {
		errorLog.Add("TPCal_Eva","H_SAT");
		return 0;
		}

	const double Qtp_wet = e_wet*ma_tp*(H_A_I- H_SAT);

	double Q_tp=0;

	if(Qtp_wet>Qtp_dry) {Q_tp=Qtp_wet;}
	else{Q_tp=Qtp_dry;}

	Q_tp = D->P.r_tp*Q_tp;

	D->HPi.H = D->HPo.H - Q_tp/D->mr;

	const double V_TP = L*D->P.A_TP;
	const double DP_TP = D->P.DP_TP*V_TP;

	D->HPi.P = D->HPo.P + DP_TP;

	D->m = V_TP*D->P.rho_TP;

	D->L = L;

	return 0;
}

//Lumped model
double Evaporator_Lumped(
				double mr,//refrigerant mass flow rate
			   HP HPo,//refrigerant outlet state
			   TP TPi,//inlet air state
			   HP *HPi,//inlet refrigerant state of the segment
			   MASS *Sm,//mass and volume in the evaporator
			   ETdim *P)//adjustment parameters
{

	const double cp_a= wair.Cp(TPi.T,TPi.P);
	if(errorLog.IsError()) {
		errorLog.Add("Evaporator_Lumped", "cp_a");
		return 0;
	}

	TXP TXP_prop={0,0,0};

	TXP_prop.P=HPo.P;
	TXP_prop.X=1;

	const double Tsat = PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(HPo.P);
	if(errorLog.IsError()) {
		errorLog.Add("Evaporator_Lumped", "Tsat");
		return 0;
	}
	
/*	const double H_SAT=wair.h(Tsat,0.999);
	if(errorLog.IsError()) {
	errorLog.Add("Evaporator_Lumped","H_SAT");
	return 0;
	}
	
	const double H_Air = wair.h(TPi.T,TPi.P);
	if(errorLog.IsError()) {
	errorLog.Add("Evaporator_Lumped","H_Air");
	return 0;
	}*/

	//dry heat transfer-----------------------------
	const double NTU = P->UA_TOT/(P->ma_TOT*cp_a);
	const double Q_dry = (1-exp(-NTU))*P->ma_TOT*cp_a*(TPi.T-Tsat);

	HPi->P = HPo.P + P->DP_TOT;
	double Q_e = 0;

	//wet heat transfer-------------------------------
	/*const double NTU_w = P->UAw_TOT/(P->ma_TOT);
	const double Q_wet = (1-exp(-NTU_w))*P->ma_TOT*(H_Air-H_SAT);

	

	if(Q_wet>Q_dry) Q_e = Q_wet;
	else */
	Q_e = Q_dry;
	
	HPi->H = HPo.H - Q_e/mr;

	Sm->m=0; Sm->V=0;
	return 0;
}