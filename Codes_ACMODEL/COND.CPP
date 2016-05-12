#include "stdafx.h"

#include <stdio.h>
#include <math.h>
#include "r22.h"
#include "ref.h"
#include "mass.h"
#include "cond.h"
#include "air.h"
#include "numeric.h"
#include "num2d.h"
#include "corr.h"
#include "volume.h"
#include "pressure.h"
#include "cmine.h"
#include "getdata.h"
#include "errorlog.h"
#include "spline1d.h"
#include "spline2d.h"

//B.S.----------------------------
#include "outelem.h"//B.S.
extern OutputDataElement gOutData;//B.S.
//#include <fstream.h>
//#include <iostream.h>
//-------------------------------B.S.

extern ErrorLog errorLog;

extern AirProp air;

int condInit=1;

#define CONDENSER_PRESSURE_DROP
#define DEBUG_COND_202_AIRFLOW

struct ACRP {
	HP HPi;
	double Go,Gi,Di,vi,si;
};
struct CLRP {
	HP HPi;
	double q,K,G,D,vi;
};
struct Cond0P {
	double mr,Tai,Tro,*Prms,Ga,Tsc;
	HP HPi;
	const char *filename;
};

// Local functions
#ifdef CONDENSER_PRESSURE_DROP
//B.S. the following is new function for calculating heat transfer and pressure drop in a tube unit, considering the pressure drop
#else
static void CondTubeLnoDP(double,HP*,double,double,double*,MASS*,CGP*); 
#endif

//B.S.--------------------------
// the original function for calculating the pressure drop and heat transfer in a tube unit is still defined here in case.
#ifdef CONDENSER_PRESSURE_DROP
static void CondTubeL(double,HP*,double,double,double*,MASS*,CGP*);// this is the original function
#endif
//------------------------------B.S.

//B.S.-----------------------------------
struct Phase_Cond
{
	CGP P;
	HP HPi;
	double mr;
	double Tai;
	HP HPo;
	double L;
	double m;
};

double GasCal_Cond(double, void *);//B.S., new
double TPCal_Cond(double, void *);//B.S., new
double LiqCal_Cond(double, void *);//B.S., new
//--------------------------------------B.S.

static double AreaChangeResid(double*,double*,void*);
static double CompLossResid(double*,double*,void*);
static int HPLimitsConst(double*,double*,void*);
static double Cond0Func(double,void*);

//B.S. changes the following function for input more geometrical parameters of micro-fin tube and specially configured fin 
/**********************************************************************
Model of an air cooled condenser coil.  The model is comprised of
small pieces of single finned tubes put together to form rows.
Rows of finned tubes are put together using elbows and manifolds
to make an entire coil. This function defines the coil's circuit
pattern.  Model parameters are stored in the file "cond.dat".

Inputs:
	mr - refrigerant mass flow rate (kg/s).
	HPi - inlet refrigerant thermodynamic state (J/kg,kPa).
	Tai - inlet air temperature (C).
	Ga - air mass flux (kg/s*m^2).
		If Ga<0, then nominal mass flux is used from data file.
	Prms - array of three multiplicative correlation correction factors.
		Prms[0] - corrects air side convection coefficient.
		Prms[1] - corrects refrigerant side convection coefficient.
		Prms[2] - corrects refrigerant side pressure drop correlation.
		Prms[3] - fouling factor (1=no fouling,0=completely fouled).
		Set these parameters to 1 (Prms={1,1,1,1}) for no corrections.
		NSeg - number of segments (-1 default - use value from file)

Outputs:
	HPo - outlet refrigerant thermodynamic state (J/kg,kPa).
	Tao - outlet air temperature (C).
	m - mass of charge in condenser (kg) and volume of charge (m^3).
***********************************************************************/
int Condenser(const char *filename,
			  double mr,
			  HP HPi,
			  double Tai,
			  double GaI,
			  HP *HPo,
			  double *Tao,
			  MASS *m,
			  CGP *Cond_struc,//new
			  double *Prms,
			  int NSeg)
{
	static int type;
	static double Ga;
	static double GaNom,hAirAdj,hRefAdj,PRefAdj,FoulFac, hRefAdj_Sub;
	static CGP P;
	static int endSign;

	if(condInit || NSeg>0) {
		FILE *fp = fopen(filename,"r");
		type = PosComponentData(fp,CONDENSER,1);
		if(!type) {
			errorLog.Add("Condenser","Parameters not found");
			return 0;
		}

		P.Di=GetDouble(fp);           // tube inside diameter (m)
		P.L=GetDouble(fp);            // tube length (m)
		P.xp=GetDouble(fp);           // pipe wall thickness (m)
		P.Df=GetDouble(fp);           // fin diameter (m)
		P.z=GetDouble(fp);            // B.S., originally it is space between fins (m), here it is used as fin pitch, which includes a fin thickness
		P.th=GetDouble(fp);           // fin thickness (m)
		P.vsp=GetDouble(fp);          // vertical tube spacing (m)
		P.P_l = GetDouble(fp);		// B.S. tube spacing along the airflow direction (m)
		P.NSeg=(int)GetDouble(fp);    // number of lengthwise segments (-)
		if(NSeg>0) P.NSeg = NSeg;
		P.Brad=GetDouble(fp);         // radius of return bend (m)
		GaNom=GetDouble(fp);          // nominal air mass flux (kg/s/m^2)

		P.Nbranchs=(int)GetDouble(fp);		// number of branchs in main section
		P.Nmaintubes=(int)GetDouble(fp);	// number of tubes in main section
		P.Nsubtubes=(int)GetDouble(fp);		// number of tubes in subsooled section
		P.Ndeep = (int)GetDouble(fp);		// B.S. number of high rows
		P.Frontal_A = GetDouble(fp);		// B.S. frontal area

		hAirAdj=GetDouble(fp);        // hAirAdj
		hRefAdj=GetDouble(fp);
		hRefAdj_Sub = GetDouble(fp);//B.S., it is a newly added adjustment factor for correcting the inaccurate ratio of two-phase heat transfer and single-phase heat transfer.
		PRefAdj=GetDouble(fp);
		FoulFac=GetDouble(fp);

		//B.S. ------------------------
		//new input parameter for condenser
		P.Microfin=(int)GetDouble(fp); //microfin type, 0=smooth tube, 1=helical, 2=cross-grooved, 3=herringbone
		P.finN = GetDouble(fp); //fin number in a micro-fin tube
		P.gama = GetDouble(fp);//fin apex angle in a micro-fin tube
		P.beta = GetDouble(fp);	//fin helix angle in a micro-fin tube
		P.finH = GetDouble(fp); //fin height in a micro-fin tube
		P.w_b= GetDouble(fp); //base width of a single fin, in the axial direction
		P.w_e= GetDouble(fp); //top width of a single fin, in the axial direction
		P.w_z= GetDouble(fp); //base distance between two neighboring fins, in the axial direction
		// airside fin type
		P.airFin = (int)GetDouble(fp);//1-plain, 2-corrugated, 3-slit, 4-louvered, 5-convex louvered, 6-smooth wavy, 7-spine
		P.sub1= GetDouble(fp);//fin substructures, for the possible use later, such as the waffle height of the wavy fin, louver number of the louverd fin
		P.sub2= GetDouble(fp);//fin substructures
		P.sub3= GetDouble(fp);//fin substructures
		P.sub4= GetDouble(fp);//fin substructures
		P.sub5= GetDouble(fp);//fin substructures
		P.K_T = GetDouble(fp);//400, this is the conductance factor of copper
		P.K_F= GetDouble(fp);//237, conductance factor of fin
		//------------------------------B.S.

		if(type==202) {
			P.cfmA=GetDouble(fp);
			P.cfmB=GetDouble(fp);
			P.cfmC=GetDouble(fp);
		}

		endSign = (int)GetDouble(fp);//-256

		fclose(fp);

		if(endSign!=-256)
		{
		errorLog.Add("Condenser","parameter input wrong");
		return -1;
		}

		condInit=0 ;
		// Calculate physical sizes needed in calculations
		// from input dimensions.
		const double pi = acos(-1.0);
		// Length of a tube segment.
		P.Ls=P.L/(double)P.NSeg;

		//B.S.-----------------------------------
		// tube inside geometry
		if(P.Microfin)
		{
		P.D_b = P.Di + 2*P.finH;//tube diameter at the base of the fin
		P.Do=P.D_b+2*P.xp; //Pipe outside diameter.
		P.D_m = (P.Di + P.D_b)/2; //mean diameter of the micro-fin tube
		P.P_H = P.finN*(P.w_z+2*pow((pow(P.finH,2.0)+pow((0.5*(P.w_b-P.w_e)),2.0)),0.5)+P.w_e);// the hydraulical circumference
		P.Acs = pi*pow((P.D_b/2.0),2.0)-P.finN*0.5*P.finH*(P.w_b+P.w_e);//cross area of the micro-fin tube, this is the actual cross-section area
		P.Dh_i = 4*P.Acs/P.P_H;//inside hydraulical diameter 
		P.Ax = pi/4.0*P.Di*P.Di;// Inside pipe cross sectional area, based on fin tips
		P.Api=pi*P.Di*P.Ls ;// Inside pipe surface area (one tube segment), based on fin tips. All the two-phase heat transfer are calculated with Cavallini correlations, which are based on the heat transfer surface area on the tip of the fins. 
		}
		else 
		{
		P.finH = 0;//without fin height
		P.D_b = P.Di;//tube diameter at the base of the fin
		P.Do=P.D_b+2*P.xp; //Pipe outside diameter.
		P.D_m = P.Di ; //averge diameter of the micro-fin tube
		P.Ax = pi/4.0*P.Di*P.Di;// Inside pipe cross sectional area
		P.Api=pi*P.Di*P.Ls ;// Inside pipe surface area (one tube segment)
		P.P_H = pi*P.D_b;// the hydraulical circumference
		P.Acs = P.Ax;//cross area of the micro-fin tube, this is the actual cross-section area
		P.Dh_i = P.Di;//inside hydraulical diameter
		}
		//-----------------------------------------B.S.
		
		// Volume inside tube segment.
		P.Vs=P.Acs*P.Ls;
		// Length of travel through bend.
		P.Blen=pi*P.Brad;
	
		
		//tube outside geometry

		//B.S. --------------------------------------------------
		//get the fin reference length, according to Schimidt equation
		const double M=P.vsp/2e0;
		const double R_O=P.Do/2;
		double L = 0;

		if (P.Ndeep>1) L=0.5e0*pow((pow(M,2e0)+pow(P.P_l,2e0)),0.5);
		else L = P.P_l/2.0;

		const double BETA=L/M;
		const double PHI=M/R_O;
		const double R_EFF=R_O*1.27e0*PHI*pow((BETA-0.3e0),0.5e0);
		P.Df = R_EFF*2;// confirm the reference fin outside diameter
		//-------------------------------------------------------B.S.

		// Distance from outside of pipe to edge of fin. 
	
		P.y=(P.Df-P.Do)/2 ;//B.S. correct the original term: P.y=(P.Df-P.Do)
		// Number of fins along one tube segment.
		P.N= P.Ls/(P.z);//B.S. correct P.Ls/(P.z+P.th); here P.z is inputted as the fin pith while not the fin spacing
		// Outside pipe surface area (one tube segment).
		P.Apo=pi*P.Do*(P.Ls-P.N*P.th);//B.S.

		// Fin wetted area (one tube segment).
		const double x = P.Df + P.th/2.0;
		P.Af=P.N*pi/2.0*(x*x-P.Do*P.Do);// the fin has double faces
		// Air flow area (one tube segment).
		P.Aflow=(P.vsp-P.Do)*(P.Ls-P.N*P.th) ;//minimum airflow cross-sectional area
		//printf("Cond Min Flow Area = %lf (m^2)\n",P.Aflow);

		// Hydrolic diameter
		P.Dh=4*P.Aflow/(P.Af+P.Apo)*P.vsp;//hydraulic diameter of the fin channel

	}

	/* Process scaling factors */
	P.hAirAdj=hAirAdj*Prms[0];
	P.hRefAdj=hRefAdj*Prms[1];
	P.PRefAdj=PRefAdj*Prms[2];
	//B.S. the following is for this is for adjusting the single-phase heat transfer
	P.hRefAdj_Sub=hRefAdj_Sub*Prms[9];//B.S. Prms[9] is for avoiding being confusing with other adjust parameters

	if(type==202) {
		// Date: Tue, 07 Mar 2000 09:36:46 -0500
		// To: Todd Rossi <rossi@mail.fielddiagnostics.com>
		// From: "Natascha S. Castro" <natascha@nist.gov>
		//
		// The polynomial that I extracted from experimental
		// data is : y = 16.677x2 - 7567.8x + 861160
		// I chose a polynomial fit.  It would be good to have the model read
		// coefficients from a file.  

		// Date: Thu, 23 Mar 2000 15:10:49 -0500
		// To: Todd Rossi <rossi@mail.fielddiagnostics.com>
		// From: "Natascha S. Castro" <natascha@nist.gov>
		//
		// y=cfm
		// x=psia

		// Date: Fri, 24 Mar 2000 09:41:39 -0500
		// To: Todd Rossi <rossi@mail.fielddiagnostics.com>
		// From: "Natascha S. Castro" <natascha@nist.gov>
      //
     	// The range for normal data is 220 psia to 245 psia, and the equation holds
		// for that range.  I think selecting extreme constant values is a good idea,
		// particularly because the data looks like a step, with a lower values level
		// and upper values level.  Based on the normal data, set discharge pressures
		// less than 225 psia and below to 2400 cfm and discharge pressures greater
		// than 245 to 3500 cfm.
		const double Pdis_kPa = HPi.P;
		const double Pdis_psia = Pdis_kPa/101.3*14.7;
		double af_cfm;
		if(Pdis_psia<225.0) {
			af_cfm = 2400.0;
		} else if(Pdis_psia>245.0) {
			af_cfm = 3500.0;
		} else {
			af_cfm = P.cfmA+Pdis_psia*(P.cfmB+Pdis_psia*P.cfmC);
		}
		const double af_m3s = af_cfm/60.0*(12.0*0.0254)*(12.0*0.0254)*(12.0*0.0254);
		const double af_kgs = af_m3s/air.v(Tai);
		const double Ntubes = P.Nsubtubes+P.Nbranchs*P.Nmaintubes;
		const double Aflow_seg = P.Aflow*P.NSeg;
		const double Aflow_tot = Aflow_seg*Ntubes;
		const double GaNom = af_kgs/Aflow_tot;

		Ga=GaNom*FoulFac*Prms[3];

		#ifdef DEBUG_COND_202_AIRFLOW
		{
		FILE*fp=fopen("DebugCond202Airflow.dbg","w");
		if(fp) {
			fprintf(fp,"Pdis_kPa=%le.\n",Pdis_kPa);
			fprintf(fp,"Pdis_psia=%le.\n",Pdis_psia);
			fprintf(fp,"A,B,C = %le,%le,%le.\n",P.cfmA,P.cfmB,P.cfmC);
			fprintf(fp,"af_cfm=%le.\n",af_cfm);
			fprintf(fp,"af_m3s=%le.\n",af_m3s);
			fprintf(fp,"af_kgs=%le.\n",af_kgs);
			fprintf(fp,"Ntubes=%le.\n",Ntubes);
			fprintf(fp,"Aflow_seg=%le.\n",Aflow_seg);
			fprintf(fp,"Aflow_tot=%le.\n",Aflow_tot);
			fprintf(fp,"GaNom=%le.\n",GaNom);
			fprintf(fp,"FoulFac=%le.\n",FoulFac);
			fprintf(fp,"Prms[3]=%le.\n",Prms[3]);
			fprintf(fp,"Ga=%le.\n",Ga);
			fclose(fp);
		}
		}
		#endif
	} else {
		if(GaI<0) {
			Ga=GaNom*FoulFac*Prms[3];
		} else {
			Ga=GaI*Prms[3];//shenbo correct the original Ga=GaI
		}
	}

	// TR Calculate the air cfm after GA is established
	const int Ntubes = P.Nsubtubes+P.Nbranchs*P.Nmaintubes;
	const double A_total = P.Aflow*P.NSeg*(double)Ntubes; //total area air flows throguh for entire coil (m^2)
	const double v_air = air.v(Tai); // specific volume of air
	const double vfra = Ga*A_total*v_air; // nominal volumetric flow rate air
	P.cfma = vfra*(60.0/((12.0*0.0254)*(12.0*0.0254)*(12.0*0.0254))); // nominal cfm of air flow
	
	//B.S., prepare to set up the simple model
	P.V_TOT =0; P.V_TP=0; P.V_Liq =0; P.V_Vap = 0;
	P.m_TOT = 0; P.m_TP =0; P.m_Liq =0; P.m_Vap =0;
	P.mr = 0; P.ma_TOT =0;
	P.UA_TOT=0; P.UA_Liq =0;P.UA_Vap = 0; P.UA_TP =0;
	P.LiqL = 0; P.TPL = 0; P.VapL = 0;
	P.count1=0; P.count2=0;
	P.HP_TP1.H=0; P.HP_TP1.P=0;
	P.HP_TP2.H=0; P.HP_TP2.P=0;
	//-----------------------------------------B.S.

	P.GetP=0;//don't calculate the airside pressure drop in the following procedure
	CondCircuit(type,mr,HPi,Tai,Ga,HPo,Tao,m,&P);
	if(errorLog.IsError()) {
		errorLog.Add("Condenser","CondCircuit");
		return 0;
	}

	const double DP_circuit = Circuit_DP(Ga,Tai,*Tao,&P);//B.S. calculate the airside pressure drop cross the circuit
	gOutData.Cond_AirDP = DP_circuit;//for output

	if(errorLog.IsError()) {//the error in calculating the airside pressure drop should not influence the main iteration loop
		errorLog.ClearError("Circuit_DP","Cond");
	}

	//B.S.---------------------------------------
	P.HP_in = HPi; P.HP_out = *HPo;
	P.mr = mr; P.Ga = Ga; P.Tai=Tai; P.ma_TOT = Ga*P.Aflow*P.NSeg*(P.Nbranchs*P.Nmaintubes+P.Nsubtubes);
	Build_Cond_Lumped(&P);

	*Cond_struc = P;

	//---------------------------------------B.S.

	return 0;
}

/********************************************************************
Solves the condenser model (Condenser()) iteratively.  It takes Tsc
(subcooling of refreigerant at outlet) as an input and solves for the
required same adjustment to the Prms[0] (air-side HT coeff.) and
Prms[1] (ref-side HT coeff.). Tsc is specified to protect against
the mistake of making Tsc<=0.
********************************************************************/
void Condenser0(double Tsc,//subcooling degree
				const char *filename,
				double mr,//refrigerant mass flow rate
				HP HPi,//inlet refrigerant state
				double Tai,//inlet air temperature
				double Ga,//air flow mass 
				HP *HPo,//outlet refrigerant state
				double *Tao,//outlet air temperature
				MASS *m,// mass and volume in  the condenser
				double *Prms)//adjustment parameters
{
	TXP TXPi = HPtoTXP(HPi);
	if(errorLog.IsError()) {
		errorLog.Add("Condenser0");
		return;
	}
	TXP	TXP_prop={0,0,0};

	TXP_prop.P=TXPi.P;
	TXP_prop.X=0;
	const double Tcond = PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(TXPi.P);
	if(errorLog.IsError()) {
		errorLog.Add("Condenser0");
		return;
	}

	double Tro = Tcond-Tsc;

	Cond0P P;
	P.filename = filename;
	P.mr=mr;
	P.HPi=HPi;
	P.Tai=Tai;
	P.Ga=Ga;
	P.Prms=Prms;
	P.Tro=Tro;
	P.Tsc=Tsc;//B.S. add for including the influence from the pressure drop

	double adj = FindZero(1.0,Cond0Func,1e-7,&P);

	if(errorLog.IsError())
	{
	errorLog.ClearError("Condenser0","toBisection1DSpline");
	const double Xlo=0.3; const double Xhi=5.0;
	adj=Bisection1DSpline(Xlo,Xhi,Cond0Func,1e-7,&P);
	}

	if(errorLog.IsError()) {
		errorLog.Add("Condenser0");
		return;
	}

	Prms[0] = adj;//B.S. adjusting the airside heat transfer coefficient
	Prms[1] = 1.0;
	Prms[2] = 1.0;

	CGP Cond_stuct;
	Condenser(filename,mr,HPi,Tai,Ga,HPo,Tao,m,&Cond_stuct,Prms);
	if(errorLog.IsError()) {
		errorLog.Add("Condenser0");
		return;
	}
}

/********************************************************************
Used by Condenser0() during its iteration.
********************************************************************/
static double Cond0Func(double Adj,
						void *Params)
{
	HP HPo;
	double Tao;
	MASS m;

	Cond0P* P = (Cond0P*)Params;
	P->Prms[0]=Adj;//B.S. adjusting the airside heat transfer coefficient
	P->Prms[1]=1.0;
	P->Prms[2]=1.0;

	CGP Cond_struc;
	Condenser(P->filename,P->mr,P->HPi,P->Tai,P->Ga,&HPo,&Tao,&m,&Cond_struc,P->Prms);
	if(errorLog.IsError()) {
		errorLog.Add("Cond0Func");
		return 0;
	}
	
	TXP TXP_prop={0,0,0};
	TXP_prop.P=HPo.P;
	if(P->Tsc>=0)
	{
	TXP_prop.X=0;
	}
	else
	{
	TXP_prop.X=-P->Tsc;
	}

	const double Tcond = PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(HPo.P);//B.S., for including the pressure drop effect, this is the saturated temperature corresponding to the outlet pressure of the condenser.
	if(errorLog.IsError()) {
		errorLog.Add("Condenser0");
		return 0;
	}
	
	P->Tro=Tcond-P->Tsc;

	if(P->Tro<P->Tai) P->Tro=P->Tai; //B.S. to remove the result that has no meaning.
	
	HP HPo2={0,0};

	if(P->Tsc>=0)
	{
	HPo2 = TXPtoHP(toTXP(P->Tro,0,HPo.P));
	}
	else//B.S.
	{
	HPo2 = TXPtoHP(toTXP(Tcond,-P->Tsc,HPo.P));
	}

	if(errorLog.IsError()) {
		errorLog.Add("Cond0Func");
		return 0;
	}

	return (HPo.H-HPo2.H)/HPo2.H;
}


#ifdef CONDENSER_PRESSURE_DROP
/*********************************************************************
Condenser tube segment model.
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
void CondTubeL_new(double Gr,//refigerant mass flux
						  HP *HPi,//refrigerant inlet state
						  double Ga,//air mass flux
						  double Tai,//air inlet temperature
						  double *hao,//air outlet enthalpy
						  MASS *m,//inner volume and the refrigerant in the segment
						  CGP *P)//condenser struct
{
	// convert inlet state into TXP format
	const double pi=atan(1.0)*4.0;
	TXP TXPi = HPtoTXP(*HPi);
	if(errorLog.IsError()) {
		errorLog.Add("CondTubeL_new","TXPi");
		return;
	}

	//air heat transfer 
	const double ho = P->hAirAdj*ConvCoeffAir_CON(Tai,Ga,P);//B.S. ConvCoeffAir_CON(Tai,Ga,P) it is a new function, which can call heat transfer calculation corresponding to different fin type
	const double phi = FinEffect_Schmidt(ho,237,P->th,P->y,P->Do);//fin efficiency, using Schmidt equation 
	P->Ro = 1/(ho*(P->Apo+phi*P->Af));//airside thermal resistance
	if(errorLog.IsError()) {
		errorLog.Add("CondTube");
		return;
	}
	
	P->airT = Tai;//B.S., this may be used as the local air inlet temperature, while Tai is used as the global inlet air temperature
	P->Ga = Ga;//B.S., keep the air state parameters for late use


	// Mass flow rates of air and refrigerant.
	double ma = Ga*P->Aflow;
	double mr = Gr*P->Ax;

	// Calculate heat transfered per unit mass of refrigerant.
	double hi = 0;

	hi=P->hRefAdj*ConvCoeffInside(TXPi,Gr,P->Di, P);//B.S., this fuction can call different heat transfer calculation, correponding to different tube type
	if(errorLog.IsError()) {
		errorLog.Add("CondTubeL_new","hi");
		return;
	}

	double Ri = 1/(hi*P->Api);//refrigerant side thermal resistance
	const double R_W=log(P->Do/(P->Do-2.0*P->xp))/(2.0*pi*P->K_T*P->Ls);
	double R = 0.0;

	//R = (P->Ro+R_W+Ri);
	
	//B.S.---------------------------------------------
	if(TXPi.X>=1)//B.S., the following caculate the whole thermal resistances separately, prepared to adjust them separately.	
	{
	P->UA_Vap=P->UA_Vap+1/R;//B.S, vapor-phase heat transfer conductance
	P->VapL = P->VapL+ P->Ls;//B.S., vapor-phase tube length
	R = (P->Ro+R_W+Ri);
	}
	else if(TXPi.X<=0)
	{
	P->LiqL = P->LiqL+ P->Ls;//B.S., liquid-phase tube length
	P->UA_Liq=P->UA_Liq+1/R;//B.S. liquid-phase heat transfer conductance
	R = (P->Ro+R_W+Ri);
	}
	else 
	{
	P->UA_TP=P->UA_TP+1/R;//B.S. two-phase heat transfer conductance
	P->TPL = P->TPL+ P->Ls;//B.S. two-phase heat transfer tube length
	R = (P->Ro+R_W+Ri);
	}
	//-----------------------------------------B.S.

	double q = (TXPi.T-Tai)/mr*CmineCrossFlow_dry(R,mr,ma,TXPi,Tai);
	if(errorLog.IsError()) {
		errorLog.Add("CondTubeL_new","q");
		return;
	}

	// Calculate outlet state.
	// q = dh (energy balance)

	HPi->H -= q;


	// calculating pressure drop
	double P_out; 
	const double D_H=P->Dh_i;//B.S., get the hydraulic diameter for calculating the pressure drop with Kedzierski correlation

	Gr = mr/P->Acs;//B.S., this mass flux is for calculating the pressure drop
	const double DP_FR = FricDP(TXPi, Gr, q, P);//calculate the fictional pressure drop, no matter single-phase or two-phase flow

	if(errorLog.IsError()) {
		errorLog.Add("CondTubeL_new", "DP_FR");
		return;
	}

	if(TXPi.X>0.05&&TXPi.X<0.99) {
	//B.S. packing up the parameters for calculating the acceleration pressure drop
	PreAcc Preacc;
	Preacc.DP_FR=DP_FR;
	Preacc.G=Gr;
	Preacc.H_OUT=HPi->H;
	Preacc.P_IN=TXPi.P;
	Preacc.X_IN=TXPi.X;
	const double DP_ACC=Zbrent(10,-10,GET_PreAcc,1e-7,&Preacc);//B.S. two-phase acceleration pressure drop
	if(errorLog.IsError()) {
		errorLog.Add("CondTubeL_new", "(Preacc)");
		return;
	}
	P_out=TXPi.P-(DP_FR+DP_ACC)*P->PRefAdj;
	}
	else 
	{P_out=TXPi.P - DP_FR*P->PRefAdj; }//assume single-phase flow has no accelerattion pressure drop

	HPi->P=P_out;
	
	
	// Determine mass of charge.  It is based on the inlet state.
	//	The specific volume is recalculated so that a different model
	//	can be used from the one used to calculate the pressure drop.

	double v = VolumeALL(TXPi,Gr,P->Di,-mr*q/P->Api);		// seperate flow model
	if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"v TXPi=(%lf,%lf,%lf) Gr=%lf Di=%lf\n",TXPi.T,TXPi.X,TXPi.P,Gr,P->Di);
		errorLog.Add("CondTubeL_new",str);
		return;
	}

	m->V = P->Ls*P->Acs;//B.S. use the estimate cross-sectional area or P->Acs more exact
	m->m = m->V/v;

	// Calculate output air state.
	double hai = air.h(Tai);
	if(errorLog.IsError()) {
		errorLog.Add("CondTubeL_new","hai");
		return;
	}
	*hao = hai+q*mr/ma;

	//B.S.----------------------------------------
	TXP TXPo;
	TXPo = HPtoTXP(*HPi);
	if(errorLog.IsError()) {
		errorLog.Add("CondTubeL_new","TXPo");
		return;
	}	

	
	if(TXPi.X>=0.999999999)//B.S., the following caculate the whole thermal resistances separately, prepared to adjust them separately.	
	{
	P->m_Vap = P->m_Vap+m->V/v;//B.S. vapor-phase mass
	P->V_Vap = P->V_Vap + m->V;//B.S. vapor-phase volume

	if(TXPo.X<0.999999999)
	{
	P->HP_TP1.H=HPi->H*mr+P->HP_TP1.H;//corresponding to inlet state of the first segment that involve with two-phase heat transfer calculation
	P->HP_TP1.P=HPi->P*mr+P->HP_TP1.P;//corresponding to inlet state of the first segment that involve with two-phase heat transfer calculation
	P->count1 = P->count1 + mr;//count the number of circuit
	}
	
	}
	else if(TXPi.X<=0.00000000001)
	{
	P->m_Liq = P->m_Liq+m->V/v;//B.S. liquid-phase mass
	P->V_Liq = P->V_Liq + m->V;//B.S. liquid-phase volume
	}
	else 
	{
	P->m_TP = P->m_TP+m->V/v;//B.S. two-phase mass
	P->V_TP = P->V_TP + m->V;//B.S. two-phae volume
	
	if(TXPo.X<=0.00000000001)
	{
	P->HP_TP2.H = P->HP_TP2.H+mr*HPi->H;//corresponding to the inlet state of the last segment that involve with two-phase heat transfer calculation
	P->HP_TP2.P= P->HP_TP2.P+mr*HPi->P;//corresponding to the inlet state of the last segment that involve with two-phase heat transfer calculation
	P->count2 = P->count2 + mr;//count the number of circuit
	}

	}
	//---------------------------------------------B.S.
}
#endif

/*********************************************************************
Condenser tube return bend model.
Inputs:
	G = refrigerant mass flux (kg/s/m^2)
	HPi = refrigerant inlet state (h,P)
Outputs:
	HPi = refrigerant outlet state (h,P)
	m = mass of charge in return bend (kg)
*********************************************************************/
void CondReturnBend(double G,HP *HPi,MASS *m,CGP *P)
{
	CLRP Q;
	TXP TXPi;
	double v;

	TXPi=HPtoTXP(*HPi);
	if(errorLog.IsError()) {
		errorLog.Add("CondReturnBend");
		return;
	}

	// Q.vi=VolumeALL(TXPi,G,P->Di);	// seperate flow model
	Q.vi = PropertyTXPth(VOL,TXPi);			// homogeneous flow model
	if(errorLog.IsError()) {
		errorLog.Add("CondReturnBend");
		return;
	}

	Q.K=LossCoeffReturnBend(TXPi,G,P->Brad,P->Di);
	if(errorLog.IsError()) {
		errorLog.Add("CondReturnBend");
		return;
	}

	Q.D=P->Di;
	Q.HPi=*HPi;
	Q.G=G;
	Q.q=0;          // adiabatic

	// Solve for outlet state.
	FindZero2DConst((double*)HPi,CompLossResid,HPLimitsConst,1e-6,&Q);
	if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"G=%lf HPi=(%lf,%lf)",G,HPi->H,HPi->P);
		errorLog.Add("CondReturnBend",str);
		return;
	}

	/* Determine mass of charge.  It is based on the inlet state.
		The specific volume is recalculated so that a different model
		can be used from the one used to calculate the pressure drop.  */

	v = VolumeALL(TXPi,G,P->Di,0);	// seperate flow model
	// v = PropertyTXP(VOL,TXPi);	// homogeneous flow model
	if(errorLog.IsError()) {
		errorLog.Add("CondReturnBend");
		return;
	}

	m->V=P->Blen*P->Ax;
	m->m=m->V/v;

	//B.S.----------------------------------------------------
	if(TXPi.X==1)//B.S., the following caculate the whole thermal resistances separately, prepared to adjust them separately.	
	{
	P->m_Vap = P->m_Vap+m->V/v;//B.S. vapor-phase mass
	P->V_Vap = P->V_Vap + m->V;//B.S. vapor-phase volume
	}
	else if(TXPi.X<1&&TXPi.X>0.00)
	{
	P->m_TP = P->m_TP+m->V/v;//B.S. two-phase mass
	P->V_TP = P->V_TP + m->V;//B.S. two-phase volume
	}
	else 
	{
	P->m_Liq = P->m_Liq+m->V/v;//B.S. liquid-phase mass
	P->V_Liq = P->V_Liq + m->V;//B.S. liquid-phase volume
	}
	//---------------------------------------------B.S.
}

/*********************************************************************
Condnenser manifold model.  Combines tubes together.  Holds no
refrigerant charge.
Inputs:
	ni = number of inlet tubes
	no = number of outlet tubes
	D = tube diameters, assumed same (m)
	G = inlet refrigerant mass flux (kg/s/m^2)
	HPi = refrigerant inlet state (h,P)
Outputs:
	HPi = refrigerant outlet state (h,P)
	G = outlet refrigerant mass flux (kg/s/m^2)
*********************************************************************/
void CondMan(int ni,int no,double D,double *G,HP *HPi)
{
	ACRP P;
	CLRP Q;
	TXP TXPi;

	/* isentropic area change */
	P.Di=D;
	P.HPi=*HPi;
	P.Gi=*G;
	*G*=(double)ni/(double)no;
	P.Go=*G;

	TXPi = HPtoTXP(*HPi);
	if(errorLog.IsError()) {
		errorLog.Add("CondMan","TXPi");
		return;
	}

	P.si = PropertyTXPth(ENTR,TXPi);
	if(errorLog.IsError()) {
		errorLog.Add("CondMan","P.si");
		return;
	}

	// specific volume
	// P.vi=VolumeALL(TXPi,P.Gi,P.Di);	// seperate flow model
	P.vi = PropertyTXPth(VOL,TXPi);			// homogeneous flow model
	if(errorLog.IsError()) {
		errorLog.Add("CondMan","P.vi");
		return;
	}

	FindZero2DConst((double*)HPi,AreaChangeResid,HPLimitsConst,1e-6,&P);
	if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"CondMan: (AreaChangeResid) Gr=%lf HPi=(%lf,%lf)\n",P.Go,HPi->H,HPi->P);
		errorLog.Add("CondMan",str);
		return;
	}

	// irreversible enterance loss
	Q.K=0.5;
	Q.D=D;
	Q.HPi=*HPi;
	Q.G=*G;
	Q.q=0;

	TXPi=HPtoTXP(*HPi);
	if(errorLog.IsError()) {
		errorLog.Add("CondMan");
		return;
	}

	// specific volume
	// Q.vi = VolumeALL(TXPi,Q.G,Q.D);	      // seperate flow model
	Q.vi = PropertyTXPth(VOL,TXPi);					// homogeneous flow model
	if(errorLog.IsError()) {
		errorLog.Add("CondMan");
		return;
	}

	FindZero2DConst((double*)HPi,CompLossResid,HPLimitsConst,1e-6,&Q);
	if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"(CompLossResid) Gr=%lf HPi=(%lf,%lf)",Q.G,HPi->H,HPi->P);
		errorLog.Add("CondMan",str);
		return;
	}
}

/********************************************************************
Generates two residuals that must be zero at the corrent outlet
state of the manifold.  Energy conservation requires that the
stagnation enthalpy across the manifold be zero and the assumption
that the process is reversible and adiabatic results in and
isentropic process that forms the second residual (s1=s2).
********************************************************************/
static double AreaChangeResid(double *X,double *Y,void *Params)
{
	HP* HPo = (HP*)X;
	ACRP* P = (ACRP*)Params;

	// outlet specific volume and entropy
	TXP TXPa = HPtoTXP(*HPo);
	if(errorLog.IsError()) {
		errorLog.Add("AreaChangeResid","TXPa");
		return 0;
	}

	const double so = PropertyTXPth(ENTR,TXPa);
	if(errorLog.IsError()) {
		errorLog.Add("AreaChangeResid","so");
		return 0;
	}

	// specific volume
	// double vo = VolumeALL(TXPa,P->Go,P->Di);	// seperate flow model
	const double vo = PropertyTXPth(VOL,TXPa);				// homogeneous flow model
	if(errorLog.IsError()) {
		errorLog.Add("AreaChangeResid","vo");
		return 0;
	}

	/* inlet and outlet stagnation enthalpy */
	const double xi = P->vi*P->Gi;
	const double h0i = P->HPi.H + 0.5*xi*xi;
	const double xo = vo*P->Go;
	const double h0o = HPo->H + 0.5*xo*xo;

	// stagnation enthalpy residual
	Y[0] = (h0o-h0i)/h0i;
	// entropy residual
	Y[1] = (so-P->si)/P->si;

	return 0;
}

/********************************************************************
Used to solve for the change of state considering the energy and
momentum conservation equations for compressible flow.
********************************************************************/
static double CompLossResid(double *X,double *Y,void *Params)
{
	HP* HPo = (HP*)X;
	CLRP* P = (CLRP*)Params;

	// outlet specific volume and entropy
	TXP TXPo = HPtoTXP(*HPo);
	if(errorLog.IsError()) {
		errorLog.Add("CompLossResid");
		return 0;
	}

	// specific volume
	// double vo = VolumeALL(TXPo,P->G,P->D);	   // seperate flow model
	const double vo = PropertyTXPth(VOL,TXPo);				// homogeneous flow model
	if(errorLog.IsError()) {
		errorLog.Add("CompLossResid","vo");
		return 0;
	}

	// inlet and outlet stagnation enthalpy
	const double xi = P->vi*P->G;
	const double h0i = P->HPi.H+0.5*xi*xi;		// with kinetic energy terms
	const double xo = vo*P->G;
	const double h0o = HPo->H+0.5*xo*xo;

	// momentum conservation
	double t1 = P->G*P->G*(P->vi + vo)*P->K/4.0;
	double t2 = 1.0e3*(HPo->P-P->HPi.P);
	double t3 = P->G*P->G*(vo - P->vi);

	// stagnation enthalpy residual
	Y[0] = (h0i-h0o-P->q)/h0i;
	Y[1] = (t1+t2+t3)/t1;

	return 0;
}

/********************************************************************
Provides constraints used in solving for the outlet state of the
manifold.  The numerical solver updates guesses of enthalpy and
pressure.  This functions states a new guess to ensure that is
within the property tables.
********************************************************************/
static int HPLimitsConst(double* X,double* /*Y*/,void* /*Params*/)
{
	HP* HPo=(HP*)X;
	TXP TXP_prop={0,0,0};

	if(HPo->P<PMINth || HPo->P>PMAXth) return 1;

	TXP_prop.P=HPo->P;
	TXP_prop.T=TMIN;
	TXP_prop.X=0;
	const double hmin = PropertyTXPth(ENTH,TXP_prop);//refscthPT.h(HPo->P,TMIN);
	if(errorLog.IsError()) {
		errorLog.Add("AreaChangeResidConst","hmin");
		return 0;
	}
	if(HPo->H<hmin) return 1;
	TXP_prop.P=HPo->P;
	TXP_prop.T=TMAX;
	TXP_prop.X=1;
	const double hmax = PropertyTXPth(ENTH,TXP_prop);//refshthPT.h(HPo->P,TMAX);
	if(errorLog.IsError()) {
		errorLog.Add("AreaChangeResidConst","hmax");
		return 0;
	}
	if(HPo->H>hmax) return 1;

	return 0;
}

//----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//the followings are some old subroutines not used currently

/*********************************************************************
Models an entire condenser tube from return bend to return bend.
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
void CondTube(double Gr,//refrigerant mass flux
			  HP *HPi,//inlet refrigerant state
			  double Ga,//air flux
			  double Tai,//air inlet temperature
			  double *hao,//air outlet enthalpy
			  MASS *m,//inner volume and the charge mass 
			  CGP *P)// condenser struct
{
	MASS mi;

	// Calculate air side resistance.
	double ho = P->hAirAdj*ConvCoeffAir_CON(Tai,Ga,P);//B.S. ConvCoeffAir_CON(Tai,Ga,P) it is a new function, which can call heat transfer calculation corresponding to different fin type
	double phi = FinEffect_Schmidt(ho,237,P->th,P->y,P->Do);//fin efficiency, using Schmidt equation 
	P->Ro = 1/(ho*(P->Apo+phi*P->Af));//airside thermal resistance
	if(errorLog.IsError()) {
		errorLog.Add("CondTube");
		return;
	}
	
	P->airT = Tai;//B.S., this may be used as the local air inlet temperature, while Tai is used as the global inlet air temperature
	P->Ga = Ga;//B.S., keep the air state parameters for late use

	double Shao=0;
	m->m=0;
	m->V=0;
	for(int i=0;i<P->NSeg;i++) {

		#ifdef CONDENSER_PRESSURE_DROP
		CondTubeL_new(Gr,HPi,Ga,Tai,hao,&mi,P);//B.S. new function for calculating the heat transfer and pressure drop in a condensing tube segment
		//CondTubeL(Gr,HPi,Ga,Tai,hao,&mi,P); B.S. the old function is not used.
		#else
		CondTubeLnoDP(Gr,HPi,Ga,Tai,hao,&mi,P);	
		#endif

		if(errorLog.IsError()) {
			errorLog.Add("CondTube");
			return;
		}

		Shao+=*hao;
		m->V+=mi.V;
		m->m+=mi.m;
	}

	/* outlet air state */
	*hao=Shao/(double)P->NSeg;
	if(errorLog.IsError()) {
		errorLog.Add("CondTube");
		return;
	}
}


#ifdef CONDENSER_PRESSURE_DROP
/*********************************************************************
Condenser tube segment model.
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
static void CondTubeL(double Gr,HP *HPi,double Ga,double Tai,double *hao,MASS *m,CGP *P)
{
	CLRP Q;

	// convert inlet state into TXP format
	TXP TXPi = HPtoTXP(*HPi);
	if(errorLog.IsError()) {
		errorLog.Add("CondTubeL","TXPi");
		return;
	}

	 Q.vi = VolumeALL(TXPi,Gr,P->Di,0);	// seperate flow model
//	Q.vi = PropertyTXP(VOL,TXPi);			// homogeneous flow model
	if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"Q.vi TXPi=(%lf,%lf,%lf)\n",TXPi.T,TXPi.X,TXPi.P);
		errorLog.Add("CondTubeL",str);
		return;
	}

	// Mass flow rates of air and refrigerant
	double ma = Ga*P->Aflow;
	double mr = Gr*P->Ax;

	// Calculate heat transfered per unit mass of refrigerant.
	double hi = P->hRefAdj*ConvCoeffInside(TXPi,Gr,P->Di, P);
	if(errorLog.IsError()) {
		errorLog.Add("CondTubeL","hi");
		return;
	}

	double Ri = 1/(hi*P->Api);
	double R = P->Ro+Ri;

	Q.q = (TXPi.T-Tai)/mr*CmineCrossFlow_dry(R,mr,ma,TXPi,Tai);
	if(errorLog.IsError()) {
		errorLog.Add("CondTubeL","Q.q");
		return;
	}

	Q.K = LossCoeffStraightPipe(TXPi,Gr,P->Ls,P->Di);
	if(errorLog.IsError()) {
		errorLog.Add("CondTubeL","Q.K");
		return;
	}

	Q.D=P->Di;
	Q.HPi=*HPi;
	Q.G=Gr;

	// Solve for outlet state
	FindZero2DConst((double*)HPi,CompLossResid,HPLimitsConst,1e-6,&Q);
	if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"G=%lf HPi=(%lf,%lf)",Gr,HPi->H,HPi->P);
		errorLog.Add("CondTubeL",str);
		return;
	}

	// Determine mass of charge.  It is based on the inlet state.
	//	The specific volume is recalculated so that a different model
	//	can be used from the one used to calculate the pressure drop.

	double v = VolumeALL(TXPi,Gr,P->Di,mr*Q.q/P->Api);		// seperate flow model
	// double v = PropertyTXP(VOL,TXPi);		// homogeneous flow model
	if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"v TXPi=(%lf,%lf,%lf) Gr=%lf Di=%lf\n",TXPi.T,TXPi.X,TXPi.P,Gr,P->Di);
		errorLog.Add("CondTubeL",str);
		return;
	}

	m->V=P->Ls*P->Ax;
	m->m=m->V/v;

	/* Calculate output air state. */
	double hai = air.h(Tai);
	if(errorLog.IsError()) {
		errorLog.Add("CondTubeL","hai");
		return;
	}

	*hao = hai+Q.q*mr/ma;
}

#else
/*********************************************************************
Condenser tube segment model.
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
static void CondTubeLnoDP(double Gr,HP *HPi,double Ga,double Tai,double *hao,MASS *m,CGP *P)
{
	// convert inlet state into TXP format
	TXP TXPi = HPtoTXP(*HPi);
	if(errorLog.IsError()) {
		errorLog.Add("CondTubeLnoDP","TXPi");
		return;
	}

	// Mass flow rates of air and refrigerant.
	double ma = Ga*P->Aflow;
	double mr = Gr*P->Ax;

	// Calculate heat transfered per unit mass of refrigerant.
	double hi = P->hRefAdj*ConvCoeffInside(TXPi,Gr,P->Di);
	if(errorLog.IsError()) {
		errorLog.Add("CondTubeLnoDP","hi");
		return;
	}

	double Ri = 1/(hi*P->Api);
	double R = P->Ro+Ri;

	double q = (TXPi.T-Tai)/mr*CmineCrossFlow(R,mr,ma,TXPi,Tai);
	if(errorLog.IsError()) {
		errorLog.Add("CondTubeLnoDP","q");
		return;
	}

	// Calculate outlet state.
	// 1. No pressure drop
	// 2. q = dh (energy balance)

	HPi->H -= q;

	// Determine mass of charge.  It is based on the inlet state.
	//	The specific volume is recalculated so that a different model
	//	can be used from the one used to calculate the pressure drop.

	double v = VolumeALL(TXPi,Gr,P->Di);		// seperate flow model
	// v = PropertyTXP(VOL,TXPi);		      // homogeneous flow model
	if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"v TXPi=(%lf,%lf,%lf) Gr=%lf Di=%lf\n",TXPi.T,TXPi.X,TXPi.P,Gr,P->Di);
		errorLog.Add("CondTubeLnoDP",str);
		return;
	}

	m->V = P->Ls*P->Ax;
	m->m = m->V/v;

	// Calculate output air state.
	double hai = air.h(Tai);
	if(errorLog.IsError()) {
		errorLog.Add("CondTubeLnoDP","hai");
		return;
	}

	*hao = hai+q*mr/ma;
}

#endif




//B.S.--------------------------------------------------
//*******************************************************************************
//*******************************************************************************
//*******************************************************************************
//the following subroutine is for correlating the moving boundary and lumped heat exchanger models
void Build_Cond_Lumped(CGP* P)
{
	//lumped model-----------------------------
	const double cp_a= air.Cp(P->Tai);
	if(errorLog.IsError()) {
		errorLog.Add("Build_Cond", "cp_a");
		return;
	}

	TXP TXP_prop = {0,0,0};
	TXP TXP1;//inlet condenser state
	TXP1= HPtoTXP(P->HP_in);
	if(errorLog.IsError()) {
		errorLog.Add("Build_Cond", "TXP1");
		return;
	}

	TXP_prop.P=TXP1.P;
	TXP_prop.X=1;
	const double Tsat1 = PropertyTXPth(TSAT,TXP_prop);//reftpvthP.Tsat(TXP1.P);//saturation refrigerant temperature at the inlet condenser state
	if(errorLog.IsError()) {
		errorLog.Add("Build_Cond", "Tsat1");
		return;
	}

	const double Q_TOT = P->mr*(P->HP_in.H-P->HP_out.H);//overall heat transfer amount
	const double epsilon_TOT = Q_TOT/(P->ma_TOT*cp_a)/(Tsat1-P->Tai);//overall effectivenss
	P->DP_TOT = (P->HP_in.P-P->HP_out.P);//overall pressure drop
	const double NTU_TOT = -log(1-epsilon_TOT);
	P->UA_TOT = NTU_TOT*P->ma_TOT*cp_a;//overall heat transfer conductance
	//lumped model ended

	return;
}

void Build_Cond_MovBound(CGP* P)
{

	TXP TXP_prop={0,0,0};
	//average density, pressure drop, cross-sectional area of each phase-------
	if(P->count1>0)
	{
	P->HP_TP1.H=P->HP_TP1.H/P->count1;//averge enthalpy at the first two-phase segment of different circuits
	P->HP_TP1.P=P->HP_TP1.P/P->count1;//averge pressure at the first two-phase segment of different circuits
	}
	else//otherwise use the inlet state of the condenser
	{
	P->HP_TP1 = P->HP_in;
	}

	if(P->count2>0)
	{
	P->HP_TP2.H=P->HP_TP2.H/P->count2;//averge enthalpy at the last two-phase segment of different circuits
	P->HP_TP2.P=P->HP_TP2.P/P->count2;//averge pressure at the last two-phase segment of different circuits
	}
	else//otherwise use the outlet state of the condenser
	{
	P->HP_TP2 = P->HP_out;
	}

	
	P->L_TOT = P->VapL+P->TPL+P->LiqL;//total tube length of the condenser
	P->V_TOT = P->V_Liq+P->V_Vap+P->V_TP;//total inner volume of the condenser
	P->A_TOT=P->V_TOT/P->L_TOT;//average heat transfer area per tube length

	//------------------------------average density, pressure drop, cross-sectional area of each phase


	//prepare heat transfer calculation---------------
	P->Ga_meanL = P->ma_TOT/P->L_TOT;//averge air flow rate across per tube length

	const double ma_v = P->Ga_meanL * P->VapL;//air flow rate across the gas-phase region
	const double ma_tp = P->Ga_meanL * P->TPL;//air flow rate across the two-phase region
	const double ma_l = P->Ga_meanL * P->LiqL;//air flow rate across the liquid-phase region
	
	const double cp_a= air.Cp(P->Tai);
	if(errorLog.IsError()) {
		errorLog.Add("Build_Cond", "cp_a");
		return;
	}

	TXP TXP1;//inlet condenser state
	TXP1= HPtoTXP(P->HP_in);
	if(errorLog.IsError()) {
		errorLog.Add("Build_Cond", "TXP1");
		return;
	}

	TXP_prop.P=TXP1.P;
	TXP_prop.X=1;
	const double Tsat1 = PropertyTXPth(TSAT,TXP_prop);//reftpvthP.Tsat(TXP1.P);//saturation refrigerant temperature at the inlet condenser state
	if(errorLog.IsError()) {
		errorLog.Add("Build_Cond", "Tsat1");
		return;
	}

	TXP_prop.P=TXP1.P;
	TXP_prop.T=Tsat1;
	TXP_prop.X=1;
	const double cp_rv = PropertyTXPtr(SPEC,TXP_prop);//refshtrPT.Cp(TXP1.P,Tsat1);
	if(errorLog.IsError()) {
		errorLog.Add("Build_Cond", "cp_rv");
		return;
	}

	TXP TXP2;//state at the first two-phase segment of the condenser
	TXP2= HPtoTXP(P->HP_TP1);
	if(errorLog.IsError()) {
		errorLog.Add("Build_Cond", "TXP2");
		return;
	}
	
	TXP_prop.P=TXP2.P;
	TXP_prop.X=1;
	const double Tsat2 = PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(TXP2.P);
	if(errorLog.IsError()) {
		errorLog.Add("Build_Cond", "Tsat2");
		return;
	}

	TXP_prop.P=P->HP_TP1.P;
	TXP_prop.T=Tsat2;
	TXP_prop.X=1;

	const double Hsat1= PropertyTXPth(ENTH,TXP_prop);//refshthPT.h(P->HP_TP1.P,Tsat2);
	if(errorLog.IsError()) {
		errorLog.Add("Build_Cond","Hsat1");
		return;
	}
	P->H1_residual = Hsat1 - P->HP_TP1.H;//for maintaining the consistency in the moving boundary model analysis

	TXP TXP3;
	TXP3= HPtoTXP(P->HP_TP2);//state at the last two-phase segment of the condenser
	if(errorLog.IsError()) {
		errorLog.Add("Build_Cond", "TXP3");
		return;
	}

	TXP_prop.P=TXP3.P;
	TXP_prop.X=1;
	const double Tsat3 = PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(TXP3.P);
	if(errorLog.IsError()) {
		errorLog.Add("Build_Cond", "Tsat3");
		return;
	}

	TXP_prop.P=P->HP_TP2.P;
	TXP_prop.T=Tsat3;
	TXP_prop.X=0;

	const double Hsat2= PropertyTXPth(ENTH,TXP_prop);//refscthPT.h(P->HP_TP2.P,Tsat3);
	if(errorLog.IsError()) {
		errorLog.Add("Build_Cond","Hsat2");
		return;
	}

	if(P->LiqL>0) P->H2_residual = Hsat2 - P->HP_TP2.H;//for maintaining the consistency in the moving boundary model analysis
	//Note: without "if(P->LiqL>0)" it may cause divengence at some cases
	
	TXP_prop.P=TXP3.P;
	TXP_prop.T=Tsat3;
	TXP_prop.X=0;
	const double cp_rl = PropertyTXPtr(SPEC,TXP_prop);//refsctrPT.Cp(TXP3.P,Tsat3);
	if(errorLog.IsError()) {
		errorLog.Add("Build_Cond", "cp_rl");
		return;
	}
	
	TXP TXP4;//condenser outlet state
	TXP4= HPtoTXP(P->HP_out);
	if(errorLog.IsError()) {
		errorLog.Add("Build_Cond", "TXP4");
		return;
	}

	//---------------------------Prepare heat transfer calculation

	//-----------------------------------------heat transfer of the vapor phase
	double cmin_v=0;
	double Cr_v = 0;

	if(ma_v>0)
	{
	const double Q_v = P->mr*(P->HP_in.H-P->HP_TP1.H);//vapor-phase heat transfer
	Cr_v = (P->mr*cp_rv)/(ma_v*cp_a);

	if(ma_v*cp_a>P->mr*cp_rv) {cmin_v = P->mr*cp_rv; Cr_v = Cr_v;}
	else {cmin_v = ma_v*cp_a;Cr_v = 1/Cr_v;}
	const double Qv_max = cmin_v*(TXP1.T - P->Tai);
	const double e_v = Q_v/Qv_max;//actual heat transfer effectiveness of the vapor-phase
	const double NTU_v = P->UA_Vap/cmin_v;
	
	double e_sup =0;//theoritical heat transfer effectivenss of the vapor-phase
	if(ma_v*cp_a>P->mr*cp_rv) {//cmax unmixed
	const double l = 1-exp(-Cr_v*NTU_v);
	e_sup = 1 - exp(-l/Cr_v);
	}
	else {//cmin unmixed
	const double l = -Cr_v*(1-exp(-NTU_v));
	e_sup = (1-exp(l))/Cr_v;
	}
	
	P->r_v = e_v/e_sup;//for ajusting the theortical heat transfer effectivenss of vapor-phase
	P->U_Vap = P->UA_Vap/(P->VapL);//average vapor-phase heat transfer conduntance per tube length

	P->rho_Vap = P->m_Vap/P->V_Vap;//average vapor-phase density
	P->DP_Vap = (P->HP_in.P-P->HP_TP1.P)/P->V_Vap;//average pressure drop gradient
	P->A_Vap=P->V_Vap/P->VapL;//average heat transfer surface area per tube length

	}
	//--------------------------------------------------------heat transfer of vapor phase
	
	//heat transfer of liquid phase--------------------
	double cmin_l=0;
	double Cr_l = 0;
		
	if(ma_l>0)
	{
	const double Q_l = P->mr*(P->HP_TP2.H-P->HP_out.H);//liquid-phase heat transfer amount
	Cr_l = (P->mr*cp_rl)/(ma_l*cp_a);

	if(ma_l*cp_a>P->mr*cp_rl) {cmin_l = P->mr*cp_rl; Cr_l=Cr_l;}
	else {
		cmin_l = ma_l*cp_a; Cr_l =1/Cr_l;}
	const double Ql_max = cmin_l*(TXP3.T - P->Tai);
	const double e_l = Q_l/Ql_max;//actual liquid-phase heat tranfer effectivenss
	const double NTU_l=P->UA_Liq/cmin_l;

	double e_sub =0;//theoritical liquid-phase heat transfer effectivenss

	if(ma_l*cp_a>P->mr*cp_rl) {//cmax unmixed
	const double l = 1-exp(-Cr_l*NTU_l);
	e_sub = 1 - exp(-l/Cr_l);
	}
	else {
	const double l = -Cr_l*(1-exp(-NTU_l));
	e_sub = (1-exp(l))/Cr_l;
	}

	P->r_l = e_l/e_sub;//parameter for adjusting the theoritical liquid-phase heat transfer effectivess
	P->U_Liq = P->UA_Liq/(P->LiqL);//averge liquid-phase heat transfer conductance

	P->DP_Liq = (P->HP_TP2.P-P->HP_out.P)/P->V_Liq;//avergage liquid-phase pressure drop gradient
	P->A_Liq=P->V_Liq/P->LiqL;//avergage liquid-phase heat transfer surface area
	P->rho_Liq = P->m_Liq/P->V_Liq;//average liquid-phase density
	}
	else//otherwise use the default parameters
	{
	static int FirCond=0;
	if(FirCond==0)
	{
	P->r_l = 1.0;
	P->U_Liq = 29;

	P->DP_Liq = 102247;
	P->A_Liq=P->A_TOT;
	P->rho_Liq = 1100;
	FirCond=1;
	}
	}
	//--------------------------------------heat transfer of liquid-phase

	//two-phase heat transfer--------------------------
	const double Q_tp = P->mr*(P->HP_TP1.H-P->HP_TP2.H);//two-phase heat transfer amount
	double cmin_tp = ma_tp*cp_a;	
	const double Qtp_max = cmin_tp*(TXP2.T - P->Tai);
	const double e_TP = Q_tp/Qtp_max;//actual two-phase heat transfer effectivenss
	const double NTU_tp = P->UA_TP/cmin_tp;
	const double e_sat = 1-exp(-NTU_tp);//theoritical two-phase heat transfer effectivenss

	P->r_tp = e_TP/e_sat;//parameter for adjusting the theoritical two-phase heat transfer effectivenss	
	P->U_TP = P->UA_TP/(P->TPL);//average two-phase heat transfer conductance per tube length

	P->rho_TP = P->m_TP/P->V_TP;//average two-phase density
	P->A_TP=P->V_TP/P->TPL;//average two-phase heat transfer surface area per tube length
	P->DP_TP = (P->HP_TP1.P-P->HP_TP2.P)/P->V_TP;//average pressure drop gradient
	//------------------------------two-phase heat transfer


	return;	
}

//****************************************************************************
//***************************************************************************
//*****************************************************************************
//the following subroutines are for simplified heat exchanger model analysis, not used now

double Condenser_Simple(
			  double mr,
			  HP HPi,
			  double Tai,
			  HP *HPo,
			  MASS *m,
			  CGP *P)
{
	double L_left =0, Lmax_TP=0;
	//gas phase
	Phase_Cond Gas;
	Gas.P = *P;
	Gas.mr = mr;
	Gas.HPi = HPi;
	Gas.Tai = Tai;

	double residual=0;
	residual = GasCal_Cond(P->L_TOT, &Gas);
	if(errorLog.IsError()) {
		errorLog.Add("Condenser_Simple", "(GasCal_Cond)");
		return 0;
	};

	if(residual >0)
	{
	*HPo = Gas.HPo;
	m->m = Gas.m;
	goto Left_Volume;
	}

	Zbrent(P->L_TOT,0.001,GasCal_Cond,1e-7,&Gas);

	if(errorLog.IsError()) {
		errorLog.Add("Condenser_Simple", "(Zbrent1)");
		return 0;
	};
 
	Phase_Cond TP;
	TP.P = *P;
	TP.mr = mr;
	TP.HPi = Gas.HPo;
	TP.Tai = Tai;
	
	Lmax_TP = P->L_TOT-Gas.L;
	residual = TPCal_Cond(Lmax_TP, &TP);
	if(residual >=0)
	{
	*HPo = TP.HPo;
	m->m = Gas.m+TP.m;
	goto Left_Volume;
	}
	
	Zbrent(Lmax_TP,0.00000000001,TPCal_Cond,1e-7,&TP);

	if(errorLog.IsError()) {
		errorLog.Add("Condenser_Simple", "(Zbrent2)");
		return 0;
	};
	
	Phase_Cond Liq;
	Liq.P = *P;
	Liq.mr = mr;
	Liq.HPi = TP.HPo;
	Liq.Tai = Tai;

	L_left = P->L_TOT - Gas.L-TP.L;
	if(L_left <=0.0000000001)
	{
	Liq.HPo=TP.HPo;
	Liq.m=0;
	}
	else{
		LiqCal_Cond(L_left,&Liq);}

	if(errorLog.IsError()) {
		errorLog.Add("Condenser_Simple", "LiqCal");
		return 0;
	};

	*HPo = Liq.HPo;
	m->m = Gas.m+TP.m+Liq.m;

	P->LiqL = L_left;//B.S. new

	//shenbo add for accouting the unmeasured volume
	Left_Volume:;
		TXP TXP_left = HPtoTXP(*HPo);//B.S., get the outlet state of the condenser
		if(errorLog.IsError()) {
			errorLog.Add("CondCircuit201","TXP_left");
			return 0;
			}
		const double V_left=0.00037;//B.S., the unaccounted volume at the subcooled site, which is got by tuning procedure
		const double Volume_left= PropertyTXPth(VOL,TXP_left);//B.S., refrigerant specific volume at the unaccounted subcooled region
		if(errorLog.IsError()) {
			errorLog.Add("CondCircuit201","Volume_left");
			return 0;
			}
		const double rho_left=1/Volume_left;

		const double m_left=V_left/Volume_left;//B.S., refrigerant mass in the unaccounted volume

	//	m->m+= m_left;

	return 0;
}


double GasCal_Cond(double L, void *Params)
{
	Phase_Cond* D=(Phase_Cond*)Params;
	
	TXP TXP_prop={0,0,0};
	double cmin_v=0;
	double Cr_v = 0;
	const double ma_v = L*D->P.Ga_meanL;
	const double cp_a= air.Cp(D->Tai);
	if(errorLog.IsError()) {
		errorLog.Add("GasCal_Cond", "cp_a");
		return 0;
	}

	TXP TXP1;
	TXP1= HPtoTXP(D->HPi);
	if(errorLog.IsError()) {
		errorLog.Add("GasCal_Cond", "TXP1");
		return 0;
	}

	if(TXP1.X<1)
	{
	D->L =0;
	D->HPo=D->HPi;
	D->m =0;
	return 0;
	}

	TXP_prop.P=TXP1.P;
	TXP_prop.X=1;
	const double Tsat1 = PropertyTXPth(TSAT,TXP_prop);//reftpvthP.Tsat(TXP1.P);
	if(errorLog.IsError()) {
		errorLog.Add("GasCal_Cond", "Tsat1");
		return 0;
	}


	TXP_prop.P=TXP1.P;
	TXP_prop.T=Tsat1;
	TXP_prop.X=1;
	const double cp_rv = PropertyTXPtr(SPEC,TXP_prop);//refshtrPT.Cp(TXP1.P,Tsat1);
	if(errorLog.IsError()) {
		errorLog.Add("GasCal_Cond", "cp_rv");
		return 0;
	}

	Cr_v = (D->mr*cp_rv)/(ma_v*cp_a);

	if(ma_v*cp_a>D->mr*cp_rv) {cmin_v = D->mr*cp_rv; Cr_v = Cr_v;}
	else {cmin_v = ma_v*cp_a;Cr_v = 1/Cr_v;}

	const double Qv_max = cmin_v*(TXP1.T - D->Tai);
	const double NTU_v = D->P.U_Vap*L/cmin_v;
	
	double e_sup =0;
	if(ma_v*cp_a>D->mr*cp_rv) {//cmax unmixed
	const double l = 1-exp(-Cr_v*NTU_v);
	e_sup = 1 - exp(-l/Cr_v);
	}
	else {//cmin unmixed
	const double l = -Cr_v*(1-exp(-NTU_v));
	e_sup = (1-exp(l))/Cr_v;
	}

	e_sup =e_sup*D->P.r_v;
	
	const double Q_gas = Qv_max*e_sup;

	const double V_Vap = D->P.A_Vap*L;//( int(L/0.2)*0.2+0.2);
	const double DP_Vap = V_Vap * D->P.DP_Vap;
	D->m = V_Vap*D->P.rho_Vap;
	
	D->HPo.P = D->HPi.P - DP_Vap;
	D->HPo.H = D->HPi.H - Q_gas/D->mr;

	TXP_prop.P=D->HPo.P;
	TXP_prop.X=1;

	const double Tsat2 = PropertyTXPth(TSAT,TXP_prop);//reftpvthP.Tsat(D->HPo.P);
	if(errorLog.IsError()) {
		errorLog.Add("GasCal_Cond", "Tsat2");
		return 0;
	}

	TXP_prop.P=D->HPo.P;
	TXP_prop.T=Tsat2;
	TXP_prop.X=1;
	const double Hsat= PropertyTXPth(ENTH,TXP_prop)-D->P.H1_residual;//refshthPT.h(D->HPo.P,Tsat2)-D->P.H1_residual;
	if(errorLog.IsError()) {
		errorLog.Add("GasCal_Cond","Hsat");
		return 0;
	}
	D->L =L;

	return (D->HPo.H- Hsat)/fabs(Hsat + D->HPo.H);
}

double TPCal_Cond(double L, void *Params)
{
	Phase_Cond* D=(Phase_Cond*)Params;

	const double ma_tp = L*D->P.Ga_meanL;

	const double cp_a= air.Cp(D->Tai);
	if(errorLog.IsError()) {
		errorLog.Add("TPCal_Cond", "cp_a");
		return 0;
	}

	TXP TXP2;
	TXP2= HPtoTXP(D->HPi);
	if(errorLog.IsError()) {
		errorLog.Add("TPCal_Cond", "TXP2");
		return 0;
	}

	const double cmin_tp = ma_tp*cp_a;


	const double Qtp_max = cmin_tp*(TXP2.T - D->Tai);
	const double NTU_tp = D->P.U_TP*L/cmin_tp;
	double e_sat = 1-exp(-NTU_tp);
	e_sat = e_sat *D->P.r_tp;

	const double Q_tp = Qtp_max*e_sat;

	const double V_TP = D->P.A_TP*L;//( int(L/0.2)*0.2+0.2);
	const double DP_TP = V_TP * D->P.DP_TP;

	D->m = V_TP*D->P.rho_TP;
	
	D->HPo.P = D->HPi.P - DP_TP;

	TXP TXP_prop={0,0,0};
	TXP_prop.P=D->HPo.P;
	TXP_prop.X=1;

	const double Tsat3 = PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(D->HPo.P);
	if(errorLog.IsError()) {
		errorLog.Add("TPCal_Cond", "Tsat3");
		return 0;
	}

	D->HPo.H= D->HPi.H - Q_tp/D->mr;
		
	TXP_prop.P=D->HPo.P;
	TXP_prop.T=Tsat3;
	TXP_prop.X=0;

	const double Hsat = PropertyTXPth(ENTH,TXP_prop)-D->P.H2_residual;//refscthPT.h(D->HPo.P,Tsat3)-D->P.H2_residual;
	if(errorLog.IsError()) {
		errorLog.Add("TPCal_Cond","Hsat");
		return 0;
	}
	D->L = L;


	return (D->HPo.H- Hsat)/fabs(D->HPo.H+ Hsat);
}

double LiqCal_Cond(double L, void *Params)
{
	Phase_Cond* D=(Phase_Cond*)Params;
	
	double cmin_l=0;
	double Cr_l = 0;
	const double ma_l = L*D->P.Ga_meanL;
	const double cp_a= air.Cp(D->Tai);
	if(errorLog.IsError()) {
		errorLog.Add("LiqCal_Cond", "cp_a");
		return 0;
	}
	TXP TXP_prop={0,0,0};

	TXP TXP3;
	TXP3= HPtoTXP(D->HPi);
	if(errorLog.IsError()) {
		errorLog.Add("LiqCal_Cond", "TXP3");
		return 0;
	}


	TXP_prop.P=TXP3.P;
	TXP_prop.X=0;

	const double Tsat3 = PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(TXP3.P);
	if(errorLog.IsError()) {
		errorLog.Add("LiqCal_Cond", "Tsat3");
		return 0;
	}

	TXP_prop.P=TXP3.P;
	TXP_prop.T=Tsat3;
	TXP_prop.X=0;

	const double cp_rl = PropertyTXPtr(SPEC,TXP_prop);//refsctrPT.Cp(TXP3.P,Tsat3);
	if(errorLog.IsError()) {
		errorLog.Add("LiqCal_Cond", "cp_rl");
		return 0;
	}

	Cr_l = (D->mr*cp_rl)/(ma_l*cp_a);

	if(ma_l*cp_a>D->mr*cp_rl) {cmin_l = D->mr*cp_rl; Cr_l=Cr_l;}
	else {
		cmin_l = ma_l*cp_a; Cr_l =1/Cr_l;}

	const double Ql_max = cmin_l*(TXP3.T - D->Tai);

	const double NTU_l=D->P.U_Liq*L/cmin_l;

	double e_sub =0;

	if(ma_l*cp_a>D->mr*cp_rl) {//cmax unmixed
	const double l = 1-exp(-Cr_l*NTU_l);
	e_sub = 1 - exp(-l/Cr_l);
	}
	else {//cmin unmixed
	const double l = -Cr_l*(1-exp(-NTU_l));
	e_sub = (1-exp(l))/Cr_l;
	}
	e_sub = e_sub*D->P.r_l;
	
	const double Q_liq = Ql_max*e_sub;

	const double V_Liq = D->P.A_Liq*L;
	const double DP_Liq = V_Liq * D->P.DP_Liq;

	D->m = V_Liq*D->P.rho_Liq;
	
	D->HPo.P = D->HPi.P - DP_Liq;
	D->HPo.H = D->HPi.H - Q_liq/D->mr;

	D->L = L;

	return 0;
}


//Lumped Condener model
double Condenser_Lumped(
			  double mr,
			  HP HPi,
			  double Tai,
			  HP *HPo,
			  MASS *m,
			  CGP *P)
{	
	const double cp_a= air.Cp(Tai);
	if(errorLog.IsError()) {
		errorLog.Add("Condenser_Lumped", "cp_a");
		return 0;
	}

	TXP TXP_prop={0,0,0};
	TXP_prop.P=HPi.P;
	TXP_prop.X=1;
	const double Tsat = PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(HPi.P);
	if(errorLog.IsError()) {
		errorLog.Add("Condenser_Lumped", "Tsat");
		return 0;
	}

	const double NTU = P->UA_TOT/(P->ma_TOT*cp_a);
	const double Q = (1-exp(-NTU))*P->ma_TOT*cp_a*(Tsat-Tai);
	HPo->H = HPi.H - Q/mr;
	HPo->P = HPi.P - P->DP_TOT;


	m->m=0; m->V=0;
	return 0;
}
