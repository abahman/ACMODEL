#ifndef _EVAP
#define _EVAP

#include "r22.h"
#include "wair.h"
#include "mass.h"

struct ETdim {
	double Di,L,xp,Df,z,th,vsp,Brad;
	int NSeg;
	double Dh,Do,y,Ls,N,Ax,Api,Apo,Af,Aflow,Vs,Blen,BVs;
	double Ro;
	double Gr,Ga;
	HP HPo;
	TXP TXPo;
	TP TPi;
	double hAirAdj,hRefAdj,PRefAdj,WAirAdj;
	int type;
	int Nrows,Ndeep;
	int NBranchs, NBraTube;

	double P_l; //spacing between tubs in the longitudual direction (m)
	//new parameters for micro-fin tubes, and specially configured fins
	int microfin;// decide whether is micro-fin tubes microfin=1
		
	double w_b;
	//base length for single micro-fin
	double w_e;
	// top length for single micro-fin
	double w_z;
	//width between the neighboring micro-fins
	double finH;
	//micro-fin height
	double beta;
	//micro-fin helical angle
	double gama;
	//micro-fin apex angle
	double finN;
	//total micro-fin number
	double Acs;
	// micro-fin cross-sectional area
	double P_H;
	//micro-fin hydraulic circumference
	double Dh_i;
	//micro-fin tube hydraulic diameter
	double K_T;
	//conductance factor of tube wall
	double K_F;
	//conductance factor of fin material
	double L_F;
	//reference fin length for Schimidt fin efficiency calculation
	double D_b,D_m;
	//fin base diameter, and mean diameter
	int airFin;
	double sub1,sub2,sub3,sub4,sub5;//for inputing sub-structures of fin surface
	double ho, wetadj;
	//airside heat transfer coefficient
	double Frontal_A;//frontal area
	double GetP;//calculate the airside pressure drop

	//variables for generating the simple evaporator function
	double V_TOT, V_TP, V_Liq, V_Vap;//inner volume of different phase
	double L_TOT,LiqL,VapL,TPL;//tube length of different phase
	double L_dry, L_wet;//tube lengh of dry and wet heat transfer

	double A_TOT,A_Liq,A_Vap,A_TP;//heat transfer surface area of different phase
	
	double m_TOT, m_TP, m_Liq, m_Vap;//airflow rate across different phase
	double rho_TOT, rho_TP, rho_Liq, rho_Vap;//average density of different phase
	
	double U_TP, U_Liq, U_Vap;//averge dry heat conductance of different phase 
	double Uw_TP, Uw_Liq, Uw_Vap;//average wet heat conductance of different phase
	
	double DP_TOT, DP_TP, DP_Liq, DP_Vap;//average pressure gradient of different phase
	
	double UA_TOT, UA_Liq,UA_Vap,UA_TP;//overall dry heat conductance of different phase
	double UAw_TOT, UAw_Liq,UAw_Vap,UAw_TP;//overall wet heat conductance of different phase

	double mr, ma_TOT;//overall refrigerant and air mass flow rate
	double Ga_meanL;//average air flow rate per tube length

	HP HP_out, HP_dry, HP_wet, HP_in,HP_TP1,HP_TP2;//state parameters at important locations
	double count1,count2;//count the state points of two-phase flow begin point and end point 
	double Qtp_dry, Qtp_wet;//two-phase dry heat transfer and wet heat transfer amount
	double r_v, r_l,r_tp;//parameters for adjusting the theoretical heat transfer effectiveness of each phase
	double H2_residual;//for the consistency of the moving boundary model analysis
	//------------------------------B.S.
	double q_flux, T_w;//segment heat flux and inside tube wall temperature
	int wet;//wet=0 to calculate dry heat transfer, wet=1 to calculate wet heat transfer
	int REV;
	double cfma;
	double Hout8[10],DPr[10];//superheat and pressure drop at each evaporator branch
};

void Build_Eva_MovBound(ETdim*);//B.S. fuction for building the simplified evaporator model, moving boundary model
void Build_Eva_Lumped(ETdim*);//B.S. fuction for building the simplified evaporator model, lumped model
double Evaporator_Simple(double,HP,TP,HP*,MASS*,ETdim*);//moving boundary model
double Evaporator_Lumped(double,HP,TP,HP*,MASS*,ETdim*);//lumped model

int EvapCircuit(int,double,HP*,double,TP,HP*,TP*,MASS*,double*,ETdim*);

int Evaporator(const char *filename,double,HP*,double,TP,HP*,TP*,MASS*,double*,ETdim*,double*,int NSeg=-1);
void Evaporator0(int,double,const char *filename,double,HP,double,TP,HP*,TP*,MASS*,double*,double*,int);
void Evaporator1(double,const char *filename,double,HP,double,TP,HP*,TP*,MASS*,double*,double*);

int EvapTube(double,HP*,double,TP,WH*,MASS*,ETdim*);
void EvapTubeBend(double,HP*,MASS*,ETdim*);
void EvapTubeBend_Fwd(double,HP*,MASS*,ETdim*);
int EvapTubeL_Rev(double,HP*,double,TP,WH*,MASS*,ETdim*);//B.S., new
int EvapTubeL_Fwd(double,HP*,double,TP,WH*,MASS*,ETdim*);//B.S., new


//
// Shell-and-tube heat exchanger with refrigerant int the tubes
// and water in the shell.
int Evaporator(const char* filename,double mr,HP HPo,double mw,double Twi,HP& HPi,double& Two,MASS& m,double* Prms);

#endif
