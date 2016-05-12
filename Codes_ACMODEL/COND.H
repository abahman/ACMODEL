#ifndef _COND
#define _COND

#include "r22.h"
#include "mass.h"

struct CGP {
	double Di,L,xp,Df,z,th,vsp,Brad;
	int NSeg;
	double Dh,Do,y,Ls,N,Ax,Api,Apo,Af,Aflow,Vs,Blen;
	double Ro;
	double hAirAdj,hRefAdj,PRefAdj;
	double hRefAdj_Sub;//B.S. parameter for adjusting the heat transfer ratio between the two-phase heat transfer and single-phase heat transfer
	int Nbranchs,Nmaintubes,Nsubtubes;
	int Ndeep;//B.S. accounting for the row number in the airflow direction
	double cfmA,cfmB,cfmC;
	int Microfin;// B.S. for judging if it is micro-fin tube
	double finN, gama, beta, finH;//B.S. micro-fin geometry
	double w_b, w_e, w_z; //B.S. micro-fin geometry
	double D_b, D_m; //B.S. micro-fin geometry
	double Acs,P_H, Dh_i;//B.S. micro-fin geometry
	int airFin; //B.S. airside fin type
	double sub1,sub2,sub3,sub4,sub5;//for inputing sub-structures of fin surface
	double P_l; //airside fin geometry
	double K_T, K_F;// conductance factor
	double Ga;//airside velocity
	double airT;//air temperature
	double Frontal_A; //frontal area
	double GetP;//calculate the airside pressure drop
	
	//variables for generating the simple condenser function
	double V_TOT, V_TP, V_Liq, V_Vap;//B.S., inner volume of different phase
	double L_TOT,LiqL,VapL,TPL;//B.S., tube length of different phase
	double A_TOT,A_Liq,A_Vap,A_TP;//heat tranfer surface area of different phase
	
	double m_TOT, m_TP, m_Liq, m_Vap;//air mass flow rate across different phase region
	double rho_TOT, rho_TP, rho_Liq, rho_Vap;//average density of different phases
	
	double U_TP, U_Liq, U_Vap;//average heat transfer conductance per tube length
	double DP_TOT, DP_TP, DP_Liq, DP_Vap;//average pressure drop gradient of different phase
	double UA_TOT, UA_Liq,UA_Vap,UA_TP;//heat transfer conductance of diffferent phase
	
	double Tai;//inlet air temperature
	double mr, ma_TOT;//refrigerant and air mass flow rate
	
	HP HP_in, HP_TP1, HP_TP2, HP_out;//state parameters at important locations
	double Ga_meanL;//averager mass flow rate per tube length

	double r_v, r_l,r_tp;//parameters for adjusting the theoritical heat transfer effectivenss of different phase
	double count1,count2;//count the state points of two-phase flow end and beginning
	double H1_residual, H2_residual;//for consistency of the moving boundary model analysis
	//------------------------------B.S.
	
	//Lumped model included above---------------------------B.S. 

	double cfma; // cfm air

};

void Build_Cond_MovBound(CGP*);//fuction for building the simplified condenser model, phase by phase and lumped model
void Build_Cond_Lumped(CGP*);//fuction for building the lumped condenser model
double Condenser_Simple(double,HP,double,HP*,MASS*,CGP*);//for conducting moving boundary model analysis
double Condenser_Lumped(double,HP,double,HP*,MASS*,CGP*);//for conducting lumped model analysis

int CondCircuit(int,double,HP,double,double,HP*,double*,MASS*,CGP*);

void CondTube(double,HP*,double,double,double*,MASS*,CGP*);
void CondMan(int,int,double,double*,HP*);
void CondReturnBend(double,HP*,MASS*,CGP*);
void CondTubeL_new(double,HP*,double,double,double*,MASS*,CGP*);

int Condenser(const char *filename, double,HP,double,double,HP*,double*,MASS*,CGP*,double*,int NSeg=-1);
void Condenser0(double,const char *flename,double,HP,double,double,HP*,double*,MASS*,double*);

#endif

