#include "stdafx.h"

#include "compstat.h"
#include "ref.h"
#include "r22.h"
#include "numeric.h"
#include "math.h"
#include "getdata.h"
#include "errorlog.h"
#include "spline1d.h"


extern ErrorLog errorLog;

//#define DEBUG_COMPSTAT
//#define DEBUG_COMPCCFUNC

int compCCLowInit = 1;
int compCCHighInit = 1;

// Local data structures
struct CompCCP {
	HP HP1;
	double v1,P2,np,npAdj;
	double H2,n;
};

// Local functions
double CompCCFunc(double v2a,void *Params);

// INPUTS
//
// H1 - inlet enthalpy
// P1 - inlet pressure
// P2 - intermediate pressure
// P3 - discharge pressure
//
// OUTPUTS
//
// H2 - outlet enthalpy
// mr - mass flow rate
// Ei - electric power consumption
// m - total rerigerant mass
void
CompressorCCLowStage(const char* filename,double H1,double P1,double P2,double P3,double& H2,double& mr,double& Ei,MASS& m,double* Prms)
{
	static CompCCP Q;

	// Compressor parameters from data file
	static double powerP[11]={1,1,1,1,1,1,1,1,1,1,1};
	static double massFlowP[11]={1,1,1,1,1,1,1,1,1,1,1};
	static double adj[3]={1,1,1};

	if(compCCLowInit) {

		#ifdef DEBUG_COMPSTAT
		{
		FILE *fp = fopen("complow.dbg","w");
		fprintf(fp,"Opening CompressorCCLowStage...\n");
		fclose(fp);
		}
		#endif

		FILE *fp = fopen(filename,"r");
		int type = PosComponentData(fp,COMPRESSOR,1);
		if(type==0 || type!=102) {
			errorLog.Add("CompressorCCLowStage","Compressor parameters not found");
			return;
		}
		Q.np = GetDouble(fp);
		for(int i=0;i<11;i++) powerP[i] = GetDouble(fp);
		for(int i=0;i<11;i++) massFlowP[i] = GetDouble(fp);
		for(int i=0;i<3;i++)  adj[i] = GetDouble(fp);
		fclose(fp);

		compCCLowInit=0;
	}

	double nvAdj = adj[0]*Prms[0];
	double neAdj = adj[1]*Prms[1];
	Q.npAdj = adj[2]*Prms[2];

	/* Set up parameters needed to solve for v2 */
	Q.HP1.H = H1;
	Q.HP1.P = P1;
	Q.P2 = P2;
	TXP TXP1 = HPtoTXP(Q.HP1);
	if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"(TXP1) H=%lf P=%lf",Q.HP1.H,Q.HP1.P);
		errorLog.Add("CompressorCCLowStage",str);
		return;
	}
	Q.v1 = PropertyTXPth(VOL,TXP1);
	if(errorLog.IsError()) {
		errorLog.Add("CompressorCCLowStage","Q.v1");
		return;
	}

	// Bound v2 between v2min and v2max
	//
	// v2max is the largest specific volume for P2.
	double v2max = PropertyTXPth(VOL,toTXP(TMAX,1,P2));
	if(errorLog.IsError()) {
		errorLog.Add("CompressorCCLowStage","v2max");
		return;
	}

	// v2min is the smallest specific volume for P2.
	double v2min = PropertyTXPth(VOL,HPtoTXP(toHP(H1,P2)));
	if(errorLog.IsError()) {
		errorLog.Add("CompressorCCLowStage","v2min");
		return;
	}

	/* Iteratively solve for v2 */
	Zbrent(v2min,v2max,CompCCFunc,1e-7,&Q);
	if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"v2 HTXP1=(%lf,%lf,%lf,%lf) P2=%lf P3=%lf",H1,TXP1.T,TXP1.X,TXP1.P,P2,P3);
		errorLog.Add("CompressorCCLowStage",str);
		return;
	}

	// output state
	H2=Q.H2;

	double Psuc = P1;
	double Pdisc = P3;

	// mass flow rate
	double* map = massFlowP;
	mr =  map[0]+
	      Psuc*(map[1]+Psuc*(map[2]+Psuc*map[3]))+
	      Pdisc*(map[4]+Pdisc*(map[5]+Pdisc*map[6]))+
	      Psuc*Pdisc*(map[7]+Pdisc*map[8]+Psuc*map[9]+Psuc*Pdisc*map[10]);
	mr *= nvAdj;

	// superheat correction
	TXP TXP_prop={0,0,0};
	TXP_prop.P = P1;
	TXP_prop.X=1;

	const double tsat_std = PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(P1);
	if(errorLog.IsError()) {
		errorLog.Add("CompressorCCHighStage","tsat_std");
		return;
	}
	TXP txp_std = toTXP(tsat_std+11.1,1,P1);   // 20F standard superheat
	double v_std = PropertyTXPth(VOL,txp_std);
	if(errorLog.IsError()) {
		errorLog.Add("CompressorCCHighStage","v_std");
		return;
	}
	mr *= v_std/Q.v1;

	// electric power consumption
	map = powerP;
	Ei =  map[0]+
	      Psuc*(map[1]+Psuc*(map[2]+Psuc*map[3]))+
	      Pdisc*(map[4]+Pdisc*(map[5]+Pdisc*map[6]))+
	      Psuc*Pdisc*(map[7]+Pdisc*map[8]+Psuc*map[9]+Psuc*Pdisc*map[10]);
	Ei *= neAdj;

	// charge mass
	const double h=0.3302,d=0.2183,fv=0.7;
	m.V = fv*h*acos((-1)/4.0*d*d);
	m.m = m.V/Q.v1;

	#ifdef DEBUG_COMPSTAT
	{
	FILE *fp = fopen("complow.dbg","a");
	fprintf(fp,"HTXP1=(%lf,%lf,%lf,%lf) P2=%lf P3=%lf v1=%lf T1_std=%lf v_std=%lf\n",1e-3*H1,TXP1.T,TXP1.X,TXP1.P,P2,P3,Q.v1,txp_std.T,v_std);
	fclose(fp);
	}
	#endif
}


void
CompressorCCHighStage(const char* filename,double H2,double P2,double P3,double& H3)
{
	static CompCCP Q;

	//
	// Compressor parameters from data file
	static double powerP[11]={1,1,1,1,1,1,1,1,1,1,1};
	static double massFlowP[11]={1,1,1,1,1,1,1,1,1,1,1};
	static double adj[3]={1,1,1};
	//
	if(compCCHighInit) {

		#ifdef DEBUG_COMPSTAT
		{
		FILE *fp = fopen("comphigh.dbg","w");
		fprintf(fp,"Opening CompressorCCHighStage...\n");
		fclose(fp);
		}
		#endif

		FILE *fp = fopen(filename,"r");
		int type = PosComponentData(fp,COMPRESSOR,1);
		if(type==0 || type!=102) {
			errorLog.Add("CompressorCCHighStage","Compressor parameters not found");
			return;
		}
		Q.np = GetDouble(fp);
		for(int i=0;i<11;i++) powerP[i] = GetDouble(fp);
		for(int i=0;i<11;i++) massFlowP[i] = GetDouble(fp);
		for(int i=0;i<3;i++)  adj[i] = GetDouble(fp);
		fclose(fp);

		compCCHighInit=0;
	}

	// Set up parameters needed to solve for v2
	Q.npAdj = adj[2];
	Q.HP1.H = H2;
	Q.HP1.P = P2;
	Q.P2 = P3;
	TXP TXP1 = HPtoTXP(Q.HP1);
	if(errorLog.IsError()) {
		char str[64];
		sprintf(str,"(TXP1) H=%lf P=%lf\n",Q.HP1.H,Q.HP1.P);
		errorLog.Add("CompressorCCHighStage",str);
		return;
	}
	Q.v1 = PropertyTXPth(VOL,TXP1);
	if(errorLog.IsError()) {
		errorLog.Add("CompressorCCHighStage","1");
		return;
	}

	// Bound v2 between v2min and v2max
	double v2max = PropertyTXPth(VOL,toTXP(TMAX,1,P3));
	if(errorLog.IsError()) {
		errorLog.Add("CompressorCCHighStage","2");
		return;
	}
	HP HPmin = toHP(H2,P3);
	double v2min = PropertyTXPth(VOL,HPtoTXP(HPmin));
	if(errorLog.IsError()) {
		errorLog.Add("CompressorCCHighStage","3");
		return;
	}

	// Iteratively solve for v2
	const double v2 = Zbrent(v2min,v2max,CompCCFunc,1e-7,&Q);
	if(errorLog.IsError()) {
		errorLog.Add("CompressorCCHighStage","4");
		return;
	}

	// outlet state
	H3 = Q.H2;

	#ifdef DEBUG_COMPSTAT
	{
	FILE *fp = fopen("comphigh.dbg","a");
	fprintf(fp,"H2=%lf P2=%lf P3=%lf H3=%lf v2min=%lf v2max=%lf v2=%lf\n",H2,P2,P3,H3,v2min,v2max,v2);
	fclose(fp);
	}
	#endif
}

double
CompCCFunc(double v2a,void *Params)
{
	CompCCP* P = (CompCCP*)Params;

	const double np = P->np*P->npAdj;

	P->n = log(P->HP1.P/P->P2)/log(v2a/P->v1);

	const double rn = (P->n-1)/P->n;
	const double rP = P->P2/P->HP1.P;
	if(rP<0) {
		errorLog.Add("CompCCFunc","rP<0.0");
		return -1;
	}

	HP HP2;
	HP2.H = P->H2 = P->HP1.H+1.0e3*P->HP1.P*P->v1*(pow(rP,rn)-1)/(rn*np);
	HP2.P = P->P2;

	TXP TXP2 = HPtoTXP(HP2);
	if(errorLog.IsError()) {
		char str[64];
		sprintf(str,"TXP2 (HP2=%lf,%lf)\n",1e-3*HP2.H,HP2.P);
		errorLog.Add("CompCCFunc",str);
		return -1;
	}

	const double v2b = PropertyTXPth(VOL,TXP2);
	if(errorLog.IsError()) {
		errorLog.Add("CompCCFunc","v2b");
		return -1;
	}

	const double dv = 2.0*(v2a-v2b)/(v2a+v2b);

	#ifdef DEBUG_COMPCCFUNC
	{
	static int init = 1;
	FILE *fp;
	if(init) {
		fp = fopen("compfunc.dbg","w");
		init = 0 ;
	} else {
		fp = fopen("compfunc.dbg","a");
	}
	if(fp) {
		fprintf(fp,"P1=%lf P2=%lf v2a=%lf dv=%lf\n",P->HP1.P,P->P2,v2a,dv);
		fclose(fp);
	}
	}
        #endif

	return dv;
}
