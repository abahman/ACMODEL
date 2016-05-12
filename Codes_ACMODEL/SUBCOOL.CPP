#include "stdafx.h"

#include "subcool.h"
#include "cmine.h"
#include "getdata.h"
#include "ref.h"
#include "errorlog.h"
#include "corr.h"

extern ErrorLog errorLog;

//#define DEBUG_SUBCOOL

int subcoolInit = 1;

// Subcooler transfers heat between two refrigerant streams.
// The inlet state and flow rate of stream 1 are given and
// the outlet state and flow rate of stream 2 are provided.
// The function provides the outlet state of stream 2, the
// inlet state of stream 1, and the mass of refrigerant in the
// device.
//
// INPUTS
// mr1 - mass flow rate of stream 1
// HP1i - inlet state of stream 1
// mr2 - mass flow rate of stream 2
// HP2o - outlet state of stream 2
void
Subcooler(const char* filename,double mr1,HP HP1i,double mr2,HP HP2o,HP& HP1o,HP& HP2i,MASS& m)
{
	// Load model parameters from parameter file
	static double UA = 300;
	static double V1total = 0.15;
	static double V2total = 0.15;
	static int nPasses = 30;

	if(subcoolInit) {

		#ifdef DEBUG_SUBCOOL
		{
		FILE *fp = fopen("subcool.dbg","w");
		fprintf(fp,"Opening...\n");
		fclose(fp);
		}
		#endif

		FILE *fp = fopen(filename,"r") ;
		int type = PosComponentData(fp,SUBCOOLER,1);
		if(type!=501) {
			errorLog.Add("Subcooler","parameters not found");
			return;
		}
		UA = GetDouble(fp);
		V1total = GetDouble(fp);
		V2total = GetDouble(fp);
		nPasses = (int)GetDouble(fp);
		fclose(fp);
		subcoolInit=0;
	}

	#ifdef DEBUG_SUBCOOL
	{
	FILE *fp = fopen("subcool.dbg","a");
	fprintf(fp,"HP2o=(%lf,%lf)\n",HP2o.H,HP2o.P);
	fclose(fp);
	}
	#endif

	const double Rseg = (double)nPasses/UA;
	const double V1pass = V1total/(double)nPasses;
	const double V2pass = V2total/(double)nPasses;

	// no pressure drop
	HP1o.P = HP1i.P;
	HP2i.P = HP2o.P;

	m.V = m.m = 0;
	for(int i=0;i<nPasses;i++) {

		// material properties
		TXP TXP1i = HPtoTXP(HP1i);
		if(errorLog.IsError()) {
			errorLog.Add("Subcooler","TXP1i");
			return;
		}

		TXP TXP2o = HPtoTXP(HP2o);
		if(errorLog.IsError()) {
			char str[128];
			sprintf(str,"TXP2o HP2o=(%lf,%lf)\n",HP2o.H,HP2o.P);
			errorLog.Add("Subcooler",str);
			return;
		}

		#ifdef DEBUG_SUBCOOL
		{
		FILE *fp = fopen("subcool.dbg","a");
		fprintf(fp,"HP2o=(%lf,%lf)\n",HP2o.H,HP2o.P);
		fclose(fp);
		}
		#endif

		// Stream 2 is computed backwards.  Therefore, if conditions allow for
		// too much heat transfer, then the inlet of stream 2 can be low enough
		// enthalpy to be off the table.  Stop going through the segments, if
		// enough segments have been computed to saturate the amount of heat the
		// stream 2 can absorb.
		if(TXP2o.T<TMIN+40) {
			errorLog.Add("Subcooler","TXP2o.T<TMIN+40");
			return;
		}

		const double Cmine = CmineCounterFlow(Rseg,mr1,mr2,TXP1i,TXP2o);
		if(errorLog.IsError()) {
			errorLog.Add("Subcooler","Cmine");
			return;
		}

		double C;
		if(IsTwoPhase(TXP2o.X)) {
			C = Cmine;
		} else {
			const double Cp2 = PropertyTXPtr(SPEC,TXP2o);
			if(errorLog.IsError()) {
				char str[128];
				sprintf(str,"Cp2=%lf TXP=(%lf,%lf,%lf)",Cp2,TXP2o.T,TXP2o.X,TXP2o.P);
				errorLog.Add("Subcooler",str);
				return;
			}
			const double C2 = mr2*Cp2;
			C = C2*Cmine/(C2-Cmine);
		}

		const double q = (TXP1i.T-TXP2o.T)*C;

		HP1i.H -= q/mr1;
		HP2o.H -= q/mr2;

		const double v1 = PropertyTXPth(VOL,TXP1i);   // homogeneous flow model
		if(errorLog.IsError()) {
			errorLog.Add("Subcooler","v1");
			return;
		}
		m.V += V1pass;
		m.m += V1pass/v1;

		const double v2 = PropertyTXPth(VOL,TXP2o);   // homogeneous flow model
		if(errorLog.IsError()) {
			errorLog.Add("Subcooler","v2");
			return;
		}
		m.V += V2pass;
		m.m += V2pass/v2;
	}

	HP1o.H = HP1i.H;
	HP2i.H = HP2o.H;
}
