#include "stdafx.h"

#include <stdio.h>
#include <math.h>
#include <string.h>

#include "r22.h"
#include "wair.h"
#include "mass.h"
#include "comp.h"
#include "lines.h"
#include "evap.h"
#include "cond.h"
#include "expdev.h"
#include "subcool.h"
#include "compstat.h"
#include "getdata.h"
#include "errorlog.h"
#include "defs.h"
//#include "dataelem.h"
#include "spline1d.h"
#include "dataST.h"
#include "conversion.h"
#include "mainfunc.h"
#include "compu.h"
#include "outelem.h"
#include "spline2d.h"

// Bo Shen added a new short tube expansion device model.  Include this #define
// to use the new model.  Otherwise, the ExpDev model will be used.
#define USE_SHORT_TUBE_MODEL

extern ErrorLog errorLog;

//extern DataElement gData;


extern OutputDataElement gOutData;

extern WetAirPropTR wair;

//#define EXTRA_REPORTING

/*********************************************************************
ACTune is a program that solves the cycle the easy way first.  

It starts with known driving conditions:
35.14   condenser inlet air temperature (C)
3.74    condenser air mass flux (kg/m^2/s)
21.7    evaporator inlet air drybulb temperature (C)
0.421   evaporator inlet air relative humidity (NU)
4.48    evaporator air mass flux (kg/m^2/s)

and known states:
1980    head pressure (kPa)
573     suction pressure (kPa)
7.95    superheat (C)
0.76    subcooling (C)

and solves for the variables forming the solution in the hard solver.

1. difference enthalpy between condenser outlet and measured
2. difference enthalpy between evaporator inlet and measured
3. pressure difference between condenser outlet and expdev inlet

Adjust the hx UAs and orifice sizes and set the charge to the specified
amount to get a known starting point for the hard solver.  Then the hard
solver, starting frm this known point, solves for how this machine will 
operate under different driving conditions and with faults (e.g. fouled
hx) imposed.

****

Main program for ACTUNE.  This program tunes model parameters so
that predictions match measured values.

ARGUMENTS
modelParamsFilename = Input file containing cycle parameters.
measurementsFilename = Input file name defining measured conditions.
outputFilename = Output file name with parameter fitting results.

RETURNS
0 = successful
1 = error

NOTES
View file [error.log] for detailed error messages.
*********************************************************************/
int ActuneRTU(const char* modelParamsFilename, const char* measurementsFilename, const char* outputFilename, int REV)
{
	//HP HPs[5];
	//MASS m[9];
	TP TPi;

	if(errorLog.IsError()) errorLog.ClearError("ActuneRTU");

	// Get input file containing data to tune cyle too.
	FILE *fp = fopen(measurementsFilename,"r");
	if(fp==NULL) {
		errorLog.Add("Can not open input file with measurements!");
		return 1;
	}

	// Driving conditions
	const int CORRECT = GetInt(fp);
	const double Tai = GetDouble(fp);			// ambient temperature
	const double Gca = GetDouble(fp);			// cond air mass flux
	TPi.T = GetDouble(fp);						// return air drybulb temperature
	TPi.P = GetDouble(fp);						// return air relative humidity
	const double Gea = GetDouble(fp);			// evap air mass flux

	// Measured cycle states
	const double P0 = GetDouble(fp);			// head pressure - must be head pressure for compressor model input
	const double P4 = GetDouble(fp)-0;			// suction pressure - suction line pressure drop B.S
	const double Tsh = GetDouble(fp);			// superheat
	const double Tsc = GetDouble(fp);			// subcooling
	fclose(fp);

	// Calculate the suction line state (4) from the input pressure and superheat
	// Tevap is the saturation temperature based on the suction line pressure

	TXP TXP_prop={0,0,0};
	TXP_prop.P=P4;
	TXP_prop.X=1;

	const double Tevap = PropertyTXPth(TSAT,TXP_prop);//reftpvthP.Tsat(P4);
	if(errorLog.IsError()) {
		errorLog.Add("ActuneRTU","Tevap");
		return 2;
	}

	const HP HP4 = TXPtoHP(toTXP(Tevap+Tsh,1,P4));
	if(errorLog.IsError()) {
		errorLog.Add("ActuneRTU","HP4");
		return 3;
	}

	// Run MainFunc1() to evaluate the refrigeration cycle
	// Run MainFunc1() to evaluate the refrigeration cycle
	// Run MainFunc1() to evaluate the refrigeration cycle

	InVars I;
	I.Tai = Tai;
	I.Goa = Gca;
	I.Gra = Gea;
	I.Charge = 2.26;
	I.TPi = TPi;
	I.Dadj = 1;
	I.nvAdj = 1;
	strcpy(I.filename,modelParamsFilename);

	double X[3],Y[3];
	X[0] = HP4.H;
	X[1] = P4;
	X[2] = P0;

	InitializeComponents();
	int errorCode;

	I.REV_CAL=REV;
	if(errorCode=MainFunc1(X,Y,&I)) {
		printf("MainFunc1: errorCode=%d\n",errorCode);

		fp = fopen(outputFilename,"a");
		fprintf(fp,"<P>Can not solve model with selected inputs. ErrorCode=%d.\n",errorCode);
		fclose(fp);

		errorLog.Add("ActuneRTU","MainFunc1");

		return 5;
	}
	
	const double mr = gOutData.mr[1];
	const double P1 = gOutData.P1;
	const double Hcondout = gOutData.H[4];
	const double Pcondout = gOutData.P[4];
	const double Hevapin = gOutData.H[7];
	HP HP0;
	HP0.H = gOutData.H[2];
	HP0.P = gOutData.P[2];

	const double Charge = gOutData.Charge;
	const double Tco = gOutData.Tco;
	const double Teo = gOutData.Teo;

#ifdef USE_SHORT_TUBE_MODEL
	const double mr_exp = gOutData.mr_exp;
#endif

	// The measured state in the liquid line is determined by the measured liquid line temperature
	// and the pressure calculated by the model.  Do not use the head pressure specified in the input file.
	// Not considering the condenser pressure drop introduces a large tuning error.
	// Tcond is the saturation temperature based on the liquid line pressure.

	TXP_prop.P=Pcondout;
	TXP_prop.X=0;
	const double Tcond = PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(Pcondout);
	if(errorLog.IsError()) {
		errorLog.Add("ActuneRTU","Tcond");
		return 4;
	}
	HP HP1 = TXPtoHP(toTXP(Tcond-Tsc,0,Pcondout));//B.S., this struct variable will be used for iterating to get the evaporator heat transfer coefficient
	if(errorLog.IsError()) {
		errorLog.Add("ActuneRTU","HP1");
		return 6;
	}

	// diff enth cond out and measured
	const double ce = (Hcondout-HP1.H)/HP1.H;
	// diff enth evap in and measured
	const double ee = (Hevapin-HP1.H)/HP1.H;
	// diff pres cond outlet expdev+restrict inlet
#ifdef USE_SHORT_TUBE_MODEL
	const double pe = 2*(mr-mr_exp)/(mr+mr_exp);
#else
	const double pe = 2*(P1-Pcondout)/(P1+Pcondout);
#endif

	// print report
	// print report
	// print report

	fp = fopen(outputFilename,"w");
	fprintf(fp,"<TABLE border=\"1\">\n");
	fprintf(fp,"<CAPTION><B>Input data</B></CAPTION>\n");
	fprintf(fp,"\t<TR>\n");
	fprintf(fp,"\t\t<TD>Tai = %.1lfC (%.1lfF)</TD><TD>Gca = %.4lf kg/s/m^2</TD>\n",Tai,CtoF(Tai),Gca);
	fprintf(fp,"\t</TR>\n");
	fprintf(fp,"\t<TR>\n");
	const double wetBulb = wair.WetBulb(TPi.T,TPi.P);
	fprintf(fp,"\t\t<TD>TPi.T = %.1lfC (%.1lfF), TPi.P = %.3lf (%.1lf%%), WetBulb = %.3lf</TD><TD>Gea= %.4lf kg/s/m^2</TD>\n",TPi.T,CtoF(TPi.T),TPi.P,100*TPi.P,CtoF(wetBulb),Gea);
	fprintf(fp,"\t</TR>\n");
	fprintf(fp,"\t<TR>\n");
	fprintf(fp,"\t\t<TD>Ph = %.0lfkPa (%.1lfpsig)</TD><TD>Tcond(ll) = %.1lfC (%.1lfF) Tctoa = %.1lfC (%.1lfF)</TD>\n",P0,KPAtoPSIG(P0),Tcond,CtoF(Tcond),Tcond-Tai,DCtoDF(Tcond-Tai));
	fprintf(fp,"\t</TR>\n");
	fprintf(fp,"\t<TR>\n");
	fprintf(fp,"\t\t<TD>Ps = %.0lfkPa (%.1lfpsig)</TD><TD>Tevap(sl) = %.1lfC (%.1lfF)</TD>\n",P4,KPAtoPSIG(P4),Tevap,CtoF(Tevap));
	fprintf(fp,"\t</TR>\n");
	fprintf(fp,"\t<TR>\n");
	fprintf(fp,"\t\t<TD>Tsh = %.1lfC (%.1lfF)</TD><TD>Tsc = %.1lfC (%.1lfF)</TD>\n",Tsh,DCtoDF(Tsh),Tsc,DCtoDF(Tsc));
	fprintf(fp,"\t</TR>\n");
	fprintf(fp,"</TABLE>\n");
	fprintf(fp,"<P>View: <A HREF=\"ErrorLog.htm\">Error Log File</A>\n");
	fprintf(fp,"<P>View: <A HREF=\"./OutputDoc/ACModelOutput.htm\">ACModel Output File</A>\n");
	fclose(fp);

	fp = fopen(outputFilename,"a");
	fprintf(fp,"<P>\n");
	fprintf(fp,"<TABLE border=\"1\">\n");
	fprintf(fp,"<CAPTION><B>Results</B></CAPTION>\n");
	fprintf(fp,"\t<TR>\n");
	fprintf(fp,"\t\t<TD COLSPAN=3><B>Air flow</B></TD>\n");
	fprintf(fp,"\t</TR>\n");
	fprintf(fp,"\t<TR>\n");
	fprintf(fp,"\t\t<TD>CFM air cond/evap</TD><TD>%.0lf (cfm)</TD><TD>%.0lf (cfm)</TD>\n",gOutData.cfmca,gOutData.cfmea);
	fprintf(fp,"\t</TR>\n");
	fprintf(fp,"\t<TR>\n");
	fprintf(fp,"\t\t<TD COLSPAN=3><B>Charge</B></TD>\n");
	fprintf(fp,"\t</TR>\n");
	fprintf(fp,"\t<TR>\n");
	fprintf(fp,"\t\t<TD>Mass</TD><TD>%lf (kg)</TD><TD>&nbsp;</TD>\n",Charge);
	fprintf(fp,"\t</TR>\n");
	fprintf(fp,"\t<TR>\n");
	fprintf(fp,"\t\t<TD COLSPAN=3><B>Liquid line enthalpy</B></TD>\n");
	fprintf(fp,"\t</TR>\n");
	fprintf(fp,"\t<TR>\n");
	fprintf(fp,"\t\t<TD>measured</TD><TD>%lf (kJ/kg)</TD><TD>&nbsp;</TD>\n",HP1.H/1000);
	fprintf(fp,"\t</TR>\n");
	fprintf(fp,"\t<TR>\n");
	fprintf(fp,"\t\t<TD>condenser model prediction</TD><TD>%lf (kJ/kg)</TD><TD>error = %.4f%%</TD>\n",Hcondout/1000,100*ce);
	fprintf(fp,"\t</TR>\n");
	fprintf(fp,"\t<TR>\n");
	fprintf(fp,"\t\t<TD>evaporator model prediction</TD><TD>%lf (kJ/kg)</TD><TD>error = %.4f%%</TD>\n",Hevapin/1000,100*ee);
	fprintf(fp,"\t</TR>\n");
#ifdef USE_SHORT_TUBE_MODEL
	fprintf(fp,"\t<TR>\n");
	fprintf(fp,"\t\t<TD COLSPAN=3><B>Refrigerant mass flow rate</B></TD>\n");
	fprintf(fp,"\t</TR>\n");
	fprintf(fp,"\t<TR>\n");
	fprintf(fp,"\t\t<TD>compressor mass flow rate</TD><TD>%lf (kg/s)</TD><TD>&nbsp;</TD>\n",mr);
	fprintf(fp,"\t</TR>\n");
	fprintf(fp,"\t<TR>\n");
	fprintf(fp,"\t\t<TD>expansion device mass flow rate</TD><TD>%lf (kg/s)</TD><TD>error = %.4f%%</TD>\n",mr_exp,100*pe);
	fprintf(fp,"\t</TR>\n");
#else
	fprintf(fp,"\t<TR>\n");
	fprintf(fp,"\t\t<TD COLSPAN=3><B>Liquid line pressure</B></TD>\n");
	fprintf(fp,"\t</TR>\n");
	fprintf(fp,"\t<TR>\n");
	fprintf(fp,"\t\t<TD>measured head pressure</TD><TD>%lf (kJ/kg)</TD><TD>&nbsp;</TD>\n",P0);
	fprintf(fp,"\t</TR>\n");
	fprintf(fp,"\t<TR>\n");
	fprintf(fp,"\t\t<TD>condenser model prediction</TD><TD>%lf (kPa)</TD><TD>&nbsp;</TD>\n",Pcondout);
	fprintf(fp,"\t</TR>\n");
	fprintf(fp,"\t<TR>\n");
	fprintf(fp,"\t\t<TD>expdev+restrition model predictions</TD><TD>%lf (kPa)</TD><TD>error = %.4f%%</TD>\n",P1,100*pe);
	fprintf(fp,"\t</TR>\n");
#endif
	fprintf(fp,"\t<TR>\n");
	fprintf(fp,"\t\t<TD COLSPAN=3><B>Air temperature differences</B></TD>\n");
	fprintf(fp,"\t</TR>\n");
	fprintf(fp,"\t<TR>\n");
	fprintf(fp,"\t\t<TD>evaporator delta T</TD><TD>%.1lfC (%.1lfF)</TD><TD>&nbsp;</TD>\n",TPi.T-Teo,DCtoDF(TPi.T-Teo));
	fprintf(fp,"\t</TR>\n");
	fprintf(fp,"\t<TR>\n");
	fprintf(fp,"\t\t<TD>condenser delta T</TD><TD>%.1lfC (%.1lfF)</TD><TD>&nbsp;</TD>\n",Tco-Tai,DCtoDF(Tco-Tai));
	fprintf(fp,"\t</TR>\n");
	fprintf(fp,"\t<TR>\n");
	fprintf(fp,"\t\t<TD COLSPAN=3><B>Capacity & Efficiency</B></TD>\n");
	fprintf(fp,"\t</TR>\n");
	fprintf(fp,"\t<TR>\n");
	fprintf(fp,"\t\t<TD>Capacity</TD><TD>%.4lf (tons)</TD><TD>&nbsp;</TD>\n",WtoBTUH(-1*gOutData.Qe)/12000);
	fprintf(fp,"\t</TR>\n");
	fprintf(fp,"\t<TR>\n");
	fprintf(fp,"\t\t<TD>EER</TD><TD>%.4lf (BTU/hr/kW)</TD><TD>&nbsp;</TD>\n",WtoBTUH((-1*gOutData.Qe/gOutData.Ecmp)));
	fprintf(fp,"\t</TR>\n");
	fprintf(fp,"</TABLE>\n");

	fprintf(fp,"<P>\n");
	fprintf(fp,"<B>Analysis</B>\n");
	fprintf(fp,"<P>\n");
	fprintf(fp,"Use Charge = %lf (kg)\n",Charge);
	fclose(fp);
	
	fp = fopen(outputFilename,"a");
	fprintf(fp,"<P>\n");
	if(ce>0)
		fprintf(fp,"Increase condenser capacity by %.4lf%%.\n",100.0*fabs(ce));
	else
		fprintf(fp,"Decrease condenser capacity by %.4lf%%.\n",100.0*fabs(ce));
	fclose(fp);

	errorLog.ClearError("ActuneRTU");

	if(CORRECT) {
		// Find the same value of hAirAdj and hRefAdj to match measurements
		MASS m;
		double Tao;
	//	HP HP1;//B.S. this struct variable should be a global variable, it can't be defined here
		double Prms[10];
		for(int i=0;i<10;i++) Prms[i]=1;

		Condenser0(Tsc,modelParamsFilename,mr,HP0,Tai,Gca,&HP1,&Tao,&m,Prms);

		fp = fopen(outputFilename,"a");
		fprintf(fp,"<P>\n");
		if(errorLog.IsError()) {
			errorLog.ClearError("ActuneRTU");
			fprintf(fp,"ERROR! Can not solve for HAirAdj in cond model.\n");
		} else {
			fprintf(fp,"To correct: multiply HAirAdj in cond model by [%.4lf].\n",Prms[0]);
		}
		fclose(fp);

		/*
		Find value of Ga to match measurements
		for(i=0;i<5;i++) Prms[i]=1; status=0;
		Condenser0(3,Tsc[0],mr,HPs[0],Tai[0],-1,&HPs[1],&Tao,&m[2],Prms,&status);
		if(status) {sprintf(s,"main: \n"); bprintf(s); return 0;}
		printf("To correct: multiply GaNom   in cond model by: %lf\n",Prms[3]);
		*/
	}

	fp = fopen(outputFilename,"a");
	fprintf(fp,"<P>\n");
	if(ee>0)
		fprintf(fp,"Increase evaporator capacity by %.4lf%%.\n",100.0*fabs(ee));
	else
		fprintf(fp,"Decrease evaporator capacity by %.4lf%%.\n",100.8*fabs(ee));
	fclose(fp);
	
	HP HP3;

	if(CORRECT) {
		// Find value of hAirAdj to match measurements
		double Aflow;
		MASS m;
		TP TPo;
	//	HP HP3;
		double Prms[10];
		for(int i=0;i<10;i++) Prms[i]=1;

		if(errorLog.IsError()) errorLog.ClearError("ActuneRTU");
		Evaporator0(0,HP1.H,modelParamsFilename,mr,HP4,Gea,TPi,&HP3,&TPo,&m,&Aflow,Prms,REV);

		fp = fopen(outputFilename,"a");
		fprintf(fp,"<P>\n");
		if(errorLog.IsError()) {
			errorLog.ClearError("ActuneRTU");
			fprintf(fp,"ERROR! Can not solve for HAirAdj in evap model.\n");
		} else {
			fprintf(fp,"To correct: multiply HAirAdj in evap model by [%.4lf].\n",Prms[0]);
		}
		fclose(fp);

		/*
		Find value of Ga to match measurements
		for(i=0;i<5;i++) Prms[i]=1; status=0;
		Evaporator0(4,HP1.H,filename,mr/4,HPs[4],-1,TPi,&HPs[3],&TPo,&m[7],&Aflow,Prms,&status);
		if(status) {sprintf(s,"main: \n"); bprintf(s); return 0;}
		printf("To correct: multiply GaNom   in evap model by: %lf\n",Prms[4]);
		*/
	}

	fp = fopen(outputFilename,"a");
	fprintf(fp,"<P>\n");
	if(pe>0)
		fprintf(fp,"Decrease expansion device diameter by %.4lf%%.\n",100.0*fabs(pe));
	else
		fprintf(fp,"Increase expansion device diameter by %.4lf%%.\n",100.0*fabs(pe));
	fclose(fp);

	if(CORRECT) {//no need to adjust the expansion device for TxV units
		// Find value of Dadj to match measurements 
		double Prms[10];
		for(int i=0;i<10;i++) Prms[i]=1;

		// This solution for the diameter assumes the heat exchanger models are perfectly
		// tuned to provide the measured subcooling.  This may not be the right diameter for the
		// actual hx values as they are now.

		if(errorLog.IsError()) errorLog.ClearError("ActuneRTU");

		HP1.P = HP1.P-35;//B.S., include pressure drop in the liquid line

#ifdef USE_SHORT_TUBE_MODEL
		Short_Tube_0(modelParamsFilename,HP1,HP3,mr/4,Prms);
		//TxV_0(modelParamsFilename,HP1,HP3,Tsh,mr,Prms);
#else
		double P2;
		MASS m;
		ExpDev0(0,HP1.P,modelParamsFilename,1,HP1.H,P4,mr/4,&P2,&m,Prms);
#endif

		fp = fopen(outputFilename,"a");
		fprintf(fp,"<P>\n");
		if(errorLog.IsError()) {
			errorLog.ClearError("ActuneRTU");
			fprintf(fp,"ERROR! Can not solve for Diameter in exp dev model.\n");
		} else {
			fprintf(fp,"To correct: multiply Diameter in exp dev model by [%.4lf].\n",Prms[0]);
		}
		fclose(fp);
	}

	// Wite input file for ACModel
	// Wite input file for ACModel
	// Wite input file for ACModel

/*	fp = fopen("./InputDoc/master.in","w");
	if(fp) {
		fprintf(fp,"1\t%lf\t%lf\t%lf\n",X[0],X[1],X[2]);
		fprintf(fp,"4\t1\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t1.0\t1.0\n",Tai,TPi.T,TPi.P,Gca,Gea,Charge);
		fprintf(fp,"4\t1\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t1.0\t1.0\n",Tai+1,TPi.T,TPi.P,Gca,Gea,Charge);
		fprintf(fp,"4\t1\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t1.0\t1.0\n",Tai+2,TPi.T,TPi.P,Gca,Gea,Charge);
		fprintf(fp,"4\t1\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t1.0\t1.0\n",Tai+3,TPi.T,TPi.P,Gca,Gea,Charge);
		fprintf(fp,"4\t1\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t1.0\t1.0\n",Tai+4,TPi.T,TPi.P,Gca,Gea,Charge);
		fprintf(fp,"4\t1\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t1.0\t1.0\n",Tai+5,TPi.T,TPi.P,Gca,Gea,Charge);
		fprintf(fp,"1\t%lf\t%lf\t%lf\n",X[0],X[1],X[2]);
		fprintf(fp,"4\t1\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t1.0\t1.0\n",Tai,TPi.T,TPi.P,Gca,Gea,Charge);
		fprintf(fp,"4\t1\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t1.0\t1.0\n",Tai-1,TPi.T,TPi.P,Gca,Gea,Charge);
		fprintf(fp,"4\t1\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t1.0\t1.0\n",Tai-2,TPi.T,TPi.P,Gca,Gea,Charge);
		fprintf(fp,"4\t1\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t1.0\t1.0\n",Tai-3,TPi.T,TPi.P,Gca,Gea,Charge);
		fprintf(fp,"4\t1\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t1.0\t1.0\n",Tai-4,TPi.T,TPi.P,Gca,Gea,Charge);
		fprintf(fp,"4\t1\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t1.0\t1.0\n",Tai-5,TPi.T,TPi.P,Gca,Gea,Charge);
		fclose(fp);
	}*/

	// TODD-Haorong
	fp = fopen("master_autogen.in","w");
	if(fp) {
		int i;
		int j;

		// vary outdoor ambient temperature only
		fprintf(fp,"1\t%lf\t%lf\t%lf\n",X[0],X[1],X[2]);
		const double Tamb_start = Tai;
		const double Tamb_step = DFtoDC(5);
		const double Twb_start = FtoC(60);
		const double Twb_step = DFtoDC(5);
		for(i=0;i<=6;i++) {
			for (j=0;j<=3;j++) {
				const double Twb = Twb_start + (double)j*Twb_step;
				const double Tdb = Twb + DFtoDC(10);
				TP tp = TZtoTP(Tdb,Twb);
				fprintf(fp,"4\t1\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t1.0\t1.0\n",Tamb_start+(double)i*Tamb_step,tp.T,tp.P,Gca,Gea,Charge);
			}
			for (j=0;j<=3;j++) {
				const double Twb = Twb_start - (double)j*Twb_step;
				const double Tdb = Twb + DFtoDC(10);
				TP tp = TZtoTP(Tdb,Twb);
				fprintf(fp,"4\t1\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t1.0\t1.0\n",Tamb_start+(double)i*Tamb_step,tp.T,tp.P,Gca,Gea,Charge);
			}

		}
		fprintf(fp,"1\t%lf\t%lf\t%lf\n",X[0],X[1],X[2]);
		for(i=1;i<=8;i++) {
			for (j=0;j<=3;j++) {
				const double Twb = Twb_start + (double)j*Twb_step;
				const double Tdb = Twb + DFtoDC(10);
				TP tp = TZtoTP(Tdb,Twb);
				fprintf(fp,"4\t1\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t1.0\t1.0\n",Tamb_start-(double)i*Tamb_step,tp.T,tp.P,Gca,Gea,Charge);
			}
			for (j=0;j<=3;j++) {
				const double Twb = Twb_start - (double)j*Twb_step;
				const double Tdb = Twb + DFtoDC(10);
				TP tp = TZtoTP(Tdb,Twb);
				fprintf(fp,"4\t1\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t1.0\t1.0\n",Tamb_start-(double)i*Tamb_step,tp.T,tp.P,Gca,Gea,Charge);
			}
		}
/*
		// vary indoor wetbulb only
		fprintf(fp,"1\t%lf\t%lf\t%lf\n",X[0],X[1],X[2]);
		const double Twb_start = FtoC(68);
		const double Twb_step = DFtoDC(2);
		for(i=0;i<=4;i++) {
			const double Twb = Twb_start + (double)i*Twb_step;
			const double Tdb = Twb + DFtoDC(13);
			TP tp = TZtoTP(Tdb,Twb);
			fprintf(fp,"4\t1\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t1.0\t1.0\n",Tai,tp.T,tp.P,Gca,Gea,Charge);
		}
		fprintf(fp,"1\t%lf\t%lf\t%lf\n",X[0],X[1],X[2]);
		for(i=1;i<=8;i++) {
			const double Twb = Twb_start - (double)i*Twb_step;
			const double Tdb = Twb + DFtoDC(13);
			TP tp = TZtoTP(Tdb,Twb);
			fprintf(fp,"4\t1\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t1.0\t1.0\n",Tai,tp.T,tp.P,Gca,Gea,Charge);
		}

		// vary charge only
		fprintf(fp,"1\t%lf\t%lf\t%lf\n",X[0],X[1],X[2]);
		const double Charge_start = 1.0; 
		const double Charge_step = 0.1;
		for(i=0;i<=5;i++) {
			const double Chg = Charge*(Charge_start + (double)i*Charge_step);
			fprintf(fp,"4\t1\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t1.0\t1.0\n",Tai,TPi.T,TPi.P,Gca,Gea,Chg);
		}
		fprintf(fp,"1\t%lf\t%lf\t%lf\n",X[0],X[1],X[2]);
		for(i=1;i<=9;i++) {
			const double Chg = Charge*(Charge_start - (double)i*Charge_step);
			fprintf(fp,"4\t1\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t1.0\t1.0\n",Tai,TPi.T,TPi.P,Gca,Gea,Chg);
		}
*/
		fclose(fp);
	}

	return 0;
}

/*********************************************************************
Main program for ACTUNE.  This program tunes model parameters so
that predictions match measured values.
*********************************************************************/
char* chillerMassNodeLabel[9] = {
	"not used",
	"compressor low stage",
	"compressor high stage",
	"oil separator",
	"condenser",
	"receiver",
	"subcooler",
	"accumulator",
	"evaporator"             // 8
};

int ActuneCondenser(double mdot_r,double DP,double DT,double Ga,double Tai,double LP,double LT,double Tao,double* mass)
{
	// inputs
	double Ones[5]={1,1,1,1,1};
	TXP TXPin = toTXP(DT,1,DP);
	HP HPin = TXPtoHP(TXPin);
	TXP TXPout = toTXP(LT,0,LP);
	HP HPout = TXPtoHP(TXPout);

	TXP TXP_prop={0,0,0};
	TXP_prop.P=LP;
	TXP_prop.X=0;
	const double Tsat_out = PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(LP);
	const double Tsc_out = Tsat_out-LT;

	// outputs
	double Tao_mod;
	MASS m;
	HP HPout_mod;

	CGP Cond_struc;
	while(Condenser("./InputDoc/acmodel.dat",mdot_r,HPin,Tai,Ga,&HPout_mod,&Tao_mod,&m,&Cond_struc,Ones)) {}
	if(errorLog.IsError()) {
		errorLog.Add("ActuneCondenser","Condenser");
		return 0;
	}

	*mass = m.m;

	TXP TXPout_mod = HPtoTXP(HPout_mod);

	TXP_prop.P=TXPout_mod.P;
	TXP_prop.X=0;

	const double Tsat_out_mod = PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(TXPout_mod.P);
	const double Tsc_out_mod = Tsat_out_mod-TXPout_mod.T;

	FILE* fp = fopen("actune.out","a");
	if(!fp) return 1;

	fprintf(fp,"Condenser tuning:\n");

	double Tsc_act = Tsc_out;
	double Tsc_mod = Tsc_out_mod;
	double Tsc_err = Tsc_mod-Tsc_act;
	double Tsc_per_err = 100.0*(Tsc_mod-Tsc_act)/Tsc_act;
	fprintf(fp,"Tsc_act=%lf Tsc_model=%lf Error=%lf (%+.2lf%%)\n",Tsc_act,Tsc_mod,Tsc_err,Tsc_per_err);
	if(Tsc_act<1.0) {
		fprintf(fp,">>> Actual Tsc is too low.  The condenser can not be tuned without adequate subcooling.  The actual refrigerant enthalpy change can not be determined.\n");
		fprintf(fp,"\n");
		fclose(fp);
		return 2;
	}

	double dTa_act = Tao-Tai;
	double dTa_mod = Tao_mod-Tai;
	double dTa_err = dTa_mod-dTa_act;
	double dTa_per_err = 100.0*(dTa_mod-dTa_act)/dTa_act;
	fprintf(fp,"dTa_act=%lf dTa_model=%lf Error=%lf (%+.2lf%%)\n",dTa_act,dTa_mod,dTa_err,dTa_per_err);
	double dHr_act = HPin.H-HPout.H;
	double dHr_mod = HPin.H-HPout_mod.H;
	double dHr_err = dHr_mod-dHr_act;
	double dHr_per_err = 100.0*(dHr_mod-dHr_act)/dHr_act;
	fprintf(fp,"dHr_act=%lf dHr_model=%lf Error=%lf (%+.2lf%%)\n",dHr_act,dHr_mod,dHr_err,dHr_per_err);
	double dPr_act = HPin.P-HPout.P;
	double dPr_mod = HPin.P-HPout_mod.P;
	double dPr_err = dPr_mod-dPr_act;
	double dPr_per_err = 100.0*(dPr_mod-dPr_act)/dPr_act;
	fprintf(fp,"dPr_act=%lf dPr_model=%lf Error=%lf (%+.2lf%%)\n",dPr_act,dPr_mod,dPr_err,dPr_per_err);
	fprintf(fp,"Refrigerant mass=%lf(kg).\n",*mass);

	if(dTa_per_err>+1.0) {
		fprintf(fp,">>> Increase condenser air flow rate.\n");
	} else if(dTa_per_err<-1.0) {
		fprintf(fp,">>> Reduce condenser air flow rate.\n");
	} else {
		fprintf(fp,">>> Condenser air flow rate is Ok.\n");
	}
	if(dHr_per_err>+1.0) {
		fprintf(fp,">>> Reduce condenser UA.\n");
	} else if(dHr_per_err<-1.0) {
		fprintf(fp,">>> Increase condenser UA.\n");
	} else {
		fprintf(fp,">>> Condenser UA Ok.\n");
	}

	/***
	double Prms[5]={1,1,1,1,1};
	Condenser0(Tsc_act,"./InputDoc/acmodel.dat",mdot_r,HPin,Tai,Ga,&HPout_mod,&Tao_mod,&m,Prms);
	if(errorLog.IsError()) {errorLog.Add("ActuneCondenser"); return 0;}
	fprintf(fp,">>> Multiply Air and Refrigerant side UA by %lf,%lf.\n",Prms[0],Prms[1]);
	***/

	fprintf(fp,"\n");
	fclose(fp);

	return 0;
}

// variable      description                   units
// mdot_r        refrigerant mass flow rate    kg/s
// SP            suction pressure              Pa
// ST            suction temperature           C
// mdot_f        air/water mass flow rate      kg/s
int ActuneEvaporator(double mdot_r,double SP,double ST,double mdot_f,double Tfi,double LP,double LT,double Tfo,double* mass)
{
	// inputs
	double Ones[5]={1,1,1,1,1};
	TXP TXPin = toTXP(LT,0,LP);
	HP HPin = TXPtoHP(TXPin);
	TXP TXPout = toTXP(ST,1,SP);
	HP HPout = TXPtoHP(TXPout);

	// outputs
	double Tfo_mod;
	MASS m;
	HP HPin_mod;

	while(Evaporator("./InputDoc/acmodel.dat",mdot_r,HPout,mdot_f,Tfi,HPin_mod,Tfo_mod,m,Ones)) {}
	if(errorLog.IsError()) {
		errorLog.Add("ActuneEvaporator");
		return 0;
	}

	*mass = m.m;

	FILE* fp = fopen("actune.out","a");
	if(fp) {
		fprintf(fp,"Evaporator tuning:\n");
		double dTa_act = Tfi-Tfo;
		double dTa_mod = Tfi-Tfo_mod;
		double dTa_err = dTa_mod-dTa_act;
		double dTa_per_err = 100.0*(dTa_mod-dTa_act)/dTa_act;
		fprintf(fp,"dTa_act=%lf dTa_model=%lf Error=%lf (%.2lf%%)\n",dTa_act,dTa_mod,dTa_err,dTa_per_err);
		double dHr_act = HPout.H-HPin.H;
		double dHr_mod = HPout.H-HPin_mod.H;
		double dHr_err = dHr_mod-dHr_act;
		double dHr_per_err = 100.0*(dHr_mod-dHr_act)/dHr_act;
		fprintf(fp,"dHr_act=%lf dHr_model=%lf Error=%lf (%.2lf%%)\n",dHr_act,dHr_mod,dHr_err,dHr_per_err);
		fprintf(fp,"Refrigerant mass=%lf(kg).\n",*mass);

		if(dTa_per_err>+1.0) {
			fprintf(fp,">>> Increase evaporator air flow rate.\n");
		} else if(dTa_per_err<-1.0) {
			fprintf(fp,">>> Reduce evaporator air flow rate.\n");
		} else {
			fprintf(fp,">>> Evaporator air flow rate is Ok.\n");
		}
		if(dHr_per_err>+1.0) {
			fprintf(fp,">>> Reduce evaporator UA.\n");
		} else if(dHr_per_err<-1.0) {
			fprintf(fp,">>> Increase evaporator UA.\n");
		} else {
			fprintf(fp,">>> Evaporator UA Ok.\n");
		}

		fprintf(fp,"\n");
		fclose(fp);
	}

	return 0;
}

int ActuneSubcooler(double mr1,double P4,double T4,double mr3,double P5,double T5,double P10_11_12,double T10,double T11,double* mass)
{
	// inputs
	TXP TXP4 = toTXP(T4,0,P4);
	HP HP4 = TXPtoHP(TXP4);
	TXP TXP5 = toTXP(T5,0,P5);
	HP HP5 = TXPtoHP(TXP5);
	TXP TXP10 = toTXP(T10,1,P10_11_12);
	HP HP10 = TXPtoHP(TXP10);
	TXP TXP11 = toTXP(T11,1,P10_11_12);
	HP HP11 = TXPtoHP(TXP11);

	const double mr2 = mr1 + mr3;
	HP HP9;
	HP9.H = (mr2*HP11.H-mr1*HP10.H)/mr3;
	HP9.P = P10_11_12;
	//TXP TXP9 = HPtoTXP(HP9);

	HP HP12;
	HP12.H = HP4.H;
	HP12.P = P10_11_12;
	//TXP TXP12 = HPtoTXP(HP12);

	// outputs
	MASS m;
	HP HP5_mod;
	HP HP12_mod;

	Subcooler("./InputDoc/acmodel.dat",mr1,HP4,mr3,HP9,HP5_mod,HP12_mod,m);
	if(errorLog.IsError()) {
		errorLog.Add("ActuneSubcooler");
		return 0;
	}

	*mass = m.m;

	FILE* fp = fopen("actune.out","a");
	if(fp) {
		fprintf(fp,"Subcooler tuning:\n");
		double dHhp_act = HP4.H-HP5.H;
		double dHhp_mod = HP4.H-HP5_mod.H;
		double dHhp_err = dHhp_mod-dHhp_act;
		double dHhp_per_err = 100.0*(dHhp_mod-dHhp_act)/dHhp_act;
		fprintf(fp,"dHhp_act=%lf dHhp_model=%lf Error=%lf (%+.2lf%%)\n",dHhp_act,dHhp_mod,dHhp_err,dHhp_per_err);
		double dHlp_act = HP9.H-HP12.H;
		double dHlp_mod = HP9.H-HP12_mod.H;
		double dHlp_err = dHlp_mod-dHlp_act;
		double dHlp_per_err = 100.0*(dHlp_mod-dHlp_act)/dHlp_act;
		fprintf(fp,"dHlp_act=%lf dHlp_model=%lf Error=%lf (%+.2lf%%)\n",dHlp_act,dHlp_mod,dHlp_err,dHlp_per_err);
		fprintf(fp,"Refrigerant mass=%lf(kg).\n",*mass);
		fprintf(fp,"\n");
		fclose(fp);
	}

	return 0;
}

// Tsh_mod (C)
int
ActuneExpansionDevice(double P,double T,double Tsh_mod)
{
	TXP TXP_prop={0,0,0};
	TXP_prop.P=P;
	TXP_prop.X=1;
	const double Tsat = PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(P);
	const double Tsh = T - Tsat;

	FILE* fp = fopen("actune.out","a");
	if(fp) {
		fprintf(fp,"Expansion Device tuning:\n");
		double Tsh_act = Tsh;
		double Tsh_err = Tsh_mod-Tsh_act;
		fprintf(fp,"Tsh_act=%lf Tsh_model=%lf Error=%lf\n",Tsh_act,Tsh_mod,Tsh_err);
		fprintf(fp,"\n");
		fclose(fp);
	}

	return 0;
}

int ActuneCCLowStage(double P8,double T8,double P10,double P1,double T10,double mr,double E,double* mass)
{
	double Ones[5]={1,1,1,1,1};
	TXP TXP8 = toTXP(T8,1,P8);
	HP HP8 = TXPtoHP(TXP8);
	TXP TXP10 = toTXP(T10,1,P10);
	HP HP10 = TXPtoHP(TXP10);

	double H10_mod;
	double mr_mod;
	double E_mod;
	MASS m;

	TXP TXP_prop={0,0,0};

	CompressorCCLowStage("./InputDoc/acmodel.dat",HP8.H,HP8.P,P10,P1,H10_mod,mr_mod,E_mod,m,Ones);

	*mass = m.m;

	FILE* fp = fopen("actune.out","a");
	if(fp) {
		//fprintf(fp,"TXP8=%lf,%lf,%lf\n",TXP8.T,TXP8.X,TXP8.P);
		//fprintf(fp,"HP8=%lf,%lf\n",HP8.H,HP8.P);

		fprintf(fp,"Low stage compressor tuning:\n");

		double dH_act = HP10.H-HP8.H;
		double dH_mod = H10_mod-HP8.H;
		double dH_err = dH_mod-dH_act;
		double dH_per_err = 100.0*(dH_mod-dH_act)/dH_act;
		fprintf(fp,"dH_act=%lf dH_model=%lf Error=%lf (%+.2lf%%)\n",dH_act,dH_mod,dH_err,dH_per_err);

		double mdot_act = mr;
		double mdot_mod = mr_mod;
		double mdot_err = mdot_mod-mdot_act;
		double mdot_per_err = 100.0*(mdot_mod-mdot_act)/mdot_act;
		fprintf(fp,"mdot_act=%lf mdot_model=%lf Error=%lf (%+.2lf%%)\n",mdot_act,mdot_mod,mdot_err,mdot_per_err);
		if(fabs(mdot_per_err)>+1.0) {
			fprintf(fp,">>> Compare model to compressor map to resolve errors.\n");
			const double sp_psig = P8/101.3*14.7-14.7;
			
			TXP_prop.P=P8;
			TXP_prop.X=0;
			const double et_f = PropertyTXPth(TSAT,TXP_prop)*1.8+32;;//reftplthP.Tsat(P8)*1.8+32;
			fprintf(fp,">>> SP=%lf (psig) ET=%lf (F)\n",sp_psig,et_f);
			const double dp_psig = P1/101.3*14.7-14.7;

			TXP_prop.P=P1;
			TXP_prop.X=0;
			const double ct_f = PropertyTXPth(TSAT,TXP_prop)*1.8+32;;//reftplthP.Tsat(P1)*1.8+32;
			fprintf(fp,">>> DP=%lf (psig) CT=%lf (F)\n",dp_psig,ct_f);
			fprintf(fp,">>> mdot_act=%lf mdot_model=%lf (lbs/hr)\n",mdot_act*3600.0*2.205,mdot_mod*3600.0*2.205);
		} else {
			fprintf(fp,">>> Mass flow rate prediction is Ok.\n");
		}

		double power_act = E;
		double power_mod = E_mod;
		double power_err = power_mod-power_act;
		double power_per_err = 100.0*(power_mod-power_act)/power_act;
		fprintf(fp,"power_act=%lf power_model=%lf Error=%lf (%+.2lf%%)\n",power_act,power_mod,power_err,power_per_err);

		fprintf(fp,"Refrigerant mass=%lf(kg).\n",*mass);
		fprintf(fp,"\n");
		fclose(fp);
	}

	return 0;
}

int ActuneCCHighStage(double P11,double T11,double P1,double T1)
{
	TXP TXP11 = toTXP(T11,1,P11);
	HP HP11 = TXPtoHP(TXP11);
	TXP TXP1 = toTXP(T1,1,P1);
	HP HP1 = TXPtoHP(TXP1);

	double H1_mod;

	CompressorCCHighStage("./InputDoc/acmodel.dat",HP11.H,HP11.P,HP1.P,H1_mod);

	FILE* fp = fopen("actune.out","a");
	if(fp) {
		fprintf(fp,"High stage compressor tuning:\n");
		double dH_act = HP1.H-HP11.H;
		double dH_mod = H1_mod-HP11.H;
		double dH_err = dH_mod-dH_act;
		double dH_per_err = 100.0*(dH_mod-dH_act)/dH_act;
		fprintf(fp,"dH_act=%lf dH_model=%lf Error=%lf (%+.2lf%%)\n",dH_act,dH_mod,dH_err,dH_per_err);
		fprintf(fp,"\n");
		fclose(fp);
	}

	return 0;
}

int ActuneNISTChillerComponents(SteadyData* steady)
{
	FILE* fp = fopen("actune.out","w");
	if(fp) {
		fprintf(fp,"Opening...\n");
		fclose(fp);
	}
 
	/***
	// Natascha's original data with no subcooling
	const double mr2 = 0.379;     //        kg/s
	const double T2 = 90.70;      //        C, compressor outlet
	const double P2 = 1540;       //        kPa
	const double T3 = 33.54;      //        C, condenser outlet
	const double P3 = 1302;       //        kPa
	const double Gca = 40;        //        kg/s/m^2
	const double Tcai = 2;        //        C
	const double Tcao = 38;       //        C
	const double mr1 = 0.28;      //        kg/s
	const double P7 = 230;        //        kPa, evaporator outlet
	const double T7 = -19.12;     //        C
	const double P5 = 1298.4;     //        kPa, subcooler outlet
	const double T5 = 15.43;      //        C
	const double mw = 3.25;       //        kg/s
	const double Tewi = -3.79;    //        C
	const double Tewo = -10.0;    //        C
	const double P10_11_12=686.5; //        kPa
	const double T10 = 42.65;     //        C
	const double T11 = 31.24;     //        C
	const double E = 63.23;       //        kW, compressor power
	***/

	/***
	// 3/9/99 data
	// sent on 3/15/99
	// Tamb = -3.28C
	const double mr2 = 0.33355;   //        kg/s
	const double T2 = 93.58;      //        C, compressor outlet
	const double P2 = 1520.83;    //        kPa
	const double T3 = 35.58;      //        C, condenser outlet
	const double P3 = 1440.423;   //        kPa
	const double Gca = 40;        //        kg/s/m^2
	const double Tcai = -3.28;    //        C
	const double Tcao = 36.87;    //        C
	const double mr1 = 0.231745;  //        kg/s
	const double P7 = 151.7178;   //        kPa, evaporator outlet
	const double T7 = -28.85;     //        C
	const double P5 = 1360.017;   //        kPa, subcooler outlet
	const double T5 = 6.41;       //        C
	const double mw = 3.200103;   //        kg/s
	const double Tewi = -15.73;   //        C
	const double Tewo = -19.86;   //        C
	const double P10_11_12=515.5787; //        kPa
	const double T10 = 32.66;     //        C
	const double T11 = 19.59;     //        C
	const double E = 63.23;       //  ???   kW, compressor power
	***/

	/***
	// 3/17/99 1:58pm
	// steady state, 65F ambient, 2 cond fans running
	// Tamb = 23.28C
	const double mr2 = 0.343721;  //        kg/s
	const double T2 = 91.21;      //        C, compressor outlet
	const double P2 = 1469.499;   //        kPa
	const double T3 = 35.17;      //        C, condenser outlet
	const double P3 = 1399.359;   //        kPa
	const double Gca = 40;        //        kg/s/m^2
	const double Tcai = 23.28;    //        C
	const double Tcao = 40.71;    //        C
	const double mr1 = 0.233258;  //        kg/s
	const double P7 = 146.068;    //        kPa, evaporator outlet
	const double T7 = -29.58;     //        C
	const double P5 = 1329.219;   //        kPa, subcooler outlet
	const double T5 = 5.38;       //        C
	const double mw = 3.15713;    //        kg/s
	const double Tewi = -15.73;   //        C
	const double Tewo = -20.53;   //        C
	const double P10_11_12=506.1394; //        kPa
	const double T10 = 32.2;      //        C
	const double T11 = 11.94;     //        C
	const double E = 63.23;       //  ???   kW, compressor power
	***/

	/***
	// data set send 3/23/99
	// Tamb = 2.74C
	const double mr2 = 0.321595;  //        kg/s
	const double T2 = 93.12;      //        C, compressor outlet
	const double P2 = 1513.526;   //        kPa
	const double T3 = 35;         //        C, condenser outlet
	const double P3 = 1447.038;   //        kPa
	//const double Gca = 40;      //        kg/s/m^2
	const double Gca = -1;        //        kg/s/m^2
	const double Tcai = 2.74;     //        C
	const double Tcao = 38.07;    //        C
	const double mr1 = 0.220017;  //        kg/s
	const double P7 = 137.8;      //        kPa, evaporator outlet
	const double T7 = -31.73;     //        C
	const double P5 = 1380.549;   //        kPa, subcooler outlet
	const double T5 = 5.58;       //        C
	const double mw = 3.121733;   //        kg/s
	const double Tewi = -17.58;   //        C
	const double Tewo = -22.05;   //        C
	const double P10_11_12=489.9479; //        kPa
	const double T10 = 31.02;     //        C
	const double T11 = 15.88;     //        C
	const double E = 63.23;       //  ???   kW, compressor power
	***/

	/***
	// data set send 5/27/99
	// normal data after flow rate correction
	const double mr2 = 0.221639;  //        kg/s
	const double T2 = 93.69;      //        C, compressor outlet
	const double P2 = 1546.645;   //        kPa
	const double T3 = 30.4;       //        C, condenser outlet
	const double P3 = 1490.284;   //        kPa
	//const double Gca = 40;      //        kg/s/m^2
	const double Gca = -1;        //        kg/s/m^2
	const double Tcai = 7.8895;   //        C
	const double Tcao = 38.12;    //        C
	const double mr1 = 0.182148;  //        kg/s
	const double P7 = 141.56;     //        kPa, evaporator outlet
	const double T7 = -30.07;     //        C
	const double P5 = 1433.923;   //        kPa, subcooler outlet
	const double T5 = 9.81;       //        C
	const double mw = 3.1355;     //        kg/s
	const double Tewi = -16.89;   //        C
	const double Tewo = -21.191;  //        C
	const double P10_11_12=471.29;//        kPa
	const double T10 = 26.41;     //        C
	const double T11 = 13.52;     //        C
	const double E = 63.23;       //  ???   kW, compressor power
	***/

	/***
	// data set send 7/22/99 as part of steady-state data set
	const double mr2 = 0.228958;     //        kg/s
	const double T2 = 89.45884;      //        C, compressor outlet
	const double P2 = 1530.491;      //        kPa
	const double T3 = 32.29615;      //        C, condenser outlet
	const double P3 = 1473.639;      //        kPa
	//const double Gca = 40;         //        kg/s/m^2
	const double Gca = -1;           //        kg/s/m^2
	const double Tcai = 11.46957;    //        C
	const double Tcao = 37.13032;    //        C
	const double mr1 = 0.190459;     //        kg/s
	const double P7 = 141.8417;      //        kPa, evaporator outlet
	const double T7 = -31.2247;      //        C
	const double P5 = 1416.788;      //        kPa, subcooler outlet
	const double T5 = 12.03681;      //        C
	const double mw = 3.131245;      //        kg/s
	const double Tewi = -17.3451;    //        C
	const double Tewo = -21.4402;    //        C
	const double P10_11_12=466.3234; //        kPa
	const double T10 = 21.66723;     //        C
	const double T11 = 11.12536;     //        C
	const double E = 20.84763;       //        kW, compressor power
	***/

	double mass_cond;
	if(ActuneCondenser(steady->mr2,steady->P2,steady->T2,steady->Gca,steady->Tcai,steady->P3,steady->T3,steady->Tcao,&mass_cond)) return 1;
	if(errorLog.IsError()) {errorLog.Add("ActuneNISTChillerComponents"); return 0;}

	double mass_evap;
	ActuneEvaporator(steady->mr1,steady->P7,steady->T7,steady->mw,steady->Tewi,steady->P5,steady->T5,steady->Tewo,&mass_evap);
	if(errorLog.IsError()) {errorLog.Add("ActuneNISTChillerComponents"); return 0;}

	double mass_subcooler;
	ActuneSubcooler(steady->mr1,steady->P3,steady->T3,steady->mr2-steady->mr1,steady->P5,steady->T5,steady->P10_11_12,steady->T10,steady->T11,&mass_subcooler);
	if(errorLog.IsError()) {errorLog.Add("ActuneNISTChillerComponents"); return 0;}

	// entering low stage compressor from evaporator
	ActuneExpansionDevice(steady->P7,steady->T7,TSH_LP);
	if(errorLog.IsError()) {errorLog.Add("ActuneNISTChillerComponents"); return 0;}

	// entering high stage compressor at intermediate pressure
	ActuneExpansionDevice(steady->P10_11_12,steady->T11,TSH_HP);
	if(errorLog.IsError()) {errorLog.Add("ActuneNISTChillerComponents"); return 0;}

	double mass_cclowstage;
	ActuneCCLowStage(steady->P7,steady->T7,steady->P10_11_12,steady->P2,steady->T10,steady->mr1,steady->E,&mass_cclowstage);
	if(errorLog.IsError()) {errorLog.Add("ActuneNISTChillerComponents"); return 0;}

	ActuneCCHighStage(steady->P10_11_12,steady->T11,steady->P2,steady->T2);
	if(errorLog.IsError()) {errorLog.Add("ActuneNISTChillerComponents"); return 0;}

	fp = fopen("actune.out","a");
	if(fp) {
		fprintf(fp,"Total mass=%lf(kg).\n",mass_cond+mass_evap+mass_subcooler+mass_cclowstage);
		fclose(fp);
	}

	return 0;
}

int ActuneNISTChiller()
{
	if(errorLog.IsError()) errorLog.ClearError("ActuneNISTChiller");

	// Get parameter file name for "prms"
	const char *filename = "./InputDoc/acmodel.dat";

	/***
	{
	// Get input file "actune.in"
	FILE* fp = fopen("actune3.in","r");
	GetInt(fp);
	// Driving conditions
	double Tcai = GetDouble(fp);    // ambient temperature
	double Gca  = GetDouble(fp);    // cond air mass flux
	double Tewi = GetDouble(fp);    // return water temperature
	double mw   = GetDouble(fp);    // evaporator water flow rate
	// Measured cycle states
	double P1   = GetDouble(fp);      // head pressure
	double Pint = GetDouble(fp);      // intermediate pressure
	double P8   = GetDouble(fp);      // suction pressure
	fclose(fp);
	}
	***/

	const double Tcai = 11.46957;
	const double Gca = -1;
	const double Tewi = -17.3451;
	const double mw = 2.0;
	const double P1 = 1530.491;
	const double Pint = 466.3234;
	const double P8 = 141.8417;
	TXP TXP_prop={0,0,0};

	/***
	const double Tcai = gData.input[3];
	const double Gca = gData.input[13];
	const double Tewi = gData.input[4];
	const double mw = gData.input[14];
	const double P1 = gData.input[0];
	const double Pint = gData.input[1];
	const double P8 = gData.input[2];
	***/

	FILE* fp = fopen("actune.out","w");
	if(fp) {
		fprintf(fp,"Input deck:\n");
		fprintf(fp,"Tcai  = %14lf  Gca = %14lf\n",Tcai,Gca);
		fprintf(fp,"Tewi  = %14lf  mw  = %14lf\n",Tewi,mw);
		fprintf(fp,"Phead = %14lf\n",P1);
		fprintf(fp,"Pint  = %14lf\n",Pint);
		fprintf(fp,"Psuc  = %14lf\n\n",P8);
		fclose(fp);
	}

	double mr1,mr2,mr3;
	double Ones[5]={1,1,1,1,1};
	MASS m[20];
	HP HPs[20];
	TXP TXPs[20];

	double Y[3];

	//
	// insert guess values into thermodynamic states
	HPs[8].P = P8;
	HPs[9].P = HPs[10].P = HPs[11].P = Pint;
	HPs[1].P = P1;

	//
	// low stage inlet enthalpy (TxV 1)
	double Tsh_8 = 0;
	{
		TXP_prop.P=HPs[8].P;
		TXP_prop.X=1;
		double T8sat = PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(HPs[8].P);
		if(errorLog.IsError()) {
			errorLog.Add("ActuneNISTChiller","T8sat");
			return 0;
		}
		double T8 = T8sat + TSH_LP;
		Tsh_8 = T8 - T8sat;

		TXPs[8] = toTXP(T8,1,HPs[8].P);
		HPs[8] = TXPtoHP(TXPs[8]);
		if(errorLog.IsError()) {
			errorLog.Add("ActuneNISTChiller","HPs[8]");
			return 0;
		}
		#ifdef EXTRA_REPORTING
		{
		FILE* fp = fopen("actune.dbg","w");
		if(fp) {
			fprintf(fp,"Low stage inlet enthalpy.\n");
			fprintf(fp,"T8sat  = %8.3lf\n",T8sat);
			fprintf(fp,"T8     = %8.3lf\n",T8);
			fprintf(fp,"HTXP8  = %8.3lf,%8.3lf,%8.3lf,%8.3lf (low stage compressor inlet)\n\n",0.001*HPs[8].H,TXPs[8].T,TXPs[8].X,TXPs[8].P);
			fclose(fp);
		}
		}
		#endif
	}
	//
	// low stage compressor
	{
		double Ei;
		CompressorCCLowStage(filename,HPs[8].H,HPs[8].P,HPs[10].P,HPs[1].P,HPs[10].H,mr1,Ei,m[1],Ones);
		if(errorLog.IsError()) {
			errorLog.Add("ActuneNISTChiller","CompressorCCLowStage");
			return 0;
		}
		#ifdef EXTRA_REPORTING
		{
		FILE* fp = fopen("actune.dbg","a");
		if(fp) {
			fprintf(fp,"Low stage compressor.\n");
			fprintf(fp,"mr     = %8.3lf\n",mr1);
			fprintf(fp,"Ei     = %8.3lf\n",Ei);
			fprintf(fp,"HP8    = %8.3lf,%8.3lf\n",HPs[8].H/1000,HPs[8].P);
			TXPs[10] = HPtoTXP(HPs[10]);
			if(errorLog.IsError()) {
				errorLog.Add("ActuneNISTChiller","TXPs[10]");
				fclose(fp);
				return 0;
			}
			fprintf(fp,"HTXP10 = %8.3lf,%8.3lf,%8.3lf,%8.3lf\n\n",0.001*HPs[10].H,TXPs[10].T,TXPs[10].X,TXPs[10].P);
			fclose(fp);
		}
		}
		#endif
	}
	//
	// high stage inlet enthalpy (TxV 2)
	double Tsh_11 = 0;
	{
		TXP_prop.P=HPs[11].P;
		TXP_prop.X=0;

		double T11sat = PropertyTXPth(TSAT,TXP_prop);//reftplthP.Tsat(HPs[11].P);
		if(errorLog.IsError()) {
			errorLog.Add("ActuneNISTChiller","T11sat");
			return 0;
		}
		double T11 = T11sat + TSH_HP;
		Tsh_11 = T11 - T11sat;

		TXPs[11] = toTXP(T11,1,HPs[11].P);
		HPs[11] = TXPtoHP(TXPs[11]);
		if(errorLog.IsError()) {
			errorLog.Add("ActuneNISTChiller","HPs[11]");
			return 0;
		}
	}
	//
	// high stage compressor
	{
		CompressorCCHighStage(filename,HPs[11].H,HPs[11].P,HPs[1].P,HPs[1].H);
		m[2].m = m[2].V = 0.0;
		if(errorLog.IsError()) {
			errorLog.Add("ActuneNISTChiller","CompressorCCHighStage");
			return 0;
		}
		#ifdef EXTRA_REPORTING
		{
			FILE* fp = fopen("actune.dbg","a");
			if(fp) {
				fprintf(fp,"High stage compressor.\n");
				fprintf(fp,"HP11   = %8.3lf,%8.3lf\n",HPs[11].H/1000,HPs[11].P);
				TXPs[1] = HPtoTXP(HPs[1]);
				if(errorLog.IsError()) {
					errorLog.Add("ActuneNISTChiller","TXPs[1]");
					fclose(fp);
					return 0;
				}
				fprintf(fp,"HTXP1  = %8.3lf,%8.3lf,%8.3lf,%8.3lf\n\n",0.001*HPs[1].H,TXPs[1].T,TXPs[1].X,TXPs[1].P);
				fclose(fp);
			}
		}
		#endif
	}
	//
	// volume flow rate of low stage is twice the high stage
	{
		double v8  = PropertyTXPth(VOL,TXPs[8]);
		if(errorLog.IsError()) {
			errorLog.Add("ActuneNISTChiller","v8");
			return 0;
		}
		double v11 = PropertyTXPth(VOL,TXPs[11]);
		if(errorLog.IsError()) {
			errorLog.Add("ActuneNISTChiller","v11");
			return 0;
		}
		mr2 = 0.5*mr1*v8/v11;
	}
	//
	// Refrigerant entering intermediate stage from subcooler.
	// Use mass and energy balance.
	{
		mr3 = mr2 - mr1;
		HPs[9].H = (mr2*HPs[11].H-mr1*HPs[10].H)/mr3;
		TXPs[9] = HPtoTXP(HPs[9]);
		if(errorLog.IsError()) {
			errorLog.Add("ActuneNISTChiller","TXPs[9]");
			fclose(fp);
			return 0;
		}
		#ifdef EXTRA_REPORTING
		{
			FILE* fp = fopen("actune.dbg","a");
			if(fp) {
				fprintf(fp,"Intermediate stage.\n");
				fprintf(fp,"mr2 = %lf\n",mr2);
				fprintf(fp,"mr3 = %lf\n",mr3);
				fprintf(fp,"HTXP9  = %8.3lf,%8.3lf,%8.3lf,%8.3lf\n\n",0.001*HPs[9].H,TXPs[9].T,TXPs[9].X,TXPs[9].P);
				fclose(fp);
			}
		}
		#endif
	}
	//
	// oil separator
	HPs[2] = HPs[1];
	m[3].m = m[3].V = 0.0;
	//
	// condenser
	{
		double Tao;
		CGP Cond_struc;
		while(Condenser(filename,mr2,HPs[2],Tcai,Gca,&HPs[3],&Tao,&m[4],&Cond_struc,Ones)) {}
		if(errorLog.IsError()) {
			errorLog.Add("ActuneNISTChiller","Condenser");
			return 0;
		}
		#ifdef EXTRA_REPORTING
		{
			FILE* fp = fopen("actune.dbg","a");
			if(fp) {
				fprintf(fp,"Condenser.\n");
				TXPs[2] = HPtoTXP(HPs[2]);
				if(errorLog.IsError()) {
					errorLog.Add("ActuneNISTChiller","TXPs[2]");
					fclose(fp);
					return 0;
				}
				fprintf(fp,"HTXP2  = %8.3lf,%8.3lf,%8.3lf,%8.3lf\n",0.001*HPs[2].H,TXPs[2].T,TXPs[2].X,TXPs[2].P);
				TXPs[3] = HPtoTXP(HPs[3]);
				if(errorLog.IsError()) {
					errorLog.Add("ActuneNISTChiller","TXPs[3]");
					fclose(fp);
					return 0;
				}
				fprintf(fp,"HTXP3  = %8.3lf,%8.3lf,%8.3lf,%8.3lf\n\n",0.001*HPs[3].H,TXPs[3].T,TXPs[3].X,TXPs[3].P);
				fclose(fp);
			}
		}
		#endif
	}
	//
	// receiver
	HPs[4] = HPs[3];
	m[5].m = m[5].V = 0.0;
	TXPs[4] = HPtoTXP(HPs[4]);
	if(errorLog.IsError()) {
		errorLog.Add("ActuneNISTChiller","TXPs[4]");
		return 0;
	}
	if(TXPs[9].T>TXPs[4].T) {
		fp = fopen("actune.out","a");
		if(fp) {
			fprintf(fp,"Results:\n");
			fprintf(fp,"T9(%lf)>T4(%lf) is unphysical.  Decrease condenser UA.\n",TXPs[9].T,TXPs[4].T);
			fclose(fp);
		}
		return 0;
	}
	//
	// subcooler
	{
		Subcooler(filename,mr1,HPs[4],mr3,HPs[9],HPs[5],HPs[12],m[6]);
		if(errorLog.IsError()) {
			errorLog.Add("ActuneNISTChiller","Subcooler");
			return 0;
		}
		TXPs[5] = HPtoTXP(HPs[5]);
		if(errorLog.IsError()) {
			errorLog.Add("ActuneNISTChiller","TXPs[5]");
			return 0;
		}
		#ifdef EXTRA_REPORTING
		{
			FILE* fp = fopen("actune.dbg","a");
			if(fp) {
				fprintf(fp,"Subcooler.\n");
				fprintf(fp,"HTXP5  = %8.3lf,%8.3lf,%8.3lf,%8.3lf\n",0.001*HPs[5].H,TXPs[5].T,TXPs[5].X,TXPs[5].P);
				TXPs[12] = HPtoTXP(HPs[12]);
				if(errorLog.IsError()) {
					errorLog.Add("ActuneNISTChiller","TXPs[12]");
					fclose(fp);
					return 0;
				}
				fprintf(fp,"HTXP12 = %8.3lf,%8.3lf,%8.3lf,%8.3lf\n\n",0.001*HPs[12].H,TXPs[12].T,TXPs[12].X,TXPs[12].P);
				fclose(fp);
			}
		}
		#endif
	}
	//
	// accumulator
	HPs[7] = HPs[8];
	m[7].m = m[7].V = 0.0;
	#ifdef EXTRA_REPORTING
	{
		FILE* fp = fopen("actune.dbg","a");
		if(fp) {
			fprintf(fp,"Accumulator.\n");
			TXPs[7] = HPtoTXP(HPs[7]);
			if(errorLog.IsError()) {
				errorLog.Add("ActuneNISTChiller","TXPs[7]");
				fclose(fp);
				return 0;
			}
			fprintf(fp,"HTXP7 = %8.3lf,%8.3lf,%8.3lf,%8.3lf\n\n",0.001*HPs[7].H,TXPs[7].T,TXPs[7].X,TXPs[7].P);
			fclose(fp);
		}
	}
	#endif
	//
	// evaporator
	{
		double Tewo;
		while(Evaporator(filename,mr1,HPs[7],mw,Tewi,HPs[6],Tewo,m[8],Ones)) {}
		if(errorLog.IsError()) {
			errorLog.Add("ActuneNISTChiller","Evaporator: (1)UA may be too large to support outlet superheat constraint.");
			return 0;
		}
		#ifdef EXTRA_REPORTING
		{
			FILE* fp = fopen("actune.dbg","a");
			if(fp) {
				fprintf(fp,"Evaporator.\n");
				TXPs[6] = HPtoTXP(HPs[6]);
				if(errorLog.IsError()) {
					errorLog.Add("ActuneNISTChiller","TXPs[6] (evaporator inlet state)");
					fclose(fp);
					return 0;
				}
				fprintf(fp,"HTXP6  = %8.3lf,%8.3lf,%8.3lf,%8.3lf (evaporator inlet state)\n",0.001*HPs[6].H,TXPs[6].T,TXPs[6].X,TXPs[6].P);
				fprintf(fp,"Tewi/o = %8.3lf,%8.3lf\n\n",Tewi,Tewo);
				fclose(fp);
			}
		}
		#endif
	}
	//
	// residuals
	// dh=0 across TxV1
	Y[0]=(HPs[5].H-HPs[6].H)/(HPs[5].H+HPs[6].H);
	// dh=0 across TxV2
	Y[1]=(HPs[4].H-HPs[12].H)/(HPs[4].H+HPs[12].H);
	// charge is conserved
	double Charge = 0;
	double Volume = 0;
	for(int i=1;i<=8;i++) {
		Charge += m[i].m;
		Volume += m[i].V;
	}
	#ifdef EXTRA_REPORTING
	{
		FILE* fp = fopen("actune.dbg","a");
		if(fp) {
			fprintf(fp,"Charge.\n");
			for(int i=1;i<=8;i++) {
				fprintf(fp,"%-25.25s charge=%lf volume=%lf.\n",chillerMassNodeLabel[i],m[i].m,m[i].V);
			}
			fprintf(fp,"%-25.25s charge=%lf volume=%lf.\n\n","total",Charge,Volume);
			fclose(fp);
		}
	}
	#endif

	// print report
	fp = fopen("actune.out","a");
	if(fp) {
		fprintf(fp,"Results:\n");
		fprintf(fp,"Tsh_8 = %lf (C,entering compressor low stage)\n",Tsh_8);
		fprintf(fp,"Tsh_ll = %lf (C,entering compressor high stage)\n",Tsh_11);
		fprintf(fp,"\n");

		const double dH1 = 1e-3*(HPs[5].H-HPs[6].H);
		const double H1_ave = 0.5*1e-3*(HPs[5].H+HPs[6].H);
		const double dH1_pe = 1e2*dH1/H1_ave;
		fprintf(fp,"H5=%lf H6=%lf dH1=%lf (kJ/kg) %lf (%%error)\n",1e-3*HPs[5].H,1e-3*HPs[6].H,dH1,dH1_pe);
		if(dH1_pe<-1.0) {
			fprintf(fp,">>> Increase evaporator UA.\n");
		} else if(dH1_pe>+1.0) {
			fprintf(fp,">>> Decrease evaporator UA.\n");
		} else {
			fprintf(fp,">>> Evaporator UA Ok.\n");
		}
		fprintf(fp,"\n");

		const double dH2 = 1e-3*(HPs[4].H-HPs[12].H);
		const double H2_ave = 0.5*1e-3*(HPs[4].H+HPs[12].H);
		const double dH2_pe = 1e2*dH2/H2_ave;
		fprintf(fp,"H4=%lf H12=%lf dH2=%lf (kJ/kg) %lf (%%error)\n",1e-3*HPs[4].H,1e-3*HPs[12].H,dH2,dH2_pe);
		if(dH2_pe<-1.0) {
			fprintf(fp,">>> Decrease condenser UA.\n");
		} else if(dH2_pe>+1.0) {
			fprintf(fp,">>> Increase condenser UA.\n");
		} else {
			fprintf(fp,">>> Condenser UA Ok.\n");
		}
		fprintf(fp,"\n");

		fprintf(fp,">>> Set charge to %lf (kg)\n\n",Charge);

		fclose(fp);
	}

	return 0;
}

