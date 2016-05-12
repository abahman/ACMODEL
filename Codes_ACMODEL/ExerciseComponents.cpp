#include "stdafx.h"

#include "ExerciseComponents.h"
#include "expdev.h"
#include "evap.h"
#include "cond.h"
#include "comp.h"
#include "r22.h"
#include "Corr.h"
#include "conversion.h"
#include "spline1d.h"
#include "errorlog.h"
#include "spline2d.h"

// these guys enable the API call that opens browers and Excel to show results
// after the calculations are complete.
#include <windows.h>
#include <shellapi.h>

extern ErrorLog errorLog;
extern RefPropTPthT reftplthT;
extern RefPropTPthP reftplthP;
extern WetAirPropTR wair;

// design conditions
const double SH = 12; // F
const double SC = 12; // F
const double ET_SL = 45; // F, at the suction line
const double ET_DT = 45+4; // F, add 4F due to evaporator coild pressure drop
const double CTOA_DL = 18+4; // F, add 4F due to pressure drop in cond coil, 18F is at the liquid line
const double CTOA_LL = 18; // F, at the liquid line
const double AMB = 95; // F
const double RA = 80; // F
const double RARH = 0.51;
const double RAWB = 67; // F

const double MR = 530; // kg/s, refrigerant mass flow rate
const double GOA = 3.82; // kg/s/m^2, outside air mass flux
const double GIA = 2.79; // kg/s/m^2, inside air mass flux
const double DSH = 72; // F, discharge superheat

//------------------------------------------------------------------------------
//------------------------ Compresssor Model -----------------------------------
//------------------------------------------------------------------------------

/***
Exercise compressor model.  Plot its performance over a range of reasonable inout values.
See if it is working well.

Returns 0 if successful.
***/

int ExerciseCompressorDesign()
{
	double mr; // mass flow rate (kg/s)
	HP HPo; // outlet state
	HP HPi; // inlet state
	double power; // compressor power (W)
	MASS m; // refrigerant mass
	double CompPrms[5] = {1,1,1,1,1};

	FILE*fp=fopen("ExerciseCompressor.htm","a");
	if(fp==NULL) return 1;

	fprintf(fp,"<TABLE BORDER=1>\n");
	fprintf(fp,"<CAPTION>Design Conditions</CAPTION>\n");
	fprintf(fp,"<TR>\n");
	fprintf(fp,"<TH colspan=4>Inputs</TH>\n");
	fprintf(fp,"<TH colspan=7>Outputs</TH>\n");
	fprintf(fp,"</TR>\n");
	fprintf(fp,"<TR>\n");
	// inputs
	fprintf(fp,"<TH>Tevap (F)</TH>\n");
	fprintf(fp,"<TH>Tcond (F)</TH>\n");
	fprintf(fp,"<TH>Tsh (F)</TH>\n");
	fprintf(fp,"<TH>Xsl (0-1)</TH>\n");
	// outputs
	fprintf(fp,"<TH>mr (lbs/hr)</TH>\n");
	fprintf(fp,"<TH>power (kW)</TH>\n");
	fprintf(fp,"<TH>T_dl (F)</TH>\n");
	fprintf(fp,"<TH>Tsh_dl (F)</TH>\n");
	fprintf(fp,"<TH>mass (kg)</TH>\n");
	fprintf(fp,"<TH>DH (kJ/kg)</TH>\n");
	fprintf(fp,"<TH>Wdot (BTU/hr)</TH>\n");
	fprintf(fp,"</TR>\n");

	if(errorLog.IsError()) errorLog.ClearError("ExerciseCompressor");

	const double Tai = FtoC(AMB);		// C

	const double Tevap = FtoC(ET_SL);
	const double P1 = reftplthT.Psat(Tevap);
	const double Tsh = DFtoDC(SH);
	const TXP txp1 = {Tevap+Tsh,1,P1};
	const HP hp1 = TXPtoHP(txp1);

	const double Tcond = FtoC(AMB+CTOA_DL);
	const double P2 = reftplthT.Psat(Tcond);

	Compressor("./InputDoc/acmodel.dat",hp1.H,P1,P2,Tai,&mr,&HPo,&HPi,&power,&m,CompPrms);
	if(errorLog.IsError()) {
		errorLog.Add("ExerciseCompressor","Design");
		return 1;
	}

	const double Tsat_dl = reftplthP.Tsat(HPo.P);

	const TXP txp2 = HPtoTXP(HPo);
	if(errorLog.IsError()) {
		errorLog.Add("ExerciseCompressor","txp2");
		return 1;
	}

	fprintf(fp,"<TR>\n");
	// inputs
	fprintf(fp,"<TD>%.1lf</TD>\n",CtoF(Tevap));
	fprintf(fp,"<TD>%.1lf</TD>\n",CtoF(Tcond));
	fprintf(fp,"<TD>%.1lf</TD>\n",DCtoDF(txp1.T-Tevap));
	fprintf(fp,"<TD>%.2lf</TD>\n",txp1.X);
	// outputs
	fprintf(fp,"<TD>%.1lf</TD>\n",KGStoLBMH(mr));
	fprintf(fp,"<TD>%.3lf</TD>\n",power/1000);
	fprintf(fp,"<TD>%.1lf</TD>\n",CtoF(txp2.T));
	fprintf(fp,"<TD>%.1lf</TD>\n",DCtoDF(txp2.T-Tsat_dl));
	fprintf(fp,"<TD>%.4lf</TD>\n",m.m);
	fprintf(fp,"<TD>%.1lf</TD>\n",(HPo.H-HPi.H)/1000.0);
	fprintf(fp,"<TD>%.1lf</TD>\n",WtoBTUH(mr*(HPo.H-HPi.H)));
	fprintf(fp,"</TR>\n");

	fprintf(fp,"</TABLE>\n");
	fprintf(fp,"<P>\n");
	fclose(fp);

	return errorLog.IsError()?1:0;
}

// Vary the evaporating temperature
int ExerciseCompressorVaryTevap()
{
	double mr; // mass flow rate (kg/s)
	HP HPo; // outlet state
	HP HPi; // inlet state
	double power; // compressor power (W)
	MASS m; // refrigerant mass
	double CompPrms[5] = {1,1,1,1,1};

	FILE*fp=fopen("ExerciseCompressor.htm","a");
	if(fp==NULL) return 1;

	fprintf(fp,"<TABLE BORDER=1>\n");
	fprintf(fp,"<CAPTION>Vary Evaporating Temperature</CAPTION>\n");
	fprintf(fp,"<TR>\n");
	fprintf(fp,"<TH colspan=4>Inputs</TH>\n");
	fprintf(fp,"<TH colspan=7>Outputs</TH>\n");
	fprintf(fp,"</TR>\n");
	fprintf(fp,"<TR>\n");
	// inputs
	fprintf(fp,"<TH>Tevap (F)</TH>\n");
	fprintf(fp,"<TH>Tcond (F)</TH>\n");
	fprintf(fp,"<TH>Tsh (F)</TH>\n");
	fprintf(fp,"<TH>Xsl (0-1)</TH>\n");
	// outputs
	fprintf(fp,"<TH>mr (lbs/hr)</TH>\n");
	fprintf(fp,"<TH>power (kW)</TH>\n");
	fprintf(fp,"<TH>T_dl (F)</TH>\n");
	fprintf(fp,"<TH>Tsh_dl (F)</TH>\n");
	fprintf(fp,"<TH>mass (kg)</TH>\n");
	fprintf(fp,"<TH>DH (kJ/kg)</TH>\n");
	fprintf(fp,"<TH>Wdot (BTU/hr)</TH>\n");
	fprintf(fp,"</TR>\n");

	if(errorLog.IsError()) errorLog.ClearError("ExerciseCompressor");

	const double Tai = FtoC(AMB);		// C
	const double Tcond = FtoC(AMB+CTOA_DL);
	const double P2 = reftplthT.Psat(Tcond);

	for(int i=0;i<=60;i+=5)	{

		const double Tevap = FtoC(i);
		const double P1 = reftplthT.Psat(Tevap);
		const double Tsh = DFtoDC(SH);
		const TXP txp1 = {Tevap+Tsh,1,P1};
		const HP hp1 = TXPtoHP(txp1);

		Compressor("./InputDoc/acmodel.dat",hp1.H,P1,P2,Tai,&mr,&HPo,&HPi,&power,&m,CompPrms);
		if(errorLog.IsError()) {
			errorLog.Add("ExerciseCompressor","Vary evaporating temperature");
			break;
		}

		const double Tsat_dl = reftplthP.Tsat(HPo.P);

		const TXP txp2 = HPtoTXP(HPo);
		if(errorLog.IsError()) {
			errorLog.Add("ExerciseCompressor","txp2");
			break;
		}

		fprintf(fp,"<TR>\n");
		// inputs
		fprintf(fp,"<TD>%.1lf</TD>\n",CtoF(Tevap));
		fprintf(fp,"<TD>%.1lf</TD>\n",CtoF(Tcond));
		fprintf(fp,"<TD>%.1lf</TD>\n",DCtoDF(txp1.T-Tevap));
		fprintf(fp,"<TD>%.2lf</TD>\n",txp1.X);
		// outputs
		fprintf(fp,"<TD>%.1lf</TD>\n",KGStoLBMH(mr));
		fprintf(fp,"<TD>%.3lf</TD>\n",power/1000);
		fprintf(fp,"<TD>%.1lf</TD>\n",CtoF(txp2.T));
		fprintf(fp,"<TD>%.1lf</TD>\n",DCtoDF(txp2.T-Tsat_dl));
		fprintf(fp,"<TD>%.4lf</TD>\n",m.m);
		fprintf(fp,"<TD>%.1lf</TD>\n",(HPo.H-HPi.H)/1000.0);
		fprintf(fp,"<TD>%.1lf</TD>\n",WtoBTUH(mr*(HPo.H-HPi.H)));
		fprintf(fp,"</TR>\n");

	}

	fprintf(fp,"</TABLE>\n");
	fprintf(fp,"<P>\n");
	fclose(fp);

	return errorLog.IsError()?1:0;
}

int ExerciseCompressorVaryTcond()
{
	double mr; // mass flow rate (kg/s)
	HP HPo; // outlet state
	HP HPi; // inlet state
	double power; // compressor power (W)
	MASS m; // refrigerant mass
	double CompPrms[5] = {1,1,1,1,1};

	FILE*fp=fopen("ExerciseCompressor.htm","a");
	if(fp==NULL) return 1;

	fprintf(fp,"<TABLE BORDER=1>\n");
	fprintf(fp,"<CAPTION>Vary Condensing Temperature</CAPTION>\n");
	fprintf(fp,"<TR>\n");
	fprintf(fp,"<TH colspan=4>Inputs</TH>\n");
	fprintf(fp,"<TH colspan=7>Outputs</TH>\n");
	fprintf(fp,"</TR>\n");
	fprintf(fp,"<TR>\n");
	// inputs
	fprintf(fp,"<TH>Tevap (F)</TH>\n");
	fprintf(fp,"<TH>Tcond (F)</TH>\n");
	fprintf(fp,"<TH>Tsh (F)</TH>\n");
	fprintf(fp,"<TH>Xsl (0-1)</TH>\n");
	// outputs
	fprintf(fp,"<TH>mr (lbs/hr)</TH>\n");
	fprintf(fp,"<TH>power (kW)</TH>\n");
	fprintf(fp,"<TH>T_dl (F)</TH>\n");
	fprintf(fp,"<TH>Tsh_dl (F)</TH>\n");
	fprintf(fp,"<TH>mass (kg)</TH>\n");
	fprintf(fp,"<TH>DH (kJ/kg)</TH>\n");
	fprintf(fp,"<TH>Wdot (BTU/hr)</TH>\n");
	fprintf(fp,"</TR>\n");

	if(errorLog.IsError()) errorLog.ClearError("ExerciseCompressor");

	const double Tai = FtoC(AMB);		// C

	const double Tevap = FtoC(ET_SL);
	const double P1 = reftplthT.Psat(Tevap);
	const double Tsh = DFtoDC(SH);
	const TXP txp1 = {Tevap+Tsh,1,P1};
	const HP hp1 = TXPtoHP(txp1);

	for(int i=50;i<=150;i+=5)	{

		const double Tcond = FtoC(i);
		const double P2 = reftplthT.Psat(Tcond);

		Compressor("./InputDoc/acmodel.dat",hp1.H,P1,P2,Tai,&mr,&HPo,&HPi,&power,&m,CompPrms);
		if(errorLog.IsError()) {
			errorLog.Add("ExerciseCompressor","Vary condensing temperature");
			break;
		}

		const double Tsat_dl = reftplthP.Tsat(HPo.P);

		const TXP txp2 = HPtoTXP(HPo);
		if(errorLog.IsError()) {
			errorLog.Add("ExerciseCompressor","txp2");
			break;
		}

		fprintf(fp,"<TR>\n");
		// inputs
		fprintf(fp,"<TD>%.1lf</TD>\n",CtoF(Tevap));
		fprintf(fp,"<TD>%.1lf</TD>\n",CtoF(Tcond));
		fprintf(fp,"<TD>%.1lf</TD>\n",DCtoDF(txp1.T-Tevap));
		fprintf(fp,"<TD>%.2lf</TD>\n",txp1.X);
		// outputs
		fprintf(fp,"<TD>%.1lf</TD>\n",KGStoLBMH(mr));
		fprintf(fp,"<TD>%.3lf</TD>\n",power/1000);
		fprintf(fp,"<TD>%.1lf</TD>\n",CtoF(txp2.T));
		fprintf(fp,"<TD>%.1lf</TD>\n",DCtoDF(txp2.T-Tsat_dl));
		fprintf(fp,"<TD>%.4lf</TD>\n",m.m);
		fprintf(fp,"<TD>%.1lf</TD>\n",(HPo.H-HPi.H)/1000.0);
		fprintf(fp,"<TD>%.1lf</TD>\n",WtoBTUH(mr*(HPo.H-HPi.H)));
		fprintf(fp,"</TR>\n");

	}

	fprintf(fp,"</TABLE>\n");
	fprintf(fp,"<P>\n");
	fclose(fp);

	return errorLog.IsError()?1:0;
}

int ExerciseCompressorVaryTsh()
{
	double mr; // mass flow rate (kg/s)
	HP HPo; // outlet state
	HP HPi; // inlet state
	double power; // compressor power (W)
	MASS m; // refrigerant mass
	double CompPrms[5] = {1,1,1,1,1};

	FILE*fp=fopen("ExerciseCompressor.htm","a");
	if(fp==NULL) return 1;

	fprintf(fp,"<TABLE BORDER=1>\n");
	fprintf(fp,"<CAPTION>Vary Inlet Superheat</CAPTION>\n");
	fprintf(fp,"<TR>\n");
	fprintf(fp,"<TH colspan=4>Inputs</TH>\n");
	fprintf(fp,"<TH colspan=7>Outputs</TH>\n");
	fprintf(fp,"</TR>\n");
	fprintf(fp,"<TR>\n");
	// inputs
	fprintf(fp,"<TH>Tevap (F)</TH>\n");
	fprintf(fp,"<TH>Tcond (F)</TH>\n");
	fprintf(fp,"<TH>Tsh (F)</TH>\n");
	fprintf(fp,"<TH>Xsl (0-1)</TH>\n");
	// outputs
	fprintf(fp,"<TH>mr (lbs/hr)</TH>\n");
	fprintf(fp,"<TH>power (kW)</TH>\n");
	fprintf(fp,"<TH>T_dl (F)</TH>\n");
	fprintf(fp,"<TH>Tsh_dl (F)</TH>\n");
	fprintf(fp,"<TH>mass (kg)</TH>\n");
	fprintf(fp,"<TH>DH (kJ/kg)</TH>\n");
	fprintf(fp,"<TH>Wdot (BTU/hr)</TH>\n");
	fprintf(fp,"</TR>\n");

	if(errorLog.IsError()) errorLog.ClearError("ExerciseCompressor");

	const double Tai = FtoC(AMB);		// C

	const double Tcond = FtoC(AMB+CTOA_DL);
	const double P2 = reftplthT.Psat(Tcond);

	// Vary quality - two-phase compressor inlet conditions
	for(int i=70;i<=100;i+=2)	{

		const double Tevap = FtoC(ET_SL);
		const double P1 = reftplthT.Psat(Tevap);
		const TXP txp1 = {Tevap,0.01*(double)i,P1};
		const HP hp1 = TXPtoHP(txp1);

		Compressor("./InputDoc/acmodel.dat",hp1.H,P1,P2,Tai,&mr,&HPo,&HPi,&power,&m,CompPrms);
		if(errorLog.IsError()) {
			errorLog.Add("ExerciseCompressor","Vary superheat");
			break;
		}

		const double Tsat_dl = reftplthP.Tsat(HPo.P);

		const TXP txp2 = HPtoTXP(HPo);
		if(errorLog.IsError()) {
			errorLog.Add("ExerciseCompressor","txp2");
			break;
		}

		fprintf(fp,"<TR>\n");
		// inputs
		fprintf(fp,"<TD>%.1lf</TD>\n",CtoF(Tevap));
		fprintf(fp,"<TD>%.1lf</TD>\n",CtoF(Tcond));
		fprintf(fp,"<TD>%.1lf</TD>\n",DCtoDF(txp1.T-Tevap));
		fprintf(fp,"<TD>%.2lf</TD>\n",txp1.X);
		// outputs
		fprintf(fp,"<TD>%.1lf</TD>\n",KGStoLBMH(mr));
		fprintf(fp,"<TD>%.3lf</TD>\n",power/1000);
		fprintf(fp,"<TD>%.1lf</TD>\n",CtoF(txp2.T));
		fprintf(fp,"<TD>%.1lf</TD>\n",DCtoDF(txp2.T-Tsat_dl));
		fprintf(fp,"<TD>%.4lf</TD>\n",m.m);
		fprintf(fp,"<TD>%.1lf</TD>\n",(HPo.H-HPi.H)/1000.0);
		fprintf(fp,"<TD>%.1lf</TD>\n",WtoBTUH(mr*(HPo.H-HPi.H)));
		fprintf(fp,"</TR>\n");

	}

	// Vary superheat
	for(int i=0;i<=100;i+=5)	{

		const double Tevap = FtoC(45);
		const double P1 = reftplthT.Psat(Tevap);
		const double Tsh = DFtoDC(i);
		const TXP txp1 = {Tevap+Tsh,1,P1};
		const HP hp1 = TXPtoHP(txp1);

		Compressor("./InputDoc/acmodel.dat",hp1.H,P1,P2,Tai,&mr,&HPo,&HPi,&power,&m,CompPrms);
		if(errorLog.IsError()) {
			errorLog.Add("ExerciseCompressor","Vary superheat");
			break;
		}

		const double Tsat_dl = reftplthP.Tsat(HPo.P);

		const TXP txp2 = HPtoTXP(HPo);
		if(errorLog.IsError()) {
			errorLog.Add("ExerciseCompressor","txp2");
			break;
		}

		fprintf(fp,"<TR>\n");
		// inputs
		fprintf(fp,"<TD>%.1lf</TD>\n",CtoF(Tevap));
		fprintf(fp,"<TD>%.1lf</TD>\n",CtoF(Tcond));
		fprintf(fp,"<TD>%.1lf</TD>\n",DCtoDF(txp1.T-Tevap));
		fprintf(fp,"<TD>%.2lf</TD>\n",txp1.X);
		// outputs
		fprintf(fp,"<TD>%.1lf</TD>\n",KGStoLBMH(mr));
		fprintf(fp,"<TD>%.3lf</TD>\n",power/1000);
		fprintf(fp,"<TD>%.1lf</TD>\n",CtoF(txp2.T));
		fprintf(fp,"<TD>%.1lf</TD>\n",DCtoDF(txp2.T-Tsat_dl));
		fprintf(fp,"<TD>%.4lf</TD>\n",m.m);
		fprintf(fp,"<TD>%.1lf</TD>\n",(HPo.H-HPi.H)/1000.0);
		fprintf(fp,"<TD>%.1lf</TD>\n",WtoBTUH(mr*(HPo.H-HPi.H)));
		fprintf(fp,"</TR>\n");

	}

	fprintf(fp,"</TABLE>\n");
	fprintf(fp,"<P>\n");
	fclose(fp);

	return errorLog.IsError()?1:0;
}

int ExerciseCompressor()
{
	FILE*fp=fopen("ExerciseCompressor.htm","w");
	if(fp==NULL) return 1;
	fprintf(fp,"<B>Exercise Compressor</B>\n");
	fclose(fp);

	ExerciseCompressorDesign();
	ExerciseCompressorVaryTevap();
	ExerciseCompressorVaryTcond();
	ExerciseCompressorVaryTsh();

	return 0;
}

//------------------------------------------------------------------------------
//--------------------------- Condenser Model ----------------------------------
//------------------------------------------------------------------------------

// Vary inlet air temp
int ExerciseCondenserVary1D(int vary)
{
	HP hpo;
	MASS m;
	double Tao;
	CGP Cond_struc;
	double CondPrms[5] = {1,1,1,1,1};

	FILE*fp=fopen("ExerciseCondenser.htm","a");
	if(fp==NULL) return 1;

	fprintf(fp,"<TABLE BORDER=1>\n");
	switch(vary) {
		case 0:
			fprintf(fp,"<CAPTION>Design conditions</CAPTION>\n");
			break;
		case 1:
			fprintf(fp,"<CAPTION>Vary inlet air temperature</CAPTION>\n");
			break;
		case 2:
			fprintf(fp,"<CAPTION>Vary inlet air flow rate</CAPTION>\n");
			break;
		case 3:
			fprintf(fp,"<CAPTION>Vary CTOA</CAPTION>\n");
			break;
		case 4:
			fprintf(fp,"<CAPTION>Vary ref mass flow rate</CAPTION>\n");
			break;
		case 5:
			fprintf(fp,"<CAPTION>Vary NSeg</CAPTION>\n");
			break;
		default:
			fprintf(fp,"<CAPTION>ERROR</CAPTION>\n");
			break;
		};

	fprintf(fp,"<TR>\n");
	fprintf(fp,"<TH colspan=6>Inputs</TH>\n");
	fprintf(fp,"<TH colspan=9>Outputs</TH>\n");
	fprintf(fp,"</TR>\n");
	fprintf(fp,"<TR>\n");
	// inputs
	fprintf(fp,"<TH>Tai (F)</TH>\n");
	fprintf(fp,"<TH>CFM (cfm)</TH>\n");
	fprintf(fp,"<TH>mr (lbm/hr)</TH>\n");
	fprintf(fp,"<TH>CTOA_dl (F)</TH>\n");
	fprintf(fp,"<TH>Tsh_dl (F)</TH>\n");
	fprintf(fp,"<TH>X_dl (0-1)</TH>\n");
	// outputs
	fprintf(fp,"<TH>CTOA_ll (F)</TH>\n");
	fprintf(fp,"<TH>Tsc_ll (F)</TH>\n");
	fprintf(fp,"<TH>X_ll (0-1)</TH>\n");
	fprintf(fp,"<TH>dP (psi)</TH>\n");
	fprintf(fp,"<TH>DTa (F)</TH>\n");
	fprintf(fp,"<TH>mass (kg)</TH>\n");
	fprintf(fp,"<TH>Qdot (BTU/hr)</TH>\n");
	fprintf(fp,"<TH>UA (BTU/hr/F)</TH>\n");
	fprintf(fp,"<TH>Hout (kJ/kg)</TH>\n");
	fprintf(fp,"</TR>\n");

	if(errorLog.IsError()) errorLog.ClearError("ExerciseCompressor");

	const int Imax = 20;
	for(int i=0;i<=Imax;i++)	{

		// inlet air state
		double Tai = FtoC(AMB);
		// air mass flux (kg/s/m^2)
		double Goa = GOA;
		// refrigerant mass flow rate
		double mr = LBMHtoKGS(MR);
		// inlet refrigerant state
		double ctoa = DFtoDC(CTOA_DL); // condensing temperature over ambient
		const double Tsh = DFtoDC(DSH);
		int NSeg = -1;

		const double frac = (double)i/(double)Imax;
		switch(vary) {
		case 0:
			if(i>0) continue;
			break;
		case 1:
			Tai = FtoC(70+80*frac);
			break;
		case 2:
			if(i==0) continue;
			Goa = GOA*2*frac;
			break;
		case 3:
			if(i<3) continue;
			ctoa = DFtoDC(50.0*frac);
			break;
		case 4:
			if(i==0) continue;
			mr = LBMHtoKGS(frac*2*MR);
			break;
		case 5:
			if(i==0) continue;
			NSeg = (int)(40.0*frac);
			break;
		default:
			return 1;
		};

		const double Tsat_dl = Tai+ctoa;
		const double P_dl = reftplthT.Psat(Tsat_dl);
		const TXP txpi = {Tsat_dl+Tsh,1,P_dl};
		const HP hpi = TXPtoHP(txpi);

		Condenser("./InputDoc/acmodel.dat",mr,hpi,Tai,Goa,&hpo,&Tao,&m,&Cond_struc,CondPrms,NSeg);
		if(errorLog.IsError()) {
			errorLog.Add("ExerciseCondenser","Vary");
			break;
		}

		const TXP txpo = HPtoTXP(hpo);
		if(errorLog.IsError()) {
			errorLog.Add("ExerciseCondenser","txpo");
			break;
		}
		const double Tsat_ll = reftplthP.Tsat(txpo.P);
		if(errorLog.IsError()) {
			errorLog.Add("ExerciseCondenser","Tsat_ll");
			break;
		}

		// Calculate the total heat transfer and then the UA = Qdot/DT
		const double Qdot = WtoBTUH(mr*(hpi.H-hpo.H)); // BTU/hr - for all 4 sections of the heat exchanger
		const double DT = DCtoDF(Tsat_ll - Tai); // F
		const double UA = Qdot/DT; // BTU/hr/F - for all 4 sections of the heat exchanger

		fprintf(fp,"<TR>\n");
		// inputs
		fprintf(fp,"<TD>%.1lf</TD>\n",CtoF(Tai));
		fprintf(fp,"<TD>%.1lf</TD>\n",Cond_struc.cfma);
		fprintf(fp,"<TD>%.1lf</TD>\n",KGStoLBMH(mr));
		fprintf(fp,"<TD>%.1lf</TD>\n",DCtoDF(Tsat_dl-Tai));
		fprintf(fp,"<TD>%.1lf</TD>\n",DCtoDF(txpi.T-Tsat_dl));
		fprintf(fp,"<TD>%.2lf</TD>\n",txpi.X);
		// outputs
		fprintf(fp,"<TD>%.1lf</TD>\n",DCtoDF(Tsat_ll-Tai));// CTOA_ll
		fprintf(fp,"<TD>%.1lf</TD>\n",DCtoDF(Tsat_ll-txpo.T));// Tsc_ll
		fprintf(fp,"<TD>%.2lf</TD>\n",txpo.X);
		fprintf(fp,"<TD>%.2lf</TD>\n",DKPAtoDPSIG(txpi.P-txpo.P));
		fprintf(fp,"<TD>%.1lf</TD>\n",DCtoDF(Tao-Tai));
		fprintf(fp,"<TD>%.4lf</TD>\n",m.m);
		fprintf(fp,"<TD>%.0lf</TD>\n",Qdot);
		fprintf(fp,"<TD>%.0lf</TD>\n",UA);
		fprintf(fp,"<TD>%.2lf</TD>\n",hpo.H/1000.0);
		fprintf(fp,"</TR>\n");

	}

	fprintf(fp,"</TABLE>\n");
	fprintf(fp,"<P>\n");
	fclose(fp);

	return errorLog.IsError()?1:0;
}

int ExerciseCondenser()
{
	FILE*fp=fopen("ExerciseCondenser.htm","w");
	if(fp==NULL) return 1;
	fprintf(fp,"<B>Exercise Condenser</B>\n");
	fclose(fp);

	for(int i=0;i<6;i++) {
		ExerciseCondenserVary1D(i);
	}

	return 0;
}

//------------------------------------------------------------------------------
//--------------------------- Evaporator Model ----------------------------------
//------------------------------------------------------------------------------

int ExerciseEvaporatorVary1D(int vary)
{
	HP hpi;
	MASS m;
	TP TPo;
	ETdim Evap_struc;
	double Aflow;
	double Evap_prms[5] = {1,1,1,1,1};

	FILE*fp=fopen("ExerciseEvaporator.htm","a");
	if(fp==NULL) return 1;

	fprintf(fp,"<TABLE BORDER=1>\n");
	switch(vary) {
		case 0:
			fprintf(fp,"<CAPTION>Design conditions</CAPTION>\n");
			break;
		case 1:
			fprintf(fp,"<CAPTION>Vary inlet air temperature</CAPTION>\n");
			break;
		case 2:
			fprintf(fp,"<CAPTION>Vary inlet air flow rate</CAPTION>\n");
			break;
		case 3:
			fprintf(fp,"<CAPTION>Vary ETBRA</CAPTION>\n");
			break;
		case 4:
			fprintf(fp,"<CAPTION>Vary ref mass flow rate</CAPTION>\n");
			break;
		case 5:
			fprintf(fp,"<CAPTION>Vary NSeg</CAPTION>\n");
			break;
		default:
			fprintf(fp,"<CAPTION>ERROR</CAPTION>\n");
			break;
		};

	fprintf(fp,"<TR>\n");
	fprintf(fp,"<TH colspan=7>Inputs</TH>\n");
	fprintf(fp,"<TH colspan=9>Outputs</TH>\n");
	fprintf(fp,"</TR>\n");
	fprintf(fp,"<TR>\n");
	// inputs
	fprintf(fp,"<TH>Tai (F), RHai (%%)</TH>\n");
	fprintf(fp,"<TH>CFM (cfm)</TH>\n");
	fprintf(fp,"<TH>mr (lbm/hr)</TH>\n");
	fprintf(fp,"<TH>ETBRA_sl (F)</TH>\n");
	fprintf(fp,"<TH>Tsh_sl (F)</TH>\n");
	fprintf(fp,"<TH>X_sl (0-1)</TH>\n");
	fprintf(fp,"<TH>NSeg</TH>\n");
	// outputs
	fprintf(fp,"<TH>ETBRA_dt (F)</TH>\n");
	fprintf(fp,"<TH>X_dt (0-1)</TH>\n");
	fprintf(fp,"<TH>dP (psi)</TH>\n");
	fprintf(fp,"<TH>DTa (F)</TH>\n");
	fprintf(fp,"<TH>DTwb (F)</TH>\n");
	fprintf(fp,"<TH>mass (kg)</TH>\n");
	fprintf(fp,"<TH>Qdot (BTU/hr)</TH>\n");
	fprintf(fp,"<TH>UA (BTU/hr/F)</TH>\n");
	fprintf(fp,"<TH>Hin (kJ/kg)</TH>\n");
	fprintf(fp,"</TR>\n");

	if(errorLog.IsError()) errorLog.ClearError("ExerciseCompressor");

	const int Imax = 20;
	for(int i=0;i<=Imax;i++)	{

		// inlet air state
		TP TPi = {FtoC(RA),RARH};
		// air mass flux (kg/s/m^2)
		double Ga = GIA;
		// refrigerant mass flow rate
		double mr = LBMHtoKGS(MR);
		// outlet refrigerant state
		double ETBRA = DFtoDC(RA-ET_SL); // evaporating temperature below return air
		const double Tsh = DFtoDC(SH); // suction line superheat
		int NSeg = -1;

		const double frac = (double)i/(double)Imax;
		switch(vary) {
		case 0:
			if(i>0) continue;
			break;
		case 1:
			TPi.T = FtoC(90-50*frac);
			break;
		case 2:
			if(i==0) continue;
			Ga = GIA*2*frac;
			break;
		case 3:
			if(i<3) continue;
			ETBRA = DFtoDC(100.0*frac);
			break;
		case 4:
			mr = LBMHtoKGS((1-frac)*2*MR);
			break;
		case 5:
			if(i==0) continue;
			NSeg = (int)(40.0*frac);
			break;
		default:
			return 1;
		};

		const double Tsat_sl = TPi.T-ETBRA;
		const double P_sl = reftplthT.Psat(Tsat_sl);
		const TXP txpo = {Tsat_sl+Tsh,1,P_sl};
		HP hpo = TXPtoHP(txpo);
		const double Twbi = wair.WetBulb(TPi.T,TPi.P);

		Evaporator("./InputDoc/acmodel.dat",mr/4,&hpo,Ga,TPi,&hpi,&TPo,&m,&Aflow,&Evap_struc,Evap_prms,NSeg);
		if(errorLog.IsError()) {
			errorLog.Add("ExerciseEvaporator","Evaporator");
			break;
		}

		const TXP txpi = HPtoTXP(hpi);
		if(errorLog.IsError()) {
			errorLog.Add("ExerciseEvaporator","txpi");
			break;
		}
		const double Tsat_dt = reftplthP.Tsat(txpi.P);
		if(errorLog.IsError()) {
			errorLog.Add("ExerciseEvaporator","Tsat_dt");
			break;
		}
		const double Twbo = wair.WetBulb(TPo.T,TPo.P);

		// Calculate the total heat transfer and then the UA = Qdot/DT
		const double Qdot = WtoBTUH(mr*(hpo.H-hpi.H)); // BTU/hr - for all 4 sections of the heat exchanger
		const double DT = DCtoDF(TPi.T - Tsat_sl); // F
		const double UA = Qdot/DT; // BTU/hr/F - for all 4 sections of the heat exchanger

		fprintf(fp,"<TR>\n");
		// inputs
		fprintf(fp,"<TD>%.1lfF, %.1lf%%</TD>\n",CtoF(TPi.T),100.0*TPi.P);
		fprintf(fp,"<TD>%.1lf</TD>\n",4*Evap_struc.cfma);
		fprintf(fp,"<TD>%.1lf</TD>\n",KGStoLBMH(mr));
		fprintf(fp,"<TD>%.1lf</TD>\n",DCtoDF(TPi.T-Tsat_sl));
		fprintf(fp,"<TD>%.1lf</TD>\n",DCtoDF(txpo.T-Tsat_sl));
		fprintf(fp,"<TD>%.2lf</TD>\n",txpo.X);
		if(NSeg>0) {
			fprintf(fp,"<TD>%d</TD>\n",NSeg);
		} else {
			fprintf(fp,"<TD>Default</TD>\n");
		}
		// outputs
		fprintf(fp,"<TD>%.1lf</TD>\n",DCtoDF(TPi.T-Tsat_dt));// ETBRA_dt
		fprintf(fp,"<TD>%.2lf</TD>\n",txpi.X);
		fprintf(fp,"<TD>%.2lf</TD>\n",DKPAtoDPSIG(txpi.P-txpo.P));
		fprintf(fp,"<TD>%.1lf</TD>\n",DCtoDF(TPi.T-TPo.T));
		fprintf(fp,"<TD>%.1lf</TD>\n",DCtoDF(Twbi-Twbo));
		fprintf(fp,"<TD>%.4lf</TD>\n",4*m.m);
		fprintf(fp,"<TD>%.0lf</TD>\n",Qdot);
		fprintf(fp,"<TD>%.0lf</TD>\n",UA);
		fprintf(fp,"<TD>%.2lf</TD>\n",hpi.H/1000.0);
		fprintf(fp,"</TR>\n");

	}

	fprintf(fp,"</TABLE>\n");
	fprintf(fp,"<P>\n");
	fclose(fp);

	return errorLog.IsError()?1:0;
}

int ExerciseEvaporator()
{
	FILE*fp=fopen("ExerciseEvaporator.htm","w");
	if(fp==NULL) return 1;
	fprintf(fp,"<B>Exercise Evaporator</B>\n");
	fclose(fp);

	for(int i=0;i<6;i++) {
		ExerciseEvaporatorVary1D(i);
	}

	return 0;
}

//------------------------------------------------------------------------------
//--------------------------- Expansion Device Model ---------------------------
//------------------------------------------------------------------------------

int ExerciseExpDevVary1D(int vary)
{
	MASS m;
	double P_up;
	double Ones[5] = {1,1,1,1,1};

	FILE*fp=fopen("ExerciseExpDev.htm","a");
	if(fp==NULL) return 1;

	fprintf(fp,"<TABLE BORDER=1>\n");
	switch(vary) {
		case 0:
			fprintf(fp,"<CAPTION>Design conditions</CAPTION>\n");
			break;
		case 1:
			fprintf(fp,"<CAPTION>Vary inlet state</CAPTION>\n");
			break;
		case 2:
			fprintf(fp,"<CAPTION>Vary ref flow rate</CAPTION>\n");
			break;
		default:
			fprintf(fp,"<CAPTION>ERROR</CAPTION>\n");
			break;
		};

	fprintf(fp,"<TR>\n");
	fprintf(fp,"<TH colspan=6>Inputs</TH>\n");
	fprintf(fp,"<TH colspan=2>Outputs</TH>\n");
	fprintf(fp,"</TR>\n");
	fprintf(fp,"<TR>\n");
	// inputs
	fprintf(fp,"<TH>Tsat_ll (F)</TH>\n");
	fprintf(fp,"<TH>Tsat_dt (F)</TH>\n");
	fprintf(fp,"<TH>Tsc (F)</TH>\n");
	fprintf(fp,"<TH>X_ll (0-1)</TH>\n");
	fprintf(fp,"<TH>mr (lbm/hr)</TH>\n");
	fprintf(fp,"<TH>P_ll-P_dt (psi)</TH>\n");
	// outputs
	fprintf(fp,"<TH>txp_ll.P-P_up (psi)</TH>\n");
	fprintf(fp,"<TH>mass (kg)</TH>\n");
	fprintf(fp,"</TR>\n");

	if(errorLog.IsError()) errorLog.ClearError("ExerciseCompressor");

	const int Imax = 20;
	for(int i=0;i<=Imax;i++)	{

		// liquid line
		double Tsc = DFtoDC(SC);
		const double Tsat_ll = FtoC(AMB+CTOA_LL);

		// low side after valve
		const double Tsat_dt = FtoC(ET_DT);

		double mr = LBMHtoKGS(MR);

		const double frac = (double)i/(double)Imax;
		switch(vary) {
		case 0:
			if(i>0) continue;
			break;
		case 1:
			Tsc = DFtoDC(20*frac);
			break;
		case 2:
			if(i==0) continue;
			mr = LBMHtoKGS(frac*2*MR);
			break;
		default:
			return 1;
		};

		const double P_ll = reftplthT.Psat(Tsat_ll);
		const TXP txp_ll = {Tsat_ll-Tsc,0,P_ll};
		const HP hp_ll = TXPtoHP(txp_ll);

		const double P_dt = reftplthT.Psat(Tsat_dt); // pressure - distribution tubes

		ExpDev("./InputDoc/acmodel.dat",1,hp_ll.H,P_dt,mr/4,&P_up,&m,Ones); 
		if(errorLog.IsError()) {
			errorLog.Add("ExerciseEvaporator","Evaporator");
			break;
		}

		fprintf(fp,"<TR>\n");
		// inputs
		fprintf(fp,"<TD>%.1lfF</TD>\n",CtoF(Tsat_ll));
		fprintf(fp,"<TD>%.1lfF</TD>\n",CtoF(Tsat_dt));
		fprintf(fp,"<TD>%.1lfF</TD>\n",DCtoDF(Tsc));
		fprintf(fp,"<TD>%.4f</TD>\n",txp_ll.X);
		fprintf(fp,"<TD>%.1lf</TD>\n",KGStoLBMH(mr));
		fprintf(fp,"<TD>%.2lf</TD>\n",DKPAtoDPSIG(P_ll-P_dt));
		// outputs
		fprintf(fp,"<TD>%.2lf</TD>\n",DKPAtoDPSIG(txp_ll.P-P_up));
		fprintf(fp,"<TD>%.4lf</TD>\n",m.m);
		fprintf(fp,"</TR>\n");

	}

	fprintf(fp,"</TABLE>\n");
	fprintf(fp,"<P>\n");
	fclose(fp);

	return errorLog.IsError()?1:0;
}

int ExerciseExpDev()
{
	FILE*fp=fopen("ExerciseExpDev.htm","w");
	if(fp==NULL) return 1;
	fprintf(fp,"<B>Exercise ExpDev</B>\n");
	fclose(fp);

	for(int i=0;i<3;i++) {
		ExerciseExpDevVary1D(i);
	}

	return 0;
}


//---------------------------------------------------------------------------
#define EXERCISE_COMPRESSOR
#define EXERCISE_CONDENSER
#define EXERCISE_EVAPORATOR
#define EXERCISE_EXPDEV
#define EXERCISE_CORR

int ExerciseComponents()
{
	int errorCode;

#ifdef EXERCISE_COMPRESSOR
	errorCode = ExerciseCompressor();
	if(errorCode>0 || errorLog.IsError()) {
		char msg[256];
		sprintf(msg,"ExerciseCompressor: errorCode=%d",errorCode);
		errorLog.Add("main",msg);
		ShellExecute(NULL, "open", "ErrorLog.htm", NULL, NULL, SW_SHOWNORMAL);
	}
	ShellExecute(NULL, "open", "ExerciseCompressor.htm", NULL, NULL, SW_SHOWNORMAL);
#endif

#ifdef EXERCISE_CONDENSER
	errorCode = ExerciseCondenser();
	if(errorCode>0 || errorLog.IsError()) {
		char msg[256];
		sprintf(msg,"ExerciseCondenser: errorCode=%d",errorCode);
		errorLog.Add("main",msg);
		ShellExecute(NULL, "open", "ErrorLog.htm", NULL, NULL, SW_SHOWNORMAL);
	}
	ShellExecute(NULL, "open", "ExerciseCondenser.htm", NULL, NULL, SW_SHOWNORMAL);
#endif

#ifdef EXERCISE_EVAPORATOR
	errorCode = ExerciseEvaporator();
	if(errorCode>0 || errorLog.IsError()) {
		char msg[256];
		sprintf(msg,"ExerciseEvaporator: errorCode=%d",errorCode);
		errorLog.Add("main",msg);
		ShellExecute(NULL, "open", "ErrorLog.htm", NULL, NULL, SW_SHOWNORMAL);
	}
	ShellExecute(NULL, "open", "ExerciseEvaporator.htm", NULL, NULL, SW_SHOWNORMAL);
#endif

#ifdef EXERCISE_EXPDEV
	errorCode = ExerciseExpDev();
	if(errorCode>0 || errorLog.IsError()) {
		char msg[256];
		sprintf(msg,"ExerciseExpDev: errorCode=%d",errorCode);
		errorLog.Add("main",msg);
		ShellExecute(NULL, "open", "ErrorLog.htm", NULL, NULL, SW_SHOWNORMAL);
	}
	ShellExecute(NULL, "open", "ExerciseExpDev.htm", NULL, NULL, SW_SHOWNORMAL);
#endif


	return 0;
}