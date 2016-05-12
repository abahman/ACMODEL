#ifndef _MAINFUNC
#define _MAINFUNC

//B.S----------------------
#include "mass.h"
#include "cond.h"
#include "evap.h"

void InitializeComponents(void);
//----------------------------B.S.

// data structure providing inputs to MainFunc()
struct InVars {
	double Tai;				// condenser air inlet temperature (C)
	double Goa;				// condenser air mass flux
	double Gra;				// evaporator air mass flux
	double Charge;			// system charge level (kg)
	TP TPi;					// evaporator air inlet state (temp + humidity)
	double Dadj;			// expansion device diameter scaler
	double nvAdj;			// compressor capacity scaler
	char filename[128];		// filename describing equipment

	//B.S--------------------
	double ConvTol;//overal residual tolerance
	CGP Cond;//condenser struct
	ETdim Evap;//evaporator struct
	double Tsup;//superheat degree at the suction
	double Tsub;//subcooling degree at the exit of liquid line
	double CompT;
	double Qual;//quality at the exit of the liquid line
	double subset;//subcooling degree set
	double ChargMass;//charge inventory of the system
	int Lumped;//when Lumped=1 conduct the lumped heat exchanger model analysis
	int EvaSign;//when EvaSign=0 conduct the detailed evaporator model analysis, otherwise conduct the moving boundary evaporator model analysis
	int CondSign;//when CondSign=0 conduct the detailed condenser model analysis, otherwise conduct the moving boundary condenser model analysis
	double CE;//residual
	HP HP5;//state upstream the orifice
	HP HP6;//state downstream the orifice
	double Tsup8;
	double mr;//mass flow rate

	//definitions
	int REV_CAL;
	int SPEC_SUB;
	int PISTON;
	double T[9];
	double P[9];
	double X[4];	//quality at inlet and outlet of the suction line, and liquid line, order:1, 4, 5, 8
	double Cadj;	// condenser airside heat transfer adjustment			
	double Eadj;	// evaporator airside heat transfer adjustment
	double ComP;
};

int MainFunc1(double*,double*,void*);
int MainFunc2(double*,double*,void*);
int MainFunc3(double*,double*,void*);
int MainFunc_Validating(void*);
int MainFunc_Tuning(void*);
void PrintCycleTimes();

#endif
