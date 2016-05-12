#include "stdafx.h"

#include "r22.h"
#include "air.h"
#include "corr.h"
#include "cmine.h"
#include <stdio.h>
#include <math.h>
#include "errorlog.h"
#include "spline1d.h"

extern ErrorLog errorLog;

extern AirProp air;

// #define EXTRA_REPORTING

// Calculates the product of minimum heat capacity (Cmin) and
// effectiveness (e) for transferring heat between refrigerant
// and dry air in crossflow.

//-----------------------------------------------------B.S.
//the Cmine is treated separately for dry condition and wet condition, 
//for there may exist difference between them in some methods.
//here it is simply to distinguish them with different air specific heats. 
//-----------------------------------------------------B.S.

//B.S.-----------------------------------------------------
//The following function is only for calculating the Cmine under dry condition

double
CmineCrossFlow_dry(double R,//overall heat resistance
				   double mr,//refrigerant mass flow rate
				   double ma,//airflow mass flow rate
				   TXP TXP,//refrigerant state
				   double Ta)//air temperature
{
	if (IsTwoPhase(TXP.X)) {
		//
		// all heat exchangers
		// refrigerant has infinite heat capacity rate
		// Cr=Cmin/Cmax=0 (Cmax=Cr=inf)
		double Cmin = ma*air.Cp(Ta);
		if(errorLog.IsError()) {
			errorLog.Add("CmineCrossFlow");
			return -1;
		}
		double Ntu = 1/(R*Cmin);
		double e = 1-exp(-1*Ntu);
		return Cmin*e;
	} else {
		double Ca = ma*air.Cp(Ta);
		if(errorLog.IsError()) {
			errorLog.Add("CmineCrossFlow");
			return -1;
		}
		double Cr = mr*PropertyTXPtr(SPEC,TXP);
		if(errorLog.IsError()) {
			errorLog.Add("CmineCrossFlow");
			return -1;
		}
		if (Ca>Cr) {
			//
			// crossflow heat exchanger
			// refrigerant has minimum heat capacity rate
			// air unmixed, refrigerant mixed
			double Cratio = Cr/Ca;
			double Ntu = 1/(R*Cr);
			double l = 1-exp(-1*Ntu*Cratio);
			double e = 1-exp(-1*l/Cratio);
			return Cr*e;
		} else {
			//
			// crossflow heat exchanger
			// air has minimum heat capacity rate
			// air unmixed, refrigerant mixed
			double Cratio = Ca/Cr;
			double Ntu = 1/(R*Ca);
			double l = -1*Cratio*(1-exp(-1*Ntu));
			double e = (1/Cratio)*(1-exp(l));
			return Ca*e;
		}
	}
}
//-----------------------------------------------------B.S.


//B.S.-----------------------------------------------------
//The following function is only for calculating the Cmine under wet condition, 
//however, it is not functioned now, and only the specific heat of the wet air is 
//added in

double
CmineCrossFlow_wet(double R,//overall thermal resistance
				   double mr,//refrigerant mass flow rate
				   double ma,//air mass flow rate
				   TXP TXP,//refrigerant state
				   double Ta,//air temperature
				   double CP_A)//wet air specific heat
{
	if (IsTwoPhase(TXP.X)) {
		//
		// all heat exchangers
		// refrigerant has infinite heat capacity rate
		// Cr=Cmin/Cmax=0 (Cmax=Cr=inf)
		double Cmin = ma*CP_A;
		if(errorLog.IsError()) {
			errorLog.Add("CmineCrossFlow");
			return -1;
		}
		double Ntu = 1/(R*Cmin);
		double e = 1-exp(-1*Ntu);
		return Cmin*e;
	} else {
		double Ca = ma*CP_A;
		if(errorLog.IsError()) {
			errorLog.Add("CmineCrossFlow");
			return -1;
		}
		double Cr = mr*PropertyTXPtr(SPEC,TXP);
		if(errorLog.IsError()) {
			errorLog.Add("CmineCrossFlow");
			return -1;
		}
		if (Ca>Cr) {
			//
			// crossflow heat exchanger
			// refrigerant has minimum heat capacity rate
			// air unmixed, refrigerant mixed
			double Cratio = Cr/Ca;
			double Ntu = 1/(R*Cr);
			double l = 1-exp(-1*Ntu*Cratio);
			double e = 1-exp(-1*l/Cratio);
			return Cr*e;
		} else {
			//
			// crossflow heat exchanger
			// air has minimum heat capacity rate
			// air unmixed, refrigerant mixed
			double Cratio = Ca/Cr;
			double Ntu = 1/(R*Ca);
			double l = -1*Cratio*(1-exp(-1*Ntu));
			double e = (1/Cratio)*(1-exp(l));
			return Ca*e;
		}
	}
}
//-----------------------------------------------------B.S.


// Calculates the product of minimum heat capacity (Cmin) and
// effectiveness (e) for transferring heat between two refrigerants
// in counterflow.
//
// UA = 1/R (W/K)
// mr1,mr2 (kg/s)
// TXP1,TXP2 (C,-,kPa)
double
CmineCounterFlow(double R,double mr1,double mr2,TXP TXP1,TXP TXP2)
{
	if(IsTwoPhase(TXP1.X)) {
		if(IsTwoPhase(TXP2.X)) {
			//
			// both fluids two phase, no temperature change
			return 1/R;
		} else {
			//
			// Cratio=Cmin/Cmax=0
			double Cmin = mr2*PropertyTXPtr(SPEC,TXP2);
			if(errorLog.IsError()) {
				errorLog.Add("CmineCounerFlow");
				return -1;
			}
			double Ntu = 1/(R*Cmin);
			double e = 1-exp(-1*Ntu);
			return Cmin*e;
		}
	} else {
		if(IsTwoPhase(TXP2.X)) {

			// Cratio=Cmin/Cmax=0;
			double Cmin = mr1*PropertyTXPtr(SPEC,TXP1);
			if(errorLog.IsError()) {
				errorLog.Add("CmineCounerFlow");
				return -1;
			}
			if(fabs(R*Cmin)<1e-20) {
				errorLog.Add("CmineCounterFlow","fabs(R*cmin)<1e-20");
				return -1;
			}
			double Ntu = 1/(R*Cmin);
			double e = 1-exp(-1*Ntu);
			return Cmin*e;
		} else {

			// counterflow
			double Cr1 = mr1*PropertyTXPtr(SPEC,TXP1);
			if(errorLog.IsError()) {
				errorLog.Add("CmineCounerFlow");
				return -1;
			}
			double Cr2 = mr2*PropertyTXPtr(SPEC,TXP2);
			if(errorLog.IsError()) {
				errorLog.Add("CmineCounerFlow");
				return -1;
			}
			double Cmin = Cr1>Cr2 ? Cr2 : Cr1;
			double Cmax = Cr1>Cr2 ? Cr1 : Cr2;
			double Cratio = Cmin/Cmax;
			double Ntu = 1/(R*Cmin);
			double l = exp(-1*Ntu*(1-Cratio));
			double e = (1-l)/(1-Cratio*l);
			return Cmin*e;
		}
	}
}

// Calculates the product of minimum heat capacity (Cmin) and
// effectiveness (e) for transferring heat in a shell-and-tube
// heat exchanger with one shell pass.  Refrigerant is flowing
// through the tubes and water/gylcol is in the shell.
//
// Refer to Incropera and DeWitt p.661 for details
//
// UA = 1/R (W/K)
// mr,mw (kg/s)
// TXPr (C,-,kPa)
// Tw (C)
double
CmineShellAndTube(double R,double mr,double mw,TXP TXPr,double /*Tw*/)
{
		//
		// heat capacity rate for liquid water
		// Cp(water) = 4.200 +/- 0.022 (kJ/kg*K)
		// between T = 273.15K and 375K (0C-100C)
		// Incropera & DeWitt p.A22
		double Cw = mw*4200.0;

		if(IsTwoPhase(TXPr.X)) {
			// all heat exchangers
			// Cratio=Cmin/Cmax=0 (Cmax=Cr=inf)
			double Ntu = 1/(R*Cw);
			double e = 1-exp(-1*Ntu);
			#ifdef EXTRA_REPORTING
				printf("Cw = %lf, Ntu = %lf, e = %lf\n",Cw,Ntu,e);
			#endif
			return Cw*e;
		} else {
			// shell-and-tube
			// one shell pass
			// 2-4-6 tube passes
			// returns e per tube pass (i.e. multiply
			const double Cr = mr*PropertyTXPtr(SPEC,TXPr);
			if(errorLog.IsError()) {
				errorLog.Add("CmineShellAndTube");
				return -1;
			}
			const double Cmin = Cr<Cw ? Cr : Cw;
			const double Cratio = Cr<Cw ? Cr/Cw : Cw/Cr;
			const double Ntu = 1.0/(R*Cmin);
			const double arg = 1.0+Cratio*Cratio;
			if(arg<0) {
				errorLog.Add("CmineShellAndTube","1.0+Cratio*Cratio<0.0");
				return -1;
			}
			const double l = pow(arg,0.5);
			const double m = exp(-1.0*Ntu*l);
			const double e = 2.0/(1+Cratio+l*(1+m)/(1-m));
			#ifdef EXTRA_REPORTING
				printf("Cr=%lf Cw=%lf\n",Cr,Cw);
			#endif
			return Cmin*e;
		}
}

// Calculates the heat transfer resistance (1/UA) when transferring
// heat between refrigerant and dry air in crossflow.
double
ResistanceCrossFlow(double CminE,double mr,double ma,TXP TXP,double Ta)
{
	if (IsTwoPhase(TXP.X)) {

		// all heat exchangers
		// Cratio=Cmin/Cmax=0 (Cmax=Cr=inf)
		double Cmin = ma*air.Cp(Ta);
		if(errorLog.IsError()) {
			errorLog.Add("ResistanceCrossFlow");
			return -1;
		}
		double e = CminE/Cmin;
		if(e>=1) {
			errorLog.Add("ResistanceCrossFlow","log domain error");
			return -1;
		}
		double Ntu = -1*log(1-e);
		return 1/(Ntu*Cmin);
	} else {
		double Ca = ma*air.Cp(Ta);
		if(errorLog.IsError()) {
			errorLog.Add("ResistanceCrossFlow");
			return -1;
		}

		double Cr = mr*PropertyTXPtr(SPEC,TXP);
		if(errorLog.IsError()) {
			errorLog.Add("ResistanceCrossFlow");
			return -1;
		}

		if (Ca>Cr) {
			//
			// crossflow
			// refirgerant mixed (Cmin), air unmixed (Cmax)
			double Cratio = Cr/Ca;
			double e = CminE/Cr;
			if(e>=1) {
				errorLog.Add("ResistanceCrossFlow","log domain error");
				return -1;
			}
			double l = Cratio*log(1-e)+1;
			if(l<=0) {
				errorLog.Add("ResistanceCrossFlow","log domain error");
				return -1;
			}
			double Ntu = -1*(1/Cratio)*log(l);
			return 1/(Ntu*Cr);
		} else {
			//
			// crossflow
			// refirgerant mixed (Cmax), air unmixed (Cmin)
			double Cratio = Ca/Cr;
			double e = CminE/Ca;
			double l = 1-e*Cratio;
			if(l<=0) {
				errorLog.Add("ResistanceCrossFlow","log domain error");
				return -1;
			}
			double m = 1+log(l)/Cratio;
			if(m<=0) {
				errorLog.Add("ResistanceCrossFlow","log domain error");
				return -1;
			}
			double Ntu = -1*log(m);
			return 1/(Ntu*Ca);
		}
	}
}
