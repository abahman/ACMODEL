#ifndef _SPLINE2D
#define _SPLINE2D

#include "spline.h"
#include "defs.h"

class Spline2d : public Spline {
public:
	Spline2d(const char* fn,int xcol=0,int ycol=1);
	~Spline2d();
	double Evaluate_y(int prop,double x,double y);
protected:
	double Evaluate(int prop,double x,double y, int return_value);
	double Evaluate_dydx1(int prop,double x,double y);
	double Evaluate_dydx2(int prop,double x,double y);
	void Evaluate_y_dy(int prop,double x1,double x2,int return_value,double* a,double* dadx1,double* dadx2);
private:
	// table dimensions
	int nProp;             // number of properties
	int nData;             // number of rows
	int nX;                // number of x values
	int nY;                // maximum number of y values

	int xcol;
	int ycol;

	double dx2;            // minimum distance between x1 values

	int is_dx1_constant;

	// data arrays
	double **TBL,***y2,***y;
	double *x1new,**x;
	int *x2n;
	int x1n;
	double x2min,x2max;

	void ScanFile();
	void AllocateMemory();
	void RestoreMemory();
	void LoadData();
};

// Refrigerant properties
// Single phase (SP) - implies 2D spline
// thermodynamic properties (th) - not transport properties
// pressure and enthalpy (PH) are independent variables

// Assumes following table format with log10 properties:
// column      property
// 0          *log10(pressure(Pa))
// 1           log10(termperature(C)+5e2)
// 2          *log10(enthalpy(J/kg)+1e6)
// 3           log10(volume(m^3/kg)
// 4           log10(entropy(J/kg/K+1e6)
// 5           log10(IntEnergy(J/kg)+1e6)
// 6           log10(density(kg/m^3))

// Assumes following table format:
// column      property
// 0          *pressure
// 1           termperature
// 2          *enthalpy
// 3           volume
// 4           entropy
// 5           int energy (J/kg)
// 6           density (kg/m^3)
class RefPropSPthPH : public Spline2d {
public:
	RefPropSPthPH(const char* fn) : Spline2d(fn,0,2) {}

	double T(double P,double h);
	double v(double P,double h) {return OUT(Evaluate_y(1,IN1(P),IN2(h)));}
	double s(double P,double h) {return OUT_HUS(Evaluate_y(2,IN1(P),IN2(h)));}
	double u(double P,double h) {return OUT_HUS(Evaluate_y(3,IN1(P),IN2(h)));}
	double rho(double P,double h) {return OUT(Evaluate_y(4,IN1(P),IN2(h)));}

	double dT_dP(double P,double h) {return dy_dx1(0,P,h);}
	double dT_dh(double P,double h) {return dy_dx2(0,P,h);}
	double dv_dP(double P,double h) {return dy_dx1(1,P,h);}
	double dv_dh(double P,double h) {return dy_dx2(1,P,h);}
	double ds_dP(double P,double h) {return dy_dx1(2,P,h);}
	double ds_dh(double P,double h) {return dy_dx2(2,P,h);}
	double drho_dP(double P,double h) {return dy_dx1(4,P,h);}
	double drho_dh(double P,double h) {return dy_dx2(4,P,h);}

	void rhoPlusDerivatives(double P,double h,double* rho,double* drho_dP,double* drho_dh);
private:
#if defined(REFPROP_LOG10_2D)
	double IN1(double P) {return log10(P)-3;}
	double IN2(double h) {return log10(h+1e6);}
	double OUT(double y) {return pow(10,y);}
	double OUT_T(double y) {return pow(10,y)-5e2;}
	double OUT_HUS(double y) {return pow(10,y)-1e6;}
#elif defined(REFPROP_ACMODEL)
	double IN1(double P) {return P;}
	double IN2(double h) {return h;}
	double OUT(double y) {return y;}
	double OUT_T(double y) {return y;}
	double OUT_HUS(double y) {return y;}
#else
	double IN1(double P) {return 1e-3*P;}
	double IN2(double h) {return h;}
	double OUT(double y) {return y;}
	double OUT_T(double y) {return y;}
	double OUT_HUS(double y) {return y;}
#endif
	double dy_dx1(int prop,double P,double h);
	double dy_dx2(int prop,double P,double h);
};

// Refrigerant properties
// Single phase (SP) - implies 2D spline
// transport properties (tr) - not thermodynamic properties
// pressure and temperature (PT) are independent variables

// Assumes following table format:
// column      property
// 0           pressure
// 1           temperature
// 2           conductivity
// 3           viscosity
// 4           specific heat
class RefPropSPtrPT : public Spline2d {
public:
	RefPropSPtrPT(const char* fn) : Spline2d(fn,0,1) {
		extrap_frac = 1.0;
	}
	enum {
		COND,
		VISC,
		SPEC
	};
#ifdef REFPROP_ACMODEL
	double k(double P,double T) {return Evaluate_y(COND,P,T);}
	double mu(double P,double T) {return Evaluate_y(VISC,P,T);}
	double Cp(double P,double T) {return Evaluate_y(SPEC,P,T);}
#else
	double k(double P,double T) {return Evaluate_y(COND,1e-3*P,T);}
	double mu(double P,double T) {return Evaluate_y(VISC,1e-3*P,T);}
	double Cp(double P,double T) {return Evaluate_y(SPEC,1e-3*P,T);}
#endif
};

// Refrigerant properties
// Single phase (SP) - implies 2D spline
// thermodynamic properties (th) - not transport properties
// pressure and temperature (PT) are independent variables

// Assumes following table format:
// column      property
// 0          *pressure
// 1          *termperature
// 2           enthalpy
// 3           volume
// 4           entropy
// 5           int energy
// 6           density (kg/m^3)
class RefPropSPthPT : public Spline2d {
public:
	RefPropSPthPT(const char* fn) : Spline2d(fn,0,1) {}
	enum {
		ENTH,
		VOL,
		ENTR,
		INTENERGY,
		DENSITY
	};
	double h(double P,double T) {
		return OUT_HUS(Evaluate_y(ENTH,IN1(P),IN2(T)));
	}
	double v(double P,double T) {return OUT(Evaluate_y(VOL,IN1(P),IN2(T)));}
	double s(double P,double T) {return OUT_HUS(Evaluate_y(ENTR,IN1(P),IN2(T)));}
	double u(double P,double T) {return OUT_HUS(Evaluate_y(INTENERGY,IN1(P),IN2(T)));}
	double rho(double P,double T) {return OUT_HUS(Evaluate_y(DENSITY,IN1(P),IN2(T)));}
private:
#ifdef REFPROP_LOG10_2D
	double IN1(double P) {return log10(P)-3;}
	double IN2(double T) {return log10(T+5e2);}
	double OUT(double y) {return pow(10,y);}
	double OUT_HUS(double y) {return pow(10,y)-1e6;}
#elif defined(REFPROP_ACMODEL)
	double IN1(double P) {return P;}
	double IN2(double T) {return T;}
	double OUT(double y) {return y;}
	double OUT_HUS(double y) {return y;}
#else
	double IN1(double P) {return 1e-3*P;}
	double IN2(double T) {return T;}
	double OUT(double y) {return y;}
	double OUT_HUS(double y) {return y;}
#endif
};

// Refrigerant properties
// two phase (TP) for refrigerant mixtures (R410A, R407C) - implies 2D spline
// thermodynamic properties (th) - not transport properties
// pressure and quality (PX) are independent variables

// Assumes following table format:
// column      property
// 0          *pressure
//1		*quality		
// 2          termperature
// 3           enthalpy
// 4           volume
// 5           entropy
// 6           int energy
// 7           density (kg/m^3)
class RefPropTPthPX : public Spline2d {//B.S.
public:
	RefPropTPthPX(const char* fn) : Spline2d(fn,0,1) {}
	enum {
		TEMP,
		ENTH,
		VOL,
		ENTR,
		INTENERGY,
		DENSITY
	};
	double T(double P,double X) {return OUT(Evaluate_y(TEMP,IN1(P),IN2(X)));}
	double h(double P,double X) {return OUT_HUS(Evaluate_y(ENTH,IN1(P),IN2(X)));}
	double v(double P,double X) {return OUT(Evaluate_y(VOL,IN1(P),IN2(X)));}
	double s(double P,double X) {return OUT_HUS(Evaluate_y(ENTR,IN1(P),IN2(X)));}
	double u(double P,double X) {return OUT_HUS(Evaluate_y(INTENERGY,IN1(P),IN2(X)));}
	double rho(double P,double X) {return OUT_HUS(Evaluate_y(DENSITY,IN1(P),IN2(X)));}
private:
#ifdef REFPROP_LOG10_2D
	double IN1(double P) {return log10(P)-3;}
	double IN2(double X) {return log10(X+5e2);}
	double OUT(double y) {return pow(10,y);}
	double OUT_HUS(double y) {return pow(10,y)-1e6;}
#elif defined(REFPROP_ACMODEL)
	double IN1(double P) {return P;}
	double IN2(double X) {return X;}
	double OUT(double y) {return y;}
	double OUT_HUS(double y) {return y;}
#else
	double IN1(double P) {return 1e-3*P;}
	double IN2(double X) {return X;}
	double OUT(double y) {return y;}
	double OUT_HUS(double y) {return y;}
#endif
};

// Refrigerant properties
// two phase (TP) for refrigerant mixtures (R410A, R407C) - implies 2D spline
// transport properties (tr) - not thermodynamic properties
// pressure and quality (PX) are independent variables

// Assumes following table format:
// column      property
// 0			*pressure
//	1			*quality		
// 2			Surface Tension

class RefPropTPtrPX : public Spline2d {//B.S.
public:
	RefPropTPtrPX(const char* fn) : Spline2d(fn,0,1) {}
	enum {
		TENSION,//the following enum are not used now.
		COND,
		VISC,
		SPEC
	};
	double Sigma(double P,double X) {return OUT(Evaluate_y(TENSION,IN1(P),IN2(X)));}
	
private:
#ifdef REFPROP_LOG10_2D
	double IN1(double P) {return log10(P)-3;}
	double IN2(double X) {return log10(X+5e2);}
	double OUT(double y) {return pow(10,y);}
	double OUT_HUS(double y) {return pow(10,y)-1e6;}
#elif defined(REFPROP_ACMODEL)
	double IN1(double P) {return P;}
	double IN2(double X) {return X;}
	double OUT(double y) {return y;}
	double OUT_HUS(double y) {return y;}
#else
	double IN1(double P) {return 1e-3*P;}
	double IN2(double X) {return X;}
	double OUT(double y) {return y;}
	double OUT_HUS(double y) {return y;}
#endif
};

// Refrigerant properties
// two phase (TP) for refrigerant mixtures (R410A, R407C) - implies 2D spline
// thermodynamic properties (th) - not transport properties
// pressure and enthalpy (PH) are independent variables

// Assumes following table format:
// column      property
// 0          *pressure
//	1			Qulity	
// 2           termperature
// 3          *enthalpy
// 4          volume
// 5           entropy
// 6           int energy
// 7           density (kg/m^3)
class RefPropTPthPH : public Spline2d {//B.S.
public:
	RefPropTPthPH(const char* fn) : Spline2d(fn,0,3) {}
	double X(double P,double h) {return OUT_T(Evaluate_y(0,IN1(P),IN2(h)));}
	double T(double P,double h) {return OUT_T(Evaluate_y(1,IN1(P),IN2(h)));}
	double v(double P,double h) {return OUT_HUS(Evaluate_y(2,IN1(P),IN2(h)));}
	double s(double P,double h) {return OUT_HUS(Evaluate_y(3,IN1(P),IN2(h)));}
	double u(double P,double h) {return OUT_HUS(Evaluate_y(4,IN1(P),IN2(h)));}
private:
#ifdef REFPROP_LOG10_2D
	double IN1(double P) {return log10(P)-3;}
	double IN2(double h) {return log10(h);}
	double OUT_T(double y) {return pow(10,y)-5e2;}
	double OUT_HUS(double y) {return pow(10,y)-1e6;}
#elif defined(REFPROP_ACMODEL)
	double IN1(double P) {return P;}
	double IN2(double h) {return h;}
	double OUT_T(double y) {return y;}
	double OUT_HUS(double y) {return y;}
#else
	double IN1(double P) {return 1e-3*P;}
	double IN2(double h) {return h;}
	double OUT_T(double y) {return y;}
	double OUT_HUS(double y) {return y;}
#endif
};


// Refrigerant properties
// Single phase (SP) - implies 2D spline
// thermodynamic properties (th) - not transport properties
// pressure and volume (PV) are independent variables

// Assumes following table format:
// column      property
// 0          *pressure
// 1           termperature
// 2           enthalpy
// 3          *volume
// 4           entropy
// 5           int energy
// 6           density (kg/m^3)
class RefPropSPthPV : public Spline2d {
public:
	RefPropSPthPV(const char* fn) : Spline2d(fn,0,3) {}
	double T(double P,double v) {return OUT_T(Evaluate_y(0,IN1(P),IN2(v)));}
	double h(double P,double v) {return OUT_HUS(Evaluate_y(1,IN1(P),IN2(v)));}
	double s(double P,double v) {return OUT_HUS(Evaluate_y(2,IN1(P),IN2(v)));}
	double u(double P,double v) {return OUT_HUS(Evaluate_y(3,IN1(P),IN2(v)));}
private:
#ifdef REFPROP_LOG10_2D
	double IN1(double P) {return log10(P)-3;}
	double IN2(double v) {return log10(v);}
	double OUT_T(double y) {return pow(10,y)-5e2;}
	double OUT_HUS(double y) {return pow(10,y)-1e6;}
#elif defined(REFPROP_ACMODEL)
	double IN1(double P) {return P;}
	double IN2(double v) {return v;}
	double OUT_T(double y) {return y;}
	double OUT_HUS(double y) {return y;}
#else
	double IN1(double P) {return 1e-3*P;}
	double IN2(double v) {return v;}
	double OUT_T(double y) {return y;}
	double OUT_HUS(double y) {return y;}
#endif
};

// Refrigerant properties
// Single phase (SP) - implies 2D spline
// thermodynamic properties (th) - not transport properties
// pressure and internal energy (PU) are independent variables

// Assumes following table format:
// column      property
// 0          *pressure
// 1           temperature
// 2           enthalpy
// 3           volume
// 4           entropy
// 5          *int energy
// 6           density (kg/m^3)
class RefPropSPthPU : public Spline2d {
public:
	RefPropSPthPU(const char* fn) : Spline2d(fn,0,5) {}
	double T(double P,double u) {return OUT_T(Evaluate_y(0,IN1(P),IN2(u)));}
	double h(double P,double u) {return OUT_HUS(Evaluate_y(1,IN1(P),IN2(u)));}
	double v(double P,double u) {return OUT(Evaluate_y(2,IN1(P),IN2(u)));}
	double s(double P,double u) {return OUT_HUS(Evaluate_y(3,IN1(P),IN2(u)));}
private:
#ifdef REFPROP_LOG10_2D
	double IN1(double P) {return log10(P)-3;}
	double IN2(double u) {return log10(u+1e6);}
	double OUT(double y) {return pow(10,y);}
	double OUT_T(double y) {return pow(10,y)-5e2;}
	double OUT_HUS(double y) {return pow(10,y)-1e6;}
#elif defined(REFPROP_ACMODEL)
	double IN1(double P) {return P;}
	double IN2(double u) {return u;}
	double OUT(double y) {return y;}
	double OUT_T(double y) {return y;}
	double OUT_HUS(double y) {return y;}
#else
	double IN1(double P) {return 1e-3*P;}
	double IN2(double u) {return u;}
	double OUT(double y) {return y;}
	double OUT_T(double y) {return y;}
	double OUT_HUS(double y) {return y;}
#endif
};

// Wet air properties
// temperature and relative humidity (TR) are independent variables
// file format
// drybulb temperature
// relative humidity
// wetbulb
// dew point
// specific heat
// humidity ratio
// specific enthalpy
class WetAirPropTR : public Spline2d {
public:
	WetAirPropTR(const char *fn) : Spline2d(fn,0,1) {}
	enum {
		WETBULB,
		DEWPNT,
		SPEC,
		HUMRAT,
		ENTH
	};
	double WetBulb(double T,double R) {return Evaluate_y(WETBULB,T,R);}
	double DewPoint(double T,double R) {return Evaluate_y(DEWPNT,T,R);}
	double Cp(double T,double R) {return Evaluate_y(SPEC,T,R);}
	double HumidityRatio(double T,double R) {return Evaluate_y(HUMRAT,T,R);}
	double h(double T,double R) {return Evaluate_y(ENTH,T,R);}
};

// Lubricant properties
//  implies 2D spline
// temperature and pressure (TP) are independent variables

// Assumes following table format:
// column      property
// 0           temperature
// 1           pressure

class LubriProp : public Spline2d {
public:
	LubriProp(const char* fn) : Spline2d(fn,0,1) {
		extrap_frac = 1.0;
	}

	enum {
		Xi
	};

#ifdef REFPROP_ACMODEL
	double Solubility(double T,double P) {return Evaluate_y(Xi,T,P);}
#endif

};

#endif
