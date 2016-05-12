#ifndef _SPLINE1D
#define _SPLINE1D

#include "spline.h"
#include "defs.h"

class Spline1d : public Spline {
public:
	Spline1d(const char* fn, int _xcol=0);
	~Spline1d();

	double Evaluate(int prop,double x);

protected:
	double Evaluate_dy(int prop,double x);
	void Evaluate_y_dy(int prop,double x,double* a,double* da);

private:
	int nProp;             // number of properties
	int nData;             // number of rows
	double **y,**y2;
	double xmin,xmax;
	int xcol;

	int is_dx_constant;

	void ScanFile();
	void AllocateMemory();
	void RestoreMemory();
	void LoadData();
};

// x        pressure (kPa)
// prop
// 0        temeprature
// 1        enthalpy
// 2        volume
// 3        entropy
// 4        int energy (J/kg)
// 5        density (kg/m^3)
class RefPropTPthP : public Spline1d {
public:
	enum {
		TSAT,
		ENTH,
		VOL,
		ENTR,
		NTENERGY,
		DENSITY
	};

	RefPropTPthP(const char* fn):Spline1d(fn) {}

	double Tsat(double P) {return OUT_T(Evaluate(0,IN(P)));}
	double h(double P) {return OUT_HUS(Evaluate(1,IN(P)));}
	double v(double P) {return OUT(Evaluate(2,IN(P)));}
	double s(double P) {return OUT_HUS(Evaluate(3,IN(P)));}
	double u(double P) {return OUT_HUS(Evaluate(4,IN(P)));}
	double rho(double P) {return OUT(Evaluate(5,IN(P)));}

	double dTsat_dP(double P) {return dy_dx(0,P);}
	double dh_dP(double P) {return dy_dx(1,P);}
	double dv_dP(double P) {return dy_dx(2,P);}
	double ds_dP(double P) {return dy_dx(3,P);}
	double du_dP(double P) {return dy_dx(4,P);}
	double drho_dP(double P) {return dy_dx(5,P);}

	void h_dhdP(double P,double* h,double* dh_dP);
	void rho_drhodP(double P,double* rho,double* drho_dP);

private:
#if defined(REFPROP_LOG10_1D)
	double IN(double P) {return log10(P)-3;}
	double OUT(double y) {return pow(10,y);}
	double OUT_T(double y) {return pow(10,y)-5e2;}
	double OUT_HUS(double y) {return pow(10,y)-1e6;}
#elif defined(REFPROP_ACMODEL)
	double IN(double P) {return P;}
	double OUT(double y) {return y;}
	double OUT_T(double y) {return y;}
	double OUT_HUS(double y) {return y;}
#else
	double IN(double P) {return 1e-3*P;}
	double OUT(double y) {return y;}
	double OUT_T(double y) {return y;}
	double OUT_HUS(double y) {return y;}
#endif

	double dy_dx(int prop,double x);
};

// x        pressure (kPa)
// prop
// 0        temeprature
// 1           conductivity
// 2           viscosity
// 3          specific heat
// 4		 liquid surface tension
class RefPropTPtrP : public Spline1d {
public:
	enum {
		COND,
		VISC,
		SPEC,
		TENSION
	};
	RefPropTPtrP(const char* fn):Spline1d(fn) {}

	double k(double P) {return OUT(Evaluate(0,IN(P)));}
	double mu(double P) {return OUT(Evaluate(1,IN(P)));}
	double Cp(double P) {return OUT(Evaluate(2,IN(P)));}
	double Tension(double P) {return OUT(Evaluate(3,IN(P)));}

	double dk_dP(double P) {return dy_dx(0,P);}
	double dmu_dP(double P) {return dy_dx(1,P);}
	double dCp_dP(double P) {return dy_dx(2,P);}
	double dTension_dP(double P) {return dy_dx(3,P);}

private:
#if defined(REFPROP_LOG10_1D)
	double IN(double P) {return log10(P)-3;}
	double OUT(double y) {return pow(10,y);}
	double OUT_T(double y) {return pow(10,y)-5e2;}
	double OUT_HUS(double y) {return pow(10,y)-1e6;}
#elif defined(REFPROP_ACMODEL)
	double IN(double P) {return P;}
	double OUT(double y) {return y;}
	double OUT_T(double y) {return y;}
	double OUT_HUS(double y) {return y;}
#else
	double IN(double P) {return 1e-3*P;}
	double OUT(double y) {return y;}
	double OUT_T(double y) {return y;}
	double OUT_HUS(double y) {return y;}
#endif

	double dy_dx(int prop,double x);
};



/********************************/

// x        temperature (C)
// prop
// 0        pressure (Pa)
// 1        enthalpy
// 2        volume
// 3        entropy
// 4        int energy
// 5        density
class RefPropTPthT : public Spline1d {
public:
	RefPropTPthT(const char* fn):Spline1d(fn,1) {}

	double Psat(double T) {return OUT(Evaluate(0,IN(T)));}
	double h(double T) {return OUT_H(Evaluate(1,IN(T)));}

private:
#ifdef REFPROP_LOG10_1D
	double IN(double T) {return log10(T+5e2);}
	double OUT(double y) {return 1e3*pow(10,y);}
	double OUT_H(double y) {return pow(10,y)-1e6;}
#elif defined(REFPROP_ACMODEL)
	double IN(double T) {return T;}
	double OUT(double y) {return y;}
	double OUT_H(double y) {return y;}
#else
	double IN(double T) {return T;}
	double OUT(double y) {return 1e3*y;}
	double OUT_H(double y) {return y;}
#endif
};

// x       temperature
// prop
// 0       conductivity
// 1       viscosity
// 2       specific heat
// 3       ???
// 4       enthalpy
// 5       specific volume
// 6       specific entropy
class AirProp : public Spline1d {
public:
	AirProp(const char* fn):Spline1d(fn) {}
	double k(double T) {return Evaluate(0,T);}
	double mu(double T) {return Evaluate(1,T);}
	double Cp(double T) {return Evaluate(2,T);}
	double h(double T) {return Evaluate(4,T);}
	double v(double T) {return Evaluate(5,T);}
	double s(double T) {return Evaluate(6,T);}
};
#endif