#ifndef _CORR
#define _CORR

struct ETdim;
struct CGP;
struct TP;

double Circuit_DP(double, TP, TP, ETdim*);//function for calculating the airside pressure drop of the evaporater
double Circuit_DP(double, double, double, CGP*);//function for calculating the airside pressure drop of the condenser

double ConvCoeffAir_EVA(TP,double, ETdim*);
double ConvCoeffAir_CON(double,double,CGP*);
double ConvCoeffInside(TXP,double,double, CGP*);

double FinEffect_Schmidt(double,double,double,double,double);
double ConvCoeffEvapTP_microfin(TXP,double, ETdim*);
double ConvCoeffCondTP_Microfin(TXP,double,CGP*);

double ConvCoeffSP(TXP,double,ETdim*);
double ConvCoeffSP(TXP,double,CGP*);

double FricDP(TXP, double, double, ETdim*);
double FricDP(TXP, double, double, CGP*);
double Coef_Hilpert(double,double,double);
double ConvCoeffSP_Smooth(TXP,double,double);//single-phase heat transfer in smooth tube
double FricDPSP_Smooth(TXP,double,double,double);//single phase pressure drop in smooth tube
double FricDPTP_Smooth(TXP,double,double,double);//two phase pressure drop in smooth tube
//B.S.--------------------------unused
int Exercise_Corr();
double ConvCoeffEvapInside(TXP,double,double,double);//shenbo change ConvCoeffEvapInside(TXP,double,double,double,double);
int IsTwoPhase(double);
double ConvCoeffAir(double,double,double,double,double,double);

#endif
