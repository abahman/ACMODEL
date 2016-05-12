#ifndef _DATAST
#define _DATAST
#include <stdio.h>

class SteadyData {
public:
	 SteadyData();
	 double mr2;
	 double T2;
	 double P2;
	 double T3;
	 double P3;
	 double Gca ;
	 double Tcai ;
	 double Tcao;
	 double mr1;
	 double P7 ;
	 double T7 ;
	 double P5 ;
	 double T5 ;
	 double mw ;
	 double Tewi ;
	 double Tewo ;
	 double P10_11_12;
	 double T10;
	 double T11 ;
	 double E ;
	 void SetSteadyDefault();
};
#endif