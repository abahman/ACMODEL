#include "stdafx.h"

#include "dataST.h"

SteadyData::SteadyData()
{

	SetSteadyDefault();
}

void SteadyData::SetSteadyDefault(){
	//struct Steady* steady=new Steady;
	mr2 =0.228958;     //        kg/s
	T2 = 89.45884;      //        C, compressor outlet
	P2 = 1530.491;      //        kPa
	T3 = 32.29615;      //        C, condenser outlet
	P3 = 1473.639;      //        kPa
	Gca = -1;           //        kg/s/m^2
	Tcai = 11.46957;    //        C
	Tcao = 37.13032;    //        C
	mr1 = 0.190459;     //        kg/s
	P7 = 141.8417;      //        kPa, evaporator outlet
	T7 = -31.2247;      //        C
	P5 = 1416.788;      //        kPa, subcooler outlet
	T5 = 12.03681;      //        C
	mw = 2.0;  //3.131245;      //        kg/s
	Tewi = -17.3451;    //        C
	Tewo = -21.4402;    //        C
	P10_11_12=466.3234; //        kPa
	T10 = 21.66723;     //        C
	T11 = 11.12536;     //        C
	E = 20.84763;       //        kW, compressor power

}

