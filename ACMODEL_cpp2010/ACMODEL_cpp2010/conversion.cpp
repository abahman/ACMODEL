#include "stdafx.h"
#include "conversion.h"

double CtoF(double C)
{
	return 1.8*C+32.0;
}

double FtoC(double F)
{
	return (F-32.0)/1.8;
}

double DCtoDF(double DC)
{
	return 1.8*DC;
}

double DFtoDC(double DF)
{
	return DF/1.8;
}

double PSIGtoKPA(double psig)
{
	return (psig+14.7)/14.7*101.3;
}

double KPAtoPSIG(double kPa)
{
	return kPa/101.3*14.7-14.7;
}

double DKPAtoDPSIG(double kPa)
{
	return kPa/101.3*14.7;
}

// mass flow rate conversion
// KGS = kg/s
// LBMH = lbm/hr
double KGStoLBMH(double kps)
{
	return kps*2.2*3600;
}
double LBMHtoKGS(double lbmh)
{
	return lbmh/(2.2*3600);
}

// power conversion
// W = watt
// BTUH = BTU/hr
double WtoBTUH(double watt)
{
	return watt*3.413;
}
double BTUHtoW(double btuh)
{
	return btuh/3.413;
}
