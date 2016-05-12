#ifndef _CONVERSION
#define _CONVERSION

double CtoF(double C);
double FtoC(double F);
double DCtoDF(double DC);
double DFtoDC(double DF);
double PSIGtoKPA(double psig);
double KPAtoPSIG(double kPa);
double DKPAtoDPSIG(double kPa);
// mas flow rate
double KGStoLBMH(double kps);
double LBMHtoKGS(double lbmh);
// power
double WtoBTUH(double watt);
double BTUHtoW(double btuh);

#endif