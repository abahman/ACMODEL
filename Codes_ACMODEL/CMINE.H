#ifndef _CMINE
#define _CMINE

// refrigerant to dry air in crossflow
double CmineCrossFlow_dry(double,double,double,TXP,double);//B.S.
double CmineCrossFlow_wet(double,double,double,TXP,double,double);//B.S.

double ResistanceCrossFlow(double,double,double,TXP,double);

// refrigerant to refrigerant in counterflow
double CmineCounterFlow(double R,double mr1,double mr2,TXP TXP1,TXP TXP2);

// refrigerant to water/gylcol in shell-and-tube
double CmineShellAndTube(double R,double mr,double mw,TXP TXPr,double Tw);

#endif
