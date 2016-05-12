#ifndef _PRESSURE
#define _PRESSURE

double dPdLfriction(TXP,double,double);
double dPmom(double,double,double);
double dPelbow(TXP,double,double,double);
double dPdLfrictionEvap(TXP,double,double,double,double);
double LossCoeffReturnBend(TXP,double,double,double);
double LossCoeffStraightPipe(TXP,double,double,double);

struct PreAcc {
		double DP_FR,P_IN, H_OUT, G, X_IN;
			};

double GET_PreAcc(double,void*);// shenbo add new for getting the acceleration pressure drop



#endif
