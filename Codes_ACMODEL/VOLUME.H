#ifndef _VOLUME
#define _VOLUME

double Xtt(TXP);
double VolumeALL(TXP,double,double,double);
double VolumeTP(TXP,double,double,double);
/* Local pressure drop calculations. */
double VolumeTPFunc(double,void*);

/* Structure containing parameters used in the function iterated
	on to get the void fraction.                                     */
struct VolParams {
	double vv,vl,ul,uv,Beta,G,D,x;
};


struct FlowPattern//output struct variable from Kattan-Thome flow-pattern-dependent heat transfer model
{
	int JudgPattern;
	double G_wavy;
	double G_strat;
	double G_mist;
	double G_bub;
	double X_lA;
	double epsilon;
	double theta_dry;
	double delta;
	double h_nb,h_cb,h_v,h_wet,h_tp;
	int Pattern;
};

double VolumeTP_Thome(TXP,double,double,double);
double VolumeTP_FlowPattern(TXP,double,double,double);
double Eva_FlowPattern(TXP,double,double,double,FlowPattern*);//Kattan-Thome flow-pattern-specific flow boiling heat transfer model
double Cond_FlowPattern(TXP,double,double,double,FlowPattern*);//Thome-Cavallini flow-pattern-specific flow condensation heat transfer model

#endif
