#ifndef _WAIR
#define _WAIR

struct TP {
	double T,P;
} ;
struct WH {
	double W,H;
} ;

TP WHtoTP(WH,TP);
TP HPtoTP(double,double);
TP WPtoTP(double,double);
TP TWtoTP(double,double);
TP TDtoTP(double,double);
TP TZtoTP(double,double);
TP THtoTP(double,double);

#define TWAIRMAX 50
#define TWAIRMIN -20

#endif
