#ifndef _SPLINE
#define _SPLINE

class Spline {
public:
	Spline(const char* fn);
protected:
	char filename[128];
	int setup_error;
	double extrap_frac;    // fraction of full x,y range that table
                           // is allowed to extrapolate

	void spline(double *x,double *y,int n,double yp1,double ypn,double *y2);
	void splint(double *xa,double *ya,double *y2a,int n,double x,double *y,double *dy,int is_constant=0);

	int IsPartOfDouble(char c);
	int ScanStringForDouble(char** str, double& value);
};

#endif