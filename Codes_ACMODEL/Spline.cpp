#include "stdafx.h"

#include <string.h>
#include "spline.h"
#include "errorlog.h"

extern ErrorLog errorLog;

Spline::Spline(const char *fn)
{
	strncpy(filename,fn,120);

	extrap_frac = 0;    // allow 1% extrapolation

	setup_error = 1;
}

/********************************************************************
Calculates second derivatives to be stored in memory and then used
for cubic spline interpolation.  nDATA3 is used to allocate memory
for local buffer.  It much be greater than or egual to the largest
n that may occur.
INPUTS
x - pointer to array of "x" values in table.
y - pointer to array of "y" values in table.
n - number of data points in table.
yp1 - first derivative of first point. Usually infinity.
ypn - first derivative of last point. Usually infinity.
OUTPUTS
y2 - pointer to array of second derivatives.
********************************************************************/
void Spline::spline(double *x,double *y,int n,double yp1,double ypn,double *y2)
{
	int i,k;
	double p,qn,sig,un,*u;

	u = new double[n];

	if (yp1 > 0.99e30)
		y2[0]=u[0]=0.0;
	else {
		y2[0] = -0.5;
		u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
	}
	for (i=1;i<n-1;i++) {
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	if (ypn > 0.99e30)
		qn=un=0.0;
	else {
		qn=0.5;
		un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
	}
	y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
	for (k=n-2;k>=0;k--) y2[k]=y2[k]*y2[k+1]+u[k];

	delete u;
}

/********************************************************************
Evaluates a cubic spline interpolation.
********************************************************************/
void
Spline::splint(double *xa,double *ya,double *y2a,int n,double x,double *y,double *dy,int is_constant)
{
	/***
	double xf = extrap_frac*(xa[n-1]-xa[0]);
	if(x>xa[n-1]+xf || x<xa[0]-xf) {
		char msg[128];
		sprintf(msg,"Data out-of-range (x=%le,xmax=%le,xmin=%le)",x,xa[n-1]+xf,xa[0]-xf);
		errorLog.Add("Spline::splint",msg);
		return;
	}
	if(x>xa[n-1] || x<xa[0]) {
		char msg[128];
		sprintf(msg,"Extrapolating within allowed %.2lf%% of table full width. (x=%le,xmax=%le,xmin=%le)",100.0*extrap_frac,x,xa[n-1],xa[0]);
		errorLog.AddWarning("Spline","splint",msg);
	}
	***/

	// Find the array positions bordering x

	int klo,khi;
	if(is_constant) {
		// Faster routine if independent variable x is evenly spaced.
		double k = (x-xa[0])/(xa[1]-xa[0]);
		klo = (int)k;
		khi = klo + 1;
	} else {
		// Search routine
		klo = 0; 
		khi = n-1;
		while (khi-klo > 1) {
			int k = (khi+klo) >> 1;
			if(xa[k] > x) khi=k; else klo=k;
		}
	}

	double h = xa[khi] - xa[klo];
	double a = (xa[khi]-x)/h;
	double aa = a*a;
	double b = (x-xa[klo])/h;
	double bb = b*b;

	*y = a*ya[klo] + b*ya[khi] + ( a*(aa-1.0)*y2a[klo] + b*(bb-1.0)*y2a[khi] )*(h*h)/6.0;
	//Calculate dy/dx
	*dy = (ya[khi]-ya[klo])/h - ( (3.0*aa-1.0)*y2a[klo] + (3.0*bb-1.0)*y2a[khi] )*h/6.0;
}

int
Spline::IsPartOfDouble(char c)
{
	if(c>='0' && c<='9') return 1;
	if(c=='-' || c=='+') return 1;
	if(c=='.') return 1;
	if(c=='E' || c=='e') return 1;
	return 0;
}

// Returns:
//    0 = all ok
//    1 = double not found
int
Spline::ScanStringForDouble(char** str, double& value)
{
	// Advance until a valid character is found or end-of-string
	// is encountered.
	while(**str!=0 && !IsPartOfDouble(**str)) { (*str)++; }

	// Scan the string for the desired value
	int found = sscanf(*str,"%lf",&value);

	// Advance until a invalid character is found or end-of-string
	// is encountered.
	while(**str!=0 && IsPartOfDouble(**str)) { (*str)++; }

	return found==1 ? 0 : 1;
}

