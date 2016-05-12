#include "stdafx.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "r22.h"
#include "numeric.h"
#include "errorlog.h"

extern ErrorLog errorLog;

#define TOL 1e-4

/********************************************************************
Finds the zero of a function.  Uses Newton-Raphson method with
forward difference derivative estimates.
********************************************************************/
double FindZero(double Xg,double(*F)(double,void*),double tol,void *P)
{
	double Xo,X,Yp,Y,m,b,dx;
	const double e=1e-3;
	const int Imax=100,Jmax=20;
	int i,j;

	X = Xg;
	Y = (*F)(X,P);
	if(errorLog.IsError()) {
		errorLog.Add("FindZero");
		return X;
	}
	for (i=1; i<=Imax; i++ ) {
		dx=e*fabs(X);
		if(dx==0) {
			errorLog.Add("FindZero","dx=0");
			return X;
		}
		Yp=(*F)(X+dx,P);
		if(errorLog.IsError()) {
			errorLog.Add("FindZero");
			return X;
		}
		m = (Yp-Y)/dx;
		while (m==0) {
			if(fabs(dx/X)>0.5) {
				errorLog.Add("FindZero","m=0");
				return X ;
			}
			dx*=10;
			Yp=(*F)(X+dx,P);
			if(errorLog.IsError()) {
				errorLog.Add("FindZero");
				return X;
			}
			m = (Yp-Y)/dx;
		}
		Xo=X;
		b=Y-m*X;
		X=-b/m ;
		Y=(*F)(X,P);
		j=0;
		while(errorLog.IsError()) {
			j++;
			if(j>Jmax) {
				char str[128];
				sprintf(str,"(%d) Max backoff iteration",j);
				errorLog.AddWarning("FindZero",str);
				break;
			}
			char str[128];
			sprintf(str,"(%d) Backing off. Restarting.\n",j);
			errorLog.ClearError("FindZero",str);
			X=(X+Xo)/2;
			Y=(*F)(X,P);
		}
		if(errorLog.IsError()) {
			errorLog.Add("FindZero","Can not resolve error by backing off. Aborting.");
			return X;
		}
		if (fabs(Y)<=tol) return X;
	}

	{
	char str[128];
	sprintf(str,"FindZero: Max iteration error. Y=%le\n",Y);
	errorLog.Add("FindZero",str);
	}

	return X ;
}

/********************************************************************
Finds the zero of a function.  Uses Newton-Raphson method with
forward difference derivative estimates.  Adds constraint that
X must always be between Xl and Xh.
********************************************************************/
int FindZeroConst(double *Xg,double Xl,double Xh,double(*F)(double,void*),double tol,void *P)
{
	static int i,j;
	static double Xo,X,Yp,Y,m,b,dx;
	const double e=1e-3;
	const int Imax=100,Jmax=20;

	X = *Xg;
	Y = (*F)(X,P);
	if(errorLog.IsError()) {
		errorLog.Add("FindZeroConst");
		return 0;
	}
	for (i=1; i<=Imax; i++ ) {
		dx=e*fabs(X);
		if(dx==0) {
			errorLog.Add("FindZeroConst","dx=0");
			return 0;
		}
		Yp=(*F)(X+dx,P);
		if(errorLog.IsError()) {
			errorLog.Add("FindZeroConst");
			return 0;
		}
		m=(Yp-Y)/dx;
		while (m==0) {
			if(fabs(dx/X)>0.5) {
				errorLog.Add("FindZeroConst","m=0");
				return 0;
			}
			dx*=10;
			Yp=(*F)(X+dx,P);
			if(errorLog.IsError()) {
				errorLog.Add("FindZeroConst");
				return 0;
			}
			m=(Yp-Y)/dx;
		}
		Xo=X;
		b=Y-m*X ;
		X=-b/m ;
		j=0; while(X<Xl||X>Xh) {
			j++;
			if(j>Jmax) break;
			X=(X+Xo)/2;
		}
		Y=(*F)(X,P);
		j=0;
		while(errorLog.IsError()) {
			j++;
			if(j>Jmax) {
				char str[128];
				sprintf(str,"(%d) Max backoff iteration",j);
				errorLog.Add("FindZeroConst",str);
				break;
			}
			char str[128];
			sprintf(str,"(%d) Backing off. Restarting.\n",j);
			errorLog.ClearError("FindZeroConst",str);
			X=(X+Xo)/2;
			Y=(*F)(X,P);
		}
		if(errorLog.IsError()) {
			errorLog.Add("FindZeroConst");
			return 0;
		}
		if(fabs(Y)<=tol) {
			*Xg=X;
			return 0;
		}
	}

	*Xg = X;
	return 0;
}

/********************************************************************
Uses the golden section method to minimize a one dimensional
function.
INPUTS:
	X1 = lower bound of free variable
	X2 = upper bound of free variable
	F = function to be minimized
	tol = convergence tolerance
	P = parameters passed on the F
	status = returns non zero to indicate error
OUTPUT:
	return = minimum value of X
********************************************************************/
double GoldSecMin(double X1,double X2,double(*F)(double,void*),double tol,void *P)
{
	double dL,L,x0,x1,x2,f1,f2;
	const double d=1-2/(1+pow(5,0.5));
	const int n=100;
	int i;

	L=X2-X1;
	x1=X1+d*L;
	f1=F(x1,P);
	x2=X2-d*L;
	f2=F(x2,P);

	i=2;
	dL=2*L/(X1+X2);
	while (i<=n && fabs(dL)>=tol) {
		if (f1>f2) {
			X1=x1;
			x1=x2;
			f1=f2;
			L=X2-X1;
			x2=X2-d*L;
			f2=F(x2,P);
			i++;
		} else if (f1<f2) {
			X2=x2;
			x2=x1;
			f2=f1;
			L=X2-X1;
			x1=X1+d*L;
			f1=F(x1,P);
			i++;
		} else {
			X1=x1;
			X2=x2;
			L=X2-X1;
			x1=X1+d*L;
			f1=F(x1,P);
			x2=X2-d*L;
			f2=F(x2,P);
			i+=2;
		}
		dL=2*L/(X1+X2);
	}
	x0=X1+L/2;

	if(i>=n) {
		errorLog.Add("GoldSecMin","Max iteration error");
		return x0;
	}

	return x0;
}

/********************************************************************
Finds the zero of a function using bisection method. Refer to
Rao's book for details.
********************************************************************/
double Bisection1D(double Xlo,double Xhi,double(*F)(double,void*),double thres,void *P)
{
	double Xbuf,X,Y;
	double Fbuf;
	const int Imax=40;

	double Flo = F(Xlo,P);
	double Fhi = F(Xhi,P);
	if(Flo*Fhi>0) {
		char str[128];
		sprintf(str,"Root not bracketed. F(%lf)=%lf F(%lf)=%lf",Xlo,Flo,Xhi,Fhi);
		errorLog.Add("Bisection1D",str);
		return -999.0;
	}
	if(Flo>Fhi) {
		Fbuf=Flo; Xbuf=Xlo;
		Flo=Fhi; Xlo=Xhi;
		Fhi=Fbuf; Xhi=Xbuf;
	}

	int i=0;
	do {
		X=(Xhi+Xlo)/2;
		i++ ;
		if(i>Imax) {
			errorLog.Add("Bisection1D","Max iteration error");
			return X;
		}
		Y = F(X,P);
		if(Y>0) Xhi=X; else Xlo=X;
	} while(fabs(Y)>thres);

	return X;
}

/********************************************************************
Finds the zero of a function using the false position method.
Refer to numerical recipes for details.
********************************************************************/
double FalsePos(double x1,double x2,double(*F)(double,void*),double xacc,void *P)
{
	int j;
	double xl,xh,swap,dx,f,rtf;
	const int MAXIT=100;

	double fl = (*F)(x1,P);
	double fh = (*F)(x2,P);
	if (fl*fh > 0.0) {
		errorLog.Add("FalsePos","Root must be bracketed in RTFLSP");
		return -1.0;
	}
	if (fl < 0.0) {
		xl=x1;
		xh=x2;
	} else {
		xl=x2;
		xh=x1;
		swap=fl;
		fl=fh;
		fh=swap;
	}
	dx = xh-xl;
	for(j=1;j<=MAXIT;j++) {
		rtf = xl+dx*fl/(fl-fh);
		f = (*F)(rtf,P);
		if (f<0.0) {
			xl=rtf;
			fl=f;
		} else {
			xh=rtf;
			fh=f;
		}
		dx = xh-xl;
		if(fabs(f)<xacc) return rtf;
	}
	errorLog.Add("FalsePos","Max iteration error");
	return rtf;
}

#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
/********************************************************************
Finds the minimum of a function using Brent's method.
Refer to numerical recipes for details.
********************************************************************/
double Brent(double ax,double bx,double cx,double(*f)(double,void*),double tol,void *P)
{
	int iter;
	double etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	double d=0.0,e=0.0;

	double a = ((ax < cx) ? ax : cx);
	double b = ((ax > cx) ? ax : cx);
	x=w=v=bx;
	fw=fv=fx=(*f)(x,P);

	if(errorLog.IsError()) {
		errorLog.Add("Brent");
		return 0;
	}

	for (iter=1;iter<=ITMAX;iter++) {
		xm=0.5*(a+b);
		tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
		if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
			return x;
		}
		if (fabs(e) > tol1) {
			r=(x-w)*(fx-fv);
			q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r;
			q=2.0*(q-r);
			if (q > 0.0) p = -p;
			q=fabs(q);
			etemp=e;
			e=d;
			if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
				d=CGOLD*(e=(x >= xm ? a-x : b-x));
			else {
				d=p/q;
				u=x+d;
				if (u-a < tol2 || b-u < tol2)
					d=SIGN(tol1,xm-x);
			}
		} else {
			d=CGOLD*(e=(x >= xm ? a-x : b-x));
		}
		u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
		fu=(*f)(u,P);
		if(errorLog.IsError()) {
			errorLog.Add("Brent");
			return 0;
		}
		if (fu <= fx) {
			if (u >= x) a=x; else b=x;
			SHFT(v,w,x,u)
			SHFT(fv,fw,fx,fu)
		} else {
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) {
				v=w;
				w=u;
				fv=fw;
				fw=fu;
			} else if (fu <= fv || v == x || v == w) {
				v=u;
				fv=fu;
			}
		}
	}
	errorLog.Add("BRENT","Max iteration error");
	return x;
}
#undef ITMAX

#define ITMAX 40//20
#define EPS 3.0e-8
/********************************************************************
Finds the zero of a function using Brent's method.
Refer to numerical recipes for details.
********************************************************************/
double Zbrent(double x1,double x2,double(*func)(double,void*),double tol,void *P)
{
	double a=x1,b=x2,c,d,e,min1,min2;
	double fc,p,q,r,u,tol1,xm;

	double fa = (*func)(a,P);
	if(errorLog.IsError()) {
		errorLog.Add("ZBrent","Initial fa");
		return 0;
	}
	double fb = (*func)(b,P);
	if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"Initial fb (fa=%lf,fb=%lf)",fa,fb);
		errorLog.Add("ZBrent",str);
		return 0;
	}
	if (fb*fa > 0.0) {
		char str[256];
		sprintf(str,"Root must be bracketed: fa(a=%.3lf)=%.3lf fb(b=%.3lf)=%.3lf",a,fa,b,fb);
		errorLog.Add("Zbrent",str);
		return b;
	}
	fc=fb;
	for (int i=1;i<=ITMAX;i++) {
		if (fb*fc > 0.0) {
			c=a;
			fc=fa;
			e=d=b-a;
		}
		if (fabs(fc) < fabs(fb)) {
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}
		tol1=2.0*EPS*fabs(b)+0.5*tol;
		xm=0.5*(c-b);
//		if (fabs(xm) <= tol1 || fb == 0.0) return b;
//		if(i>90) {sprintf(s,"%d, %12.10le, %12.10le, %12.10le, %12.10le\n",i,xm,fb,b,tol); bprintf(s);}
		if(fabs(fb) <= tol) {
			//errorLog.ClearError("ZBrent","fabs(fb)<=tol",0);
			return b;
		}
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
			u=fb/fa;
			if (a == c) {
				p=2.0*xm*u;
				q=1.0-u;
			} else {
				q=fa/fc;
				r=fb/fc;
				p=u*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(u-1.0);
			}
			if (p > 0.0)  q = -q;
			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e=d;
				d=p/q;
			} else {
				d=xm;
				e=d;
			}
		} else {
			d=xm;
			e=d;
		}
		a=b;
		fa=fb;
		if (fabs(d) > tol1)
			b += d;
		else
			b += (xm > 0.0 ? fabs(tol1) : -fabs(tol1));
		fb = (*func)(b,P);
		if(errorLog.IsError()) {
			errorLog.Add("ZBrent","fb");
			return b;
		}
	}
	if(fabs(fb)>1e-5) {
		char str[128];
		sprintf(str,"Max iterations. tol=%le/%le",fabs(fb),tol);
		errorLog.Add("ZBrent",str);
	}
	return b;
}
#undef ITMAX

// Plot function that is solved for by ZBrent.  It should have a zero
// between x1 and x2.
int ZbrentPlot(double x1,double x2,double(*func)(double,void*),double tol,void *P)
{
	FILE*fp=fopen("Zbrent.txt","w");
	if(fp==NULL) return 1;

	const double width = x2 - x1;
	for(int i=0;i<=10000;i++) {
		const double frac = (double)i/10.0;
		const double x = x1 + width*frac;
		errorLog.ClearError("ZbrentPlot");
		const double error = (*func)(x,P);
		fprintf(fp,"%lf,%le\n",x,error);
	}

	fclose(fp);
	return 0;
}

//One Demension spline

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
void spline(double *x,double *y,int n,double yp1,double ypn,double *y2)
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
splint(double *xa,double *ya,double *y2a,int n,double x,double *y,double *dy,int is_constant)
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

/***********************************************************************
get the output data from the spline interpolation
nData--number of data points
nProp--number of property
*Data_points--input data array
*Array--output data array
x--independent variable
********************************************************************/
void outdataSpline(int nData, int nProp, dataSpline Data_in,double x,double *Array)
{
	//make the data for spline
	int k=0,i=0,j=0,h=0,L=0;
	double min=10000.0,min2=-10000.0;
	double Data_list[58][20];
	double aa,da;
	double y2[58][20];
	int is_dx_constant=0;

	if(nData>19||nProp>57)
	{
	for(i=0;i<nProp;i++) Array[i+1]=0.0;
	return;
	}

	//put the data in uprising order
	k=nData-1;
	for(i=0;i<=k;i++)
	{
	for(h=0;h<=k;h++)
	{
	if(Data_in.Data_points[0][h]<min&&Data_in.Data_points[0][h]>min2) 
	{
	L=h;
	min=Data_in.Data_points[0][h];
	}
	}
	min2=min;
	min=10000;
	for(j=0;j<=nProp;j++)
	{
	Data_list[j][i]=Data_in.Data_points[j][L];
	}
	}
	
	//Spline
	for(i=0;i<=nProp;i++) {
	spline(Data_list[0],Data_list[i+1],k,1.0e30,1.0e30,y2[i]);
	if(errorLog.IsError()) {
	errorLog.ClearError("outdataSpline","Spline");
	}
	}
	
	for(i=0;i<=nProp;i++)
	{
	splint(Data_list[0],Data_list[i+1],y2[i],k,x,&aa,&da,is_dx_constant);
	Array[i+1] = aa;
	if(errorLog.IsError()) {
	errorLog.ClearError("outdataSpline","Spline");
	Array[i+1]=0.0;
	}
	}
	return;
}

/********************************************************************
Finds the zero of a function using bisection method and spline technique. Refer to
Rao's book for details.
********************************************************************/
double Bisection1DSpline(double Xlo,double Xhi,double(*F)(double,void*),double thres,void *P)
{
	double Xbuf,X,Y;
	double Fbuf;
	const int Imax=25;
	double Flo=0, Fhi=0;
	dataSpline Data_out;
	int j=0,h=0,k=0;
	double a=0,b=0;

	do{
	
	do{
	if(errorLog.IsError()) 
	{
	errorLog.ClearError("Bisection1DSpline", "Increase the low restriction");
	Xlo=Xlo+0.1*(Xhi-Xlo);
	}
	Flo = F(Xlo,P);
	j=j+1;
	}while(errorLog.IsError()&&j<10);

	if(errorLog.IsError()) 
	{
	errorLog.Add("Bisection1DSpline", "Increase the low restriction");
	return 0;
	}
	
	j=0;

	do{
	if(errorLog.IsError()) 
	{
	errorLog.ClearError("Bisection1DSpline", "Decreas the upper restriction");
	Xhi=Xhi-0.1*(Xhi-Xlo);
	}
	Fhi = F(Xhi,P);
	j=j+1;
	}while(errorLog.IsError()&&j<10);

	if(errorLog.IsError()) 
	{
	errorLog.Add("Bisection1DSpline", "Decrease the upper restriction");
	return 0;
	}

	j=0;

	a=Xlo; b=Xhi;

	if(Flo*Fhi>0) {
	Xlo=Xlo-0.05*(b-a);
	Xhi=Xhi+0.05*(b-a);		
	}

	h=h+1;
	}while(Flo*Fhi>0&&h<20);

	if(Flo*Fhi>0) 
	{
	char str[128];
	sprintf(str,"Root not bracketed. F(%lf)=%lf F(%lf)=%lf",Xlo,Flo,Xhi,Fhi);
	errorLog.Add("Bisection1D",str);
	return -999.0;
	}
	
	Data_out.Data_points[0][0]=Flo;
	Data_out.Data_points[0][1]=Fhi;
	Data_out.Data_points[1][0]=Xlo;
	Data_out.Data_points[1][1]=Xhi;
	j=2;

	if(Flo>Fhi) {
		Fbuf=Flo; Xbuf=Xlo;
		Flo=Fhi; Xlo=Xhi;
		Fhi=Fbuf; Xhi=Xbuf;
	}

	int i=0;
	k=0;
	do {
		X=(Xhi+Xlo)/2;
		i++ ;
		if(i>Imax) {

		if(k<20) k=j;

		if(k>0)
		{
		outdataSpline(k,1,Data_out,0,&X);
		return X;
		}
		else
		{
		errorLog.Add("Bisection1DSpline", "can't spline");
		return 0;
		}

		}

		Y = F(X,P);
		if(Y>0) Xhi=X; else Xlo=X;
		if(errorLog.IsError())
		{
		errorLog.ClearError("Bisection1DSpline", "Result by spline");
		/*if(k<20) k=j;
		outdataSpline(k,1,Data_out,0,&X);
		return X;*/
		Xhi=X;
		}
	
		//backup the data for spline
		{	
		Data_out.Data_points[0][j]=Y;
		Data_out.Data_points[1][j]=X;
		j=j+1;
		if(j==20) 
		{j=0;
		k=20;//overflow
		}
	}
	
	} while(fabs(Y)>thres);

	return X;
}
