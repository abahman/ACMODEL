#include "stdafx.h"

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "spline2d.h"
#include "errorlog.h"

extern ErrorLog errorLog;

LubriProp lubriPropComp("./Lubriprop/Lubriprop.tbl");

#ifdef _R22
RefPropSPthPT refshthPT("./R22prop/r22shth.tbl");
RefPropSPthPT refscthPT("./R22prop/r22scth.tbl");
RefPropSPthPH refshthPH("./R22prop/r22shth.tbl");
RefPropSPthPH refscthPH("./R22prop/r22scth.tbl");
RefPropSPtrPT refshtrPT("./R22prop/r22shtr.tbl");
RefPropSPtrPT refsctrPT("./R22prop/r22sctr.tbl");
//RefPropSPthPU refshthPU("r22shth.tbl");
//RefPropSPthPU refscthPU("r22scth.tbl");
//RefPropSPthPV refshthPV("r22shth.tbl");
#endif

#ifdef _R410A
RefPropSPthPT refshthPT("./R410Aprop/r410ashth.tbl");
RefPropSPthPT refscthPT("./R410Aprop/r410ascth.tbl");
RefPropSPthPH refshthPH("./R410Aprop/r410ashth.tbl");
RefPropSPthPH refscthPH("./R410Aprop/r410ascth.tbl");
RefPropSPtrPT refshtrPT("./R410Aprop/r410ashtr.tbl");
RefPropSPtrPT refsctrPT("./R410Aprop/r410asctr.tbl");
RefPropTPthPH refTPthPH("./R410Aprop/r410atpthx.tbl");
RefPropTPthPX refTPthPX("./R410Aprop/r410atpthx.tbl");
RefPropTPtrPX refTPtrPX("./R410Aprop/r410atptrx.tbl");
#endif

#ifdef _R407C
RefPropSPthPT refshthPT("./R407Cprop/r407cshth.tbl");
RefPropSPthPT refscthPT("./R407Cprop/r407cscth.tbl");
RefPropSPthPH refshthPH("./R407Cprop/r407cshth.tbl");
RefPropSPthPH refscthPH("./R407Cprop/r407cscth.tbl");
RefPropSPtrPT refshtrPT("./R407Cprop/r407cshtr.tbl");
RefPropSPtrPT refsctrPT("./R407Cprop/r407csctr.tbl");
RefPropTPthPH refTPthPH("./R407Cprop/r407ctpthx.tbl");
RefPropTPthPX refTPthPX("./R407Cprop/r407ctpthx.tbl");
RefPropTPtrPX refTPtrPX("./R407Cprop/r407ctptrx.tbl");
#endif


WetAirPropTR wair("./Airprop/wair.tbl");

Spline2d::Spline2d(const char* fn,int _xcol,int _ycol) : Spline(fn)
{
	TBL = NULL;
	y2 = NULL;
	y = NULL;
	x1new = NULL;
	x = NULL;
	x2n = NULL;

	xcol = _xcol;
	ycol = _ycol;

	is_dx1_constant = 0;

	ScanFile();
	if(errorLog.IsError()) {
		errorLog.Add("Spline","Spline");
		return;
	}
	AllocateMemory();
	if(errorLog.IsError()) {
		errorLog.Add("Spline","Spline");
		return;
	}
	LoadData();
	if(errorLog.IsError()) {
		errorLog.Add("Spline","Spline");
		return;
	}
	setup_error = 0;
}

Spline2d::~Spline2d()
{
	RestoreMemory();
}

// Opens file. Determines memory requirements.  Allocates needed
// memory. Loads file into memory.  Closes file.
//
// Assumes that first column is fixed and the second column
// changes first, such as:
// 1.0 1.0 ...
// 1.0 2.0 ...
// 1.0 3.0 ...
// 2.0 1.0 ...
// 2.0 2.0 ...
// ...
//
// Returns:
//    0 = all ok
//    1 = can not open file
//    2 = number of columns in each row is not identical
void Spline2d::ScanFile()
{
	FILE *fp = fopen(filename,"r") ;
	if(fp==NULL) {
		char msg[128];
		sprintf(msg,"Can not open file %s",filename);
		errorLog.Add("Spline2d::ScanFile",msg);
		return;
	}

	// Count the number of columns in the first row.  Assume
	// that the same number of columns are in ever row.  Save
	// as "nProp" which excludes the first two columns containing
	// the independent variables.
	int cnt=0;
	char s[256];
	char *pnt=s;
	double value;
	fgets(s,250,fp);
	while(!ScanStringForDouble(&pnt,value)) { cnt++; }
	nProp = cnt-2;

	// Load x (column 1) and y (column 2) from first row
	int x_fnd=0,y_fnd=0;
	double x_save,y_save;
	pnt=s;
	for(int i=0;i<nProp+2;i++) {
		if(i==xcol) {
			ScanStringForDouble(&pnt,x_save);
			x_fnd = 1;
		}
		if(i==ycol) {
			ScanStringForDouble(&pnt,y_save);
			y_fnd = 1;
		}
	}

	// Return with error is specified xcol or ycol not found in table.
	if(!x_fnd || !y_fnd) {
		errorLog.Add("Spline2d::ScanFile","xcol or ycol not found");
		return;
	}

	// Count the number of X groups (nX) and the maximum
	// number of Y groups for any particular X (nY).  Also
	// count the total number of rows (nData).
	double x,y;
	int nYcnt=1;
	nX=2;       // one extra
	nY=1;       // one extra
	nData=2;    // one extra
	while(fgets(s,250,fp)!=NULL) {

		pnt=s;
		cnt=0;
		while(!ScanStringForDouble(&pnt,value)) {
			if(cnt==xcol) x = value;
			if(cnt==ycol) y = value;
			cnt++;
		}

		if(cnt!=nProp+2) {
			char msg[128];
			sprintf(msg,"cnt(%d) != nProp(%d)",cnt,nProp);
			errorLog.Add("Spline2d::ScanFile",msg);
			return;
		}

		nData++;
		if(x!=x_save) {
			nX++;
			x_save=x;
			if(nYcnt>nY) nY=nYcnt;
			nYcnt=1;
		} else {
			nYcnt++;
		}
	}

	fclose(fp);
}

void Spline2d::AllocateMemory()
{
	int i,j;

	TBL = new double*[nProp+2];
	for(i=0;i<nProp+2;i++) TBL[i] = new double[nData];

	y2 = new double **[nProp];
	for(i=0;i<nProp;i++) {
		y2[i] = new double*[nX];
		for(j=0;j<nX;j++) y2[i][j] = new double [nY];
	}

	y = new double**[nProp];
	for(i=0;i<nProp;i++) y[i] = new double*[nX];

	x1new = new double[nX];
	x = new double*[nX];
	x2n = new int[nX];
}

void Spline2d::RestoreMemory()
{
	int i,j;

	if(TBL) {
		for(i=0;i<nProp+2;i++) delete TBL[i];
		delete TBL;
	}

	if(y2) {
		for(i=0;i<nProp;i++) {
			for(j=0;j<nX;j++) delete y2[i][j];
			delete y2[i];
		}
		delete y2;
	}

	if(y) {
		for(i=0;i<nProp;i++) delete y[i];
		delete y;
	}

	if(x1new) delete x1new;
	if(x) delete x;
	if(x2n) delete x2n;
}

// Returns:
//     0 = all ok
//     1 = can not open file
//     2 = too many x values for memory allocated
//     3 = too many y values for memory allocated
//     4 = too many data points for memory allocated
//
// Required file data format:
//
// Fix x1, change x2 in following lines, then change x1
// and hold fixed again while changing x2, etc.
//
// The data must have increasing values of x1 and x2
//
// x1   x2   y1   y2   y3   y4  ...
// 10   10   2    3    2    2
// 10   20
// 10   30
// 20   15
// 20   25
void Spline2d::LoadData()
{
	x2min = +1e20;
	x2max = -1e20;
	dx2 = 1e20;

	FILE *fp=fopen(filename,"r");
	if(fp==NULL) {
		char msg[128];
		sprintf(msg,"Can not open file %s",filename);
		errorLog.Add("Spline2d::LoadData",msg);
		return;
	}

	// x1 is the "xcol" column in the file and the first independent variable
	// x2 is the "ycol" column in the file and the second independent variable
	double *x1 = TBL[xcol];
	double *x2 = TBL[ycol];

	// load the first line in the file
	int cnt=0;
	char s[256];
	char *pnt=s;
	double value;
	fgets(s,250,fp);
	while(!ScanStringForDouble(&pnt,value)) {
		TBL[cnt][0] = value;
		cnt++;
	}

	int i = 0;         // "i" increments for each line in the file
	int k = 0;         // "k" increments for each new value of x1
	int i_x1new = 0;   // line in file (i) for first new x1 value
	x1new[k] = x1[0];  // first value of x1new is in first line of file
	x[k] = x2;         // x[k] is a pointer to an array of x2 values for a fixed value of x1

	// y[j][k] is an array of y values for each property and each
	// new value of x1.  The properties are numbered in the order
	// in which they apprear in the input table, skipping over the
	// independent variables.
	cnt = 0;
	for(int j=0;j<nProp+2;j++) {
		if(j!=xcol && j!=ycol) y[cnt++][k]=&TBL[j][0];
	}

	is_dx1_constant = 1;

	double dx1;
	while(fgets(s,250,fp)!=NULL) {

		// Test if "i" is incrementing beyond limits of table.
		// If so, generate an error.
		if (++i>=nData) {
			char msg[128];
			sprintf(msg,"i(%d) >= nData(%d)",i,nData);
			errorLog.Add("Spline2d::LoadData",msg);
			return;
		}

		// Load new line from file.
		cnt = 0;
		pnt = s;
		while(!ScanStringForDouble(&pnt,value)) {
			TBL[cnt][i] = value;
			cnt++;
		}

		// Test if x1 value has changed.
		if(x1[i]!=x1new[k]) {

			// Test if dx1 is constant throughout the file.  If so, a more
			// efficient search algorithm can be used.
			if(is_dx1_constant) {
				if(k>0) {
					const double this_dx1 = x1[i] - x1new[k];
					const double error = fabs(dx1-this_dx1);
					if(error>1e-6) is_dx1_constant=0;
				} else {
					// first available dx1
					dx1 = x1[i] - x1new[k];
				}
			}

			// x2n[k] is the number of x2 values for the previous value
			// of x1 that just changed.
			x2n[k] = i - i_x1new;

			// Test to see if too many x2 values are being used.
			// If so, generate an error.
			if(x2n[k]>nY) {
				char msg[128];
				sprintf(msg,"x2n[k](%d) > nY(%d)",x2n[k],nY);
				errorLog.Add("Spline2d::LoadData",msg);
				return;
			}

			// Increment "k", corresponding to a new value of x1.
			// Test if too many values of "k" are being used.
			// Generate an error if this is true.
			if(++k>=nX) {
				char msg[128];
				sprintf(msg,"k(%d) >= nX(%d)",k,nX);
				errorLog.Add("Spline2d::LoadData",msg);
				return;
			}

			// Update the x and y arrays for the new value of x1.
			x[k] = &x2[i];
			cnt = 0;
			for(int j=0;j<nProp+2;j++) {
				if(j!=xcol && j!=ycol) y[cnt++][k]=&TBL[j][i];
			}

			// "x1new[k]" is the x1 value for this new group of x1 values.
			x1new[k] = x1[i];

			// "i_x1new" is the line in the file cooresponding to a new
			// value of x1.  It is use to calculate the number of x2 values
			// in each group of x1 values.
			i_x1new = i;

			// Assuming that x2 increases for each value of x1, the largest
			// value of x2 for the previous value of x1 is its last value.
			// Therefore, compare this value to the current value of the
			// overall greatest x2 (x2max) to see if x2max needs updating.
			//if(x2[i-1]>x2max) x2max=x2[i-1];

		} else {

			// This line's x1 value is the same as the previous line

			// "dx2" is the minimum distance between two x2 values
			// for a fixed x1 value
			if(fabs(x2[i]-x2[i-1])<dx2) dx2 = fabs(x2[i]-x2[i-1]);

		}

		// find largest and smallest x2.
		if(x2[i]<x2min) x2min=x2[i];
		if(x2[i]>x2max) x2max=x2[i];

	}

	//printf("filename=%s is_dx1_constant=%d dx1=%lf\n",filename,is_dx1_constant,dx1);
	//printf("filename=%s x2max=%lf\n",filename,x2max);

	// Number of x2 values included in last x1 value in table.
	x2n[k] = i - i_x1new;

	// Total number of x1 values.
	x1n = k+1;

	fclose(fp);

	// Calculate the array of second derivatives for each property
	// and each value of x1.
	for(k=0;k<nProp;k++) {
		for(int j=0;j<x1n;j++) {
			spline(x[j],y[k][j],x2n[j],1.0e30,1.0e30,y2[k][j]);
		}
	}
}

// return_value = 1 or 2  -> fills both a and dadx1
//                3       -> fills a, dadx1, and dadx2
void Spline2d::Evaluate_y_dy(int prop,double x1,double x2,int return_value,double* a,double* dadx1,double* dadx2)
{
	int i,j,k;
	double ytmp[4],yytmp[4],w2[4];

	if(setup_error) {
		errorLog.Add("Spline2d::Evaluate","Spline2d never initialized properly.");
		return;
	}

	if(prop<0 || prop>=nProp) {
		char msg[128];
		sprintf(msg,"Property out-of-range (0 < %d < %d)",prop,nProp);
		errorLog.Add("Spline2d::Evaluate",msg);
		return;
	}

	// check if x2 is out of range.
	const double wf2 = extrap_frac*(x2max-x2min);
	if(x2<x2min-wf2 || x2>x2max+wf2) {
		char msg[128];
		sprintf(msg,"Data x2 out-of-range (%le<%le<%le)",x2min,x2,x2max);
		errorLog.Add("Spline2d::Evaluate",msg);
		return;
	}

	const double wf1 = extrap_frac*(x1new[x1n-1]-x1new[0]);
	if(x1<x1new[0]-wf1 || x1>x1new[x1n-1]+wf1) {
		char msg[128];
		sprintf(msg,"Data x1 out-of-range (%le<%le<%le) in file [%s]",x1new[0],x1,x1new[x1n-1],filename);
		errorLog.Add("Spline2d::Evaluate",msg);
		return;
	}
	if(x1<x1new[0] || x1>x1new[x1n-1]) {
		char msg[256];
		sprintf(msg,"Extrapolating x1 within allowed %lf%% of table full width (%le<%le<%le).",100.0*extrap_frac,x1new[0],x1,x1new[x1n-1]);
		errorLog.AddWarning("Spline2d::Evaluate",msg);
	}

	// Find place in table bracketing x1.
	int klo,khi;
	if(is_dx1_constant) {
		// Assumes x1new[] are equally spaced.  This is more
		// efficient than direct search.
		const double k = (x1-x1new[0])/(x1new[1]-x1new[0]);
		klo = (int)k;
		if(klo<0) klo=0; else if(klo>x1n-2) klo=x1n-2;
		khi = klo + 1;
	} else {
		// Direct search
		klo=0;
		khi=x1n-1;	
		while (khi-klo > 1) {
			k=(khi+klo) >> 1;
			if (x1new[k]>x1) khi=k; else klo=k;
		}
	}

	// Deterimine which x1 values are needed.
	if(klo>0) klo-=1;
	if(khi<x1n-1) khi+=1;

	// Evaluate splines for needed x1 values.
	i = 0;
	for(j=klo;j<=khi;j++) {
		if(i>3) {
			errorLog.Add("Spline2d::Evaluate","Bracketing error 2");
			return;
		}
		double dy;
		splint(x[j],y[prop][j],y2[prop][j],x2n[j],x2,&yytmp[i],&dy);
		if(errorLog.IsError()) {
			errorLog.Add("Spline2d::Evaluate","splint");
			return;
		}
		w2[i]=x1new[j];
		i++ ;
	}

	// Generate and evaluate spline for x1 dependence,
	// and calculate da/dx1.
	spline(w2,yytmp,i,1.0e30,1.0e30,ytmp);
	splint(w2,yytmp,ytmp,i,x1,a,dadx1);
	if(errorLog.IsError()) {
		errorLog.Add("Spline2d::Evaluate");
		return;
	}

	//Calculate da/dx2
	if (return_value == 3) {

		//Get the 4 table values in the x2 direction
		double y_x2[4];
		double x0[5];
		x0[0] = x2 - 1.5*dx2;
		for (i=0; i<4; i++) {
			y_x2[i] = Evaluate(prop,x1,x0[i],1);
			//printf("%lf  %lf  %lf\n",x1,x0[i],y_x2[i]);
			x0[i+1]= x0[i]+dx2;
		}

		//Calculate the spline for the new tabulated values
		spline(x0,y_x2,4,1.0e30,1.0e30,ytmp);

		splint(x0,y_x2,ytmp,4,x2,a,dadx2);
		if(errorLog.IsError()) {
			errorLog.Add("Spline2d","Evaluate");
			return;
		}
	}
}

// Evaluates the 2D spline table.
//
// return_value = 1  -> returns property value a(x1,x2)
//                2  -> returns da_dx1(x1,x2)
//                3  -> returns da_dx2(x1,x2)
double Spline2d::Evaluate(int prop,double x1,double x2,int return_value)
{
	double a,dadx1,dadx2;
	Evaluate_y_dy(prop,x1,x2,return_value,&a,&dadx1,&dadx2);

	switch(return_value) {
		case 2: return dadx1;
		case 3: return dadx2;
	}
	return a;
}

double Spline2d::Evaluate_y(int prop,double x1,double x2)
{
	double a,dadx1,dadx2;
	Evaluate_y_dy(prop,x1,x2,1,&a,&dadx1,&dadx2);
	return a;
}

double Spline2d::Evaluate_dydx1(int prop,double x1,double x2)
{
	double a,dadx1,dadx2;
	Evaluate_y_dy(prop,x1,x2,2,&a,&dadx1,&dadx2);
	return dadx1;
}

double Spline2d::Evaluate_dydx2(int prop,double x1,double x2)
{
	double a,dadx1,dadx2;
	Evaluate_y_dy(prop,x1,x2,3,&a,&dadx1,&dadx2);
	return dadx2;
}

//----------------------- RefPropSPthPH -----------------------------
// dv_dP = v/P * dlog10(v)_dlog10(P)
// A constant added to y does not propagate through the derivative.
double RefPropSPthPH::dy_dx1(int prop,double P,double h) 
{
#ifdef REFPROP_LOG10_2D
	double y,dydx1,dydx2;
	Evaluate_y_dy(prop,IN1(P),IN2(h),2,&y,&dydx1,&dydx2);
	return pow(10,y)/P*dydx1;
#else
	return 1e-3*Evaluate_dydx1(prop,IN1(P),IN2(h));
#endif
}

// dv_dh = v/(h+1e6) * dlog10(v)_dlog10(h+1e6)
// A constant added to y does not propagate through the derivative.
double RefPropSPthPH::dy_dx2(int prop,double P,double h) 
{
#ifdef REFPROP_LOG10_2D
	double y,dydx1,dydx2;
	Evaluate_y_dy(prop,IN1(P),IN2(h),3,&y,&dydx1,&dydx2);
	return pow(10,y)/(h+1e6)*dydx2;
#else
	return Evaluate_dydx2(prop,IN1(P),IN2(h));
#endif
}

double RefPropSPthPH::T(double P,double h)
{
	double T = OUT_T(Evaluate_y(0,IN1(P),IN2(h)));
	if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"T(P=%lf,h=%lf)=%lf (%s)",P,h,T,filename);
		errorLog.Add("RefPropSPthPH::T",str);
		return 0;
	}
	return T;
}

void RefPropSPthPH::rhoPlusDerivatives(double P,double h,double* rho,double* drho_dP,double* drho_dh)
{
	Evaluate_y_dy(4,IN1(P),IN2(h),3,rho,drho_dP,drho_dh);
	if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"rhoPlusDerivatives(P=%lf,h=%lf) (%s)",P,h,filename);
		errorLog.Add("RefPropSPthPH::rhoPlusDerivatives",str);
	}
	printf("drho_dP=%le,drho_dh=%le\n",*drho_dP,*drho_dh);
#ifdef REFPROP_LOG10_2D
	(*drho_dP) *= pow(10,*rho)/P;
	(*drho_dh) *= pow(10,*rho)/(h+1e6);
	(*rho) = OUT(*rho);
#else
	(*drho_dP) *= 1e-3;
#endif
}


