#include "stdafx.h"

#include <stdio.h>
#include <math.h>
#include "spline1d.h"
#include "errorlog.h"

extern ErrorLog errorLog;

#ifdef _R22
RefPropTPthP reftpvthP("./R22prop/r22tpvth.tbl");
RefPropTPthP reftplthP("./R22prop/r22tplth.tbl");
RefPropTPthT reftpvthT("./R22prop/r22tpvth.tbl");
RefPropTPthT reftplthT("./R22prop/r22tplth.tbl");
RefPropTPtrP reftpltrP("./R22prop/r22tpltr.tbl");
RefPropTPtrP reftpvtrP("./R22prop/r22tpvtr.tbl");
#endif

#ifdef _R410A
RefPropTPthP reftpvthP("./R410Aprop/r410atpvth.tbl");
RefPropTPthP reftplthP("./R410Aprop/r410atplth.tbl");
RefPropTPthT reftpvthT("./R410Aprop/r410atpvth.tbl");
RefPropTPthT reftplthT("./R410Aprop/r410atplth.tbl");
RefPropTPtrP reftpltrP("./R410Aprop/r410atpltr.tbl");
RefPropTPtrP reftpvtrP("./R410Aprop/r410atpvtr.tbl");
#endif

#ifdef _R407C
RefPropTPthP reftpvthP("./R407Cprop/r407ctpvth.tbl");
RefPropTPthP reftplthP("./R407Cprop/r407ctplth.tbl");
RefPropTPthT reftpvthT("./R407Cprop/r407ctpvth.tbl");
RefPropTPthT reftplthT("./R407Cprop/r407ctplth.tbl");
RefPropTPtrP reftpltrP("./R407Cprop/r407ctpltr.tbl");
RefPropTPtrP reftpvtrP("./R407Cprop/r407ctpvtr.tbl");
#endif


AirProp air("./Airprop/air.tbl");

double RefPropTPthP::dy_dx(int prop,double x)
{
#ifdef REFPROP_LOG10_1D
	double y,dy_dx;
	Evaluate_y_dy(prop,IN(x),&y,&dy_dx);
	return pow(10,y)/x*dy_dx;
#else
	return 1e-3*Evaluate_dy(prop,IN(x));
#endif
}

void 
RefPropTPthP::h_dhdP(double P,double* h,double* dh_dP) 
{
	Evaluate_y_dy(1,IN(P),h,dh_dP);
#ifdef REFPROP_LOG10_1D
	(*dh_dP)*=pow(10,*h)/P;
	(*h)=OUT_HUS(*h);
#else
	(*dh_dP)*=1e-3;
#endif
}

void 
RefPropTPthP::rho_drhodP(double P,double* rho,double* drho_dP) 
{
	Evaluate_y_dy(5,IN(P),rho,drho_dP); 
#ifdef REFPROP_LOG10_1D
	(*drho_dP)*=pow(10,*rho)/P;
	(*rho)=OUT(*rho);
#else
	(*drho_dP)*=1e-3;
#endif
}

//----------------------- Spline1d -------------------------------
Spline1d::Spline1d(const char* fn,int _xcol) : Spline(fn)
{
	y = NULL;
	y2 = NULL;

	xcol = _xcol;

	is_dx_constant = 0;

	ScanFile();
	if(errorLog.IsError()) {
		errorLog.Add("Spline1d::Spline1d","ScanFile");
		return;
	}
	AllocateMemory();
	if(errorLog.IsError()) {
		errorLog.Add("Spline1d::Spline1d","AllocateMemory");
		return;
	}
	LoadData();
	if(errorLog.IsError()) {
		errorLog.Add("Spline1d::Spline1d","LoadData");
		return;
	}
	setup_error = 0;
}

Spline1d::~Spline1d()
{
	RestoreMemory();
}

void
Spline1d::AllocateMemory()
{
	int i;

	y = new double*[nProp+1];
	for(i=0;i<nProp+1;i++) y[i] = new double[nData];
	y2 = new double*[nProp];
	for(i=0;i<nProp;i++) y2[i] = new double[nData];
}

void
Spline1d::RestoreMemory()
{
	int i;

	if(y) {
		for(i=0;i<nProp+1;i++) delete y[i];
		delete y;
	}
	if(y2) {
		for(i=0;i<nProp;i++) delete y2[i];
		delete y2;
	}
}

void
Spline1d::ScanFile()
{
	FILE *fp = fopen(filename,"r");
	if(fp==NULL) {
		errorLog.Add("Spline1d::ScanFile","File not found");
		return;
	}

	// Count the number of columns in the first row.  Assume
	// that the same number of columns are in ever row.  Save
	// as "nProp" which excludes the first column containing
	// the independent variable.
	int cnt=0;
	char s[256];
	char *pnt=s;
	double value;
	fgets(s,250,fp);
	while(!ScanStringForDouble(&pnt,value)) { cnt++; }
	nProp = cnt-1;

	// Count the total number of rows (nData).
	double x;
	nData=2;    // one extra
	while(fgets(s,250,fp)!=NULL) {
		pnt=s;
		if(ScanStringForDouble(&pnt,x)) continue;

		cnt=0;
		while(!ScanStringForDouble(&pnt,value)) { cnt++; }
		if(cnt!=nProp) {
			char msg[128];
			sprintf(msg,"cnt(%d) != nProp(%d)",cnt,nProp);
			errorLog.Add("Spline1d::ScanFile",msg);
			return;
		}

		nData++;
	}

	fclose(fp);
}

void
Spline1d::LoadData()
{
	char s[256];
	char *pnt=s;

	FILE *fp = fopen(filename,"r");
	if(fp==NULL) {
		errorLog.Add("Spline1d::LoadData","File not found");
		return;
	}

	xmax = -1e20;
	xmin = +1e20;
	is_dx_constant = 1;

	if(xcol<0 || xcol>nProp) {
		errorLog.Add("Spline1d::LoadData","Invalid xcol.");
		return;
	}

	int i=0;
	double x0,dx;
	while(fgets(s,250,fp)!=NULL) {
		int col = 0;
		pnt = s;

		int j = 0;
		double a;
		while(!ScanStringForDouble(&pnt,a)) {
			if(j==xcol) {
				y[0][i] = a;
			} else {
				y[++col][i] = a;
			}
			j++;
		}
		if(y[0][i]>xmax) xmax = y[0][i];
		if(y[0][i]<xmin) xmin = y[0][i];

		// Test here to see if dx is constant.  If it is, then is_dx_constant
		// will be true and dx will be the constant value.  This will allow
		// for a more efficient search algorithm.
		if(is_dx_constant) {
			if(i>1) {
				double this_dx = y[0][i]-y[0][i-1];
				double error = fabs(dx-this_dx);
				if(error>1e-6) is_dx_constant=0;
			} else if(i>0) {
				dx = y[0][i] - y[0][i-1];
			} else {
				x0 = y[0][i];
			}
		}

		if(i>=nData) {
			errorLog.Add("Spline1d::LoadData","Too many data points");
			return;
		}
		i++;
	}

	fclose(fp);

	//printf("filename=%s is_dx_constant=%d x0=%lf dx=%lf\n",filename,is_dx_constant,x0,dx);

	for(i=0;i<nProp;i++) {
		spline(y[0],y[i+1],nData-1,1.0e30,1.0e30,y2[i]);
	}
}

void
Spline1d::Evaluate_y_dy(int prop,double x,double* a,double* da)
{
	if(setup_error) {
		errorLog.Add("Spline1d::Evaluate","Spline1d never initialized properly.");
		return;
	}

	if(prop<0 || prop>=nProp) {
		errorLog.Add("Spline1d::Evaluate","Property out of range");
		return;
	}

	double wf = extrap_frac*(xmax-xmin);
	if(x<xmin-wf || x>xmax+wf) {
	//if(x<xmin || x>xmax) {
		char str[128];
		sprintf(str,"Input out of range (x=%lf >xmax=%lf or <xmin=%lf)",x,xmax,xmin);
		errorLog.Add("Spline1d::Evaluate",str);
		return;
	}

	splint(y[0],y[prop+1],y2[prop],nData-1,x,a,da);
	if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"splint: prop=%d x=%lf a=%lf da=%lf",prop,x,*a,*da,is_dx_constant);
		errorLog.Add("Spline1d::Evaluate",str);
	}
}

double
Spline1d::Evaluate(int prop,double x)
{
	double a,da;
	Evaluate_y_dy(prop,x,&a,&da);
	return a;
}

double
Spline1d::Evaluate_dy(int prop,double x)
{
	double a,da;
	Evaluate_y_dy(prop,x,&a,&da);
	return da;
}


