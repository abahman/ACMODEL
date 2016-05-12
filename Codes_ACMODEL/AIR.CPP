#include "stdafx.h"

#include <stdio.h>
#include "air.h"
#include "spline.h"
#include "numeric.h"
#include "errorlog.h"
#include "spline1d.h"

extern ErrorLog errorLog;

extern AirProp air;

/* Local functions. */
double HaTaFunc(double,void*);

/********************************************************************
Converts the enthalpy of dry air to temperature of dry air.
********************************************************************/
double
HatoTa(double ha)
{
	double dh = HaTaFunc(TAIRMAX,&ha);
	if(errorLog.IsError()) {
		errorLog.Add("HatoTa");
		return 0;
	}

	if(dh>0) return 1e20;

	double T = Zbrent(TAIRMIN,TAIRMAX,HaTaFunc,1e-8,&ha);
	if(errorLog.IsError()) {
		errorLog.Add("HatoTa");
		return 0;
	}

	return T;
}

/********************************************************************
Function called by HatoTa().  It returns a residual that is zero
when the temperature corresponding to known enthalpy is selected.
********************************************************************/
double
HaTaFunc(double T,void *Params)
{
	double *ha = (double *)Params;
	double h = air.h(T);
	double dh = (*ha-h)/(*ha);
	return dh;
}
