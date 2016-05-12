#include "stdafx.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "wair.h"
#include "r22.h"
#include "num3d.h"
#include "numeric.h"
#include "mainfunc.h"
#include "errorlog.h"

double CEo = 0.0;			// convergence error

static double CError(double *x) {
	return sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
}

extern ErrorLog errorLog;

/********************************************************************
Solves three equations and three unknowns using a modified Newton's
method that optimizes the shotting distance if the convergence
error increases.
Inputs:  Xg[3] (vector containing initial guess)
			F(*x,*y,*P) (function that returns a vector of 3
				residuals.  The solution is defined when all residuals
				are zero.)
				x[3] (input vector)
				y[3] (residual vector)
				P (parameters)
			tol (solution tolerance.  Defines when to stop iteration)
			P (parameters passed to F())
Outputs: Xg[3] (solution vector)
			returns 0 -ok, >1 if error
********************************************************************/
int FindZero3DGolden(double *Xg,int(*F)(double*,double*,void*),double tol,void *P)
{
	static double ConvError,X0[3],X[4][3],Y[4][3];
	static double b,dx[3],a[3][3],c[3],d[2],e[2][2];
	static int i,j,k,l,Z,GStries;
	static double XG[3],YG[3],x0,x1,x2,X1,X2,f1,f2,D,L,dL;
	const int n=3,Imax=10;//10;20;30;

	for(i=0;i<3;i++) X[0][i]=Xg[i];

	/********
	// test
	{
	double XX = X[0][0];
	FILE *fp = fopen("num3d.dbg","w");
	for(int g=-10;g<=10;g++) {
		X[0][0] = XX *(1.0+(double)g/50.0);
		while(F(X[0],Y[0],P)) {}
		fprintf(fp,"%lf,%lf,%lf : %lf,%lf,%lf : %d\n",X[0][0],X[0][1],X[0][2],Y[0][0],Y[0][1],Y[0][2]);
	}
	fclose(fp);
	}
	// test
	*********/


	// Evaluate F at guess point.
	F(X[0],Y[0],P);
	ConvError=CError(Y[0]);
	if(errorLog.IsError()) {
		errorLog.Add("FindZero3DGolden","Evaluating initial guess");
		return 1;
	}
	CEo=CError(Y[0]);
	for(l=0;l<n;l++) Xg[l]=X[0][l];
	if(CEo<tol) {return 0;}

	for(i=1;i<=Imax;i++) {

		printf("FindZero3DGolden: i=%d\n",i);

	if(i==25)
	{
		int shenb=0;
	}

		// Function evaluations for partial derivatives.
		b=1e-3*CEo;
		for(j=1;j<=n;j++) {
			for(k=0;k<n;k++) X[j][k]=X[0][k];
			dx[j-1]=b*fabs(X[0][j-1]); X[j][j-1]+=dx[j-1];
			F(X[j],Y[j],P);
			ConvError=CError(Y[j]);
			// If moving a little dx in the positive X direction causes an
			// error, then try moving in the other direction.
			if(errorLog.IsError()) {
				for(k=0;k<n;k++) X[j][k]=X[0][k];
				dx[j-1]=-b*fabs(X[0][j-1]); X[j][j-1]+=dx[j-1] ;

				errorLog.ClearError("FindZero3DGolden","resetting");

				F(X[j],Y[j],P);
				ConvError = CError(Y[j]);

				// If this still doesn't work, then give up.
				if(errorLog.IsError()) {
					errorLog.Add("FindZero3DGolden","status<>0 finding dF/dx");
					return 2;
				}
			}

			if(ConvError<CEo) {
				CEo=ConvError;
				for(l=0;l<n;l++) Xg[l]=X[j][l];
			}
			if(CEo<=tol)
				return 0;

			// Test if dx is big enough to make a change in Y.
			// If not, then increas dx.
			Z=0; for(k=0;k<n;k++) if(Y[j][k]==Y[0][k]) Z=k+1 ;
			while(Z) {
				dx[j-1]*=2 ;
				if(fabs(dx[j-1]/X[0][j-1])>0.05) {
					for(k=0;k<n;k++) if(Y[j][k]==Y[0][k]) Y[j][k]=Y[0][k]+1e-20 ;
					Z=0 ;
				} else {
					X[j][j-1]=X[0][j-1]+dx[j-1] ;
					F(X[j],Y[j],P);
					ConvError=CError(Y[j]);
					if(errorLog.IsError()) {
						errorLog.Add("FindZero3DGolden","status<>0 increasing dx");
						return 3;
					}
					if(ConvError<CEo) {
						CEo=ConvError;
						for(l=0;l<n;l++) Xg[l]=X[j][l];
					}
					if (CEo<=tol) 0;
					Z=0; for(k=0;k<n;k++) if(Y[j][k]==Y[0][k]) Z=k+1 ;
				}
			}
		}

		// Calculate partial derivatives.
		for(j=0;j<n;j++) {
			c[j]=Y[0][j] ;
			for(k=0;k<n;k++) {
				a[j][k]=(Y[k+1][j]-Y[0][j])/dx[k] ;
				c[j]-=a[j][k]*X[0][k];
			}
		}
		for(l=0;l<n;l++) X0[l]=X[0][l];

		// Solve for new X using 1st order Taylor series approx.
		for(j=0;j<2;j++) d[j]=c[0]-a[0][0]/a[j+1][0]*c[j+1];
		for(j=0;j<2;j++) for(k=0;k<2;k++) e[j][k]=a[0][k+1]-a[0][0]/a[j+1][0]*a[j+1][k+1];
		X[0][2]=-(d[0]-e[0][0]/e[1][0]*d[1])/(e[0][1]-e[0][0]/e[1][0]*e[1][1]);
		X[0][1]=-(e[0][1]*X[0][2]+d[0])/e[0][0];
		X[0][0]=-(a[0][2]*X[0][2]+a[0][1]*X[0][1]+c[0])/a[0][0];

		// Test for possible increase in error.
		// If so, then find minimum using golden section method.
		F(X[0],Y[0],P);
		ConvError=CError(Y[0]);
		Z=0;
		while(errorLog.IsError()) {
			if(++Z>3) {
				errorLog.Add("FindZero3DGolden","Backing off from new guess");
				return 4;
			}
			for(j=0;j<n;j++) X[0][j]=X[0][j]-(X[0][j]-X0[j])/2;
			errorLog.AddWarning("FindZero3DGolden","Backing off from new step");

			errorLog.ClearError("FindZero3DGolden");

			F(X[0],Y[0],P);
			ConvError=CError(Y[0]);
		}
		if(ConvError<CEo) {
			CEo=ConvError;
			for(l=0;l<n;l++) Xg[l]=X[0][l];
		} else {
			D = 1-2/(1+sqrt(5.0));
			X1=0;
			X2=1;
			L=X2-X1;
			x1=X1+D*L;
			for(j=0;j<n;j++) XG[j]=X[0][j]-x1*(X[0][j]-X0[j]);
			F(XG,YG,P);
			f1=CError(YG);
			if(errorLog.IsError()) {
				errorLog.Add("FindZero3DGolden","Evaluating golden section");
				return 5;
			}

			// save this solution if it is the best so far.
			if(f1<CEo) {
				CEo=f1;
				for(l=0;l<n;l++) Xg[l]=XG[l];
			}

			// quit if solution is within specified tolerance
			if(CEo<=tol) return 0;

			x2=X2-D*L;
			for(j=0;j<n;j++) XG[j]=X[0][j]-x2*(X[0][j]-X0[j]);
			F(XG,YG,P);
			f2=CError(YG);
			if(errorLog.IsError()) {
				errorLog.Add("FindZero3DGolden","Evaluating golden section");
				return 6;
			}
			if(f2<CEo) {CEo=f2; for(l=0;l<n;l++) Xg[l]=XG[l];}
			if(CEo<=tol)
				return 0;
			dL=2*L/(X1+X2);

			GStries = 0;
			while (++GStries<5 && fabs(dL)>=0.1) {
				if (f1>f2) {
					X1=x1;
					x1=x2;
					f1=f2;
					L=X2-X1;
					x2=X2-D*L;
					dL=2*L/(X1+X2);
					for(j=0;j<n;j++) XG[j]=X[0][j]-x2*(X[0][j]-X0[j]);
					F(XG,YG,P);
					f2=CError(YG);
					if(errorLog.IsError()) {
						errorLog.Add("FindZero3DGolden","Evaluating golden section");
						return 7;
					}
					if(f2<CEo) {CEo=f2; for(l=0;l<n;l++) Xg[l]=XG[l];}
					if(CEo<=tol) return 0;
				} else if (f1<f2) {
					X2=x2;
					x2=x1;
					f2=f1;
					L=X2-X1;
					x1=X1+D*L;
					dL=2*L/(X1+X2) ;
					for(j=0;j<n;j++) XG[j]=X[0][j]-x1*(X[0][j]-X0[j]);
					F(XG,YG,P);
					f1=CError(YG);
					if(errorLog.IsError()) {
						errorLog.Add("FindZero3DGolden","Evaluating golden section");
						return 8;
					}
					if(f1<CEo) {CEo=f1; for(l=0;l<n;l++) Xg[l]=XG[l];}
					if(CEo<=tol)
						return 0;
				} else {
					X1=x1;
					X2=x2;
					L=X2-X1;
					dL=2*L/(X1+X2) ;
					x1=X1+D*L;
					for(j=0;j<n;j++) XG[j]=X[0][j]-x1*(X[0][j]-X0[j]);
					F(XG,YG,P);
					f1=CError(YG);
					if(errorLog.IsError()) {
						errorLog.Add("FindZero3DGolden","Evaluating golden section");
						return 9;
					}
					if(f1<CEo) {CEo=f1; for(l=0;l<n;l++) Xg[l]=XG[l];}
					if(CEo<=tol)
						return 0;
					x2=X2-D*L;
					for(j=0;j<n;j++) XG[j]=X[0][j]-x2*(X[0][j]-X0[j]);
					F(XG,YG,P);
					f2=CError(YG);
					if(errorLog.IsError()) {
						errorLog.Add("FindZero3DGolden","Evaluating golden section");
						return 10;
					}
					if(f2<CEo) {CEo=f2; for(l=0;l<n;l++) Xg[l]=XG[l];}
					if(CEo<=tol)
						return 0;
				}
			}
			x0=X1+L/2;
			for(j=0;j<n;j++) X[0][j]-=x0*(X[0][j]-X0[j]);
			F(XG,YG,P);
			ConvError=CError(YG);
			if(errorLog.IsError()) {
				errorLog.Add("FindZero3DGolden","Evaluating golden section");
				return 11;
			}
			if(ConvError<CEo) {CEo=ConvError; for(l=0;l<n;l++) Xg[l]=X[0][l];}
			if(CEo<=tol) return 0;
		}
	}

	if(i>=Imax) {
		errorLog.Add("FindZero3DGolden","Max iteration");
		return 12;
	}

	return 0;
}

/***
FindZero3DGuess

Solve a 3D nonlinear problem.  Start with known solution and work towards new solution.
Xg = guess cooresponding to solution for Pold parameters.  Set Pold=NULL if no previous solution.
Pold = parameters for which Xold is the solution
Pnew = paramters for desired solution.  Use Xold as initial guess.
F = residual function
tol = solution tolerance

Returns 0 if successful, and >0 if error.


TR 11/12/03 - THIS CODE IS BROKEN, BUT I DECIDED TO PUT OFF FIXING IT
***/
int FindZero3DGuess(double *Xg,int(*F)(double*,double*,void*),double tol,void *Pold,void *Pnew)
{
	int j;
	double X[3];
	InVars *Po,*Pn,R;
	const int Imax=20;

	Po = (InVars*)Pold;
	Pn = (InVars*)Pnew;

	for(j=0;j<3;j++) X[j]=Xg[j];

	FindZero3DGolden(X,F,tol,Pnew);
	if(!errorLog.IsError()) {
		// There is no error. Done.
		for(j=0;j<3;j++) Xg[j]=X[j];
		return 0;
	}

	double a = 0;
	double l = 0.5; // start with first back off, because if we are here there has already been a problem.
	double b0 = a + l;
	double b = b0;
	int i = 0;
	do {
		if(i>Imax) break;

		errorLog.ClearError("FindZero3DGuess");

		R.Tai = Po->Tai + b*(Pn->Tai - Po->Tai);
		R.TPi.T = Po->TPi.T + b*(Pn->TPi.T - Po->TPi.T);
		R.TPi.P = Po->TPi.P + b*(Pn->TPi.P - Po->TPi.P);
		R.Goa = Po->Goa + b*(Pn->Goa - Po->Goa);
		R.Gra = Po->Gra + b*(Pn->Gra - Po->Gra);
		R.Charge = Po->Charge + b*(Pn->Charge - Po->Charge);
		R.Dadj = Po->Dadj + b*(Pn->Dadj - Po->Dadj);
		R.nvAdj = Po->nvAdj + b*(Pn->nvAdj - Po->nvAdj);

		for(j=0;j<3;j++) X[j]=Xg[j];

		printf("++++++ FindZero3DGuess: i=%d, a=%lf, b=%lf, l=%lf\n",i,a,b,l);

		FindZero3DGolden(X,F,tol,&R);

		if(errorLog.IsError()) {
			i++; 
			l/=2;
		} else {
			a=b;
			for(j=0;j<3;j++) Xg[j]=X[j];
		}
		b=a+l;

	} while(errorLog.IsError() || b<=b0);

	if(errorLog.IsError()) {
		errorLog.Add("FindZero3DGuess");
		return 1;
	}

	return 0;
}


