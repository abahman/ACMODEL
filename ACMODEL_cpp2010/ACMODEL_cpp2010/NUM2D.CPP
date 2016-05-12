#include "stdafx.h"

#include "r22.h"
#include "numeric.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#define TOL 1e-8
#include "errorlog.h"

extern ErrorLog errorLog;
extern double CEo;

static double CError(double *x) 
{
	return sqrt(x[0]*x[0]+x[1]*x[1]);
}

/********************************************************************
Solves two simultaneous equations using Newton-Raphson method.
********************************************************************/
double FindZero2D(double *Xg,double(*F)(double*,double*,void*,int*),
	double tol,void *P,int *status) {

	double ConvError = 0.0;
	double X[3][2];
	double Y[3][2];
	double dx[2], a[2][2], c[2] ;
	const double b=1e-3 ;
	const int n=2, Imax=20 ;
	int i, j, k, Z ;

	for(int i=0;i<n;i++) X[0][i]=Xg[i];
	F(X[0],Y[0],P,status) ;
	if(*status) {printf("FindZero2D: \n"); return ConvError;}
	for(int i=1;i<=Imax;i++) {
		// Function evaluations.
		for(j=1;j<=n;j++) {
			for(k=0;k<n;k++) X[j][k]=X[0][k] ;
			dx[j-1]=b*fabs(X[0][j-1]); X[j][j-1]+=dx[j-1] ;
			F(X[j],Y[j],P,status) ;
			if(*status) {printf("FindZero2D: \n"); return ConvError;}
			Z=0; for(k=0;k<n;k++) if(Y[j][k]==Y[0][k]) Z=k+1 ;
			while(Z) {
				dx[j-1]*=2 ;
				if(fabs(dx[j-1]/X[0][j-1])>0.05) {
					for(k=0;k<n;k++) if(Y[j][k]==Y[0][k]) Y[j][k]=Y[0][k]+1e-20 ;
					Z=0 ;
				} else {
					X[j][j-1]=X[0][j-1]+dx[j-1] ;
					F(X[j],Y[j],P,status) ;
					Z=0; for(k=0;k<n;k++) if(Y[j][k]==Y[0][k]) Z=k+1 ;
				}
				if(*status) {printf("FindZero2D: \n"); return ConvError;}
			}
		}
		// Calculate partial derivatives.
		for(int j=0;j<n;j++) {
			c[j]=Y[0][j] ;
			for(k=0;k<n;k++) {
				a[j][k]=(Y[k+1][j]-Y[0][j])/dx[k] ;
				c[j]-=a[j][k]*X[0][k];
			}
		}
		// Solve for new X using 1st order Taylor series approx.
		X[0][1]=-(c[0]-a[0][0]/a[1][0]*c[1])/(a[0][1]-a[0][0]/a[1][0]*a[1][1]) ;
		X[0][0]=-(a[0][1]*X[0][1]+c[0])/a[0][0] ;
		// Evaluate new X.
//		printf(">>>>>> X[0]=%lf, X[1]=%lf\n",X[0][0],X[0][1]) ;
		F(X[0],Y[0],P,status) ;
		ConvError=sqrt(Y[0][0]*Y[0][0]+Y[0][1]*Y[0][1]) ;
		if(*status) {printf("FindZero2D: \n"); return ConvError;}
		if (ConvError<=tol) break ;
	}
	if (i>=Imax) {
		printf("FindZero2D: Max iteration error. ConvError=%le\n",ConvError);
		*status=1 ;
		return ConvError ;
	}
	for(int i=0;i<n;i++) Xg[i]=X[0][i];

	return ConvError ;
}

/********************************************************************
Solves to simultaneous equations using Newton-Raphson method, 
but in addition it forces the solution to remain within a feasible
space defined by a constrain function.  For example, this can be used
to avoid regions where the function is not defined.
********************************************************************/
double FindZero2DConst(double *Xg,double (*F)(double*,double*,void*),int(*G)(double*,double*,void*),double tol,void *P)
{
	double ConvError = 0.0;
	double X[3][2];
	double Y[3][2];
	double X0[2];
	double dx[2],a[2][2],c[2];
	const double b=1e-4 ;
	const int n=2,Imax=200;
	int i,j,k,Z;

	for(i=0;i<n;i++) X[0][i]=Xg[i];

	F(X[0],Y[0],P);
	if(errorLog.IsError()) {
		errorLog.Add("FindZero2DConst");
		return ConvError;
	}

	for(i=1;i<=Imax;i++) {
		// Function evaluations.
		for(j=1;j<=n;j++) {
			for(k=0;k<n;k++) X[j][k]=X[0][k] ;
			dx[j-1]=b*fabs(X[0][j-1]); X[j][j-1]+=dx[j-1] ;
			// Test for constraint violation.
			// If so, then move the other way.
			if(G(X[j],Y[j],P)) {
				for(k=0;k<n;k++) X[j][k]=X[0][k] ;
				dx[j-1]=-b*fabs(X[0][j-1]); X[j][j-1]+=dx[j-1] ;
			}
			if(errorLog.IsError()) {
				errorLog.Add("FindZero2DConst");
				return ConvError;
			}
			F(X[j],Y[j],P);
			if(errorLog.IsError()) {
				errorLog.Add("FindZero2DConst");
				return ConvError;
			}
			Z=0; for(k=0;k<n;k++) if(Y[j][k]==Y[0][k]) Z=k+1;
			while(Z) {
				dx[j-1]*=2;
				if(fabs(dx[j-1]/X[0][j-1])>0.05) {
					for(k=0;k<n;k++) if(Y[j][k]==Y[0][k]) Y[j][k]=Y[0][k]+1e-20;
					Z=0;
				} else {
					X[j][j-1]=X[0][j-1]+dx[j-1];
					F(X[j],Y[j],P);
					Z=0; for(k=0;k<n;k++) if(Y[j][k]==Y[0][k]) Z=k+1;
				}
				if(errorLog.IsError()) {
					errorLog.Add("FindZero2DConst");
					return ConvError;
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

		// Solve for new X using 1st order Taylor series approx.
		X0[0]=X[0][0]; X0[1]=X[0][1];
		X[0][1]=-(c[0]-a[0][0]/a[1][0]*c[1])/(a[0][1]-a[0][0]/a[1][0]*a[1][1]) ;
		X[0][0]=-(a[0][1]*X[0][1]+c[0])/a[0][0] ;

		// Test for constraint violation.
		Z=0 ;
		while(G(X[0],Y[0],P)) {
			if(errorLog.IsError()) {
				errorLog.Add("FindZero2DConst");
				return ConvError;
			}
			Z++ ;
			if(Z>10) {
				errorLog.Add("FindZero2DConst","Can not find feasible space");
				for(int i=0;i<n;i++) Xg[i]=X[0][i];
				return -1 ;
			}
			X[0][0]=X[0][0]-(X[0][0]-X0[0])/2;
			X[0][1]=X[0][1]-(X[0][1]-X0[1])/2;
		}

		// Evaluate new X.
		F(X[0],Y[0],P);
		ConvError=sqrt(Y[0][0]*Y[0][0]+Y[0][1]*Y[0][1]);
		if(errorLog.IsError()) {
			errorLog.Add("FindZero2DConst");
			return ConvError;
		}
		if(ConvError<=tol) break;
	}
	if(i>=Imax) {
		errorLog.Add("FindZero2DConst","Max iteration error");
		for(int i=0;i<n;i++) Xg[i]=X[0][i];
		return -1 ;
	}
	for(i=0;i<n;i++) Xg[i]=X[0][i];

	return ConvError;
}

/********************************************************************
Solve two simultaneous equations using a two-dimensional bisection
method.  This technique is stable, but extremely inefficient.
********************************************************************/
double Bisection2D(double *Xlo,double *Xhi,double *X,
	double(*F)(double*,double*,void*,int*),double thres,
	void *P,int *status) {

	double Xh,Xl,Y[2],dx[2],CE;
	int j=0;
	const int Jmax=500 ;

	j=0;
	do {
		X[0]=(Xhi[0]+Xlo[0])/2 ;
		Xh=Xhi[1]; Xl=Xlo[1];
		do {
			j++ ;
			X[1]=(Xh+Xl)/2 ;
			F(X,Y,P,status) ;
			if(*status) {printf("Bisection2D: \n"); return 0;}
			if(j>Jmax) {
				printf("Bisection2D: Max iteration error.\n") ;
				*status=1 ;
				CE = sqrt(Y[0]*Y[0]+Y[1]*Y[1]);
				return CE ;
			}
			if(Y[1]<0) Xh=X[1]; else Xl=X[1];
			dx[1]=2*(Xh-Xl)/(Xh+Xl) ;
		} while(fabs(Y[1])>thres) ;
		// } while(fabs(Y[1])>thres && fabs(dx[1])>1e-2*thres) ;
		if(Y[0]<0) Xhi[0]=X[0]; else Xlo[0]=X[0];
		dx[0]=2*(Xhi[0]-Xlo[0])/(Xhi[0]+Xlo[0]) ;
		//	printf("%le\t%le\n",dx[0],Y[0]) ;
	} while(fabs(Y[0])>thres) ;
	// } while(fabs(Y[0])>thres && fabs(dx[0])>1e-2*thres) ;

	CE = sqrt(Y[0]*Y[0]+Y[1]*Y[1]);
	return CE;
}

/********************************************************************
Solves two simultaneous equations using Newton-Raphson method.
Uses a golden section minimization between new and old guesses
if the convergence error increases.
********************************************************************/
int FindZero2DGolden(double *Xg,int(*F)(double*,double*,void*),double tol,void *P)
{
	static double ConvError, X0[3], X[4][3];
	static double Y[4][3]={{0,0,0},{0,0,0},{0,0,0},{0,0,0}};
	static double dx[3],a[3][3],c[3];
	static double b=1e-4;//1e-4
	static int i,j,k,l,Z;
	static double XG[3],YG[3],x0,x1,x2,X1,X2,f1,f2,D,L,dL;
	const int n=2,Imax=10;//20

	// Determine dx.
	for(int i=0;i<n;i++) {
		X[0][i]=Xg[i] ;
		dx[i]=b*fabs(Xg[i]) ;
		if(dx[i]==0) {
			errorLog.Add("FindZero2DGolden","dx=0");
			return 1;
		}
	}

	// Evaluate F at guess point.
	F(X[0],Y[0],P);
	ConvError=CError(Y[0]);
	if(errorLog.IsError()) {
		errorLog.Add("FindZero2DGolden","Evaluating initial guess");
		return 2;
	}
	CEo=ConvError; for(l=0;l<n;l++) Xg[l]=X[0][l];
	if(CEo<tol) return 0;

	for(int i=1;i<=Imax;i++) {
		// Function evaluations.
		for(j=1;j<=n;j++) {
			for(k=0;k<n;k++) X[j][k]=X[0][k] ;
			dx[j-1]=b*fabs(X[0][j-1]); X[j][j-1]+=dx[j-1] ;
			F(X[j],Y[j],P);
			ConvError=CError(Y[j]);
			if(errorLog.IsError()) {
				for(k=0;k<n;k++) X[j][k]=X[0][k] ;
				dx[j-1]=-b*fabs(X[0][j-1]); X[j][j-1]+=dx[j-1] ;
				errorLog.ClearError("FindZero2DGolden");
				F(X[j],Y[j],P);
				ConvError=CError(Y[j]);
				if(errorLog.IsError()) {
					errorLog.ClearError("FindZero2DGolden","status<>0 finding dF/dx");
					return 0;
				}
			}
			if(ConvError<CEo) {CEo=ConvError; for(l=0;l<n;l++) Xg[l]=X[j][l];}
			if (CEo<=tol) return 0;

			// Test if dx is big enough to make a change in Y.
			// If not, then increas dx.
			Z=0; for(k=0;k<n;k++) if(Y[j][k]==Y[0][k]) Z=k+1 ;
			while(Z) {
				dx[j-1]*=2 ;
				if(fabs(dx[j-1]/X[0][j-1])>0.05) {
					for(k=0;k<n;k++) if(Y[j][k]==Y[0][k]) Y[j][k]=Y[0][k]+1e-20;
					Z=0 ;
				} else {
					X[j][j-1]=X[0][j-1]+dx[j-1];
					F(X[j],Y[j],P);
					ConvError=CError(Y[j]);
					if(errorLog.IsError()) {
						errorLog.ClearError("FindZero2DGolden","status<>0 increasing dx");
						return 0;
					}
					if(ConvError<CEo) {CEo=ConvError; for(l=0;l<n;l++) Xg[l]=X[j][l];}
					if (CEo<=tol) return 0;
					Z=0; for(k=0;k<n;k++) if(Y[j][k]==Y[0][k]) Z=k+1 ;
				}
			}
		}

		// Calculate partial derivatives.
		for(int j=0;j<n;j++) {
			c[j]=Y[0][j] ;
			for(k=0;k<n;k++) {
				a[j][k]=(Y[k+1][j]-Y[0][j])/dx[k] ;
				c[j]-=a[j][k]*X[0][k];
			}
		}
		for(l=0;l<n;l++) X0[l]=X[0][l];

		// Solve for new X using 1st order Taylor series approx.
		X[0][1]=-(c[0]-a[0][0]/a[1][0]*c[1])/(a[0][1]-a[0][0]/a[1][0]*a[1][1]);
		X[0][0]=-(a[0][1]*X[0][1]+c[0])/a[0][0];

		// Test for possible increase in error.
		// If so, then find minimum using golden section method.
		F(X[0],Y[0],P);
		ConvError=CError(Y[0]);
		Z=0;
		while(errorLog.IsError()) {
			Z++ ;
			if(Z>30) {//10--30B.S. change
				errorLog.ClearError("FindZero2DGolden","Backing off from new guess");//B.S. sometimes, wrong exit
				return 0;
			}
			for(int j=0;j<n;j++) X[0][j]=X[0][j]-(X[0][j]-X0[j])/2;
			errorLog.ClearError("FindZero2DNew","Backing off from new step");
			F(X[0],Y[0],P);
			ConvError=CError(Y[0]);
		}

		if(ConvError<CEo) {
			CEo=ConvError;
			for(l=0;l<n;l++) Xg[l]=X[0][l];
			if(CEo<=tol) return 0;//B.S. new
		} else {
			D = 1-2/(1+sqrt(5.0));
			X2 = 1; X1=0;
			L = X2-X1;
			x1 = X1+D*L;
			for(int j=0;j<n;j++) XG[j]=X[0][j]-x1*(X[0][j]-X0[j]);
			f1 = F(XG,YG,P);
			F(XG,YG,P);
			f1 = CError(YG);
			if(errorLog.IsError()) {
				errorLog.ClearError("FindZero2DGolden","Evaluating golden section");
				return 0;
			}
			if(f1<CEo) {CEo=f1; for(l=0;l<n;l++) Xg[l]=XG[l];}
			if(CEo<=tol) return 0;
			x2 = X2-D*L;
			for(int j=0;j<n;j++) XG[j]=X[0][j]-x2*(X[0][j]-X0[j]);
			F(XG,YG,P);
			f2 = CError(YG);
			if(errorLog.IsError()) {
				errorLog.ClearError("FindZero2DGolden","Evaluating golden section");
				return 0;
			}
			if(f2<CEo) {CEo=f2; for(l=0;l<n;l++) Xg[l]=XG[l];}
			if(CEo<=tol) return 0;
			dL=2*L/(X1+X2) ;
			while (fabs(dL)>=0.01) {//B.S.0.01--0.001
				if (f1>f2) {
					X1=x1;
					x1=x2;
					f1=f2;
					L=X2-X1;
					x2=X2-D*L;
					dL=2*L/(X1+X2) ;
					for(int j=0;j<n;j++) XG[j]=X[0][j]-x2*(X[0][j]-X0[j]);
					F(XG,YG,P);
					f2=CError(YG);
					if(errorLog.IsError()) {
						errorLog.ClearError("FindZero2DGolden","Evaluating golden section");//wrong exit
						return 0;
					}
					if(f2<CEo) {CEo=f2; for(l=0;l<n;l++) Xg[l]=XG[l];}
					if(CEo<=tol) 0;
				} else if (f1<f2) {
					X2=x2;
					x2=x1;
					f2=f1;
					L=X2-X1;
					x1=X1+D*L;
					dL=2*L/(X1+X2) ;
					for(int j=0;j<n;j++) XG[j]=X[0][j]-x1*(X[0][j]-X0[j]);
					F(XG,YG,P);
					f1=CError(YG);
					if(errorLog.IsError()) {
						errorLog.ClearError("FindZero2DGolden","Evaluating golden section");
						return 0;
					}
					if(f1<CEo) {CEo=f1; for(l=0;l<n;l++) Xg[l]=XG[l];}
					if(CEo<=tol) 0;
				} else {
					X1=x1;
					X2=x2;
					L=X2-X1;
					dL=2*L/(X1+X2) ;
					x1=X1+D*L;
					for(int j=0;j<n;j++) XG[j]=X[0][j]-x1*(X[0][j]-X0[j]);
					F(XG,YG,P);
					f1=CError(YG);
					if(errorLog.IsError()) {
						errorLog.ClearError("FindZero2DGolden","Evaluating golden section");
						return 0;
					}
					if(f1<CEo) {CEo=f1; for(l=0;l<n;l++) Xg[l]=XG[l];}
					if(CEo<=tol) 0;
					x2 = X2-D*L;
					for(int j=0;j<n;j++) XG[j]=X[0][j]-x2*(X[0][j]-X0[j]);
					F(XG,YG,P);
					f2 = CError(YG);
					if(errorLog.IsError()) {
						errorLog.ClearError("FindZero2DGolden","Evaluating golden section");
						return 0;
					}
					if(f2<CEo) {CEo=f2; for(l=0;l<n;l++) Xg[l]=XG[l];}
					if(CEo<=tol) 0;
				}
			}
			x0=X1+L/2;
			for(int j=0;j<n;j++) X[0][j]-=x0*(X[0][j]-X0[j]);
			F(X[0],Y[0],P);
			ConvError=CError(Y[0]);
			if(errorLog.IsError()) {
				errorLog.ClearError("FindZero2DGolden","Evaluating golden section");
				return 0;
			}
			if(ConvError<CEo) {CEo=ConvError; for(l=0;l<n;l++) Xg[l]=X[0][l];}
			if(CEo<=tol) 0;
		}
	}
	errorLog.Add("FindZero2DGolden","Max iteration");
	return 5;
}

/********************************************************************
Solves two simultaneous equations using Newton-Raphson method.
Tests to see if guesses are in feasible space and automatically
backs off from non-feasible guesses.
********************************************************************/
double FindZero2DTest(double *Xg,double(*F)(double*,double*,void*),double tol,void *P)
{
	double CEo,ConvError,X0[3],X[4][3];
	double Y[4][3]={{0,0,0},{0,0,0},{0,0,0},{0,0,0}};
	double dx[3],a[3][3],c[3];
	double b0=1e-3,b=1e-4;
	const int n=2,Imax=20;
	int i,j,k,l,Z;

	// Determine dx.
	for(int i=0;i<n;i++) {
		X[0][i]=Xg[i];
		dx[i]=b*fabs(Xg[i]);
		if(dx[i]==0) {
			errorLog.Add("FindZero2DTest","dx=0");
			return -1 ;
		}
	}

	// Evaluate F at guess point.
	ConvError = F(X[0],Y[0],P);
	if(errorLog.IsError()) {
		errorLog.Add("FindZero2DTest","Evaluating initial guess");
		return -1;
	}
	CEo=ConvError; for(l=0;l<n;l++) Xg[l]=X[0][l];
	if(CEo<tol) return CEo ;

	for(int i=0;i<Imax;i++) {
		b=b0/pow(2,(double)i);
		// Function evaluations.
		for(j=1;j<=n;j++) {

			for(k=0;k<n;k++) X[j][k]=X[0][k];
			dx[j-1]=b*fabs(X[0][j-1]); X[j][j-1]+=dx[j-1];

			ConvError=F(X[j],Y[j],P);
			if(errorLog.IsError()) {
				for(k=0;k<n;k++) X[j][k]=X[0][k];
				dx[j-1]=-b*fabs(X[0][j-1]); X[j][j-1]+=dx[j-1];
				errorLog.ClearError("FindZero2DTest");
				ConvError=F(X[j],Y[j],P);
				if(errorLog.IsError()) {
					errorLog.ClearError("FindZero2DTest","status<>0 finding dF/dx");
					return CEo;
				}
			}
			if(ConvError<CEo) {CEo=ConvError; for(l=0;l<n;l++) Xg[l]=X[j][l];}
			if (CEo<=tol) return CEo;

			// Test if dx is big enough to make a change in Y.
			// If not, then increas dx.
			Z=0; for(k=0;k<n;k++) if(Y[j][k]==Y[0][k]) Z=k+1 ;
			while(Z) {
				dx[j-1]*=2 ;
				if(fabs(dx[j-1]/X[0][j-1])>0.05) {
					for(k=0;k<n;k++) if(Y[j][k]==Y[0][k]) Y[j][k]=Y[0][k]+1e-20;
					Z=0;
				} else {
					X[j][j-1] = X[0][j-1]+dx[j-1];
					ConvError = F(X[j],Y[j],P);
					if(errorLog.IsError()) {
						errorLog.ClearError("FindZero2DTest","status<>0 increasing dx");
						return CEo;
					}
					if(ConvError<CEo) {CEo=ConvError; for(l=0;l<n;l++) Xg[l]=X[j][l];}
					if (CEo<=tol) return CEo;
					Z=0; for(k=0;k<n;k++) if(Y[j][k]==Y[0][k]) Z=k+1 ;
				}
			}
		}

		// Calculate partial derivatives.
		for(int j=0;j<n;j++) {
			c[j]=Y[0][j] ;
			for(k=0;k<n;k++) {
				a[j][k]=(Y[k+1][j]-Y[0][j])/dx[k] ;
				c[j]-=a[j][k]*X[0][k];
			}
		}

		for(l=0;l<n;l++) X0[l]=X[0][l];

		// Solve for new X using 1st order Taylor series approx.
		X[0][1]=-(c[0]-a[0][0]/a[1][0]*c[1])/(a[0][1]-a[0][0]/a[1][0]*a[1][1]);
		X[0][0]=-(a[0][1]*X[0][1]+c[0])/a[0][0] ;

		ConvError=F(X[0],Y[0],P);
		Z=0;
		while(errorLog.IsError()) {
			Z++ ;
			if(Z>10) {
				errorLog.ClearError("FindZero2DTest","Backing off from new guess");
				return CEo;
			}
			for(int j=0;j<n;j++) X[0][j]=X[0][j]-(X[0][j]-X0[j])/2;
			errorLog.ClearError("FindZero2DTest","Backing off from new step");
			ConvError = F(X[0],Y[0],P);
		}
		if(ConvError<CEo) {CEo=ConvError; for(l=0;l<n;l++) Xg[l]=X[0][l];}
		printf("CEo=%lf\n",CEo);
		if (CEo<=tol) return CEo;
	}
	if(i>=Imax) {
		errorLog.Add("FindZero2DTest","Max iteration");
	}

	return CEo ;
}



