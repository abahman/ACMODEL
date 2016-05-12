#include "stdafx.h"

#include <math.h>
//#include "acmodel.h"
#include "compu.h"
#include "dataelem.h"
#include "errorlog.h"
#include "define.h"
#include "spline1d.h"
#include "parse.h"
#include "outelem.h"
#include "getdata.h"
#include "conversion.h"

//B.S., new head files
//#include <fstream.h>
//#include <iostream.h>
#include "mass.h"
#include "expdev.h"
#include "r22.h"
#include "numeric.h"
//----------------------------------------B.S.

extern ErrorLog errorLog;
extern DataElement gData;                // model parameters
extern OutputDataElement gOutData;

extern RefPropTPthT reftpvthT;
extern RefPropTPthT reftplthT;

extern double CEo;//important

//B.S subfunctions--------------------------------------------
double InitialGuess_2D(double*, InVars*);//provide the initial guesses for two-dimension search
double InitialGuess_3D(double*, InVars*);//provide the initial guesses for three-dimension search
//---------------------------------------------------B.S.

//****************************************************************************
//****************************************************************************
//****************************************************************************

extern int compInit;
extern int compCCLowInit;
extern int compCCHighInit;
extern int condInit;
extern int evapInit;
extern int subcoolInit;
extern int expDevN;

// Intializes the program states.
//
// state[] is used to enable this program to run in a single threaded GUI like old
// Windows.  It could run in the GUI's IdleAction function where the computational engine
// needs to return frequently to give the UI time to interact.
//
// When the components start up, they read an input file containing there quasi-static
// variables.  This call re-iniitalizes them so they will read the file data the next time
// the component function calls are made.
void InitializeComponents(void) 
{
	// Reinitialize components
	compInit = 1;
	compCCLowInit = 1;
	compCCHighInit = 1;
	condInit = 1;
	evapInit = 1;
	subcoolInit = 1;
	expDevN = 0;
}

//****************************************************************************
//****************************************************************************
//****************************************************************************

//************************************************************************
// Main function for ACMODEL.  This function solves multiple operating
// points in batch mode.  All data input/output is contained here.
//************************************************************************

ENGINE::ENGINE()
{
	char scanstr[2048];
	int i=0,j=0;
	int endSign;
	error=0;
	strncpy(filename,"./InputDoc/define.in",120);
	FILE *fp = fopen(filename,"r");
	if(fp==NULL) {
		errorLog.Add("ENGINE::ENGINE","File not found");
		error=1;
		return;}

	fgets(scanstr,2048,fp);
	systemType=(int)GetDouble(fp);//get system type, 0--validate the models of components with measured data,1--short-tube model, 2--TxV model
	ACTUNE = (int)GetDouble(fp);  // to see if it needs to run the tuning procedure
	I.REV_CAL=(int)GetDouble(fp);   // define the iteration loop
	REF=(int)GetDouble(fp);    //what is the refrigerant
	I.SPEC_SUB=(int)GetDouble(fp);  //if specify the subcooling degree
	I.PISTON = (int)GetDouble(fp);//specify the  expansion device

	endSign=(int)GetDouble(fp);

	fclose(fp);
	if(endSign!=-256)
	{
	errorLog.Add("ENGINE::ENGINE","File Input Wrong");
	error=1;
	return;}
}

ENGINE::~ENGINE()
{
}

int ENGINE::compu(void)
{
	char scanstr[2048],*ptr;
	double Y[3];
	TXP TXP_prop={0,0,0};

	//gData.WriteParams();

	// Get parameters from file.
	strcpy(I.filename,"./InputDoc/acmodel.dat");

	// Get system type
	/*systemType = GetComponentType(I.filename,SYSTEM);
	if(!systemType) {
		errorLog.Add("ENGINE::compu","Invalid system type");
		return 0;
	}*///Rossi input

	time(&tbegin);
	Nrun0=Nrun1=Crun0=Crun1=Cerror=Cnoconv=0;
	Trun0=Trun1=0.0;

	// Count the number of operating conditions to solve for that are
	// 1) after a new initial guess is specified and (Nrun0)
	// 2) not after a new initial guess is specified. (Nrun1)
	// The idea is that these two conditions have very different solution
	// times.  Therefore, the program keeps track of the average time for
	// each type of solution to estiamte the time required to complete the
	// calculation.
	{
		//gData.WriteInput(0);

		FILE *out = NULL;
		out = fopen("./InputDoc/master.in","r");
		if(out==NULL) {
			errorLog.Add("No input file");
			return 1;
		}

		char slast = '1';
		while(1) {
			char str[128];
			fgets(str,127,out);
			if(feof(out)) break;
			if(str[0]=='3'||str[0]=='4') {
				if(slast=='1'||slast=='2') Nrun0++; else Nrun1++;
			}
			slast = str[0];
		};
		fclose(out);
	}

	// Get an initial guess from the input file
	if(systemType){
	in = fopen("./InputDoc/master.in","r");
	if(in==NULL) {
		errorLog.Add("COMPU","master.in File not found");
		error=1;
		return -1;}

	states=NULL;
	}
	else {
	states = fopen("./InputDoc/states.in","r");//input measured data at eight points to validate the components model
	if(states==NULL) {
		errorLog.Add("COMPU","states.in File not found");
		error=1;
		return -1;}
	fgets(scanstr,2048,states);

	in = fopen("./InputDoc/master.in","r");
	if(in==NULL) {
		errorLog.Add("COMPU","master.in File not found");
		error=1;
		return -1;}
	}


	switch(systemType) {
		case 0:
			{
			fscanf(states,"%d",&i);
			if(feof(states)) {
			fclose(states); 
			fclose(in); 
			return 2;}

			fgets(scanstr,2048,states);
			ptr = scanstr;
			sget(&ptr,&I.Tai);
			sget(&ptr,&I.TPi.T);
			sget(&ptr,&I.TPi.P);
			sget(&ptr,&I.Goa);
			sget(&ptr,&I.Gra);
			sget(&ptr,&I.Charge);
			sget(&ptr,&I.Tsup);
			sget(&ptr,&I.Tsub);
			sget(&ptr,&I.mr);
			sget(&ptr,&I.P[1]);
			sget(&ptr,&I.P[2]);
			sget(&ptr,&I.P[3]);
			sget(&ptr,&I.P[4]);
			sget(&ptr,&I.P[5]);
			sget(&ptr,&I.P[6]);
			sget(&ptr,&I.P[7]);
			sget(&ptr,&I.P[8]);
			sget(&ptr,&I.T[1]);
			sget(&ptr,&I.T[2]);
			sget(&ptr,&I.T[3]);
			sget(&ptr,&I.T[4]);
			sget(&ptr,&I.T[5]);
			sget(&ptr,&I.T[6]);
			sget(&ptr,&I.T[7]);
			sget(&ptr,&I.T[8]);
			sget(&ptr,&I.X[0]);
			sget(&ptr,&I.X[1]);
			sget(&ptr,&I.X[2]);
			sget(&ptr,&I.X[3]);
			sget(&ptr,&I.CompT);
			sget(&ptr,&I.Cadj);
			sget(&ptr,&I.Eadj);
			sget(&ptr,&I.ComP);
			}
		break;

		case 1:
		case 2: {
			int i;
			fscanf(in,"%d",&i);
			while(i==1) {
				fgets(scanstr,2048,in);
				sscanf(scanstr,"%lf%lf%lf",&Xg[0],&Xg[1],&Xg[2]);
				fscanf(in,"%d",&i);

				printf("New initial guess.\n");
				errorLog.ClearError("ENGINE::compu");
			}
			fgets(scanstr,2048,in);
			ptr = scanstr;
			sget(&ptr,&PRT);
			sget(&ptr,&I.Tai);
			sget(&ptr,&I.TPi.T);
			sget(&ptr,&I.TPi.P);
			sget(&ptr,&I.Goa);
			sget(&ptr,&I.Gra);
			sget(&ptr,&I.Charge);
			//----------------------------
			sget(&ptr,&I.Tsup);//B.S. new
			sget(&ptr,&I.Tsub);//B.S. new
			sget(&ptr,&I.CompT);//B.S. new
			//-----------------------------
			sget(&ptr,&I.Dadj);
			sget(&ptr,&I.nvAdj);
			
			printf(">>>>>> New input Tai=%.1fC (%.1lfF).\n",I.Tai,CtoF(I.Tai));

			gOutData.auxDataCnt = 0;
			while(gOutData.auxDataCnt<100 && sget(&ptr,&gOutData.auxData[gOutData.auxDataCnt])) {gOutData.auxDataCnt++;}
		}
		break;
		case 3: {
			int i;
			fscanf(in,"%d",&i);
			while(i==1) {
				double pi = acos(-1.0);
				double Pevap,Pint,Pcond;
				fgets(scanstr,2048,in);
				sscanf(scanstr,"%lf%lf%lf",&Pcond,&Pint,&Pevap);
				// atan() function maps -inf to +inf to +pi/2 to -pi/2.  Therefore,
				// it is used to allow the numerical method to propose valid guesses
				// for -inf to +inf.  The contraint is Pevap < Pint < Pcond.
				Xg[0] = tan(pi*(Pevap-PMINth)/(PMAXth-PMINth)-pi/2);
				Xg[1] = tan(pi*(Pint-Pevap)/(PMAXth-Pevap)-pi/2);
				Xg[2] = tan(pi*(Pcond-Pint)/(PMAXth-Pint)-pi/2);
				fscanf(in,"%d",&i);
				printf("New initial guess from input file.\n");
			}
			fgets(scanstr,2048,in);
			ptr = scanstr;
			sget(&ptr,&PRT);
			sget(&ptr,&I.Tai);
			sget(&ptr,&I.TPi.T);
			sget(&ptr,&I.Goa);
			sget(&ptr,&I.Gra);
			sget(&ptr,&I.Charge);
			sget(&ptr,&I.Tsup);//B.S. new
			sget(&ptr,&I.Tsub);//B.S. new

			printf("New input request from input file.\n");

			gOutData.auxDataCnt = 0;
			while(gOutData.auxDataCnt<100 && sget(&ptr,&gOutData.auxData[gOutData.auxDataCnt])) {gOutData.auxDataCnt++;}
		}
		break;
	}

//	gOutData.StartMasterInputFile();

	// Begin main loop
	time(&tcur);
	ConvTol=1e-6;
	I.ConvTol = ConvTol;
	do  {
		for(i=0;i<3;i++) X[i]=Xg[i];
		if(I.Tai==0) I.Tai=1e-10;
		if(I.TPi.T==0) I.TPi.T=1e-10;

		gOutData.ClearAll();
		errorLog.ClearError("begin");

		switch(systemType) {
			case 0://validate or tune the components using measured data			
			I.EvaSign =0;//use the detailed evaporator model
			I.CondSign =0;//use the detailed condenser model
			MainFunc_Validating(&I);
			//MainFunc_Tuning(&I);
			X[0]=0;X[1]=0;X[2]=0;			
			break;

			case 1:// short orifice			
				I.EvaSign =0;//use the detailed evaporator model
				I.CondSign =0;//use the detailed condenser model
				
				if(I.EvaSign)
				{
				if(I.SPEC_SUB) BuildSystem_SubC_3D(X, &I);
				else BuildSystem_Charg(X, &I);
				}
				else{				
				FindZero3DGolden(X,MainFunc1,ConvTol,&I);
				errorLog.ClearError("FindZero3DGolden");
				}
				
				MainFunc1(X,Y,&I);

				if(errorLog.IsError()||sqrt((Y[0]*Y[0]+Y[1]*Y[1]+Y[2]*Y[2])/3)>5e-6)
				{
				errorLog.Add("fixed orifice system","not converged");
				}

			
				break;
				
			case 2:   // constant superheat TxV
				
				//intialized
				X[0]=X[1];//input the evaporating pressure from master.in.
				X[1]=X[2];//input the condensing pressure from master.in			
				I.EvaSign =0;//use the detailed evaporator model for analysis
				I.CondSign =0;//use the detailed condenser model for analysis
				I.Tsup =I.Tsup;
				
				if(I.EvaSign)
				{
				if(I.SPEC_SUB) BuildSystem_SubC_2D(X, &I);
				else BuildSystem_Charg(X, &I);
				}
				else{
				FindZero2DGolden(X,MainFunc2,ConvTol,&I);

				errorLog.ClearError("FindZero2DGolden");
				}
			
				
				MainFunc2(X,Y,&I);

				if(errorLog.IsError()||sqrt((Y[0]*Y[0]+Y[1]*Y[1])/3)>5e-6)
				{
				errorLog.Add("TXV system","not converged");
				}

				X[2]=X[1];//output the condensing pressure
				X[1]=X[0];//output the evaporating pressure

				{
				TXP_prop.P=X[1];
				TXP_prop.X=1;
				TXP TXP1 = toTXP(PropertyTXPth(TSAT,TXP_prop)+I.Tsup,1,X[1]);
				X[0] = PropertyTXPth(ENTH,TXP1);//output the successful suction enthalpy
				}
				break;

			case 3:   // NIST chiller
			//	FindZero3DGuess(X,MainFunc3,ConvTol,&I);
				break;
		}

		//calculate the time used
		tlast=tcur;
		time(&tcur);
		if(fix) {
			Crun1++;
			Trun1+=difftime(tcur,tlast);
			gOutData.Time_dev = difftime(tcur,tlast);//B.S. record the time interval for the current run
		} else {
			Crun0++;
			Trun0+=difftime(tcur,tlast);
			gOutData.Time_dev = difftime(tcur,tlast);//B.S. record the time interval for the current run
		}

		gOutData.isError = errorLog.IsError() ? 1 : 0;
		gOutData.isConv = ConvError>ConvTol ? 0 : 1;

		// Print a report even if error or no convergence.
		// Put before gOutData.output because it clears the data.
		//gOutData.AppendOutputToFile();

		if(errorLog.IsError()) {

			// Increment error counter
			Cerror++;
			printf("Solution error!\n");

			gOutData.output('f');//B.S., print the states at different points under this condition, 'f' represents failure 

			// Quit if can not solve for initial guess
			if(!fix) break;

		} else {
			// Record all runs that do not converge in "noconv".

			if(ConvError>ConvTol) {
				Cnoconv++;
				gOutData.output('u');//B.S., print the states at different points under this condition, 'u' represents nonconverged 

				printf("Solution did not converge!\n");
			} else {
				for(i=0;i<3;i++) Xg[i]=X[i];

				// Only create inout file to expand in another dimension if the solution converged
				// in this dimension.
				// Must be before gOutData.output() - it clears the output data
				//gOutData.AppendMasterInputFile();

				printf("X=(%lf,%lf,%lf)",Xg[0],Xg[1],Xg[2]);// B.S., print out the initial guess for this point
				gOutData.output('s');//B.S., print the states at different points under this condition, 's' represents success

				printf("Solution successful!\n");
			}
		}

		// update status, i.e input the new working condition
		
		if(systemType){
		fscanf(in,"%d",&i);
		if(feof(in)) {fclose(in); return 2;}
		}
		else
		{
		fscanf(states,"%d",&i);
		if(feof(states)) {
		fclose(states); 
		fclose(in); 
		return 2;}
		}

		switch(systemType) {
		case 0: 
			{
			fgets(scanstr,2048,states);
			ptr = scanstr;
			sget(&ptr,&I.Tai);
			sget(&ptr,&I.TPi.T);
			sget(&ptr,&I.TPi.P);
			sget(&ptr,&I.Goa);
			sget(&ptr,&I.Gra);
			sget(&ptr,&I.Charge);
			sget(&ptr,&I.Tsup);
			sget(&ptr,&I.Tsub);
			sget(&ptr,&I.mr);
			sget(&ptr,&I.P[1]);
			sget(&ptr,&I.P[2]);
			sget(&ptr,&I.P[3]);
			sget(&ptr,&I.P[4]);
			sget(&ptr,&I.P[5]);
			sget(&ptr,&I.P[6]);
			sget(&ptr,&I.P[7]);
			sget(&ptr,&I.P[8]);
			sget(&ptr,&I.T[1]);
			sget(&ptr,&I.T[2]);
			sget(&ptr,&I.T[3]);
			sget(&ptr,&I.T[4]);
			sget(&ptr,&I.T[5]);
			sget(&ptr,&I.T[6]);
			sget(&ptr,&I.T[7]);
			sget(&ptr,&I.T[8]);
			sget(&ptr,&I.X[0]);
			sget(&ptr,&I.X[1]);
			sget(&ptr,&I.X[2]);
			sget(&ptr,&I.X[3]);
			sget(&ptr,&I.CompT);
			sget(&ptr,&I.Cadj);
			sget(&ptr,&I.Eadj);
			sget(&ptr,&I.ComP);
			}
		break;

		case 1:
		case 2:
			while(i==1) {
				fgets(scanstr,2048,in);
				sscanf(scanstr,"%lf%lf%lf",&Xg[0],&Xg[1],&Xg[2]);
				fscanf(in,"%d",&i);

				printf("New initial guess.\n");
				errorLog.ClearError("ENGINE::compu");
			}
			fgets(scanstr,2048,in);
			ptr = scanstr;
			sget(&ptr,&PRT);
			sget(&ptr,&I.Tai);
			sget(&ptr,&I.TPi.T);
			sget(&ptr,&I.TPi.P);
			sget(&ptr,&I.Goa);
			sget(&ptr,&I.Gra);
			sget(&ptr,&I.Charge);
			sget(&ptr,&I.Tsup);//B.S. new
			sget(&ptr,&I.Tsub);//B.S. new
			sget(&ptr,&I.CompT);//B.S. new
			sget(&ptr,&I.Dadj);
			sget(&ptr,&I.nvAdj);

			printf(">>>>>> New input Tai=%.1fC (%.1lfF).\n",I.Tai,CtoF(I.Tai));

			gOutData.auxDataCnt = 0;
			while(gOutData.auxDataCnt<100 && sget(&ptr,&gOutData.auxData[gOutData.auxDataCnt])) {gOutData.auxDataCnt++;}
			break;
		case 3:
			while(i==1) {
				double pi = acos(-1.0);
				double Pevap,Pint,Pcond;
				fgets(scanstr,2048,in);
				sscanf(scanstr,"%lf%lf%lf",&Pcond,&Pint,&Pevap);
				// atan() function maps -inf to +inf to +pi/2 to -pi/2.  Therefore,
				// it is used to allow the numerical method to propose valid guesses
				// for -inf to +inf.  The contraint is Pevap < Pint < Pcond.
				Xg[0] = tan(pi*(Pevap-PMINth)/(PMAXth-PMINth)-pi/2);
				Xg[1] = tan(pi*(Pint-Pevap)/(PMAXth-Pevap)-pi/2);
				Xg[2] = tan(pi*(Pcond-Pint)/(PMAXth-Pint)-pi/2);
				fscanf(in,"%d",&i);
				printf("New initial guess from input file.\n");
			}
			fgets(scanstr,2048,in);
			ptr = scanstr;
			sget(&ptr,&PRT);
			sget(&ptr,&I.Tai);
			sget(&ptr,&I.TPi.T);
			sget(&ptr,&I.Goa);
			sget(&ptr,&I.Gra);
			sget(&ptr,&I.Charge);
			sget(&ptr,&I.Tsup);//B.S. new
			sget(&ptr,&I.Tsub);//B.S. new

			printf("New input request from input file.\n");

			gOutData.auxDataCnt = 0;
			while(gOutData.auxDataCnt<100 && sget(&ptr,&gOutData.auxData[gOutData.auxDataCnt])) {gOutData.auxDataCnt++;}
			break;
		}
	} while(1);
	return 0;
}

//*****************************************************************************
//*****************************************************************************
//*****************************************************************************

//B.S. provide the reasonable initial guess for 2-D search with reversed evaporator construct
double InitialGuess_2D(double *X, InVars* P)
{
	double T_diff=0;//temperature deviation between the suction temperature and the indoor temperature
	double Y[3];//residuals
	int sign =1;
	int Num =0;
	double dev=0;
	double X_bak[2];//backup the evaporating pressure guesses
	P->EvaSign =0;//use the detailed evaporator model for analysis
	P->CondSign = 0;//use the detailed condenser model for analysis
	static int start=1,end=1;
	
	
InitialGuess_2D_start:
	if(start)
	{
	//temperature difference between the suction temperture and the indoor temperture should change with the superheat degree
	T_diff =15+(P->Tsup - 0.1)/(30-0.1)*(0.02-15);
	if(T_diff<0)
	{
	X[0] = reftpvthT.Psat(P->TPi.T-P->Tsup-0.05);
	}
	else
	{
	X[0] = reftpvthT.Psat(P->TPi.T-P->Tsup-T_diff);//initial guess for the evaporating pressure
	}
	X[1] = reftpvthT.Psat(P->Tai+15);
	start=0;
	}
	
	if(errorLog.IsError()) {
	errorLog.Add("InitialGuess_2D","T_diff");
	return 0;
	}
							
	sign = 1;//signal variable
	Num = 0;//variable for controlling the number of iteration
	dev=0.5;//offset value to adjust the evaporating pressure
	X_bak[1]=1;
	X_bak[0]=0;
	
	while(sign)
	{
	
	//if the procedure fails to adjust the evaporting pressure, change the offset value dev=dev/2.
	if(fabs(X_bak[1]-X_bak[0])<1e-5) 
	{
	dev=dev/2;
	X_bak[1]=1;
	X_bak[0]=0;
	}

	while(MainFunc2(X,Y,P)) {if(errorLog.IsError()) break;}//state initialized
	
	if(errorLog.IsError()||P->Evap.LiqL>0)//this means the evaporating pressure is too low
	{
	errorLog.ClearError("InitialGuess_2D","state initialized");
	X_bak[0] = X[0];
	X[0]=X[0]+dev;
	}
	else if((P->Evap.TPL/P->Evap.L_TOT)<0.05)//this means that the evaporating pressure is too high
	{
	X[0]=X[0]-dev;
	X_bak[1] = X[0];
	}
	else sign = 0;
	Num = Num+1;
	if(Num>100) {sign = 0; X[0]=X[0]-dev/2;}//control the iteration number
	}

	if(Num>100) {
	if(end)
	{
	start =1;
	end=0;
	goto InitialGuess_2D_start;
	}
	else
	{
//	errorLog.Add("InitialGuess_2D","fail to get a good initial guess");
	end=1;
	return 0;
	}
	}
	
	end=1;
	return 0;
}

//B.S. provide the initial guesses for 3-D search with reversed evaporator construct

double InitialGuess_3D(double *X, InVars* P)
{
	double Tmax =30.0, Tmin = -29.0;//0.1;//maximum and minimum superheat degree 30.0; 29.0
	double a =-1, b=-1, c=0;//variable for controlling the iteration
	double mr=0;//mass flow rate through the orifice
	double c_min=0, c_max=0;
 	double ConvTol=P->ConvTol;
	double Ones[5]={1,1,1,1,1};
	int Revise=0,NumRev=0;
	TXP TXP_prop={0,0,0};
	dataSpline Data_out;

	int i=0,j=0,k=0;
	
	P->Tsup = (Tmax+Tmin)/2;//set the superheat degree
	int init =1;

	X[0]=X[1];
	X[1]=X[2];

	while(init)
	{			
	InitialGuess_2D(X, P);//provide the initial guesses for 2-D search

	if(errorLog.IsError()) {
	errorLog.Add("InitialGuess_3D","InitialGuess_2D");
	return 0;
	}

	P->EvaSign =0;//use the detailed evaporator model
	P->CondSign =0;//use the detailed condenser model
	
	Revise=1;
	NumRev=0;

   while(Revise&&NumRev<3)
   {
      while(FindZero2DGolden(X,MainFunc2,ConvTol,P)) 
      {if(errorLog.IsError()) break;}//conduct 2-D search
      Revise=0;
      if(errorLog.IsError()) 
      {
         //	const double D = 1-2/(1+sqrt(5));
         errorLog.ClearError("InitialGuess_3D","FindZero2DGolden");
         P->Tsup=P->Tsup+(P->Tsup-Tmin)*0.1;
         Revise=1;
         NumRev=NumRev+1;
      }
   }
	
	errorLog.ClearError("InitialGuess_3D","FindZero2DGolden");

	mr = Short_Tube_R410A("./InputDoc/ACMODEL.dat",P->HP5,P->HP6,Ones);				
//	mr = Short_Tube_FlowDirec("./InputDoc/ACMODEL.dat",P->HP5,P->HP6,Ones);//mass flow rate through the orifice, based on the outputs from the 2-D search
//	mr = TxV_Nonlinear("./InputDoc/ACMODEL.dat",P->HP5,P->HP6,P->Tsup8,Ones);
	c = (mr-P->mr);//residual

	
	if(!errorLog.IsError()&&gOutData.CE<=ConvTol)
	{	
	gOutData.outputTOarray();
	gOutData.Array[0]=c;
	for(j=0;j<58;j++)
	{
	Data_out.Data_points[j][i]=gOutData.Array[j];
	}
	Data_out.Data_points[51][i]=c;

	i=i+1;
	if(i==20) 
	{i=0;
	k=20;//overflow
	}
	}

	if(errorLog.IsError()) {
	Tmin = P->Tsup; 
	a=-1;
	errorLog.ClearError("InitialGuess_3D","Short_Tube");
	}
	else if(c<0){Tmin = P->Tsup;a=P->Tsub;c_min = c;}
	else {Tmax = P->Tsup;b=P->Tsub;c_max = c;}

	
	//adjust the guess of superheat degree
	if(fabs(c)<5e-7){P->Tsup = P->Tsup;}
	else if(a<0||b<0){
	P->Tsup = (Tmax+Tmin)/2;
	}
	else
	{
	P->Tsup = Tmin +(Tmax-Tmin)/(c_max-c_min)*(0-c_min);
	}
			
	if(fabs(Tmax-Tmin)<0.1||fabs(c)<5e-7) init=0;//out from the iteration

	}


	X[2]=X[1];//output the sucessful condensing pressure
	X[1]=X[0];//output the sucessful evaporating pressure

	{
	TXP_prop.P=X[1];
	TXP_prop.X=1;
	if(P->Tsup<0)
	{
	TXP_prop.T=PropertyTXPth(TSAT,TXP_prop);
	X[0]=PropertyTXPth(ENTH,TXP_prop)+PropertyTXPtr(SPEC,TXP_prop)*P->Tsup;
	}
	else
	{
	TXP TXP1 = toTXP(PropertyTXPth(TSAT,TXP_prop)+P->Tsup,1,X[1]);//toTXP(reftplthP.Tsat(X[1])+P->Tsup,1,X[1]);
	X[0] = PropertyTXPth(ENTH,TXP1);//output the successful suction enthalpy
	}
	}
	if(errorLog.IsError()) {
	errorLog.Add("InitialGuess_3D","get X[0]");
	return 0;
	}

	//get the data for spline
	if(k<20) k=i;

	if(k>0)
	{
	gOutData.Array[0]=0;
	outdataSpline(k,57,Data_out,0,gOutData.Array);
	}
	return 0;
}

