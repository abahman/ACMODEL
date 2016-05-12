#include "stdafx.h"

#include <stdio.h>
#include "dataelem.h"
#include "define.h"
#include "parse.h"

void DataElement::Default()
{
	DefaultSystem(1);
}

DataElement::DataElement()
{
	Default();

	FILE *fp = NULL;
	fp = fopen("master.ini","r");
	if(fp) {
		fscanf(fp,"%s",filename);
		fclose(fp);
	}
	Import();
}

DataElement::~DataElement()
{
	WriteINIFile();
}

void DataElement::WriteINIFile()
{
	Export();

	FILE *fp = fopen("master.ini","w");
	if(fp) {
		fprintf(fp,"%s\n",filename);
		fclose(fp);
	}
}

// Read model parameters from *.acm file.
int DataElement::Import(void)
{
	FILE *fp=fopen(filename,"r");
	if(!fp) {
		strcpy(filename,"default.acm");
		return 1;
	}

	char line[256];
	while(fgets(line,256,fp)!=NULL) {
		int i;
		if(sscanf(line,"batch_filename%s",operationData.batchFilename)==1) continue;
		if(sscanf(line,"output_filename%s",operationData.outputFilename)==1) continue;

		if(sscanf(line,"system_type%d",&system_type)==1) continue;
		if(sscanf(line,"system_type%d",&system_type)==1) continue;
		if(strstr(line,"input")!=NULL) {
			if(sscanf(line,"input_type%d",&input_type)==1) continue;
			for(i=0;i<20;i++) {
				char format[128];
				sprintf(format,"input[%2d]%%lf",i);
				if(sscanf(line,format,&input[i])==1) continue;
			}
		}
		if(strstr(line,"comp")!=NULL) {
			if(sscanf(line,"comp_type%d",&comp_type)==1) continue;
			if(sscanf(line,"comp_nprms%d",&comp_nprms)==1) continue;
			for(i=0;i<40;i++) {
				char format[128];
				sprintf(format,"comp[%2d]%%lf",i);
				if(sscanf(line,format,&comp[i])==1) continue;
			}
		}
		if(strstr(line,"cond")!=NULL) {
			if(sscanf(line,"cond_type%d",&cond_type)==1) continue;
			if(sscanf(line,"cond_nprms%d",&cond_nprms)==1) continue;
			for(i=0;i<20;i++) {
				char format[128];
				sprintf(format,"cond[%2d]%%lf",i);
				if(sscanf(line,format,&cond[i])==1) continue;
			}
		}
		if(strstr(line,"evap")!=NULL) {
			if(sscanf(line,"evap_type%d",&evap_type)==1) continue;
			if(sscanf(line,"evap_nprms%d",&evap_nprms)==1) continue;
			for(i=0;i<20;i++) {
				char format[128];
				sprintf(format,"evap[%2d]%%lf",i);
				if(sscanf(line,format,&evap[i])==1) continue;
			}
		}
		if(strstr(line,"subcool")!=NULL) {
			if(sscanf(line,"subcool_type%d",&subcool_type)==1) continue;
			if(sscanf(line,"subcool_nprms%d",&subcool_nprms)==1) continue;
			for(i=0;i<10;i++) {
				char format[128];
				sprintf(format,"subcool[%2d]%%lf",i);
				if(sscanf(line,format,&subcool[i])==1) continue;
			}
		}
		if(strstr(line,"expdev")!=NULL) {
			if(sscanf(line,"expdev_type%d",&expdev_type)==1) continue;
			if(sscanf(line,"expdev_nprms%d",&expdev_nprms)==1) continue;
			for(i=0;i<10;i++) {
				char format[128];
				sprintf(format,"expdev[%2d]%%lf",i);
				if(sscanf(line,format,&expdev[i])==1) continue;
			}
		}
		if(strstr(line,"restrict")!=NULL) {
			if(sscanf(line,"restrict_type%d",&restrict_type)==1) continue;
			if(sscanf(line,"restrict_nprms%d",&restrict_nprms)==1) continue;
			for(i=0;i<10;i++) {
				char format[128];
				sprintf(format,"restrict[%2d]%%lf",i);
				if(sscanf(line,format,&restrict[i])==1) continue;
			}
		}
	}

	fclose(fp);
	return 0;
}

// Write model parameters to *.acm file.
int DataElement::Export(void)
{
	FILE *fp=fopen(filename,"w");
	if(!fp) return 1;

	int i;
	fprintf(fp,"batch_filename %s\n",operationData.batchFilename);
	fprintf(fp,"output_filename %s\n",operationData.outputFilename);
	fprintf(fp,"\n");

	fprintf(fp,"system_type %d\n",system_type);
	fprintf(fp,"\n");

	fprintf(fp,"input_type %d\n",input_type);
	for(i=0;i<20;i++) {
		fprintf(fp,"input[%2d] %+le\n",i,input[i]);
	}
	fprintf(fp,"\n");

	fprintf(fp,"comp_type %d\n",comp_type);
	fprintf(fp,"comp_nprms %d\n",comp_nprms);
	for(i=0;i<40;i++) {
		fprintf(fp,"comp[%2d] %+le\n",i,comp[i]);
	}
	fprintf(fp,"cond_type %d\n",cond_type);
	fprintf(fp,"cond_nprms %d\n",cond_nprms);
	for(i=0;i<20;i++) {
		fprintf(fp,"cond[%2d] %+le\n",i,cond[i]);
	}
	fprintf(fp,"evap_type %d\n",evap_type);
	fprintf(fp,"evap_nprms %d\n",evap_nprms);
	for(i=0;i<20;i++) {
		fprintf(fp,"evap[%2d] %+le\n",i,evap[i]);
	}
	fprintf(fp,"subcool_type %d\n",subcool_type);
	fprintf(fp,"subcool_nprms %d\n",subcool_nprms);
	for(i=0;i<10;i++) {
		fprintf(fp,"subcool[%2d] %+le\n",i,subcool[i]);
	}
	fprintf(fp,"expdev_type %d\n",expdev_type);
	fprintf(fp,"expdev_nprms %d\n",expdev_nprms);
	for(i=0;i<10;i++) {
		fprintf(fp,"expdev[%2d] %+le\n",i,expdev[i]);
	}
	fprintf(fp,"restrict_type %d\n",restrict_type);
	fprintf(fp,"restrict_nprms %d\n",restrict_nprms);
	for(i=0;i<10;i++) {
		fprintf(fp,"restrict[%2d] %+le\n",i,restrict[i]);
	}
	fprintf(fp,"\n");

	fclose(fp);
	return 0;
}

int DataElement::WriteParams(int /*restart*/)
{
	FILE *fp = fopen("./InputDoc/acmodel.dat","w");
	if(!fp) return 1;

	int i;
	switch(system_type) {
	case 1:
		fprintf(fp,"#1 RTU fixed orifice\n");
		fprintf(fp,"#%d,1 compressor\n",comp_type);
		for(i=0;i<comp_nprms;i++) fprintf(fp,"%+.6le\n",comp[i]);
		fprintf(fp,"#%d,1 condenser\n",cond_type);
		for(i=0;i<cond_nprms;i++) fprintf(fp,"%+.6le\n",cond[i]);
		fprintf(fp,"#%d,1 evaporator\n",evap_type);
		for(i=0;i<evap_nprms;i++) fprintf(fp,"%+.6le\n",evap[i]);
		fprintf(fp,"#%d,1 expdev\n",expdev_type);
		for(i=0;i<expdev_nprms;i++) fprintf(fp,"%+.6le\n",expdev[i]);
		fprintf(fp,"#%d,2 restriction\n",restrict_type);
		for(i=0;i<restrict_nprms;i++) fprintf(fp,"%+.6le\n",restrict[i]);
		break;
	case 2:
		fprintf(fp,"#2 RTU TxV\n");
		fprintf(fp,"#%d,1 compressor\n",comp_type);
		for(i=0;i<comp_nprms;i++) fprintf(fp,"%+.6le\n",comp[i]);
		fprintf(fp,"#%d,1 condenser\n",cond_type);
		for(i=0;i<cond_nprms;i++) fprintf(fp,"%+.6le\n",cond[i]);
		fprintf(fp,"#%d,1 evaporator\n",evap_type);
		for(i=0;i<evap_nprms;i++) fprintf(fp,"%+.6le\n",evap[i]);
		fprintf(fp,"#%d,1 expdev\n",expdev_type);
		for(i=0;i<expdev_nprms;i++) fprintf(fp,"%+.6le\n",expdev[i]);
		fprintf(fp,"#%d,2 restriction\n",restrict_type);
		for(i=0;i<restrict_nprms;i++) fprintf(fp,"%+.6le\n",restrict[i]);
		break;
	case 3:
		fprintf(fp,"#3 NIST chiller\n");
		fprintf(fp,"#%d,1 compressor\n",comp_type);
		for(i=0;i<comp_nprms;i++) fprintf(fp,"%+.6le\n",comp[i]);
		fprintf(fp,"#%d,1 condenser\n",cond_type);
		for(i=0;i<cond_nprms;i++) fprintf(fp,"%+.6le\n",cond[i]);
		fprintf(fp,"#%d,1 evaporator\n",evap_type);
		for(i=0;i<evap_nprms;i++) fprintf(fp,"%+.6le\n",evap[i]);
		fprintf(fp,"#%d,1 subcooler\n",subcool_type);
		for(i=0;i<subcool_nprms;i++) fprintf(fp,"%+.6le\n",subcool[i]);
		break;
	}

	fclose(fp);

	//if(restart) gClient->CmStart();

	return 0;
}

// return 0, correct file ok
// return 1, not found
int DataElement::FileBasic(char* str, FILE* fp) 
{
	char *ptr = str;
	int cnt =0;
	double d[50];
	// basic input file
	// column      description       units
	// 1           Tcai              C
	// 2           Tewi              C
	// 3           Gca_scaler        unitless      (default 1.0)
	// 4           mw_scaler         unitless      (default 1.0)
	// 5           charge_scaler     unitless      (default 1.0)
	while(cnt<50 && sget(&ptr,&d[cnt])) {
		cnt++;
	}
	if(cnt!=5) return 1;
	fprintf(fp,"4 1 %9.4lf %9.4lf %9.4lf %9.4lf %9.4lf\n",d[0],d[1],input[16]*d[2],input[17]*d[3],input[18]*d[4]);
	return 0;
}

// return 0, correct file ok
// return 1, not found
int DataElement::FileNIST(char* str, FILE* fp) 
{
	char *ptr = str;
	int cnt = 0;
	double d[50];
	char c[128];
	const int n = 40;

	if(!sget(&ptr,c)) return 3;     // date
	if(!sget(&ptr,c)) return 4;     // time
	while(cnt<n && sget(&ptr,&d[cnt])) {
		cnt++;
	}
	if(cnt!=n) return 5;

	fprintf(fp,"4 1 %9.4lf %9.4lf %9.4lf %9.4lf %9.4lf ",d[1],d[17],d[36]*input[16],d[37]*input[17],d[38]*input[18]);
	for(cnt=0;cnt<n;cnt++) fprintf(fp,"%+le ",d[cnt]);
	fprintf(fp,"\n");

	return 0;
}


int DataElement::WriteInput(int /*restart*/)
{
	FILE *fp = fopen("./InputDoc/master.in","w");
	if(!fp) return 1;

	switch(system_type) {
	case 1:
	case 2:
		fprintf(fp,"1 %le %le %le\n",input[0],input[1],input[2]);
		fprintf(fp,"4 1 ");
		fprintf(fp,"%9.4lf ",input[3]);
		fprintf(fp,"%9.4lf ",input[4]);
		fprintf(fp,"%9.4lf ",1e-2*input[5]);
		fprintf(fp,"%9.4lf ",4.8*1e-2*input[6]);
		fprintf(fp,"%9.4lf ",2.85*1e-2*input[7]);
		fprintf(fp,"%9.4lf ",1.35*1e-2*input[8]);
		fprintf(fp,"%9.4lf ",1e-2*input[9]);
		fprintf(fp,"%9.4lf ",1e-2*input[10]);
		fprintf(fp,"\n");
		fprintf(fp,"4 1 ");
		fprintf(fp,"%9.4lf ",input[11]);
		fprintf(fp,"%9.4lf ",input[12]);
		fprintf(fp,"%9.4lf ",1e-2*input[13]);
		fprintf(fp,"%9.4lf ",4.8*1e-2*input[14]);
		fprintf(fp,"%9.4lf ",2.85*1e-2*input[15]);
		fprintf(fp,"%9.4lf ",1.35*1e-2*input[16]);
		fprintf(fp,"%9.4lf ",1e-2*input[17]);
		fprintf(fp,"%9.4lf ",1e-2*input[18]);
		fprintf(fp,"\n");

		// insert batch file data
		{
		FILE *in = fopen(operationData.batchFilename,"r");
		if(in) {
			char str[256];
			while (fgets(str,255,in)!=NULL) {
				fprintf(fp,"4 1 %s",str);
			}
			fclose(in);
		}
		}

		break;
	case 3:
		fprintf(fp,"1 %le %le %le\n",input[0],input[1],input[2]);
		fprintf(fp,"4 1 ");
		fprintf(fp,"%9.4lf ",input[3]);
		fprintf(fp,"%9.4lf ",input[4]);
		fprintf(fp,"%9.4lf ",input[13]*1e-2*input[5]);
		fprintf(fp,"%9.4lf ",input[14]*1e-2*input[6]);
		fprintf(fp,"%9.4lf ",input[15]*1e-2*input[7]);
		fprintf(fp,"\n");
		fprintf(fp,"4 1 ");
		fprintf(fp,"%9.4lf ",input[8]);
		fprintf(fp,"%9.4lf ",input[9]);
		fprintf(fp,"%9.4lf ",input[16]*1e-2*input[10]);
		fprintf(fp,"%9.4lf ",input[17]*1e-2*input[11]);
		fprintf(fp,"%9.4lf ",input[18]*1e-2*input[12]);
		fprintf(fp,"\n");

		// insert batch file data
		FILE *in = fopen(operationData.batchFilename,"r");
		if(in) {
			char str[1024];
			while (fgets(str,1024,in)!=NULL) {
				if(!FileNIST(str,fp)) continue;
				FileBasic(str,fp);
			}
			fclose(in);
		}

		break;
	};

	fclose(fp);

	//if(restart) gClient->CmStart();

	return 0;
}

//------------------------ Defaults ---------------------------------
void DataElement::DefaultSystem(int type)
{
	switch(type) {
	case 1: // RTU with fixed orifice
		strcpy(filename,"rtufixed.acm");
		system_type = 1;
		DefaultInput(1);
		DefaultCompressor(101);
		DefaultCondenser(201);
		DefaultEvaporator(301);
		DefaultExpDev(401);
		DefaultRestriction(401);
		break;
	case 2: // RTU with txv
		strcpy(filename,"rtutxv.acm");
		system_type = 2;
		DefaultInput(2);
		DefaultCompressor(101);
		DefaultCondenser(201);
		DefaultEvaporator(301);
		DefaultRestriction(401);
		break;
	case 3: // NIST chiller
		strcpy(filename,"nistchlr.acm");
		system_type = 3;
		DefaultInput(3);
		DefaultCompressor(102);
		DefaultCondenser(202);
		DefaultEvaporator(303);
		DefaultSubcooler(501);
		break;
	}
}

void DataElement::DefaultInput(int system_type)
{
	switch(system_type) {
	case 1:
	case 2:
		input[0] = 2.591243e+05;
		input[1] = 6.224963e+02;
		input[2] = 2.081208e+03;
		// driving conditions corresponding to initial guess
		input[3] = 35.0000;
		input[4] = 24.00;
		input[5] = 60.00;
		input[6] = 100;
		input[7] = 100;
		input[8] = 100;
		input[9] = 100;
		input[10] = 100;
		// driving conditions corresponding to desired solution
		input[11] = 35.0000;
		input[12] = 24.00;
		input[13] = 60.00;
		input[14] = 100;
		input[15] = 100;
		input[16] = 100;
		input[17] = 100;
		input[18] = 100;
		break;
	case 3:
		// initial guesses
		input[0]=1530.5;           // compressor outlet pressure
		input[1]=466.3;            // compressor intermediate pressure
		input[2]=141.8;	           // compressor inlet pressure
		// driving conditions (corresponding to initial guess)
		input[3]=21.0000;          // condenser inlet air temp
		input[4]=-5.0;             // evaporator inlet water temp
		input[5]=100;              // condenser air flow rate
		input[6]=100;              // evaporator water flow rate
		input[7]=100;              // charge mass multi
		// driving conditions (to solve for)
		input[8]=11.4;             // condenser inlet air temp
		input[9]=-17.3;            // evaporator inlet water temp
		input[10]=100;             // condenser air flow rate multiplier
		input[11]=100;             // evaporator water flow rate multiplier
		input[12]=100;             // charge mass multiplier
		// base flow rates and charge (corresponding to initial guess)
		input[13]=-1;              // condenser air mass flow rate
		input[14]=2.0;             // evaporator water mass flow rate
		input[15]=6.0;             // charge
		// base flow rates and charge (to solve for)
		input[16]=40.0;            // condenser air mass flow rate
		input[17]=1.8;             // evaporator water mass flow rate
		input[18]=6.0;             // charge
		break;
	}
}

void DataElement::DefaultCompressor(int type)
{
	switch(type) {
	case 101:
		comp_type = 101;
		comp_nprms = 15;
		comp[0] = 3400;		// motor rpm
		comp[1] = 3.89;		// displacement volume (in^3)
		comp[2] = 0.684;	//	polytropic efficiency
		comp[3] = 5.8463e+2;	// power coefficients
		comp[4] = 7.6857e-1;
		comp[5] = 4.1250e-5;
		comp[6] = 3.8936e+3;
		comp[7] = -4.9978;
		comp[8] = 6.9486e-1;	// volumetric efficiency
		comp[9] = 2.1534e-1;
		comp[10] = 2.4364e-9;
		comp[11] = 3.8598;
		comp[12] = 1.000;       // adjusting factor for mass flow rate
		comp[13] = 1.000;       // adjusting factor for power consumption
		comp[14] = 1.000;       // adjusting factor for polytropic efficiency
		break;
	case 102:
		comp_type = 102;
		comp_nprms = 26;
		// polytropic efficiency
		comp[0] = 0.76;
		// compressor power coefficients
		comp[1] = -7.4103E-01;
		comp[2] = +0.06380395;  //*Psuc
		comp[3] = -5.0352E-05;  //*Psuc^2
		comp[4] = +6.1147E-08;  //*Psuc^3
		comp[5] = +0.00682251;  //*Pdisc
		comp[6] = -5.4431E-07;  //*Pdisc^2
		comp[7] = -9.4343E-11;  //*Pdisc^3
		comp[8] = +6.4755E-06;  //*Psuc*Pdisc
		comp[9] = +3.2070E-09;  //*Psuc*Pdisc^2
		comp[10]= -5.7209E-10;  //*Psuc^2*Pdisc
		comp[11]= -6.0781E-12;  //*Psuc^2*Pdisc^2
		// compressor mass flow coefficients
		comp[12]= -5.7896E-04;
		comp[13]= +9.2462E-04;  //*Psuc
		comp[14]= +7.8869E-07;  //*Psuc^2
		comp[15]= -1.0301E-09;  //*Psuc^3
		comp[16]= +3.4781E-07;  //*Pdisc
		comp[17]= -4.2062E-10;  //*Pdisc^2
		comp[18]= +2.5841E-13;  //*Pdisc^3
		comp[19]= -1.0067E-07;  //*Psuc*Pdisc
		comp[20]= +1.0695E-11;  //*Psuc*Pdisc^2
		comp[21]= +9.2959E-11;  //*Psuc^2*Pdisc
		comp[22]= 0.0;          //*Psuc^2*Pdisc^2
		// adjusting factors
		comp[23]= 1.00;
		comp[24]= 1.00;
		comp[25]= 1.00;
		break;
	}
}

void DataElement::DefaultCondenser(int type)
{
	switch(type) {
		case 201:
		cond_type = 201;
		cond_nprms = 17;

		cond[0] = 8.250e-3;    // inside diameter of pipe (m)
		cond[1] = 1.035;       // length of pipe (m)
		cond[2] = 0.650e-3;    // wall thickness of pipe (m)
		cond[3] = 28.7e-3;     // outside fin diameter (m)
		cond[4] = 1.304e-3;    // space between fins (m)
		cond[5] = 0.190e-3;    // fin thickness (m)
		cond[6] = 25.4e-3;     // spacing between tubes in bank (m)
		cond[7] = 10;          // number of segments to break up each pipe
		cond[8] = 0.013;       // radius of return bend (m)
		cond[9] = 4.80;        // nominal air mass flux (kg/s/m^2)
		// Circuit pattern parameters
		cond[10] = 2;			  // number of branchs in main section
		cond[11] = 10; 		  // number of tubes in main section
		cond[12] = 8;	 		  // number of tubes in subcooling section
		// Performance scaling factors.
		cond[13] = 1.0;        // air side convection coefficient
		cond[14] = 1.0;        // refrigerant side convection coefficient
		cond[15] = 1.0;        // refrigerant side pressure drop coefficient
		cond[16] = 1.0;        // fouling factor (1.0-no, 0.0-completely)

		/***
		// NIST chiller
		cond[0]=16.5e-3;
		cond[1]=45.0;
		cond[2]=0.65e-3;
		cond[3]=28.7e-3;
		cond[4]=1.304e-3;
		cond[5]=0.19e-3;
		cond[6]=25.4e-3;
		cond[7]=30;
		cond[8]=13.0e-3;
		cond[9]=4.8;
		cond[10]=1;
		cond[11]=0;
		cond[12]=1;
		cond[13]=1.00;
		cond[14]=1.00;
		cond[15]=1.00;
		cond[16]=1.00;
		***/

		break;
		case 202:
		cond_type = 202;
		cond_nprms = 20;

		cond[0] = 8.250e-3;    // inside diameter of pipe (m)
		cond[1] = 1.035;       // length of pipe (m)
		cond[2] = 0.650e-3;    // wall thickness of pipe (m)
		cond[3] = 28.7e-3;     // outside fin diameter (m)
		cond[4] = 1.304e-3;    // space between fins (m)
		cond[5] = 0.190e-3;    // fin thickness (m)
		cond[6] = 25.4e-3;     // spacing between tubes in bank (m)
		cond[7] = 10;          // number of segments to break up each pipe
		cond[8] = 0.013;       // radius of return bend (m)
		cond[9] = 4.80;        // nominal air mass flux (kg/s/m^2)
		// Circuit pattern parameters
		cond[10] = 2;			  // number of branchs in main section
		cond[11] = 10; 		  // number of tubes in main section
		cond[12] = 8;	 		  // number of tubes in subcooling section
		// Performance scaling factors.
		cond[13] = 1.0;        // air side convection coefficient
		cond[14] = 1.0;        // refrigerant side convection coefficient
		cond[15] = 1.0;        // refrigerant side pressure drop coefficient
		cond[16] = 1.0;        // fouling factor (1.0-no, 0.0-completely)
		// CFM=A+B*Pdis+C*Pdis^2
		cond[17] = 861160;        // A
		cond[18] = 7567.8;        // B
		cond[19] = 16.677;        // C
		
		/***
		// NIST chiller
		cond[0]=16.5e-3;
		cond[1]=45.0;
		cond[2]=0.65e-3;
		cond[3]=28.7e-3;
		cond[4]=1.304e-3;
		cond[5]=0.19e-3;
		cond[6]=25.4e-3;
		cond[7]=30;
		cond[8]=13.0e-3;
		cond[9]=4.8;
		cond[10]=1;
		cond[11]=0;
		cond[12]=1;
		cond[13]=1.00;
		cond[14]=1.00;
		cond[15]=1.00;
		cond[16]=1.00;
		// CFM=A+B*Pdis+C*Pdis^2
		cond[17] = +861160;        // A
		cond[18] = -7567.8;        // B
		cond[19] = +16.677;        // C
		***/

		break;
	}
}

void DataElement::DefaultEvaporator(int type)
{
	switch(type) {
	case 301:
		evap_type = 301;
		evap_nprms = 17;
		evap[0] = 8.250e-3;    // inside diameter of pipe (m)
		evap[1] = 0.654;       // length of pipe (m)
		evap[2] = 0.650e-3;    // wall thickness of pipe (m)
		evap[3] = 21.5e-3;     // outside fin diameter (m)
		evap[4] = 1.22e-3;     // space between fins (m)
		evap[5] = 0.190e-3;    // fin thickness (m)
		evap[6] = 25.4e-3;     // spacing between tubes in bank (m)
		evap[7] = 10;          // number of segments to break up each pipe
		evap[8] = 0.013;       // radius of return bend (m)
		evap[9] = 2.85;        // nominal air mass flux (kg/s/m^2)
		// Circuit pattern parameters
		evap[10] = 6;			  // number of rows high
		evap[11] = 2;			  // number of rows deep
		// Performance scaling factors.
		evap[12] = 1.0;        // air side convection coefficient (dry coil)
		evap[13] = 1.0;        // refrigerant side convection coefficient
		evap[14] = 1.0;        // refrigerant side pressure drop coefficient
		evap[15] = 1.0;        // air side convestion coefficient (wet coil)
		evap[16] = 1.0;        // fouling factor (1.0-no, 0.0-completely)
		break;
	case 303:
		// Shell-in-tube evaporator parameters (type 303)
		evap_type = 303;
		evap_nprms = 3;
		evap[0]=1100;
		evap[1]=0.15;
		evap[2]=30;
		break;
	}
}

void DataElement::DefaultExpDev(int type)
{
	switch(type) {
	case 401:
		expdev_type = 401;
		expdev_nprms = 3;
		expdev[0] = 0.9398e-3;     // orifice diameter (m).
		expdev[1] = 12e-3;         // expansion device length(m).
		//expdev[2] = 1.135412;      // adjusting factor for orifice diameter.
		expdev[2]=1.0;
		break;
	}
}

void DataElement::DefaultRestriction(int type)
{
	switch(type) {
	case 401:
		restrict_type = 401;
		restrict_nprms = 3;
		restrict[0] = 8.25e-3;       // orifice diameter (m).
		restrict[1] = 12e-3;         // expansion device length(m).
		restrict[2] = 1.0;           // adjusting factor for orifice diameter.
		break;
	}
}

void DataElement::DefaultSubcooler(int type)
{
	switch(type) {
		case 501:
		// Subcooler (type 501)
		subcool_type = 501;
		subcool_nprms = 4;
		subcool[0]=300;
		subcool[1]=0.0015;
		subcool[2]=0.015;
		subcool[3]=30;
		break;
	}
}

// assumes filename is less than 256 chars
void DataElement::strip_path(char* filename)
{
	char fn_tmp[256];
	strcpy(fn_tmp,filename);

	char *ptr_tmp = fn_tmp - 1;
	char *ptr;
	do {
		ptr = ptr_tmp + 1;
		ptr_tmp = strchr(ptr,'\\');
	} while(ptr_tmp!=NULL);

	strcpy(filename,ptr);
}

