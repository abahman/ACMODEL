#include "stdafx.h"
#include "outelem.h"
#include "wair.h"

// these guys enable the API call that opens browers and Excel to show results
// after the calculations are complete.
#include <windows.h>
#include <shellapi.h>

	// functions
OutputDataElement::OutputDataElement() 
	{
		ClearAll();
		outputFilename[0]=0;
	}
	OutputDataElement::~OutputDataElement() 
	{
		CloseOutputFile();
	}

	//**************************************************************************
	//**************************************************************************
	//**************************************************************************
	// Create HTML output file

	void OutputDataElement::StartOutputFile(const char* fn)
	{
		strcpy(outputFilename,fn);

		FILE* fp = fopen(outputFilename,"w");
		fprintf(fp,"<TABLE border=\"1\">\n");
		fprintf(fp,"<CAPTION><B>ACModel Output Results</B></CAPTION>\n");
		fprintf(fp,"<TR>\n");
		fprintf(fp,"<TH colspan=3>Convergance</TH>\n");
		fprintf(fp,"<TH colspan=5>Equipment</TH>\n");
		fprintf(fp,"<TH colspan=4>Inputs</TH>\n");
		fprintf(fp,"<TH colspan=6>Outputs</TH>\n");
		fprintf(fp,"</TR>\n");
		fprintf(fp,"\t<TR>\n");
		fprintf(fp,"\t\t<TH>Error?</TH>\n");
		fprintf(fp,"\t\t<TH>Conv?</TH>\n");
		fprintf(fp,"\t\t<TH>CE</TH>\n");
		fprintf(fp,"\t\t<TH>ID_CFM (cfm)</TH>\n");
		fprintf(fp,"\t\t<TH>OD_CFM (cfm)</TH>\n");
		fprintf(fp,"\t\t<TH>Charge (kg)</TH>\n");
		fprintf(fp,"\t\t<TH>Dadj</TH>\n");
		fprintf(fp,"\t\t<TH>nvAdj</TH>\n");
		fprintf(fp,"\t\t<TH>Tcai (F)</TH>\n");
		fprintf(fp,"\t\t<TH>Teai (F)</TH>\n");
		fprintf(fp,"\t\t<TH>Reai (%%)</TH>\n");
		fprintf(fp,"\t\t<TH>Twb_eai (F)</TH>\n");
		fprintf(fp,"\t\t<TH>ET (F)</TH>\n");
		fprintf(fp,"\t\t<TH>CTOA (F)</TH>\n");
		fprintf(fp,"\t\t<TH>SH (F)</TH>\n");
		fprintf(fp,"\t\t<TH>SC (F)</TH>\n");
		fprintf(fp,"\t\t<TH>Qevap (tons)</TH>\n");
		fprintf(fp,"\t\t<TH>EER (BTU/hr/kW)</TH>\n");
		fprintf(fp,"\t</TR>\n");
		fclose(fp);
	}
	void OutputDataElement::AppendOutputToFile()
	{
		if(strlen(outputFilename)<=0) {
			StartOutputFile("./OutputDoc/ACModelOutput.htm");
		}

		FILE* fp = fopen(outputFilename,"a");
		fprintf(fp,"\t<TR>\n");
		fprintf(fp,"\t\t<TD>%s</TD>\n",isError?"ERROR":"OK");
		fprintf(fp,"\t\t<TD>%s</TD>\n",isError?"N/A":isConv?"Y":"N");
		fprintf(fp,"\t\t<TD>%.4le</TD>\n",CE);
		fprintf(fp,"\t\t<TD>%.1lf</TD>\n",cfmea);
		fprintf(fp,"\t\t<TD>%.1lf</TD>\n",cfmca);
		fprintf(fp,"\t\t<TD>%.4lf</TD>\n",Charge/*Charge_target*/);
		fprintf(fp,"\t\t<TD>%.4lf</TD>\n",Dadj);
		fprintf(fp,"\t\t<TD>%.4lf</TD>\n",nvAdj);
		fprintf(fp,"\t\t<TD>%.1lf</TD>\n",CtoF(Tci));
		fprintf(fp,"\t\t<TD>%.1lf</TD>\n",CtoF(Tei));
		fprintf(fp,"\t\t<TD>%.0lf</TD>\n",100*HumPei);
		fprintf(fp,"\t\t<TD>%.1lf</TD>\n",CtoF(Twb_ei));
		if(0/*isError*/) {
			fprintf(fp,"\t\t<TD>&nbsp;</TD>\n");
			fprintf(fp,"\t\t<TD>&nbsp;</TD>\n");
			fprintf(fp,"\t\t<TD>&nbsp;</TD>\n");
			fprintf(fp,"\t\t<TD>&nbsp;</TD>\n");
			fprintf(fp,"\t\t<TD>&nbsp;</TD>\n");
			fprintf(fp,"\t\t<TD>&nbsp;</TD>\n");
		} else {
			fprintf(fp,"\t\t<TD>%.1lf</TD>\n",CtoF(Tsat[1]));
			fprintf(fp,"\t\t<TD>%.1lf</TD>\n",DCtoDF(Tsat[4]-Tci));
			fprintf(fp,"\t\t<TD>%.1lf</TD>\n",DCtoDF(T[8]-Tsat[8]));
			fprintf(fp,"\t\t<TD>%.1lf</TD>\n",DCtoDF(Tsat[4]-T[4]));
			fprintf(fp,"\t\t<TD>%.4lf</TD>\n",WtoBTUH(-1*Qe)/12000);
			fprintf(fp,"\t\t<TD>%.4lf</TD>\n",WtoBTUH((-1*Qe/Ecmp)));
		}
		fprintf(fp,"\t</TR>\n");
		fclose(fp);

		//ShellExecute(NULL, "open", "./OutputDoc/ACModelOutput.htm", NULL, NULL, SW_SHOWNORMAL);
	}

	//**************************************************************************
	//**************************************************************************
	//**************************************************************************
	// Create in input file that expands run in anothter dimension

	void OutputDataElement::StartMasterInputFile()
	{
		{
		FILE* fp = fopen("master_rwb.in","w");
		fprintf(fp,"\n");
		fclose(fp);
		}
		{
		FILE* fp = fopen("master_charge.in","w");
		fprintf(fp,"\n");
		fclose(fp);
		}
	}

	void OutputDataElement::AppendMasterInputFile()
	{
		//int i;

		// vary charge only
		{
		FILE* fp = fopen("master_charge.in","a");
		fprintf(fp,"1\t%lf\t%lf\t%lf\n",H[1],P[1],P[2]);
		const double Charge_start = 1.0; 
		const double Charge_step = 0.1;
		for(int i=0;i<=5;i++) {
			const double Chg = Charge_target*(Charge_start + (double)i*Charge_step);
			fprintf(fp,"4\t1\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t1.0\t1.0\n",Tci,Tei,HumPei,Gca,Gea,Chg);
		}
		fprintf(fp,"1\t%lf\t%lf\t%lf\n",H[1],P[1],P[2]);
		for(int i=1;i<=9;i++) {
			const double Chg = Charge_target*(Charge_start - (double)i*Charge_step);
			fprintf(fp,"4\t1\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t1.0\t1.0\n",Tci,Tei,HumPei,Gca,Gea,Chg);
		}
		fclose(fp);
		}

		// vary charge only
		{
		FILE* fp = fopen("master_rwb.in","a");
		fprintf(fp,"1\t%lf\t%lf\t%lf\n",H[1],P[1],P[2]);
		const double Twb_start = FtoC(68);
		const double Twb_step = DFtoDC(2);
		for(int i=0;i<=4;i++) {
			const double Twb = Twb_start + (double)i*Twb_step;
			const double Tdb = Twb + DFtoDC(13);
			TP tp = TZtoTP(Tdb,Twb);
			fprintf(fp,"4\t1\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t1.0\t1.0\n",Tci,tp.T,tp.P,Gca,Gea,Charge_target);
		}
		fprintf(fp,"1\t%lf\t%lf\t%lf\n",H[1],P[1],P[2]);
		for(int i=1;i<=8;i++) {
			const double Twb = Twb_start - (double)i*Twb_step;
			const double Tdb = Twb + DFtoDC(13);
			TP tp = TZtoTP(Tdb,Twb);
			fprintf(fp,"4\t1\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t1.0\t1.0\n",Tci,tp.T,tp.P,Gca,Gea,Charge_target);
		}
		fclose(fp);
		}
	}

	void OutputDataElement::CloseOutputFile()
	{
		if(strlen(outputFilename)<=0) return;

		FILE* fp = fopen(outputFilename,"a");
		fprintf(fp,"</TABLE>\n");
		fclose(fp);

		outputFilename[0]=0;
	}

	void OutputDataElement::ClearAll() //set zero to all the output variables
	{
		//int i;
		time = 0;
		for(int i=0;i<15;i++) {
			T[i]=0;
			X[i]=0;
			P[i]=0;
			H[i]=0;
			Tsat[i]=0;
		}
		for(int i=0;i<5;i++) {
			mr[i]=0;
		}

		for(int i=0;i<128;i++)
		{ 
			Array[i]=0;
		}
		for(int i=0;i<8;i++)
		{ 
			m[i]=0;
		}
		mc=0;
		Tci=0;
		Tco=0;
		me=0;
		Tei=0;
		Teo=0;
		Ecmp=0;
		Qc=0;
		Qe=0;
		Charge=0;
		Charge_target=0;
		P1=0;
		P2=0;
		mr_TP=0;
		CE=1e20;
		HumPei=0;
		HumPeo=0;
		devMr=0;
		devH=0;
		devCharge=0;
		inputP1=inputP2=inputH1=0;
		Cond_AirDP = 0;
		Eva_AirDP = 0;
		cfmca = 0;
		cfmea = 0;
		Tsub4=0;
		Tsub5=0;
		Tsup1=0;
		LiqLCond=0;
		TPLCond=0;
		VapLCond=0;
		LiqLEvap=0;
		TPLEvap=0;
		VapLEvap=0;

		auxDataCnt=0;
		for(int i=0;i<128;i++) auxData[i]=0;
	}

	void OutputDataElement::output(char str_res)//shenboadd
	{
	static int FirstTime=0;
	std::ofstream outstrm; 
	outstrm.open("./outputDoc/output.xls",std::ios::app); // for outputing the data of all the state points
	if(FirstTime==0)
	{
		outstrm<<"State"<<"	"<<"Ecmp"<<"	"<<"Qc"<<"	"<<"Qe"<<"	"
		<<"mca"<<"	"<<"mea"<<"	"<<"Tci"<<"	"<<"Tco"<<"	"<<"Tei"<<"	"<<"HumPei"<<"	"
		<<"Teo"<<"	"<<"HumPeo"<<"	"<<"Charge"<<"	"<<"mr_TP"<<"	"
		<<"H[1]"<<"	"<<"T[1]"<<"	"<<"X[1]"<<"	"<<"P[1]"<<"	"<<"Ts[1]"<<"	"<<"Tsup"<<"	"
		<<"H[2]"<<"	"<<"T[2]"<<"	"<<"X[2]"<<"	"<<"P[2]"<<"	"<<"Ts[2]"<<"	"
		<<"H[3]"<<"	"<<"P[3]"<<"	"
		<<"H[4]"<<"	"<<"T[4]"<<"	"<<"X[4]"<<"	"<<"P[4]"<<"	"<<"Ts[4]"<<"	"<<"Tsub4"<<"	"
		<<"H[5]"<<"	"<<"T[5]"<<"	"<<"X[5]"<<"	"<<"P[5]"<<"	"<<"Ts[5]"<<"	"<<"Tsub5"<<"	"
		<<"H[6]"<<"	"<<"P[6]"<<"	"
		<<"T[7]"<<"	"<<"X[7]"<<"	"<<"P[7]"<<"	"<<"Ts[7]"<<"	"
		<<"H[8]"<<"	"<<"P[8]"<<"	"
		<<"Cond_AirDP"<<"	"<<"Eva_AirDP"<<"	"
		<<"cfmca"<<"	"<<"cfmea"<<"	"
		<<"devMr"<<"	"<<"devH"<<"	"<<"devCharge"<<"	"<<"CE"<<"	"
		<<"inputH1"<<"	"<<"inputP1"<<"	"<<"inputP2"<<"	"
		<<"LiqLCond"<<"	"<<"TPLCond"<<"	"<<"VapLCond"<<"	"
		<<"LiqLEvap"<<"	"<<"TPLEvap"<<"	"<<"VapLEvap"<<"	"
		<<"m[0]"<<"	"<<"m[1]"<<"	"<<"m[2]"<<"	"<<"m[3]"<<"	"<<"m[4]"<<"	"<<"m[5]"<<"	"<<"m[6]"<<"	"<<"m[7]"<<"	"
		<<"Time_dev"
		<<std::endl;

		FirstTime=1;
	}

	outstrm<<str_res<<"	"<<Ecmp<<"	"<<Qc<<"	"<<-1*Qe<<"	"
		<<mc<<"	"<<me<<"	"<<Tci<<"	"<<Tco<<"	"<<Tei<<"	"<<HumPei<<"	"
		<<Teo<<"	"<<HumPeo<<"	"<<Charge<<"	"<<mr_TP<<"	"
		<<H[1]<<"	"<<T[1]<<"	"<<X[1]<<"	"<<P[1]<<"	"<<Tsat[1]<<"	"<<Tsup1<<"	"
		<<H[2]<<"	"<<T[2]<<"	"<<X[2]<<"	"<<P[2]<<"	"<<Tsat[2]<<"	"
		<<H[3]<<"	"<<P[3]	<<"	"
		<<H[4]<<"	"<<T[4]<<"	"<<X[4]<<"	"<<P[4]<<"	"<<Tsat[4]<<"	"<<Tsub4<<"	"
		<<H[5]<<"	"<<T[5]<<"	"<<X[5]<<"	"<<P[5]<<"	"<<Tsat[5]<<"	"<<Tsub5<<"	"
		<<H[6]<<"	"<<P[6]<<"	"
		<<T[7]<<"	"<<X[7]<<"	"<<P[7]<<"	"<<Tsat[7]<<"	"
		<<H[8]<<"	"<<P[8]<<"	"
		<<Cond_AirDP<<"	"<<Eva_AirDP<<"	"
		<<cfmca<<"	"<<cfmea<<"	"
		<<devMr<<"	"<<devH<<"	"<<devCharge<<"	"<<CE<<"	"
		<<inputH1<<"	"<<inputP1<<"	"<<inputP2<<"	"
		<<LiqLCond<<"	"<<TPLCond<<"	"<<VapLCond<<"	"
		<<LiqLEvap<<"	"<<TPLEvap<<"	"<<VapLEvap<<"	"
		<<m[0]<<"	"<<m[1]<<"	"<<m[2]<<"	"<<m[3]<<"	"<<m[4]<<"	"<<m[5]<<"	"<<m[6]<<"	"<<m[7]<<"	"
		<<Time_dev
		<<std::endl;

/*	outstrm<<"spline";
	for(int i=1;i<=58;i++)
	{outstrm<<"	"<<Array[i];}
	outstrm<<std::endl;*/

	outstrm.close();
	ClearAll();
	}

	void OutputDataElement::outputTOarray()
	{
	Array[1]=Ecmp;
	Array[2]=Qc;
	Array[3]=-1*Qe;
	Array[4]=mc;
	Array[5]=me;
	Array[6]=Tci;
	Array[7]=Tco;
	Array[8]=Tei;
	Array[9]=HumPei;
	Array[10]=Teo;
	Array[11]=HumPeo;
	Array[12]=Charge;
	Array[13]=mr_TP;
	Array[14]=H[1];
	Array[15]=T[1];
	Array[16]=X[1];
	Array[17]=P[1];
	Array[18]=Tsat[1];
	Array[19]=T[1]-Tsat[1];
	Array[20]=H[2];
	Array[21]=T[2];
	Array[22]=X[2];
	Array[23]=P[2];
	Array[24]=Tsat[2];
	Array[25]=H[3];
	Array[26]=P[3];
	Array[27]=H[4];
	Array[28]=T[4];
	Array[29]=X[4];
	Array[30]=P[4];
	Array[31]=Tsat[4];
	Array[32]=-1*(T[4]-Tsat[4]);
	Array[33]=H[5];
	Array[34]=T[5];
	Array[35]=X[5];
	Array[36]=P[5];
	Array[37]=Tsat[5];
	Array[38]=-1*(T[5]-Tsat[5]);
	Array[39]=H[6];
	Array[40]=P[6];
	Array[41]=T[7];
	Array[42]=X[7];
	Array[43]=P[7];
	Array[44]=Tsat[7];
	Array[45]=H[8];
	Array[46]=P[8];
	Array[47]=Cond_AirDP;
	Array[48]=Eva_AirDP;
	Array[49]=cfmca;
	Array[50]=cfmea;
	Array[51]=devMr;
	Array[52]=devH;
	Array[53]=devCharge;
	Array[54]=CE;
	Array[55]=inputH1;
	Array[56]=inputP1;
	Array[57]=inputP2;
	}

	void OutputDataElement::arrayback()
	{
	Ecmp=Array[1];
	Qc=Array[2];
	Qe=-1*Array[3];
	mc=Array[4];
	me=Array[5];
	Tci=Array[6];
	Tco=Array[7];
	Tei=Array[8];
	HumPei=Array[9];
	Teo=Array[10];
	HumPeo=Array[11];
	Charge=Array[12];
	mr_TP=Array[13];
	H[1]=Array[14];
	T[1]=Array[15];
	X[1]=Array[16];
	P[1]=Array[17];
	Tsat[1]=Array[18];
	Tsup1=Array[19];
	H[2]=Array[20];
	T[2]=Array[21];
	X[2]=Array[22];
	P[2]=Array[23];
	Tsat[2]=Array[24];
	H[3]=Array[25];
	P[3]=Array[26];
	H[4]=Array[27];
	T[4]=Array[28];
	X[4]=Array[29];
	P[4]=Array[30];
	Tsat[4]=Array[31];
	Tsub4=Array[32];
	H[5]=Array[33];
	T[5]=Array[34];
	X[5]=Array[35];
	P[5]=Array[36];
	Tsat[5]=Array[37];
	Tsub5=Array[38];
	H[6]=Array[39];
	P[6]=Array[40];
	T[7]=Array[41];
	X[7]=Array[42];
	P[7]=Array[43];
	Tsat[7]=Array[44];
	H[8]=Array[45];
	P[8]=Array[46];
	Cond_AirDP=Array[47];
	Eva_AirDP=Array[48];
	cfmca=Array[49];
	cfmea=Array[50];
	devMr=Array[51];
	devH=Array[52];
	devCharge=Array[53];
	CE=Array[54];
	inputH1=Array[55];
	inputP1=Array[56];
	inputP2=Array[57];
	}
