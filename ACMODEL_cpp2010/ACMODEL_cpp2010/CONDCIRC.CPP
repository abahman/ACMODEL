#include "stdafx.h"

#include <stdio.h>
#include "r22.h"
#include "air.h"
#include "mass.h"
#include "cond.h"
#include "errorlog.h"
#include "spline1d.h"//B.S.
#include <math.h>//B.S.
#include <string.h>//B.S.
#include "getdata.h"//B.S.
#include "parse.h"//B.S.
#include "cond.h"//B.S

extern ErrorLog errorLog;
extern AirProp air;

struct CondNode
{
	int NodNo;//node number
	int InNum;//number of tubes flowing into the node
	int OutNum;//number of the tubes flowing out of the node
	int *BranIN;//index of the tube branch flowing in
	int *BranOUT;//index of the tube branch flowing out
};

struct TubCndSeg
{
	double Tai;
	double hao;
};

struct TubeCond
{
	int TubNo;//tube No
	int RowNo;//the row no where the tube located
	int RefUpstream;//Upstream tube No at the refrigerant side (flow direction)
	int AirUpstreamUpper;//Upstream tube No at the air upper side
	int AirUpstreamLower;//Upstream tube No at the air lower side
	double GaFac;//air flow distribution factor
	double Ga;//maximum air flow rate in the segment
	int even;//tube flow direction
	TubCndSeg *Seg;
	HP HPi;//refrigerant inlet state
	HP HPo;//refrigerant outlet state
	MASS m;//mass and volume of the tube
};

struct CondBranch
{
	int BranNo;//branch number
	int EqulNo;//the equivalent branch number
	int Ini;//singal variable
	double GrFac;//mass flow distribution factor of the branch
	double Gr;//mass flux of the branch
	int TubNum;//total tube numbers in the branch
	int *TubNo;//index of the tubes in the branch
	HP HPi;//inlet enthalpy and pressure 
	HP HPo;//outlet enthalpy and pressure
	MASS m;//mass and volume
	double Para_Struc[12];//output parameters of this condenser branch
};

class StructCond {
public:
	StructCond(const char* fn);//reading the data file in
	~StructCond();//delete the memory

	int _CondCircuit(double mr,HP HPi,double Tai,double Ga,HP *HPo,double *Tao,MASS *m,CGP *P);//condenser fuction
	int	Cal_HP(CGP *P, int Pos,HP *HPo);//function to calculate the heat transfer and pressure drop
	char filename[128];//data file name
	int error;//for keeping track of the errorLog
	int SegNum;//tube flow direction
	int TubeNum;	//number of all tubes
	int BranNum; // number of branch
	int NodNum;  // number of node
	int RowNum;	 //number of Row
	double AreaFront;//front area
	double Volum;//air flow volume
	double *GaRow;//average air flow distribution factor of each row
	CondBranch *Bra;//array of branch objects
	CondNode *Nod;//array of node objects
	TubeCond *Tub;//array of tube objects
	};


// local functions
int CondCircuitMultipleTube(double,HP,double,double,HP*,double*,MASS*,CGP*);
int CondCircuitMultipleTube_204(double,HP,double,double,HP*,double*,MASS*,CGP*);
int CondCircuitSingleTube(double,HP,double,double,HP*,double*,MASS*,CGP*);

/***********************************************************************
Condenser circuitry pattern.
***********************************************************************/
int CondCircuit(int type,double mr,HP HPi,double Tai,double Ga,HP *HPo,double *Tao,MASS *m,CGP *P)
{
	static StructCond CondInfo("./InputDoc/CondStruc.in");
	switch(type) {
		case 201: 
			CondInfo._CondCircuit(mr,HPi,Tai,Ga,HPo,Tao,m,P);
			break;
		case 202:// Carrier RTU at Herrick
			CondCircuitMultipleTube(mr,HPi,Tai,Ga,HPo,Tao,m,P);
			break;
		case 203: // Single finned tube condenser coil
			CondCircuitSingleTube(mr,HPi,Tai,Ga,HPo,Tao,m,P);
			break;
		case 204://user defined
			CondCircuitMultipleTube_204(mr,HPi,Tai,Ga,HPo,Tao,m,P);
			break;
		default:
			errorLog.Add("CondCircuit","Model not found");
			break;
	}
	if(errorLog.IsError()) {
		char str[128];
		sprintf(str,"type=%d",type);
		errorLog.Add("CondCircuit",str);
		return 0;
	}

	return 0;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------

//user defined
int StructCond::_CondCircuit(double mr,HP HPi,double Tai,double Ga,HP *HPo,double *Tao,MASS *m,CGP *P)
{
	int i=0,j=0,k=0;
	double hai=0,hao=0;//air inlet and outlet enthalpy
	double Gr=0;//mass flux

	if(error) 
	{
	errorLog.Add("StructCond::_CondCircuit","Class fails to be constructed");
	return -1;
	}
	

	//air side initializing
	hai = air.h(Tai);
	if(errorLog.IsError()) {
	errorLog.Add("StructCond::_CondCircuit","hai");
	return -1;
	}
	const double rho_air=1/air.v(Tai);
	if(errorLog.IsError()) {
	errorLog.Add("StructCond::_CondCircuit","Tai");
	return -1;
	}

	const double Ma =Ga*rho_air*4.719e-4; 
	Ga=Ma/AreaFront*P->vsp*P->Ls/((P->vsp-P->Do)*(P->Ls-P->N*P->th));

	for(int i=0;i<RowNum;i++)
	{
	GaRow[i]=0;
	}
	
	for(int j=0;j<RowNum;j++)
	{
		int N=0;
		for(int i=0;i<TubeNum;i++)
		{
		if(Tub[i].RowNo==j)
		{
			for(k=0;k<SegNum;k++)
			{
			Tub[i].Seg[k].Tai=Tai;
			
			Tub[i].Seg[k].hao = air.h(Tai);
			if(errorLog.IsError()) {
				errorLog.Add("CondCircuit201","hao");
				return -1;}
			}

			if(Tub[i].RowNo>0)
			{
			const int Upper= Tub[i].AirUpstreamUpper;
			const int Lower = Tub[i].AirUpstreamLower;
			if(Upper>=0&&Lower>=0)
			{
			Tub[i].GaFac=(Tub[Upper].GaFac+Tub[Lower].GaFac)/2;
			}
			else if(Upper>=0)
			{
			Tub[i].GaFac=Tub[Upper].GaFac;
			}
			else
			{Tub[i].GaFac=Tub[Lower].GaFac;}
			}
		GaRow[j]=GaRow[j]+Tub[i].GaFac;
		N=N+1;
		}//ifend	
	}//i circle end
		GaRow[j]=GaRow[j]/N;
	}//j circle end
	
		for(int i=0;i<TubeNum;i++)
		{
		const int RowN=Tub[i].RowNo;
		Tub[i].Ga=Tub[i].GaFac/GaRow[RowN]*Ga;
		}
	

	//initialize
	double H_out=0;
	double Res=0;
	int IterN=0;

	do{

	P->VapL =0;
	P->TPL = 0;
	P->LiqL = 0;
	P->V_Vap =0;
	P->V_TP =0;
	P->V_Liq = 0;
	P->m_Vap = 0;
	P->m_TP = 0;
	P->m_Liq = 0;	
	P->UA_Vap = 0;
	P->UA_TP = 0;
	P->UA_Liq = 0;

	IterN=IterN+1;
	Gr=mr/P->Ax;
	*HPo=HPi;//HPo is the intermediate variable in the calculation

	for(int i=0;i<BranNum;i++) {Bra[i].Ini=0;}

	for(int i=0;i<NodNum;i++)//inlet nodes
	{
	if(Nod[i].BranIN[0]<0)//no inlet branch, only from the hot gas line
	{
	CondMan(Nod[i].InNum,Nod[i].OutNum,P->Di,&Gr,HPo);
	
	if(errorLog.IsError()) {
	errorLog.Add("StructCond::_CondCircuit","CondMan0");
	return -1;}

	for(int j=0;j<Nod[i].OutNum;j++)//states flowing out from the node
	{
	const int jj = Nod[i].BranOUT[j];//index of the outlet branches
	Bra[jj].HPi=*HPo;
	Bra[jj].Gr=Gr*Bra[jj].GrFac;
	Bra[jj].Ini=1;
	}//end j circle

	}//endif
	}//end i circle
	
	Cal_HP(P,0,HPo);//heat transfer and pressure drop calculation

	if(errorLog.IsError()) {
	errorLog.Add("StructCond::_CondCircuit","Cal_HP(P,0)");
	return -1;}
	
	for(int i=0;i<BranNum;i++) {Bra[i].Ini=0;}

	Gr=0;
	HPo->H=0;HPo->P=0;

	//nodes in the middle
	for(int i=0;i<NodNum;i++)
	{
	if(Nod[i].BranIN[0]>=0&&Nod[i].BranOUT[0]>=0)//nodes in the middle
	{

	for(int j=0;j<Nod[i].InNum;j++)//node inlet state
	{
	const int jj= Nod[i].BranIN[j];
	Gr=Gr+Bra[jj].Gr;
	HPo->H=Bra[jj].HPo.H*Bra[jj].Gr+HPo->H;
	HPo->P=Bra[jj].HPo.P*Bra[jj].Gr+HPo->P;
	}

	HPo->H=HPo->H/Gr;
	HPo->P=HPo->P/Gr;
	Gr=Gr/Nod[i].InNum;

	CondMan(Nod[i].InNum,Nod[i].OutNum,P->Di,&Gr,HPo);
	if(errorLog.IsError()) {
	errorLog.Add("StructCond::_CondCircuit","CondMan1");
	return -1;}

	for(int j=0;j<Nod[i].OutNum;j++)
	{
	const int jj = Nod[i].BranOUT[j];//index of outlet branches
	Bra[jj].HPi=*HPo;
	Bra[jj].Gr=Gr*Bra[jj].GrFac;
	Bra[jj].Ini=1;
	}//end j circle

	}//endif
	}//end i circle

	Cal_HP(P,1,HPo);//heat transfer and pressure drop calculation

	if(errorLog.IsError()) {
	errorLog.Add("StructCond::_CondCircuit","Cal_HP(P,1)");
	return -1;}
	
	Gr=0;
	HPo->H=0;HPo->P=0;

	for(int i=0;i<BranNum;i++)	{Bra[i].Ini=0;}
	
	//end nodes
	for(int i=0;i<NodNum;i++)
	{
	if(Nod[i].BranOUT[0]<0)//no outlet branch except the liquid line
	{
	for(int j=0;j<Nod[i].InNum;j++)
	{
	const int jj= Nod[i].BranIN[j];
	Gr=Gr+Bra[jj].Gr;
	HPo->H=Bra[jj].HPo.H*Bra[jj].Gr+HPo->H;
	HPo->P=Bra[jj].HPo.P*Bra[jj].Gr+HPo->P;
	}
	HPo->H=HPo->H/Gr;
	HPo->P=HPo->P/Gr;
	Gr=Gr/Nod[i].InNum;

	CondMan(Nod[i].InNum,Nod[i].OutNum,P->Di,&Gr,HPo);
	if(errorLog.IsError()) {
	errorLog.Add("StructCond::_CondCircuit","CondMan2");
	return -1;}

	}//endif

	}//end i circle
	
	if(RowNum==1) Res=0;
	else
	{
	Res=2*(HPo->H-H_out)/(HPo->H+H_out);}

	H_out=HPo->H;
	
	}while(fabs(Res)>1e-7&&IterN<20);

	if(fabs(Res)>1e-5) 
	{
	errorLog.Add("StructCond::_CondCircuit","can't reach the required tolerance");
	return -1;
	}

	m->m=0;m->V=0;

	for(int i=0;i<BranNum;i++)
	{
	m->m=m->m+Bra[i].m.m;
	m->V=m->V+Bra[i].m.V;
	}
	
	hao=mr*(HPi.H-HPo->H)/Ma+hai;
	*Tao = HatoTa(hao);

	if(errorLog.IsError()) {
	errorLog.Add("StructCond::_CondCircuit","*Tao");
	return -1;}
	
	m->m = m->m+(2.86-2.444)+(2.87-2.785);//one point charge tuning for varied charge
	//m->m=m->m+(2.858-2.533)+(2.953-2.8419)+(3.28-3.211)/(13.981-7.667)*(P->LiqL-7.667);//corresponding to three points charge tuning
	//m->m = m->m+(2.86-2.444)+(2.87-2.785)+(3.10-3.045)/(5.863-3.608)*(P->LiqL-3.608);//For two-point charge tuning
	//m->m = m->m+(2.86-2.444); one point charge tuning for fixed charge
	return 0;
};

int StructCond::Cal_HP(CGP *P, int Pos, HP *HPo)
{
	//heat transfer and pressure drop of the two-phase region
	int i=0,j=0,k=0;
	int Logic=0;
	double hao=0;
	MASS Sm;
	MASS mi;
	Sm.m=0; Sm.V=0;
	CGP Bak;//keeping the information of the condenser

	Bak=*P;
	
	if(Pos<0||Pos>1)
	{
	errorLog.Add("StructCond::Cal_HP","Wrong position");
	return -1;
	}

	for(int i=0;i<BranNum;i++)
	{
	//branch parameter clear zero
	P->VapL =0;
	P->TPL = 0;
	P->LiqL = 0;
	P->V_Vap =0;
	P->V_TP =0;
	P->V_Liq = 0;
	P->m_Vap = 0;
	P->m_TP = 0;
	P->m_Liq = 0;	
	P->UA_Vap = 0;
	P->UA_TP = 0;
	P->UA_Liq = 0;


	switch(Pos) 
	{
		case 0: 
			if(Bra[i].Ini==1)//this branch has been initialized
			{Logic=1;}
			else
			{Logic=0;}
			break;
		
		case 1:
			if(Bra[i].Ini==1)//this branch has been initialized
			{Logic=1;}
			else
			{Logic=0;}
			break;
			
		default:
			{errorLog.Add("StructCond::Cal_HP","position not found");
			return -1;}
			break;
		}
	
	if(Logic)
	{
	*HPo=Bra[i].HPi;

	if(Bra[i].EqulNo<0)//no equivalent branch
	{
		for(int j=0;j<Bra[i].TubNum;j++)
		{
		const int TubeN=Bra[i].TubNo[j];

		Tub[TubeN].HPi=*HPo;
		Tub[TubeN].m.m=0;
		Tub[TubeN].m.V=0;

				if(Tub[TubeN].RowNo>0)//not the first row, to get the air state
				{
		
					const int Upper = Tub[TubeN].AirUpstreamUpper;
					const int Lower = Tub[TubeN].AirUpstreamLower;
					
					for(k=0;k<SegNum;k++)
					{
						if(Upper>=0&&Lower>=0)
						{
						hao=(Tub[Upper].Seg[k].hao*Tub[Upper].Ga+Tub[Lower].Seg[k].hao*Tub[Lower].Ga)
						/(Tub[Upper].Ga+Tub[Lower].Ga);}
						else if(Upper>=0)
						{
						hao=Tub[Upper].Seg[k].hao;}
						else
						{hao=Tub[Lower].Seg[k].hao;}

						Tub[TubeN].Seg[k].Tai = HatoTa(hao);
						if(errorLog.IsError()) {
						errorLog.Add("StructCond::Cal_HP","Tai");
						return -1;}
					}

				}//endif

		for(k=0;k<SegNum;k++)
		{
		
		{if(k==9&&TubeN==13) 
		int shenb=0;}//for debugging

		int realk=0;
		if(Tub[TubeN].even) realk=(SegNum-1-k);
		else realk=k;
		
		CondTubeL_new(Bra[i].Gr,HPo,Tub[TubeN].Ga,Tub[TubeN].Seg[realk].Tai,&Tub[TubeN].Seg[realk].hao,&mi,P); 
		
		if(errorLog.IsError()) {
		errorLog.Add("StructCond::Cal_HP","CondTube");
		return -1;}
		
		Tub[TubeN].m.m=Tub[TubeN].m.m+mi.m;
		Tub[TubeN].m.V=Tub[TubeN].m.V+mi.V;
		Sm.m=Sm.m+mi.m;
		Sm.V=Sm.V+mi.V;
		}//end k circle

		CondReturnBend(Bra[i].Gr,HPo,&mi,P);
		if(errorLog.IsError()) {
		errorLog.Add("StructCond::Cal_HP","CondTubeBend");
		return -1;}
		
		Tub[TubeN].m.m=Tub[TubeN].m.m+mi.m;
		Tub[TubeN].m.V=Tub[TubeN].m.V+mi.V;
		Sm.m=Sm.m+mi.m;
		Sm.V=Sm.V+mi.V;
		Tub[TubeN].HPo=*HPo;
		}//end j circle

		//output of this branch
		Bra[i].HPo=*HPo;
		Bra[i].m=Sm;
		Sm.m=0;
		Sm.V=0;
		Bra[i].Para_Struc[0]=P->VapL;
		Bra[i].Para_Struc[1]=P->TPL;
		Bra[i].Para_Struc[2]=P->LiqL;
		Bra[i].Para_Struc[3]=P->V_Vap;
		Bra[i].Para_Struc[4]=P->V_TP;
		Bra[i].Para_Struc[5]=P->V_Liq;
		Bra[i].Para_Struc[6]=P->m_Vap;
		Bra[i].Para_Struc[7]=P->m_TP;
		Bra[i].Para_Struc[8]=P->m_Liq;	
		Bra[i].Para_Struc[9]=P->UA_Vap;
		Bra[i].Para_Struc[10]=P->UA_TP;
		Bra[i].Para_Struc[11]=P->UA_Liq;
		}
	else
	{
		const int NoBra=Bra[i].EqulNo;
		Bra[i].HPo=Bra[NoBra].HPo;
		Bra[i].m=Bra[NoBra].m;	

		int NN=0;
		for(NN=0;NN<12;NN++)
		{
		Bra[i].Para_Struc[0]=Bra[NoBra].Para_Struc[0];
		Bra[i].Para_Struc[1]=Bra[NoBra].Para_Struc[1];
		Bra[i].Para_Struc[2]=Bra[NoBra].Para_Struc[2];
		Bra[i].Para_Struc[3]=Bra[NoBra].Para_Struc[3];
		Bra[i].Para_Struc[4]=Bra[NoBra].Para_Struc[4];
		Bra[i].Para_Struc[5]=Bra[NoBra].Para_Struc[5];
		Bra[i].Para_Struc[6]=Bra[NoBra].Para_Struc[6];
		Bra[i].Para_Struc[7]=Bra[NoBra].Para_Struc[7];
		Bra[i].Para_Struc[8]=Bra[NoBra].Para_Struc[8];
		Bra[i].Para_Struc[9]=Bra[NoBra].Para_Struc[9];
		Bra[i].Para_Struc[10]=Bra[NoBra].Para_Struc[10];
		Bra[i].Para_Struc[11]=Bra[NoBra].Para_Struc[11];
		}

	}//end else
	Bak.VapL=Bak.VapL+Bra[i].Para_Struc[0];
	Bak.TPL=Bak.TPL+Bra[i].Para_Struc[1];
	Bak.LiqL=Bak.LiqL+Bra[i].Para_Struc[2];
	Bak.V_Vap=Bak.V_Vap+Bra[i].Para_Struc[3];
	Bak.V_TP = Bak.V_TP+Bra[i].Para_Struc[4];
	Bak.V_Liq = Bak.V_Liq+Bra[i].Para_Struc[5];
	Bak.m_Vap = Bak.m_Vap+ Bra[i].Para_Struc[6];
	Bak.m_TP = Bak.m_TP+Bra[i].Para_Struc[7];
	Bak.m_Liq = Bak.m_Liq+Bra[i].Para_Struc[8];	
	Bak.UA_Vap = Bak.UA_Vap+Bra[i].Para_Struc[9];
	Bak.UA_TP = Bak.UA_TP+Bra[i].Para_Struc[10];
	Bak.UA_Liq = Bak.UA_Liq+Bra[i].Para_Struc[11];

	}//endif
	
	}//end i circle

	*P=Bak;
	
	return 0;
}

StructCond::StructCond(const char* fn)
{
	char scanstr[2048],*ptr;
	int i=0,j=0;
	int endSign;
	error=0;
	strncpy(filename,fn,120);
	FILE *fp = fopen(filename,"r");
	if(fp==NULL) {
		errorLog.Add("StructCond::StructCond","File not found");
		error=1;
		return;}
	
	fgets(scanstr,2048,fp);
	NodNum = (int)GetDouble(fp);    //Node Number
	BranNum=(int)GetDouble(fp);   //Branch Number
	RowNum=(int)GetDouble(fp);    //Row Number
	TubeNum=(int)GetDouble(fp);    //tube Number
	SegNum = (int)GetDouble(fp);   //Segment Number of per tube
	AreaFront=GetDouble(fp);    //frontal area
	Volum=GetDouble(fp);		//total nominal air mass flow rate [m^3/s]

	Nod = new CondNode[NodNum];
	Bra = new CondBranch[BranNum];
	GaRow = new double[RowNum];
	Tub=new TubeCond[TubeNum];

	fgets(scanstr,2048,fp);
	for(int i=0;i<NodNum;i++)
	{
	fgets(scanstr,2048,fp);
	ptr = scanstr;
	sget(&ptr,&Nod[i].NodNo);
	sget(&ptr,&Nod[i].InNum);
	sget(&ptr,&Nod[i].OutNum);
	Nod[i].BranIN=new int[Nod[i].InNum];
	Nod[i].BranOUT=new int[Nod[i].OutNum];
	for(int j=0;j<Nod[i].InNum;j++)
	{
	sget(&ptr,&Nod[i].BranIN[j]);
	}
	for(int j=0;j<Nod[i].OutNum;j++)
	{
	sget(&ptr,&Nod[i].BranOUT[j]);
	}
	}

	fgets(scanstr,2048,fp);
	for(int i=0;i<BranNum;i++)
	{
	fgets(scanstr,2048,fp);
	ptr = scanstr;
	sget(&ptr,&Bra[i].BranNo);
	sget(&ptr,&Bra[i].EqulNo);
	sget(&ptr,&Bra[i].GrFac);
	sget(&ptr,&Bra[i].TubNum);
	Bra[i].TubNo=new int[Bra[i].TubNum];
	for(int j=0;j<Bra[i].TubNum;j++)
	{
	sget(&ptr,&Bra[i].TubNo[j]);
	}
	}

	fgets(scanstr,2048,fp);
	fgets(scanstr,2048,fp);
	for(int i=0;i<TubeNum;i++)
	{
	Tub[i].Seg = new TubCndSeg[SegNum];
	fgets(scanstr,2048,fp);
	ptr = scanstr;
	sget(&ptr,&Tub[i].TubNo);
	sget(&ptr,&Tub[i].RowNo);
	sget(&ptr,&Tub[i].RefUpstream);
	sget(&ptr,&Tub[i].AirUpstreamUpper);
	sget(&ptr,&Tub[i].AirUpstreamLower);
	sget(&ptr,&Tub[i].GaFac);
	sget(&ptr,&Tub[i].even);
	}

	endSign=(int)GetDouble(fp);

	fclose(fp);
	if(endSign!=-256)
	{
	errorLog.Add("StructCond::StructCond","File Input Wrong");
	error=1;
	return;}

};


StructCond::~StructCond()
{
	int i=0;

	for(int i=0;i<TubeNum;i++)
	{
	delete Tub[i].Seg;
	}

	delete Tub;
	
	for(int i=0;i<BranNum;i++)
	{
	delete Bra[i].TubNo;
	}

	delete Bra;

	for(int i=0;i<NodNum;i++)
	{
	delete Nod[i].BranIN;
	delete Nod[i].BranOUT;
	}
	delete Nod;
};
//----------------------------------------------------------------------
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//the followings are old models
/***********************************************************************
Multi tube condenser circuitry pattern.  This is the pattern used
in the Carrier unit at Herrick.  It is one row deep, has a manifold
at the beginning that divides the hot-gas-line into "n1" parallel
paths each having "n2" passes across the face of the coil, then the
"n1" paths are combined into one for another "n3" passes across the
face of the coil.
***********************************************************************/
int CondCircuitMultipleTube(double mr,//refrigerant mass flow rate
							HP HPi,//refrigerant inlet state
							double Tai,//air inlet temperature
							double Ga,//air mass flux
							HP *HPo,//refrigerant outlet state
							double *Tao,//air outlet temperature
							MASS *m,//mass and volume in the condenser
							CGP *P)//condenser struct
{
	static int i,n,nhao;
	static double Gr;
	static double Shao,hao;
	static MASS mi;

	// initialize
	*HPo=HPi;


	Gr=mr/P->Ax;
	m->m=0;
	m->V=0;
	Shao=0;
	nhao=0;

	double SubL =0;//B.S. for keeping the liquid length in the condenser

	CondMan(1,P->Nbranchs,P->Di,&Gr,HPo);
	if(errorLog.IsError()) {
		errorLog.Add("CondCircuit201","CondMan(1)");
		return 0;
	}


	n=P->Nbranchs; /* number of sets */
	for(int i=0;i<P->Nmaintubes;i++) {

		CondTube(Gr,HPo,Ga,Tai,&hao,&mi,P);
		if(errorLog.IsError()) {
			errorLog.Add("CondCircuit201","CondTube(1)");
			return 0;
		}

		Shao+=n*hao;
		nhao+=(int)n;
		m->m+=n*mi.m;
		m->V+=n*mi.V;

		CondReturnBend(Gr,HPo,&mi,P);
		if(errorLog.IsError()) {
			errorLog.Add("CondCircuit201","CondReturnBend(1)");
			return 0;
		}

		m->m+=n*mi.m;
		m->V+=n*mi.V;
	}

	/* Combine both sets into one. */
	CondMan(P->Nbranchs,1,P->Di,&Gr,HPo);
	if(errorLog.IsError()) {
		errorLog.Add("CondCircuit201","CondMan(2)");
		return 0;
	}

	//B.S.-------------------------
	P->VapL = P->VapL*n;//B.S. vapor length in the paralelled circuits
	P->TPL = P->TPL*n;//B.S. two-phase length in the paralelled circuits
	P->LiqL = P->LiqL*n;// B.S. liquid length in the paralelled circuits
	

	P->V_Vap = P->V_Vap*n;//B.S. vapor volume in the paralelled circuits
	P->V_TP = P->V_TP*n;//B.S. two-phase volume in the paralelled circuits
	P->V_Liq = P->V_Liq*n;//B.S. liquid volume in in the paralelled circuits

	P->m_Vap = P->m_Vap*n;//vapor mass in the paralelled circuits
	P->m_TP = P->m_TP*n;//two-phase mass in the paralelled circuits
	P->m_Liq = P->m_Liq*n;//liquid mass in the paralelled circuits
	
	P->UA_Vap = P->UA_Vap*n;//vapor heat conductance in the paralelled circuits
	P->UA_TP = P->UA_TP*n;//two-phase heat conductance in the paralelled circuits
	P->UA_Liq = P->UA_Liq*n;//liquid heat conductance in the paralelled circuits
	//----------------------B.S.



	n=1; /* number of sets */
	for(int i=0;i<P->Nsubtubes;i++) {

		CondTube(Gr,HPo,Ga,Tai,&hao,&mi,P);
		if(errorLog.IsError()) {
			errorLog.Add("CondCircuit201","CondTube(2)");
			return 0;
		}

		Shao+=n*hao;
		nhao+=(int)n;
		m->m+=n*mi.m;
		m->V+=n*mi.V;

		CondReturnBend(Gr,HPo,&mi,P);
		if(errorLog.IsError()) {
			errorLog.Add("CondCircuit201","CondReturnBend(2)");
			return 0;
		}

		m->m+=n*mi.m;
		m->V+=n*mi.V;
	}

	// air outlet state
	hao = Shao/(double)nhao;
	*Tao = HatoTa(hao);
	if(errorLog.IsError()) {
		errorLog.Add("CondCircuit201","Tao");
		return 0;
	}

	SubL= P->LiqL;//B.S. the liquid length in the subcooled site

//B.S. add for accouting the unmeasured volume

		TXP TXP_left = HPtoTXP(*HPo);//B.S., get the outlet state of the condenser
		if(errorLog.IsError()) {
			errorLog.Add("CondCircuit201","TXP_left");
			return 0;
			}

		const double V_left=-0.00045;//B.S., the unaccounted volume at the subcooled site, which is got by tuning procedure
		const double Volume_left= PropertyTXPth(VOL,TXP_left);//B.S., refrigerant specific volume at the unaccounted subcooled region
		if(errorLog.IsError()) {
			errorLog.Add("CondCircuit201","Volume_left");
			return 0;
			}
		const double rho_left=1/Volume_left;

		const double m_left=V_left/Volume_left;//B.S., refrigerant mass in the unaccounted volume

		m->m+=0;//m_left;//unaccounted charge mass in the unknown volume -0.5

	return 0;
}


/***********************************************************************
Single tube condenser circuitry pattern.  This is the simpliest pattern
the only has one tube.
***********************************************************************/
int CondCircuitSingleTube(double mr,HP HPi,double Tai,double Ga,HP *HPo,double *Tao,MASS *m,CGP *P)
{
	*HPo = HPi;
	double Gr = mr/P->Ax;

	double hao;
	CondTube(Gr,HPo,Ga,Tai,&hao,m,P);
	if(errorLog.IsError()) {
		errorLog.Add("CondCircuit202");
		return 0;
	}

	// air outlet state
	*Tao = HatoTa(hao);
	if(errorLog.IsError()) {
		errorLog.Add("CondCircuit202");
		return 0;
	}

	return 0;
}


//B.S.--------------------------------------------------------
//for user defined---Harms 5-ton

int CondCircuitMultipleTube_204(double mr,HP HPi,double Tai,double Ga,HP *HPo,double *Tao,MASS *m,CGP *P)
{

	static int I;
	static double Gr;
	static double hao;
	static MASS mi;
	double P_def[65+1],H[65+1],T_I[64+1],T_O[64+1];
	HP HP_61;// point 61 is where the tubes converge
	TXP TXP_61;
   
	// initialize
	*HPo=HPi;
	m->m=0;
	m->V=0;
	hao=0;
	P_def[1]=HPi.P;
	H[1]=HPi.H;

	double R_1 = 0;
	double R_2 = 0;
	double M_DOT_RT = mr;
	
	const double ROWS = 32; //row number
	const double AA=0.0254e0;//tube space
	const double D_O = 9.91e-3;// tube outside diameter
	const double FIN_S = 1.15e-3;//fin spacing
	const double FIN_T = 0.114e-3;//fin thickness
	const double L_COND = 1.62865e0;   //  !tube length
	const double M_DOT_AT=2.060;
	const double  FACE_A= (ROWS*AA)*L_COND;//facing area
//	Ga = M_DOT_AT/(FACE_A)*(AA*FIN_S/((FIN_S-FIN_T)*(AA-D_O)));//shenbo guess, related with input
	
	const double T_AMB = Tai;//saving the inlet air temperature

	for(int i=1;I<=64;I++){
		if (R_1 < 2){
			T_I[I]=T_AMB;
			}
		else{
			R_2=R_2+1;
			T_I[I]=T_O[I-2];
			}

		if (I < 61){
			mr=M_DOT_RT/6e0;}
		else{
			mr=M_DOT_RT/1e0;}

			//heat transfer and pressure drop calculation
			{
			Gr=mr/P->Ax;
			HPo->H = H[I]; HPo->P = P_def[I];
			Tai = T_I[I];

			CondTube(Gr,HPo,Ga,Tai,&hao,&mi,P);
			if(errorLog.IsError()) {
			errorLog.Add("CondCircuit201","CondTube(1)");
			return 0;
			}

			m->m+=mi.m;
			m->V+=mi.V;

			CondReturnBend(Gr,HPo,&mi,P);
			if(errorLog.IsError()) {
			errorLog.Add("CondCircuit201","CondReturnBend(1)");
			return 0;
			}

			m->m+=mi.m;
			m->V+=mi.V;	

			hao = hao;
			*Tao = HatoTa(hao);
			if(errorLog.IsError()) {
			errorLog.Add("CondCircuit201","Tao");
			return 0;
			}
			H[I] = HPo->H; P_def[I] = HPo->P;
			T_O[I] = *Tao;
			}

		
			if (I/10e0==int(I/10e0)){
				P_def[I+1]=P_def[1];
				H[I+1]=H[1];}
			else{
				P_def[I+1]=P_def[I];
				H[I+1]=H[I];}

			if (I==60){
			P_def[61]=(P_def[10]+P_def[20]+P_def[30]+P_def[40]+P_def[50]+P_def[60])/6e0;   // !Not Physical
			H[61]=(H[10]+H[20]+H[30]+H[40]+H[50]+H[60])/6e0;
		
			HP_61.H = H[61];
			HP_61.P = P_def[61];

			}

			R_1=R_1+1;
			if (R_2==2){
			R_1=0;
			R_2=0;
			}
	}
	/*Global Balance*/

	const double P_OUT=P_def[64];
	const double H_OUT=H[64];
	
	/*Mass Balance*/
	const double PI=4e0*atan(1e0);
	const double V_A = PI*(pow(0.00864e0,2e0)*0.667e0+pow(0.00395,2e0)*0.933e0);
	
	TXP_61 = HPtoTXP(HP_61);
	if(errorLog.IsError()) {
		errorLog.Add("CondCircuit201","HP_61->TXP_61");
		return 0;
		}
	const double Volume_61= PropertyTXPth(VOL,TXP_61);//VolumeALL(TXP_left,Gr,P->Di);
	
	if(errorLog.IsError()) {
		errorLog.Add("CondCircuit201","Volume_61");
		return 0;
		}

	const double MA = V_A/(Volume_61);
	const double MT = m->m+MA;
	
	m->m = MT;
	m->V = m->V;
	HPo->H=H_OUT;
	HPo->P=P_OUT;
	*Tao =T_O[63];

//B.S. add for accouting the unmeasured volume

		TXP TXP_left = HPtoTXP(*HPo);//B.S., get the outlet state of the condenser
		if(errorLog.IsError()) {
			errorLog.Add("CondCircuit201","TXP_left");
			return 0;
			}

		const double V_left=-0.00045;//B.S., the unaccounted volume at the subcooled site, which is got by tuning procedure
		const double Volume_left= PropertyTXPth(VOL,TXP_left);//B.S., refrigerant specific volume at the unaccounted subcooled region
		if(errorLog.IsError()) {
			errorLog.Add("CondCircuit201","Volume_left");
			return 0;
			}
		const double rho_left=1/Volume_left;

		const double m_left=V_left/Volume_left;//B.S., refrigerant mass in the unaccounted volume

		m->m+=0.64;//m_left;

	return 0;
}

