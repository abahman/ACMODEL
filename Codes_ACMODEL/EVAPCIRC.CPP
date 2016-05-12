#include "stdafx.h"

#include <stdio.h>
#include "r22.h"
#include "wair.h"
#include "mass.h"
#include "evap.h"
#include "volume.h"
#include "errorlog.h"

#include <string.h>//B.S.
#include "getdata.h"//B.S.
#include "parse.h"//B.S.
#include "spline2d.h"//B.S.
#include "spline1d.h"//B.S.
#include <math.h>

extern ErrorLog errorLog;
extern WetAirPropTR wair;
extern AirProp air;

struct EvapNode
{
	int NodNo;//node number
	int InNum;//number of tubes flowing into the node
	int OutNum;//number of the tubes flowing out of the node
	int *BranIN;//index of the tube branch flowing in
	int *BranOUT;//index of the tube branch flowing out
};

struct TubEvpSeg
{
	TP TPi;
	WH WHo;
};

struct TubeEvap
{
	int TubNo;//tube No
	int RowNo;//the row no where the tube located
	int Refdownstream;//downstream tube No at the refrigerant side (flow direction)
	int AirUpstreamUpper;//Upstream tube No at the air upper side
	int AirUpstreamLower;//Upstream tube No at the air lower side
	int even;//tube flow direction
	double GaFac;//air flow distribution factor
	double Ga;//maximum air flow rate in the segment
	HP HPi;//refrigerant inlet state, flowing direction
	HP HPo;//refrigerant outlet state, flowing direction
	MASS m;//mass and volume of the tube
	TubEvpSeg *Seg;
};

struct EvapBranch
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

class StructEvap {
public:
	StructEvap(const char* fn);//reading the data file in
	~StructEvap();//delete the memory

	int _EvapCircuit_Fwd(double mr,HP *HPo,double Ga,TP TPi,HP *HPi,TP *TPo,MASS *Sm,double *Aflow,ETdim *P);
	int _EvapCircuit_Rev(double mr,HP *HPo,double Ga,TP TPi,HP *HPi,TP *TPo,MASS *Sm,double *Aflow,ETdim *P);//reversed direction
	int	Cal_HP(ETdim *P, int Pos,HP *HPi);//function to calculate the heat transfer and pressure drop
	char filename[128];//data file name
	int error;//for keeping track of the errorLog
	int TubeNum;	//number of all tubes
	int BranNum; // number of branch
	int NodNum;  // number of node
	int RowNum;	 //number of Row
	int SegNum; //number of segment
	int Rev;
	double AreaFront;//front area
	double Volum;//air flow volume
	double *GaRow;//average air flow distribution factor of each row
	EvapBranch *Bra;//array of branch objects
	EvapNode *Nod;//array of node objects
	TubeEvap *Tub;//array of tube objects
	double ChargTunFix;//for fixed amount charge tuning, not used for evaporator, it will be used if the HX is used as condenser
	double Vap_Ref;// reference vapor length for two-point charge tuning, not used for evaporator, it will be used if the HX is used as condenser
	double K_ChargTun;//coeffcient relate charge modeling errors to vapor length, K_ChargTun = 0 means fixed amount charge tuning; not used for evaporator, it will be used if the HX is used as condenser
	};



// local functions
int EvapCircuit301(double,HP*,double,TP,HP*,TP*,MASS*,double*,ETdim*);
int EvapCircuit302(double,HP*,double,TP,HP*,TP*,MASS*,double*,ETdim*);

int EvapCircuit303(double,HP*,double,TP,HP*,TP*,MASS*,double*,ETdim*);//Shenbo add for modeling CS30TR

/********************************************************************
********************************************************************/
int EvapCircuit(int type,double mr,HP *HPo,double GaI,TP TPi,HP *HPi,TP *TPo,MASS *Sm,double *Aflow,ETdim *D)
{
	static StructEvap EvapInfo("./InputDoc/EvapStruc.in");
	EvapInfo.Rev=D->REV;

	switch(type) {
		case 301://user defined
		EvapInfo._EvapCircuit_Fwd(mr,HPo,GaI,TPi,HPi,TPo,Sm,Aflow,D);
		if(errorLog.IsError()) {
			errorLog.Add("EvapCircuit","EvapCircuit301");
			return 0;
			}
		break;
		case 302: // Carrier RTU at Rurdue
			EvapCircuit301(mr,HPo,GaI,TPi,HPi,TPo,Sm,Aflow,D); 
			
			if(errorLog.IsError()) {
				errorLog.Add("EvapCircuit","EvapCircuit301");
				return 0;
			}
			break;
		case 303: // Single finned tube
			EvapCircuit302(mr,HPo,GaI,TPi,HPi,TPo,Sm,Aflow,D);

			if(errorLog.IsError()) {
				errorLog.Add("EvapCircuit","EvapCircuit302");
				return 0;
			}
			break;
		case 304:
			EvapCircuit303(mr,HPo,GaI,TPi,HPi,TPo,Sm,Aflow,D); // this function is for CS30TR

			if(errorLog.IsError()) {
				errorLog.Add("EvapCircuit","EvapCircuit301");
				return 0;
			}
			break;

		default:
			{
			char str[128];
			sprintf(str,"Evaporator model %d not found",type);
			errorLog.Add("EvapCircuit",str);
		}
		break;
	}

	return 0;
}

//------------------------------------------------------------------------
//-------------------------------------------------------------------------
//------------------------------------------------------------------------
//user defined
int StructEvap::_EvapCircuit_Fwd(double mr,HP *HPo,double Ga,TP TPi,HP *HPi,TP *TPo,MASS *Sm,double *Aflow,ETdim *P)
{
	int i=0,j=0,k=0;
	double Gr=0;//mass flux
	double H_in=0;
	
	if(error)
	{
	errorLog.Add("StructEvap::_EvapCircuit","Class fails to be constructed");
	return -1;}
	
	if(0/*Rev*/)//reversed calculation
	{_EvapCircuit_Rev(mr,HPo,Ga,TPi,HPi,TPo,Sm,Aflow,P);
	return 0;}

	*TPo=TPi;//inlet air state

	H_in=HPi->H;//inlet refrigerant enthalpy
	
	//air side initializing
	const double rho_air=1/air.v(TPi.T);
	if(errorLog.IsError()) {
	errorLog.Add("StructEvap::_EvapCircuit","rho_air");
	return -1;}

	//const double Ma =Ga*rho_air*4.719e-4; 
	Ga=Ga/AreaFront*P->vsp*P->Ls/((P->vsp-P->Do)*(P->Ls-P->N*P->th));


	for(i=0;i<RowNum;i++)
	{GaRow[i]=0;}
		
	for(j=0;j<RowNum;j++)
	{
		int N=0;
		for(i=0;i<TubeNum;i++)
		{

		if(Tub[i].RowNo==j)
		{
			for(k=0;k<SegNum;k++)
			{
			Tub[i].Seg[k].TPi=TPi;
					
			Tub[i].Seg[k].WHo.W = wair.HumidityRatio(TPi.T,TPi.P);
			Tub[i].Seg[k].WHo.H = wair.h(TPi.T,TPi.P);
			if(errorLog.IsError()) {
				errorLog.Add("StructEvap::_EvapCircuit","WHo");
				return -1;}

			}

			if(Tub[i].RowNo>0)//not the first row
			{
				const int Upper= Tub[i].AirUpstreamUpper;
				const int Lower = Tub[i].AirUpstreamLower;
				if(Upper>=0&&Lower>=0)
				{Tub[i].GaFac=(Tub[Upper].GaFac+Tub[Lower].GaFac)/2;}
				else if(Upper>=0)
				{Tub[i].GaFac=Tub[Upper].GaFac;}
				else
				{Tub[i].GaFac=Tub[Lower].GaFac;}
			}
		GaRow[j]=GaRow[j]+Tub[i].GaFac;
		N=N+1;
		}//ifend	
		}//i circle end
		GaRow[j]=GaRow[j]/N;
	}//j circle end
	
		for(i=0;i<TubeNum;i++)
		{const int RowN=Tub[i].RowNo;
		Tub[i].Ga=Tub[i].GaFac/GaRow[RowN]*Ga;}
	
		double H_out=0;//intermidiate evaporator exit enthalpy
		double Res=0;
		double DP = 20;
		double P_in=0,P_out=0;
		HPi->P=HPo->P+DP;//temperary inlet refrigerant pressure
		HPi->H=H_in;//evaporator inlet refrigerant enthalpy
		//*HPi is the intermediate variable
		P_in = HPo->P+DP;//temperary inlet refrigerant pressure
		P_out = HPo->P;//evaporator outlet refrigerant pressure
		int IterN=0;


		//refrigerant state intialize first time
	for(i=0;i<NodNum;i++)//inlet nodes
	{
		if(Nod[i].BranIN[0]<0)//no inlet branch, only from the distributor
		{
			Gr = mr/(P->Ax*Nod[i].OutNum);	
			if(errorLog.IsError()) {
			errorLog.Add("StructCond::_CondCircuit","CondMan0");
			return -1;}

			for(j=0;j<Nod[i].OutNum;j++)//states flowing out from the node
			{
				const int jj = Nod[i].BranOUT[j];//index of the outlet branches
				Bra[jj].HPi=*HPi;
				Bra[jj].HPo=*HPo;
				Bra[jj].Gr=Gr*Bra[jj].GrFac;
				Bra[jj].Ini=1;
			}//end j circle

		}//endif
	}//end i circle

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
		HPi->H=H_in;//evaporator inlet refrigerant enthalpy
		HPi->P = P_in;

		for(i=0;i<BranNum;i++)
		{Bra[i].Ini=0;}
	
		for(i=0;i<NodNum;i++)//inlet nodes
		{
		if(Nod[i].BranIN[0]<0)//from the distributor
		{
			Gr = mr/(P->Ax*Nod[i].OutNum);	
	
			for(j=0;j<Nod[i].OutNum;j++)//states flowing out from the node
			{
				const int jj = Nod[i].BranOUT[j];//index of the outlet branches
				Bra[jj].HPi.H=HPi->H;
	
				if(Nod[i].OutNum==BranNum)//no middle nodes
				{
				const double DDP = (Bra[jj].HPi.P-Bra[jj].HPo.P);
				Bra[jj].HPi.P=HPo->P+DDP;
				P->DPr[jj]=DDP;//output the pressure drop of each branch
				}
				else
				{Bra[jj].HPi.P=HPi->P;}
				Bra[jj].Gr=Gr*Bra[jj].GrFac;
				Bra[jj].Ini=1;
			}//end j circle

		}//endif
		}//end i circle
	
	Cal_HP(P,0,HPi);//heat transfer and pressure drop calculation

	if(errorLog.IsError()) {
	errorLog.Add("StructEvap::_EvapCircuit","Cal_HP(P,0)");
	return -1;}
	
	for(i=0;i<BranNum;i++)
	{Bra[i].Ini=0;}
	

	Gr=0;
	HPi->H=0;HPi->P=0;


	//nodes in the middle
	for(i=0;i<NodNum;i++)
	{
		if(Nod[i].BranIN[0]>=0&&Nod[i].BranOUT[0]>=0)//nodes in the middle
		{

			for(j=0;j<Nod[i].InNum;j++)//node inlet state
			{
			const int jj= Nod[i].BranIN[j];
			Gr=Gr+Bra[jj].Gr;
			HPi->H=Bra[jj].HPo.H*Bra[jj].Gr+HPi->H;
			HPi->P=Bra[jj].HPo.P*Bra[jj].Gr+HPi->P;
			}

			HPi->H=HPi->H/Gr;
			HPi->P=HPi->P/Gr;
			Gr=Gr/Nod[i].InNum;
			Gr=Gr*Nod[i].InNum/Nod[i].OutNum;

			for(j=0;j<Nod[i].OutNum;j++)
			{
			const int jj = Nod[i].BranOUT[j];//index of outlet branches
			Bra[jj].HPi=*HPi;
			Bra[jj].Gr=Gr*Bra[jj].GrFac;
			Bra[jj].Ini=1;
			}//end j circle

			}//endif
	}//end i circle


	Cal_HP(P,1,HPi);//heat transfer and pressure drop calculation
	if(errorLog.IsError()) {
	errorLog.Add("StructEvap::_EvapCircuit","Cal_HP(P,1)");
	return -1;}
	
	Gr=0;
	HPi->H=0;HPi->P=0;

	//exit nodes
	for(i=0;i<NodNum;i++)
	{
		if(Nod[i].BranOUT[0]<0)//no outlet branch except the suction line
		{
			for(j=0;j<Nod[i].InNum;j++)
			{
			const int jj= Nod[i].BranIN[j];
			Gr=Gr+Bra[jj].Gr;
			HPi->H=Bra[jj].HPo.H*Bra[jj].Gr+HPi->H;
			HPi->P=Bra[jj].HPo.P*Bra[jj].Gr+HPi->P;
			P->Hout8[jj] = Bra[jj].HPo.H;//output the exit enthalpy of each branch
			}
			HPi->H=HPi->H/Gr;
			HPi->P=HPi->P/Gr;
			Gr=Gr/Nod[i].InNum;
			Gr=Gr*Nod[i].InNum/Nod[i].OutNum;

		}//endif

	}//end i circle
	
	
	Res=2*(HPi->H-H_out)/(HPi->H+H_out);
	H_out=HPi->H;
	P_out=HPi->P;
	P_in=HPo->P+(P_in-P_out);

	}while(fabs(Res)>1e-7&&IterN<20);

	if(fabs(Res)>1e-4) 
	{
	errorLog.Add("StructEvap::_EvapCircuit","can't reach the required tolerance");
	return -1;}

	Sm->m=0;Sm->V=0;

	for(i=0;i<BranNum;i++)
	{
	Sm->m=Sm->m+Bra[i].m.m;
	Sm->V=Sm->V+Bra[i].m.V;}

	WH WH_out={0,0};
	double Ma_out=0;

	//air outputs
	for(i=0;i<BranNum;i++)
	{
		if(Bra[i].EqulNo<0)//no equivalent circuits
		{
			for(j=0;j<Bra[i].TubNum;j++)
			{
				int Tubj=0;
				if(Rev)//counted from the point connected to compressor suction
				{Tubj=j;}
				else
				{Tubj=Bra[i].TubNum-1-j;}//counted from the point of evaporator entrance
				const int TubeN=Bra[i].TubNo[Tubj];

				if(Tub[TubeN].RowNo==RowNum-1)//outlet row tube at air side
				{
					for(k=0;k<SegNum;k++)
					{
					WH_out.W=Tub[TubeN].Ga*Tub[TubeN].Seg[k].WHo.W+WH_out.W;
					WH_out.H=Tub[TubeN].Ga*Tub[TubeN].Seg[k].WHo.H+WH_out.H;
					Ma_out=Ma_out+Tub[TubeN].Ga;
					}
				}//if end
			}
		}
	}

	WH_out.W = WH_out.W/Ma_out;
	WH_out.H = WH_out.H/Ma_out;

	TPo->T=TPi.T-5;
	TPo->P=0.8;

	*TPo=WHtoTP(WH_out,*TPo);
	if(errorLog.IsError()) {
	errorLog.Add("StructEvap::_EvapCircuit","WH_out to TPo");
	return -1;}

	//refrigerant outputs
	HPo->H=H_out;
	HPo->P=HPo->P;
	HPi->H=H_in;
	HPi->P = P_in;

	return 0;
};

//user defined
int StructEvap::_EvapCircuit_Rev(double mr,HP *HPo,double Ga,TP TPi,HP *HPi,TP *TPo,MASS *Sm,double *Aflow,ETdim *P)
{
	int i=0,j=0,k=0;
	double Gr=0;//mass flux

	if(error)
	{
	errorLog.Add("StructEvap::_EvapCircuit","Class fails to be constructed");
	return -1;}
	
	Rev=1;
	//air side initializing
	const double rho_air=1/air.v(TPi.T);
	if(errorLog.IsError()) {
	errorLog.Add("StructEvap::_EvapCircuit","rho_air");
	return -1;}

	const double Ma =Volum*rho_air; 
	const double GaNom=Ma/AreaFront*P->vsp*P->Ls/((P->vsp-P->Do)*(P->Ls-P->N*P->th));
	if(Ga<0) Ga= GaNom;

	for(i=0;i<RowNum;i++)
	{GaRow[i]=0;}
	
	
	for(j=0;j<RowNum;j++)
	{
		int N=0;
		for(i=0;i<TubeNum;i++)
		{

			if(Tub[i].RowNo==j)
			{
				for(k=0;k<SegNum;k++)
				{
				Tub[i].Seg[k].TPi=TPi;
				Tub[i].Seg[k].WHo.W = wair.HumidityRatio(TPi.T,TPi.P);
				Tub[i].Seg[k].WHo.H = wair.h(TPi.T,TPi.P);
				if(errorLog.IsError()) {
				errorLog.Add("StructEvap::_EvapCircuit","WHo");
				return -1;}
				}

				if(Tub[i].RowNo>0)//not the first row
				{
					const int Upper= Tub[i].AirUpstreamUpper;
					const int Lower = Tub[i].AirUpstreamLower;
					if(Upper>=0&&Lower>=0)
					{
					Tub[i].GaFac=(Tub[Upper].GaFac+Tub[Lower].GaFac)/2;}
					else if(Upper>=0)
					{Tub[i].GaFac=Tub[Upper].GaFac;}
					else
					{Tub[i].GaFac=Tub[Lower].GaFac;}
				}
				GaRow[j]=GaRow[j]+Tub[i].GaFac;
				N=N+1;
			}//ifend	
		}//i circle end
		GaRow[j]=GaRow[j]/N;
	}//j circle end
	
	for(i=0;i<TubeNum;i++)
	{
	const int RowN=Tub[i].RowNo;
	Tub[i].Ga=Tub[i].GaFac/GaRow[RowN]*Ga;}
	

	// Refrigerant side initialize

	double H_in=0;
	double Res=0;
	double P_in=0;
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
		*HPi=*HPo;
	
		for(i=0;i<BranNum;i++)
		{Bra[i].Ini=0;}

		for(i=0;i<NodNum;i++)//Outlet nodes
		{
			if(Nod[i].BranOUT[0]<0)//to compressor suction
			{
			Gr = mr/(P->Ax*Nod[i].InNum);
	
			for(j=0;j<Nod[i].InNum;j++)//states flowing into the node
			{
				const int jj = Nod[i].BranIN[j];//index of the inlet branches
				Bra[jj].HPo=*HPi;
				Bra[jj].Gr=Gr*Bra[jj].GrFac;
				Bra[jj].Ini=1;
			}//end j circle

			}//endif
		}//end i circle

	Cal_HP(P,0,HPi);//heat transfer and pressure drop calculation

	if(errorLog.IsError()) {
	errorLog.Add("StructEvap::_EvapCircuit","Cal_HP(P,0)");
	return -1;}
	
	for(i=0;i<BranNum;i++)
	{Bra[i].Ini=0;}
	

	Gr=0;
	HPi->H=0;HPi->P=0;


	//nodes in the middle
	for(i=0;i<NodNum;i++)
	{
		if(Nod[i].BranOUT[0]>=0&&Nod[i].BranIN[0]>=0)
		{

		for(j=0;j<Nod[i].OutNum;j++)//node outlet state
		{
			const int jj= Nod[i].BranOUT[j];
			Gr=Gr+Bra[jj].Gr;
			HPi->H=Bra[jj].HPi.H*Bra[jj].Gr+HPi->H;
			HPi->P=Bra[jj].HPi.P*Bra[jj].Gr+HPi->P;
		}

		HPi->H=HPi->H/Gr;
		HPi->P=HPi->P/Gr;
		Gr=Gr/Nod[i].OutNum;
	
		Gr=Gr*Nod[i].OutNum/Nod[i].InNum;
	
		for(j=0;j<Nod[i].InNum;j++)
		{
			const int jj = Nod[i].BranIN[j];//index of outlet branches
			Bra[jj].HPo=*HPi;
			Bra[jj].Gr=Gr*Bra[jj].GrFac;
			Bra[jj].Ini=1;
		}//end j circle

		}//endif
	}//end i circle

	Cal_HP(P,1,HPi);//heat transfer and pressure drop calculation

	if(errorLog.IsError()) {
	errorLog.Add("StructEvap::_EvapCircuit","Cal_HP(P,1)");
	return -1;}
	
	Gr=0;
	HPi->H=0;HPi->P=0;

	//inlet nodes
	for(i=0;i<NodNum;i++)
	{
		if(Nod[i].BranIN[0]<0)//no inlet branch except the distribution line
		{
			for(j=0;j<Nod[i].OutNum;j++)
			{
			const int jj= Nod[i].BranOUT[j];
			Gr=Gr+Bra[jj].Gr;
			HPi->H=Bra[jj].HPi.H*Bra[jj].Gr+HPi->H;
			HPi->P=Bra[jj].HPi.P*Bra[jj].Gr+HPi->P;
			}
		HPi->H=HPi->H/Gr;
		HPi->P=HPi->P/Gr;
		Gr=Gr/Nod[i].OutNum;
		Gr=Gr*Nod[i].OutNum/Nod[i].InNum;
		}//endif

	}//end i circle
		
	if(RowNum==1) Res=0;
	else
	{Res=2*(HPi->H-H_in)/(HPi->H+H_in);}
	H_in=HPi->H;
	
	}while(fabs(Res)>1e-7&&IterN<20);

	if(fabs(Res)>1e-5) 
	{
	errorLog.Add("StructEvap::_EvapCircuit","can't reach the required tolerance");
	return -1;}

	Sm->m=0;Sm->V=0;
	for(i=0;i<BranNum;i++)
	{
	Sm->m=Sm->m+Bra[i].m.m;
	Sm->V=Sm->V+Bra[i].m.V;}

	WH WH_out={0,0};
	double Ma_out=0;

	//air outputs
	for(i=0;i<TubeNum;i++)
	{
	if(Tub[i].RowNo==RowNum-1)//airside outlet row
	{
		for(j=0;j<SegNum;j++)
		{
			WH_out.W=Tub[i].Ga*Tub[i].Seg[j].WHo.W+WH_out.W;
			WH_out.H=Tub[i].Ga*Tub[i].Seg[j].WHo.H+WH_out.H;
			Ma_out=Ma_out+Tub[i].Ga;
		}
	}
	}

	WH_out.W = WH_out.W/Ma_out;
	WH_out.H = WH_out.H/Ma_out;

	TPo->T=TPi.T-5;
	TPo->P=0.8;
	*TPo=WHtoTP(WH_out,*TPo);
	if(errorLog.IsError()) {
	errorLog.Add("StructEvap::_EvapCircuit","WH_out to TPo");
	return -1;}

	return 0;
};



int StructEvap::Cal_HP(ETdim *P, int Pos,HP *HPi)
{
	//heat transfer and pressure drop of the two-phase region
	int i=0,j=0,k=0;
	int Logic=0;
	WH WHo;
	MASS Sm;
	MASS mi;
	ETdim Bak;//keeping the information of the condenser

	Sm.m=0; Sm.V=0;
	Bak=*P;
	
	if(Pos<0||Pos>1)
	{
	errorLog.Add("StructEvap::Cal_HP","Wrong position");
	return -1;}

	for(i=0;i<BranNum;i++)
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

	switch(Pos) {
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
			{errorLog.Add("StructEvap::Cal_HP","Mode not found");
			return -1;}
			break;
		}//switch end
	
	if(Logic)
	{
	if(Rev)	*HPi=Bra[i].HPo;//opposite to flow direction
	else *HPi=Bra[i].HPi;//paralell to flow direction

	if(Bra[i].EqulNo<0)//no equivalent branch
	{
		for(j=0;j<Bra[i].TubNum;j++)
		{
		
			int Tubj=0;
			if(Rev)//counted from the point connected to compressor suction
			{Tubj=j;}
			else
			{Tubj=Bra[i].TubNum-1-j;}//counted from the point of evaporator entrance

		const int TubeN=Bra[i].TubNo[Tubj];

		Tub[TubeN].HPo=*HPi;
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
					WHo.W=(Tub[Upper].Seg[k].WHo.W*Tub[Upper].Ga+Tub[Lower].Seg[k].WHo.W*Tub[Lower].Ga)
						/(Tub[Upper].Ga+Tub[Lower].Ga);
					WHo.H=(Tub[Upper].Seg[k].WHo.H*Tub[Upper].Ga+Tub[Lower].Seg[k].WHo.H*Tub[Lower].Ga)
						/(Tub[Upper].Ga+Tub[Lower].Ga);
					}
					else if(Upper>=0)
					{
					WHo=Tub[Upper].Seg[k].WHo;}
					else
					{WHo=Tub[Lower].Seg[k].WHo;}

					Tub[TubeN].Seg[k].TPi = WHtoTP(WHo,Tub[TubeN].Seg[k].TPi);
					if(errorLog.IsError()) {
					errorLog.Add("StructEvap::Cal_HP","Tub[TubeN].TPi");
					return -1;}
				}//end k circle
		}//endif
 
		for(k=0;k<SegNum;k++)
		{
			{if(k==9&&TubeN==13)	
			int shenb=0;}//for debugging

			int realk=0;
			if(Tub[TubeN].even) realk=(SegNum-1-k);
			else realk=k;
		
			if(Rev)
			{EvapTubeL_Rev(Bra[i].Gr,HPi,Tub[TubeN].Ga,Tub[TubeN].Seg[realk].TPi,&WHo,&mi,P);}
			else
			{EvapTubeL_Fwd(Bra[i].Gr,HPi,Tub[TubeN].Ga,Tub[TubeN].Seg[realk].TPi,&WHo,&mi,P);}

			if(errorLog.IsError()) {
			errorLog.Add("StructEvap::Cal_HP","EvapTube");
			return -1;}
		Tub[TubeN].Seg[realk].WHo=WHo;
		Tub[TubeN].m.m=Tub[TubeN].m.m+mi.m;
		Tub[TubeN].m.V=Tub[TubeN].m.V+mi.V;
		Sm.m=Sm.m+mi.m;
		Sm.V=Sm.V+mi.V;
		}//end k circle
		
		
		
		if(Rev)
		{EvapTubeBend(Bra[i].Gr,HPi,&mi,P);}
		else
		{EvapTubeBend_Fwd(Bra[i].Gr,HPi,&mi,P);}
		if(errorLog.IsError()) {
		errorLog.Add("StructCond::Cal_HP","CondTubeBend");
		return -1;}

		Tub[TubeN].m.m=Tub[TubeN].m.m+mi.m;
		Tub[TubeN].m.V=Tub[TubeN].m.V+mi.V;
		Sm.m=Sm.m+mi.m;
		Sm.V=Sm.V+mi.V;
		if(Rev)	Tub[TubeN].HPi=*HPi;
		else Tub[TubeN].HPo=*HPi;
		}//end j circle

		//output of this branch
		if(Rev) Bra[i].HPi=*HPi;
		else Bra[i].HPo=*HPi;
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
	else//there is an equivalent branch for it
	{
		const int NoBra=Bra[i].EqulNo;
		
		if(Rev) Bra[i].HPi=Bra[NoBra].HPi;
		else Bra[i].HPo=Bra[NoBra].HPo;

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


StructEvap::StructEvap(const char* fn)
{
	char scanstr[2048],*ptr;
	int i=0,j=0;
	int endSign;
	error=0;
	strncpy(filename,fn,120);
	FILE *fp = fopen(filename,"r");
	if(fp==NULL) {
		errorLog.Add("StructEvap::StructEvap","File not found");
		error=1;
		return;}
	
	fgets(scanstr,2048,fp);
	NodNum = (int)GetDouble(fp);    //in Node Number
	BranNum=(int)GetDouble(fp);   //Branch Number
	RowNum=(int)GetDouble(fp);    //Row Number
	TubeNum=(int)GetDouble(fp);    //overall tube Number
	SegNum = (int)GetDouble(fp);    //Segment Number of per tube
	AreaFront=GetDouble(fp);    //frontal area
	Volum=GetDouble(fp);		//total nominal air mass flow rate [m^3/s]

	Nod = new EvapNode[NodNum];
	Bra = new EvapBranch[BranNum];
	GaRow = new double[RowNum];
	Tub=new TubeEvap[TubeNum];

	fgets(scanstr,2048,fp);
	for(i=0;i<NodNum;i++)
	{
	fgets(scanstr,2048,fp);
	ptr = scanstr;
	sget(&ptr,&Nod[i].NodNo);
	sget(&ptr,&Nod[i].OutNum);//get outlet branches at first
	sget(&ptr,&Nod[i].InNum);//get inlet branches second
	Nod[i].BranIN=new int[Nod[i].InNum];
	Nod[i].BranOUT=new int[Nod[i].OutNum];
	for(j=0;j<Nod[i].OutNum;j++)//get outlet branches first
	{
	sget(&ptr,&Nod[i].BranOUT[j]);}
	for(j=0;j<Nod[i].InNum;j++)//get inlet branches second
	{
	sget(&ptr,&Nod[i].BranIN[j]);}
	}

	fgets(scanstr,2048,fp);
	for(i=0;i<BranNum;i++)
	{
	fgets(scanstr,2048,fp);
	ptr = scanstr;
	sget(&ptr,&Bra[i].BranNo);
	sget(&ptr,&Bra[i].EqulNo);
	sget(&ptr,&Bra[i].GrFac);
	sget(&ptr,&Bra[i].TubNum);
	Bra[i].TubNo=new int[Bra[i].TubNum];
	for(j=0;j<Bra[i].TubNum;j++)
	{
	sget(&ptr,&Bra[i].TubNo[j]);}//important, the tube number in branch, inputted first from the compressor suction

	}

	fgets(scanstr,2048,fp);
	fgets(scanstr,2048,fp);
	for(i=0;i<TubeNum;i++)
	{
	Tub[i].Seg = new TubEvpSeg[SegNum];
	fgets(scanstr,2048,fp);
	ptr = scanstr;
	sget(&ptr,&Tub[i].TubNo);
	sget(&ptr,&Tub[i].RowNo);
	sget(&ptr,&Tub[i].Refdownstream);
	sget(&ptr,&Tub[i].AirUpstreamUpper);
	sget(&ptr,&Tub[i].AirUpstreamLower);
	sget(&ptr,&Tub[i].GaFac);
	sget(&ptr,&Tub[i].even);
	}

	endSign=(int)GetDouble(fp);

	if(endSign!=-256)
	{
	fclose(fp);
	errorLog.Add("StructEvap::StructEvap","File Input Wrong");
	error=1;
	return;}

	ChargTunFix = GetDouble(fp);
	Vap_Ref = GetDouble(fp);
	K_ChargTun = GetDouble(fp);

	endSign=(int)GetDouble(fp);
	
	if(endSign!=-256)
	{
	fclose(fp);
	errorLog.Add("StructEvap::StructEvap","File Input Wrong");
	error=1;
	return;}

	fclose(fp);
	
}


StructEvap::~StructEvap()
{
	int i=0;

	for(i=0;i<TubeNum;i++)
	{
	delete Tub[i].Seg;
	}

	delete Tub;
	
	for(i=0;i<BranNum;i++)
	{
	delete Bra[i].TubNo;
	}

	delete Bra;

	for(i=0;i<NodNum;i++)
	{
	delete Nod[i].BranIN;
	delete Nod[i].BranOUT;
	}
	delete Nod;
};
//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//-------------------------------------------------------------------------------
//the follows are old models

/********************************************************************
Circuit pattern used in the Carrier rooftop unit at Purdue.
********************************************************************/
int EvapCircuit301(double mr,HP *HPo,double Ga,TP TPi,HP *HPi,TP *TPo,MASS *Sm,double *Aflow,ETdim *D)
{
	int i,j;
	WH WH,SWH;
	MASS m;
	const int n1 = D->Nrows;
	const int n2 = D->Ndeep;

	*Aflow=n1*D->NSeg*D->Aflow;
	const double Gr = mr/D->Ax ;

	*HPi = *HPo;
	*TPo = TPi;

	SWH.W=0;
	SWH.H=0;
	Sm->m=0;
	Sm->V=0;


	for(j=0;j<n2;j++) {
		for(i=0;i<n1;i++) {

			EvapTube(Gr,HPi,Ga,*TPo,&WH,&m,D);
			if(errorLog.IsError()) {
				errorLog.Add("EvapCircuit301","EvapTube");
				return 0;
			}
			Sm->m+=m.m; Sm->V+=m.V; SWH.W+=WH.W; SWH.H+=WH.H;

			EvapTubeBend(Gr,HPi,&m,D);
			if(errorLog.IsError()) {

				TXP TXPi = HPtoTXP(*HPi);
				if(errorLog.IsError()) {
					errorLog.Add("EvapCircuit301","TXPi");
					return 0;
				}

				char str[128];
				sprintf(str,"EvapTubeBend TXPi:(%lf,%lf,%lf)",TXPi.T,TXPi.X,TXPi.P);
				errorLog.Add("EvapCircuit301",str);
				return 0;
			}
			Sm->m+=m.m; Sm->V+=m.V;
		}

		// Calculate the mix mean state of the air between the two rows
		WH.W = SWH.W/n1;
		WH.H = SWH.H/n1;
		SWH.W = 0;
		SWH.H = 0;
		*TPo = WHtoTP(WH,TPi);
		if(errorLog.IsError()) {
			char str[64];
			sprintf(str,"TPo: WH = (%lf,%lf)",WH.W,WH.H);
			errorLog.Add("EvapCircuit301",str);
			return 0;
		}
	}



	return 0;
}

/********************************************************************
Single tube circuit.  This is the simplest pattern used for
generic modeling.
********************************************************************/
int EvapCircuit302(double mr,HP *HPo,double Ga,TP TPi,HP *HPi,TP *TPo,MASS *Sm,double *Aflow,ETdim *D)
{
	WH WH;

	*Aflow = D->NSeg*D->Aflow;
	const double Gr = mr/D->Ax;

	*HPi = *HPo;
	*TPo = TPi;

	EvapTube(Gr,HPi,Ga,*TPo,&WH,Sm,D);
	if(errorLog.IsError()) {
		errorLog.Add("EvapCircuit302");
		return 0;
	}

	// Calculate outlet state
	*TPo = WHtoTP(WH,TPi);
	if(errorLog.IsError()) {
		errorLog.Add("EvapCircuit302");
		return 0;
	}

	return 0;
}


/********************************************************************
evaporator layout for Harms 5 ton, this is to assume air flow from tube to tube
without assuming mixed after each row, the definite layout is from Todd Harms
 (2002) phD thesis
********************************************************************/

int EvapCircuit303(double mr,HP *HPo,double Ga,TP TPi,HP *HPi,TP *TPo,MASS *Sm,double *Aflow,ETdim *D)
{
	int I;
	double P[65+1],H[65+1],T_I[64+1],T_O[64+1];
    double AirP_I[64+1],AirP_O[64+1];
	WH WH,SWH;
	MASS m;
	mr=mr/8;

	int n1=16+16;// row number at the front face including upper side and lower side
	
	*Aflow=n1*D->NSeg*D->Aflow;
	const double Gr = mr/D->Ax;

	SWH.W=0;
	SWH.H=0;
	Sm->m=0;
	Sm->V=0;
	*TPo = TPi;

	const double ROWS = 32; //row number
	const double AA=0.0254e0;//tube space
	const double D_O = 9.91e-3;// tube outside diameter
	const double FIN_S = 0.001579e0;//fin spacing
	const double FIN_T = 0.114e-3;//fin thickness
	const double L_COND = 0.6172e0;   //  !tube length
	const double M_DOT_AT=1.10;
	const double  FACE_A= (ROWS*AA)*L_COND;//facing area
//	Ga = M_DOT_AT/(FACE_A)*(AA*FIN_S/((FIN_S-FIN_T)*(AA-D_O)));//shenbo guess, related with input
	
	//begin to calculation the heat transfer in the evaporato, from the bottom to upper
	for(I=1;I<=64;I++){//air state
			T_I[I]=TPi.T;
			AirP_I[I]=TPi.P;
			AirP_O[I]=TPi.P;
	}

	
	for(I=1;I<=4;I++){//refrigerant state
		P[I]=HPo->P;
		H[I]=HPo->H;
	}


	for(I=1;I<=64;I++){

		if ((I-1)/4e0 != int((I-1)/4e0)){//air state
				T_I[I]=T_O[I-1];
				AirP_I[I]=AirP_O[I-1];
		}
			
		//heat transfer and pressure drop calculation
		{
			HPi->P = P[I]; HPi->H = H[I];
			TPo->P = AirP_I[I]; TPo->T = T_I[I];

			EvapTube(Gr,HPi,Ga,*TPo,&WH,&m,D);
			if(errorLog.IsError()) {
				errorLog.Add("EvapCircuit301","EvapTube");
				return 0;
			}
			Sm->m+=m.m; Sm->V+=m.V; 

			EvapTubeBend(Gr,HPi,&m,D);//
			if(errorLog.IsError()) {

				TXP TXPi = HPtoTXP(*HPi);
				if(errorLog.IsError()) {
					errorLog.Add("EvapCircuit301","TXPi");
					return 0;
				}

				char str[128];
				sprintf(str,"EvapTubeBend TXPi:(%lf,%lf,%lf)",TXPi.T,TXPi.X,TXPi.P);
				errorLog.Add("EvapCircuit301",str);
				return 0;
			}
			Sm->m+=m.m; Sm->V+=m.V;

			*TPo = WHtoTP(WH,TPi);
			if(errorLog.IsError()) {
			char str[64];
			sprintf(str,"TPo: WH = (%lf,%lf)",WH.W,WH.H);
			errorLog.Add("EvapCircuit303",str);
			return 0;
			}
			T_O[I] = TPo->T;//air output
			AirP_O[I] = TPo->P;
			P[I] = HPi->P;//refrigerant output
			H[I] = HPi->H;
		}
		

		if (I >= 4){
			P[I+1]=P[I-3];
			H[I+1]=H[I-3];
			
			if (I/8e0==int(I/8e0)){
				P[I+1]=P[I-1];
				H[I+1]=H[I-1];
			}
			if ((I-1)/8e0==int((I-1)/8e0)){
				P[I+1]=P[I-4];
				H[I+1]=H[I-4];
			}
			if ((I-2)/8e0==int((I-2)/8e0)){
				P[I+1]=P[I-2];
				H[I+1]=H[I-2];
			}
			if ((I-3)/8e0==int((I-3)/8e0)){
				P[I+1]=P[I-5];
				H[I+1]=H[I-5];
			}
		}
	}



	const double P_OUT=(P[61]+P[62]+P[63]+P[63])/4e0;     //!Not Physical
	const double H_OUT=(H[61]+H[62]+H[63]+H[63])/4e0;

	TPo->T =T_O[64];
	TPo->P = AirP_O[64];	
	
	Sm->m = Sm->m*2.0;
	Sm->V = Sm->V*2.0;

	HPi->H=H_OUT;
	HPi->P=P_OUT;
	return 0;
}

