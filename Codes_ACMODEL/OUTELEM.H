#ifndef _OUTELEM
#define _OUTELEM

#include <string.h>
#include "conversion.h"

//B.S.----------------------------
//#include <fstream.h>
//#include <iostream.h>
//-------------------------------B.S.

class OutputDataElement {
public:
	// data
	double time;
	double T[15];//temperatures of different state points
	double X[15];//qualities of different state points
	double P[15];//pressures of different state points
	double H[15];//enthalpies of different state points
	double Tsat[15];//saturation temperatures of different state points
	double mr[5];//mass flow rates of different components
	double mc;//air mass flow rate at the condenser side
	double Tci;//air inlet temperature of the condenser
	double Tco;//air outlet temperture of the condenser
	double me;//air mass flow rate at the evaporator side
	double Tei;//air inlet temperature of the evaporator
	double Twb_ei; // air inlet wetbulb temperature - evaporator
	double Teo;//air outlet temperture of the evaporator
	double HumPei;//air inlet relative humidity at the evaporator side
	double HumPeo;//air outlet relative humidity at the evaporator side
	double Ecmp;//compressor power consumption
	double Qc;//total heat transfer in the condenser
	double Qe;//total heat transfer of the evaporator
	double Charge, Charge_target;//system charge mass
	double P1;			// pressure before restriction tube
	double P2;			// pressure before expansion device
	double devMr,devH,devCharge;// residual of the iteration
	double inputP1,inputP2,inputH1;//initial guesses
	double mr_TP;//compressor mass flow rate
	double CE;//overall residual of the iteration
	double Cond_AirDP, Eva_AirDP;
	double Time_dev;
	double cfmca, cfmea;// cfm air cond and evap
	int isError, isConv; // is there an error on this run and did thir run converge to a solution within the required tolerance
	double mr_exp; // expansion device mass flow rate for short tube model
	double Gea,Gca; // air mass fluxes
	double nvAdj; // input - compressor volumeteric efficiency scaler
	double Dadj; // input - expansion device diameter scaler
	double Array[128];//backup all the output data
	double Tsub4,Tsub5;
	double Tsup1;
	double LiqLCond,TPLCond,VapLCond;
	double LiqLEvap,TPLEvap,VapLEvap;
	double m[8];

	// auxillary data, experimental data passed through for comparison
	// to model predictions
	int auxDataCnt;
	double auxData[128];

	char outputFilename[128];

	// functions
	OutputDataElement();
	~OutputDataElement(); 
	void StartOutputFile(const char* fn);
	void AppendOutputToFile();
	void CloseOutputFile();
	void ClearAll(); //set zero to all the output variables
	void output(char str_res);//shenboadd
	void outputTOarray();//shenboadd output all the data to a data array
	void arrayback();//Boshen return the data from a data array

	// creates output file to develop multidimensional runs
	void StartMasterInputFile();
	void AppendMasterInputFile();

};

#endif