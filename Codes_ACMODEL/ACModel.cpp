// ACModel.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

// internal ACModel include files
#include "outelem.h"
#include "dataelem.h"
#include "compu.h"
#include "actune.h"
#include "errorlog.h"

// these guys enable the API call that opens browers and Excel to show results
// after the calculations are complete.
#include <windows.h>
#include <shellapi.h>

DataElement gData;                // model parameters
OutputDataElement gOutData;

extern int cyclePosition;
extern ErrorLog errorLog;

//# define RUN_EXERCISE
//#ifdef  RUN_EXERCISE
//#include "ExerciseComponents.h"
//#include "Corr.h"
//#include "R22.h"
//#define RUN_EXERCISE_COMPONENTS
//#define RUN_EXERCISE_Corr
//#define RUN_EXERCISE_PROP
//#endif

#define RUN_ACMODEL

int main(int argc, char* argv[])
{
	int errorCode;

/*#ifdef RUN_EXERCISE_COMPONENTS
	ExerciseComponents();
#endif

#ifdef RUN_EXERCISE_Corr
	Exercise_Corr();
#endif

#ifdef RUN_EXERCISE_PROP
	Exercise_Prop();
#endif*/

	// Initializes program state variables.  Look in [compu.cpp] for more details.
	InitializeComponents();

	ENGINE engine;
	if(engine.error) 
	{
	errorLog.Add("ENGINE::ENGINE","Class fails to be constructed");
	return -1;
	}

	if(engine.ACTUNE)//actune starts
	{
	// Solve for a starting point at nominal conditions using ACTune first.  Look in
	// [actune.cpp] for details.
	if((errorCode=ActuneRTU("./InputDoc/acmodel.dat","./InputDoc/actune1.in","./OutputDoc/ACTuneOutput.htm", engine.I.REV_CAL))==0) {
		printf("Actune successful!\n");
	} else {
		printf("Actune error! errorCode=%d\n",errorCode);
		//errorLog.Close();
		//ShellExecute(NULL, "open", "ErrorLog.htm", NULL, NULL, SW_SHOWNORMAL);
	}

	//ShellExecute(NULL, "open", ".//OutputDoc//ACTuneOutput.htm", NULL, NULL, SW_SHOWNORMAL);
	}//actune end
	

#ifdef RUN_ACMODEL
	// ACModel
	// ACModel
	// ACModel

	// Initializes program state variables.  Look in [compu.cpp] for more details.
	InitializeComponents();

	//ENGINE engine;//original place
	int cnt=0;
	if((errorCode=engine.compu())==0) {
		printf("Acmodel successful!\n");
	} else {
		//printf("Acmodel error! errorCode=%d\n",errorCode);
	}

	//ShellExecute(NULL, "open", ".//OutputDoc//ACModelOutput.htm", NULL, NULL, SW_SHOWNORMAL);
#endif

	printf("Done.\n");
	return 0;
}
