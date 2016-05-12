#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include "r22.h"
#include "wair.h"
#include "air.h"
#include "mass.h"
#include "comp.h"
#include "lines.h"
#include "evap.h"
#include "cond.h"
#include "expdev.h"
#include "num2d.h"
#include "num3d.h"
#include "numeric.h"
//#include "print.h"
#include "ref.h"
//#include "range.h"
#include "mainfunc.h"
//#include "getdata.h"

#include "ref.h"
#include "Spline.h"

void InitializeComponents(void);

class ENGINE {
public:
	// data
	FILE *in;
	FILE *states;
	InVars I;
	int systemType;
	int PRT;
	int fix;
	int i;
	int Cerror,Cnoconv,Nrun0,Nrun1,Crun0,Crun1;
	double Trun0,Trun1;
	double X[3];
	double Xg[3];
	double ConvError;
	double ConvTol;
	time_t tlast;
	time_t tcur;
	time_t tbegin;

	//new variables for definition
	char filename[128];//data file name
	int error;
	int ACTUNE;
	int REF;
	// functions
	ENGINE();
	~ENGINE();
	int compu(void);
	int BuildSystem_Charg(double*,void*);//B.S.
	int BuildSystem_SubC_2D(double*,void*);//B.S.
	int BuildSystem_SubC_3D(double*,void*);//B.S.
};

