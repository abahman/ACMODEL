#ifndef _ACTUNE
#define _ACTUNE

#include "datast.h"

int ActuneNISTChiller();
int ActuneNISTChillerComponents(SteadyData* steady);
int ActuneRTU(const char* modelParamsFilename, const char* measurementsFilename, const char* outputFilename, int);

#endif