#ifndef _COMPSTAT
#define _COMPSTAT

#include "r22.h"
#include "mass.h"

// Compound cooling compressor
// State 1 = compressor inlet
// State 2 = intermediate state
// State 3 = compressor outlet
//
void CompressorCCLowStage(const char* filename,double H1,double P1,double P2,double P3,double& H2,double& mr,double& Ei,MASS& m,double* Prms);
//
void CompressorCCHighStage(const char* filename,double H2,double P2,double P3,double& H3);

#endif