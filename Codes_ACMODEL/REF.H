#ifndef _Ref
#define _Ref

#include "defs.h"

#ifdef _R22
#define REFRIGERANT "r22"
#define TMAX 150
#define TMIN -100
#define PMAXth 4950
#define PMINth 10
#define PMAXtr 4500
#define PMINtr 100
#define PMAX 4500
#define PMIN 100
#endif

#ifdef _R410A
#define REFRIGERANT "r410a"
#define TMAX 160
#define TMIN -50
#define PMAXth 4900
#define PMINth 55
#define PMAXtr 4600
#define PMINtr 100
#define PMAX 4600
#define PMIN 100
#endif

#ifdef _R407C
#define REFRIGERANT "r407c"
#define TMAX 160
#define TMIN -50
#define PMAXth 4200
#define PMINth 20
#define PMAXtr 3800
#define PMINtr 100
#define PMAX 3800
#define PMIN 100
#endif

#endif 

