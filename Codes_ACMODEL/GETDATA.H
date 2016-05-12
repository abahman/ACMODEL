#ifndef _GETDATA
#define _GETDATA

#define SYSTEM 0
#define COMPRESSOR 100
#define CONDENSER 200
#define EVAPORATOR 300
#define EXPDEV 400
#define SUBCOOLER 500
#define LIQUIDLINE 600
#define FILTER 700
#define SUCTIONLINE 800
#define DISCHARGELINE 900
#define DISTRILINE 1000

#include <stdio.h>

int GetComponentType(char*,int);
int PosComponentData(FILE*,int,int);

char *GetTime(void);
double GetDouble(FILE*);
void GetString(FILE*,char*);
int GetInt(FILE*);

#endif
