#include "stdafx.h"

#include <time.h>
#include <stdio.h>
#include "getdata.h"

/* Local functions */
int GetType(FILE*,int*,int*);

/********************************************************************
Gets the first occurance of a component type of the specified
base class.

base 0 - system types
base 100 - compressor models
base 200 - condenser models
base 300 - evaporator models
base 400 - expansion device models
********************************************************************/
int GetComponentType(char *fn, int base) {
	FILE *fp;
	int type,nOccur;

	fp = fopen(fn,"r");
	if(fp==NULL) return 0;
	do {
		GetType(fp,&type,&nOccur);
	} while(type<=base || type>base+99);
	fclose(fp);

	return type;
}

/********************************************************************
********************************************************************/
int PosComponentData(FILE *fp, int base, int nOccur) {
	int t,n;

	do {
		GetType(fp,&t,&n);
	} while(!(t>base && t<base+100 && n==nOccur) && t!=0);
	return t;
}

/********************************************************************
Gets the next componet type number in the specified file.
Returns 0 if an error occurs.
********************************************************************/
int GetType(FILE *fp, int *type, int *nOccur) {
	char s[200];

	do {
		fgets(s,199,fp);
	} while(s[0]!='#' && !feof(fp));
	if(feof(fp)) {
		return 0;
	} else {
		if(sscanf(s+1,"%d,%d",type,nOccur)==2) 
			return *type;
		else
			return 0;
	}
}

/********************************************************************
Returns a pointer to a string containing the current time.
********************************************************************/
char *GetTime(void) {
	time_t t;
	time(&t);
	return ctime(&t);
}

/********************************************************************
Gets a int from the beginning of a line.
********************************************************************/
int GetInt(FILE *fp) {
	char s[200];
	int d;
	int n,i=1000;

	do {
		fgets(s,199,fp);
		n=sscanf(s,"%d",&d);
		i--;
	} while(i>0 && n!=1);

	return d;
}

/********************************************************************
Gets a double from the beginning of a line.
********************************************************************/
double GetDouble(FILE *fp) {
	char s[200];
	double d;
	int n,i=1000;

	do {
		fgets(s,199,fp);
		n=sscanf(s,"%lf",&d);
		i--;
	} while(i>0 && n!=1);

	return d;
}

/********************************************************************
Gets a string from the beginning of a line.
********************************************************************/
void GetString(FILE *fp, char *str) {
	char s[200];
	int n,i=1000;

	do {
		fgets(s,199,fp);
		n=sscanf(s,"%s",str);
		i--;
	} while(i>0 && n!=1);
}



