#include "stdafx.h"

#include "Parse.h"
#include <stdio.h>
#include <string.h>

int
IsNum(char c)
{
	if(c>='0' && c<='9')
		return 1;
	else
		return 0;
}

int
IsWhiteSpace(char c)
{
	if(c==' ' || c=='\t' || c=='\n' || c=='\r')
		return 1;
	else
		return 0;
}

int
IsFloatChar(char c)
{
	if(IsNum(c) || c=='+' || c=='-' || c=='e' || c=='E' || c=='.')
		return 1;
	else
		return 0;
}

/* Scans a string for a floating point number.  Returns 1 if
	successful and 0 otherwise.  "data" will contain the number
	if successful and the pointer "str" will advance to just after
	the parsed data.  0 will be returned when the last data in
	the file is found. */

int
sget(char **str, char *fnd)
{
	int j(0);
	char buf[256];

	while(**str!=0) {
		if(!IsWhiteSpace(**str)) {
			if(j>=256) return 0;
			buf[j++]=**str;
		} else {
			if(j) {
				buf[j]=0;
				strcpy(fnd,buf);
				return 1;
			}
		}
		(*str)++;
	}

	return 0;
}

int
sget(char **str, double *data)
{
	int j(0);
	char buf[256];

	while(**str!=0) {
		if(IsFloatChar(**str)) {
			if(j>=200) return 0;
			buf[j++]=**str;
		} else {
			if(j) {
				buf[j]=0;
				return sscanf(buf,"%lf",data);
			}
		}
		(*str)++;
	}

	return 0;
}

int
sget(char **str, int *data)
{
	int j(0);
	char buf[256];

	while(**str!=0) {
		//if(IsNum(**str)) {//B.S. change
		if(IsFloatChar(**str)) {
			if(j>=200) return 0;
			buf[j++]=**str;
		} else {
			if(j) {
				buf[j]=0;
				return sscanf(buf,"%d",data);
			}
		}
		(*str)++;
	}

	return 0;
}
