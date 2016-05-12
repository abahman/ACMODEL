#ifndef _ERRORLOG
#define _ERRORLOG

#include <stdio.h>

class ErrorLog {
public:
	ErrorLog();
	~ErrorLog();
	void Add(const char* func,const char* msg = "");
	void AddInitialError(const char* func,const char* msg = "");
	void AddWarning(const char* func,const char* msg = "");
	int IsError() {return isError;}
	int IsInitialError() {return isInitialError;}
	void ClearError(const char* func,const char* msg = "",int publish = 1);
	void Open();
	void Close();
private:
	int isError;
	int isInitialError;
	FILE* fp;
};

#endif