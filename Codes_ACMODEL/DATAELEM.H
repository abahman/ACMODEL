#ifndef _DATAELEM
#define _DATAELEM

#include <stdio.h>
#include <string.h>
//#include "operwnd.h"

struct OperationData {
	char batchFilename[128];
	char outputFilename[128];
	OperationData() {
		strcpy(batchFilename,"batch.in");
		strcpy(outputFilename,"batch.out");
	}
};

class DataElement {
private:
	void strip_path(char* filename);
	int FileBasic(char* str, FILE* fp);
	int FileNIST(char* str, FILE* fp);

protected:
	char filename[128];
public:
	OperationData operationData;

	int system_type;

	int input_type;
	double input[20];

	int comp_type;
	int comp_nprms;
	double comp[40];
	int cond_type;
	int cond_nprms;
	double cond[20];
	int evap_type;
	int evap_nprms;
	double evap[20];
	int subcool_type;
	int subcool_nprms;
	double subcool[10];
	int expdev_type;
	int expdev_nprms;
	double expdev[10];
	int restrict_type;
	int restrict_nprms;
	double restrict[10];

	DataElement();
	~DataElement();

	void Default();
	void DefaultInput(int system_type);
	void DefaultSystem(int type);
	void DefaultCompressor(int type);
	void DefaultCondenser(int type);
	void DefaultEvaporator(int type);
	void DefaultExpDev(int type);
	void DefaultRestriction(int type);
	void DefaultSubcooler(int type);

	int Import(void);
	int Export(void);
	int WriteParams(int restart=0);
	int WriteInput(int restart=1);
	void Filename(const char *fn) {
		strncpy(filename,fn,99);
		strip_path(filename);
	}
	const char* Filename(void) {return filename;}

	void WriteINIFile();
};

#endif