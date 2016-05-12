#include "stdafx.h"

#include <time.h>
#include <string.h>
#include "errorlog.h"

ErrorLog errorLog;

ErrorLog::ErrorLog()
{
	isError = 0;
	isInitialError = 0;
	fp = NULL;

	Open();
}

ErrorLog::~ErrorLog()
{
	Close();
}

void ErrorLog::Open()
{
	if(fp) return;

	fp = fopen("ErrorLog.htm","w");
	if(fp) {
		char da[32],ti[32];
		_strdate(da);
		_strtime(ti);
		fprintf(fp,"Starting error log %s %s.\n\n",da,ti);
		fprintf(fp,"<P>");
		fprintf(fp,"<TABLE BORDER=\"1\">\n");
		fprintf(fp,"\t<TR>");
		fprintf(fp,"\t\t<TH>Event</TH>\n");
		fprintf(fp,"\t\t<TH>Function</TH>\n");
		fprintf(fp,"\t\t<TH>Description</TH>\n");
		fprintf(fp,"\t</TR>");
	}
}

void ErrorLog::Close()
{
	if(!fp) return;

	fprintf(fp,"</TABLE>\n");
	fclose(fp);
	fp = NULL;
}

void ErrorLog::ClearError(const char* func, const char* msg, int /*publish*/)
{
	isError = 0;
	isInitialError = 0;

	//if(!publish) return;

	Open();
	if(!fp) return;

	fprintf(fp,"\t<TR>");
	fprintf(fp,"\t\t<TD><B>Clearing!</B></TD>\n");
	fprintf(fp,"\t\t<TD>%s</TD>\n",strlen(func)>0?func:"&nbsp;");
	fprintf(fp,"\t\t<TD>%s</TD>\n",strlen(msg)>0?msg:"&nbsp;");
	fprintf(fp,"\t</TR>");
}

void ErrorLog::AddWarning(const char* /*func*/, const char* /*msg*/)
{
	//if(fp) fprintf(fp,"Warning! %s >>> %s\n",func,msg);
}

void ErrorLog::AddInitialError(const char* func, const char* msg)
{
	isInitialError = 1;

	Add(func,msg);
}

void ErrorLog::Add(const char* func, const char* msg)
{
	isError = 1;

	Open();

	Open();
	if(!fp) return;

	fprintf(fp,"\t<TR>");
	fprintf(fp,"\t\t<TD>Error!</TD>\n");
	fprintf(fp,"\t\t<TD>%s</TD>\n",strlen(func)>0?func:"&nbsp;");
	fprintf(fp,"\t\t<TD>%s</TD>\n",strlen(msg)>0?msg:"&nbsp;");
	fprintf(fp,"\t</TR>");
}

