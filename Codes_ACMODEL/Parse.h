#ifndef _PARSE
#define _PARSE

int IsNum(char c);
int IsFloatChar(char c);
int IsWhiteSpace(char c);
int sget(char **str, char *fnd);
int sget(char **str, double *data);
int sget(char **str, int *data);

#endif