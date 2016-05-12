#include "stdafx.h"

#include <sys\timeb.h>
#include <time.h>
#include "mtime.h"

unsigned long int
MilliTime(void)
{
	timeb t;
	ftime(&t);
	return 1000*(unsigned long int)t.time + (unsigned long int)t.millitm;
}

