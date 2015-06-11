/***********************************************************
 Typical file
 verbose error debug options
***********************************************************/
#include <stdio.h>

extern	FILE	*FDOUT;
long		VERBOSE;
long		DEBUG;

xdebugverbose_c(s)
/**********************************************************
 VERBOSE function
***********************************************************/
char	*s;
{
	if (s == 0L) return(0);
	if ((*s == 'y') || (*s == 'Y')) DEBUG = 1; else DEBUG = 0;
	if (DEBUG) VERBOSE = 1;
	return(0);
	}

xverbose_c(s)
/**********************************************************
 VERBOSE function
***********************************************************/
char	*s;
{
	if (s == 0L) return(0);
	if ((*s == 'y') || (*s == 'Y')) VERBOSE = 1; else VERBOSE = 0;
	return(0);
	}

xerr_c(s)
/**********************************************************
 ERROR function
***********************************************************/
char	*s;
{
	printf("ERROR : %s\n",s);
	if (FDOUT) fprintf(FDOUT,"ERROR : %s\n",s);
	exit(0);
	}
