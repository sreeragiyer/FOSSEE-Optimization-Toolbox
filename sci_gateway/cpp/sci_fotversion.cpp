// Author: Sreerag Iyer
// Organization: FOSSEE, IIT Bombay
// Email: toolbox@scilab.in

#include "CbcConfig.h"
#include "ClpConfig.h"
#include "OsiConfig.h"
#include "SymConfig.h"
#include "IpoptConfig.h"
#include "BonminConfig.h"
#include "stdio.h"

#ifdef _MSC_VER
#define popen _popen
#endif

extern "C" {
#include "api_scilab.h"
#include "sciprint.h"


int sci_fotversion(char* fname,unsigned long int fname_len)
{
	//checking number of arguments
	CheckInputArgument(pvApiCtx,0,0);
	CheckOutputArgument(pvApiCtx,1,1);

	//FOT Version
	char fotver[]="0.2";
	sciprint("FOSSEE Optimization Toolbox: Version %s",fotver);

	//Latest Git id commit	
	char gitid[128] = "sdsd2452";		//default value of Git Id
	FILE *fp = popen("git rev-parse HEAD", "r");
	if (fp != NULL)
		while (fgets(gitid, 128, fp) != NULL);
	fclose(fp);
	sciprint("\n Latest Git Commit ID: ");
	for (int i = 0; i < 7; i++)
		sciprint("%c", gitid[i]);

	//Library versions
	char cbcver[]=CBC_VERSION;
	char clpver[]=CLP_VERSION;
	char osiver[]=OSI_VERSION;
	char symver[]=SYMPHONY_VERSION;
	char ipover[]=IPOPT_VERSION;
	char bonver[]=BONMIN_VERSION;
	
	sciprint("\n\nLibraries used in toolbox:\n");	
	sciprint(" CLP: %s\n",clpver);
	sciprint(" Symphony: %s\n",symver);
	sciprint(" IPOPT (with Mumps): %s\n",ipover);
	sciprint(" OSI: %s\n",osiver);
	sciprint(" CBC: %s\n",cbcver); 
	sciprint(" Bonmin: %s\n",bonver);
	return 0;
}
}
