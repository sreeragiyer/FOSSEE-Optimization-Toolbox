


#include "CbcConfig.h"
#include "ClpConfig.h"
#include "OsiConfig.h"
#include "SymConfig.h"
#include "IpoptConfig.h"
#include "BonminConfig.h"
#include "fotConfig.h"

#include "api_scilab.h"

int sci_fotversion(char* fname,unsigned long int fname_len)
{
	//checking number of arguments
	CheckInputArgument(pvApiCtx,0,0);
	CheckOutputArgument(pvApiCtx,1,1);

	 
	//FOT Version
	char fotver[]=FOT_VERSION;
	sciprint("FOSSEE Optimization Toolbox: Version %s\n",fotver);

	//Latest Git id commit		
	char gitid[]=GIT_ID;
    	sciprint(" Latest Git Commit ID: ");
    	for(int i=0;i<7;i++)
        	sciprint("%c",gitid[i]);
	
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
	sciprint(" IPOPT: %s\n",ipover);
	sciprint(" OSI: %s\n",osiver);
	sciprint(" CBC: %s\n",cbcver); 
	sciprint(" Bonmin: %s\n",bonver);
	return 0;
}

