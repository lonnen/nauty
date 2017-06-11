/******************************************************************************
 *                                                                            *
 * This is the header file for traces() version 2.0, which is included into   *
 *   nauty() version 2.5.                                                     *
 *                                                                            *
 *   nauty is Copyright (1984-2014) Brendan McKay. All rights reserved.       *
 *   Subject to the waivers and disclaimers in nauty.h.                       *
 *   Traces is Copyright (2008-2014) Adolfo Piperno. All rights reserved.     *
 *                                                                            *
 *   CHANGE HISTORY                                                           *
 *       28-Dec-12 : final changes for version 2.0                            *
 *****************************************************************************/

#include "gtools.h"
#include "schreier.h" 

typedef struct TracesOptions {
	boolean getcanon;
	boolean writeautoms;
	boolean cartesian;
	boolean digraph;
	boolean defaultptn;
	int linelength;
	FILE* outfile;
	int strategy;
	int verbosity;
	permnode **generators;
	void (*userautomproc)(int,int*,int);
} TracesOptions;

#define DEFAULTOPTIONS_TRACES(opts) TracesOptions opts \
= { FALSE, FALSE, FALSE, FALSE, TRUE, 0, NULL, 0, 0, NULL, NULL }

typedef struct TracesStats {
	double grpsize1;
	int grpsize2;
	int numgenerators;
	int numorbits;
	int treedepth;
	int canupdates;
	unsigned long numnodes;
	unsigned long interrupted;
	unsigned long peaknodes;
} TracesStats;

#ifdef __cplusplus
extern "C" {
#endif

extern void Traces(sparsegraph*,int*,int*,int*,TracesOptions*,
				   TracesStats*,sparsegraph*);									
extern void refine_tr(sparsegraph*,int*,int*,int*,int*,TracesOptions*);		
extern void traces_freedyn(void);

#ifdef __cplusplus
}
#endif
