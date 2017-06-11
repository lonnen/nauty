/* labelg.c version 1.6; B D McKay, Jan 2011 */

#define USAGE "labelg [-qsg] [-fxxx] [-S|-t] [-i# -I#:# -K#] [infile [outfile]]"

#define HELPTEXT \
" Canonically label a file of graphs.\n\
\n\
    -s  force output to sparse6 format\n\
    -g  force output to graph6 format\n\
        If neither -s or -g are given, the output format is\n\
        determined by the header or, if there is none, by the\n\
        format of the first input graph. Also see -S. \n\
    -S  Use sparse representation internally.\n\
         Note that this changes the canonical labelling.\n\
         Multiple edges are not supported.  One loop per vertex is ok.\n\
    -t  Use Traces.\n\
         Note that this changes the canonical labelling.\n\
         Multiple edges and loops are not supported, nor invariants.\n\
\n\
    The output file will have a header if and only if the input file does.\n\
\n\
    -fxxx  Specify a partition of the point set.  xxx is any\n\
        string of ASCII characters except nul.  This string is\n\
        considered extended to infinity on the right with the\n\
        character 'z'.  One character is associated with each point,\n\
        in the order given.  The labelling used obeys these rules:\n\
         (1) the new order of the points is such that the associated\n\
        characters are in ASCII ascending order\n\
         (2) if two graphs are labelled using the same string xxx,\n\
        the output graphs are identical iff there is an\n\
        associated-character-preserving isomorphism between them.\n\
        No option can be concatenated to the right of -f.\n\
\n\
    -i#  select an invariant (1 = twopaths, 2 = adjtriang(K), 3 = triples,\n\
        4 = quadruples, 5 = celltrips, 6 = cellquads, 7 = cellquins,\n\
        8 = distances(K), 9 = indsets(K), 10 = cliques(K), 11 = cellcliq(K),\n\
       12 = cellind(K), 13 = adjacencies, 14 = cellfano, 15 = cellfano2,\n\
       16 = refinvar(K))\n\
    -I#:#  select mininvarlevel and maxinvarlevel (default 1:1)\n\
    -K#   select invararg (default 3)\n\
\n\
    -q  suppress auxiliary information\n"


/*************************************************************************/

#include "gtools.h"
#include "nautinv.h"
#include "gutils.h"
#include "traces.h"

static struct invarrec
{
    void (*entrypoint)(graph*,int*,int*,int,int,int,int*,
                      int,boolean,int,int);
    char *name;
} invarproc[]
    = {{NULL, "none"},
       {twopaths,    "twopaths"},
       {adjtriang,   "adjtriang"},
       {triples,     "triples"},
       {quadruples,  "quadruples"},
       {celltrips,   "celltrips"},
       {cellquads,   "cellquads"},
       {cellquins,   "cellquins"},
       {distances,   "distances"},
       {indsets,     "indsets"},
       {cliques,     "cliques"},
       {cellcliq,    "cellcliq"},
       {cellind,     "cellind"},
       {adjacencies, "adjacencies"},
       {cellfano,    "cellfano"},
       {cellfano2,   "cellfano2"},
       {refinvar,    "refinvar"}
      };

#define NUMINVARS ((int)(sizeof(invarproc)/sizeof(struct invarrec)))

static nauty_counter orbtotal;
static double unorbtotal;
extern int gt_numorbits;

/**************************************************************************/

int
main(int argc, char *argv[])
{
	graph *g;
	sparsegraph sg,sh;
	int m,n,codetype;
	int argnum,j,outcode;
	char *arg,sw,*fmt;
	boolean badargs;
	boolean sswitch,gswitch,qswitch,fswitch,Oswitch;
	boolean iswitch,Iswitch,Kswitch,Mswitch,Sswitch;
	boolean uswitch,tswitch;
	int inv,mininvarlevel,maxinvarlevel,invararg;
	long minil,maxil;
	double t;
	char *infilename,*outfilename;
	FILE *infile,*outfile;
	nauty_counter nin;
	int ii,secret,loops;
	DEFAULTOPTIONS_TRACES(traces_opts);
	TracesStats traces_stats;
#if MAXN
	graph h[MAXN*MAXM];
	int lab[MAXN],ptn[MAXN],orbits[MAXN];
#else
	DYNALLSTAT(graph,h,h_sz);
	DYNALLSTAT(int,lab,lab_sz);
	DYNALLSTAT(int,ptn,ptn_sz);
	DYNALLSTAT(int,orbits,orbits_sz);
#endif

	HELP;

	nauty_check(WORDSIZE,1,1,NAUTYVERSIONID);

	sswitch = gswitch = qswitch = FALSE;
	fswitch = Oswitch = Mswitch = FALSE;
	iswitch = Iswitch = Kswitch = FALSE;
	uswitch = Sswitch = tswitch = FALSE;
	infilename = outfilename = NULL;
	inv = 0;

	argnum = 0;
	badargs = FALSE;
	for (j = 1; !badargs && j < argc; ++j)
	{
	    arg = argv[j];
	    if (arg[0] == '-' && arg[1] != '\0')
	    {
		++arg;
		while (*arg != '\0')
		{
		    sw = *arg++;
		         SWBOOLEAN('s',sswitch)
		    else SWBOOLEAN('g',gswitch)
		    else SWBOOLEAN('u',uswitch)
		    else SWBOOLEAN('q',qswitch)
		    else SWBOOLEAN('O',Oswitch)
		    else SWBOOLEAN('S',Sswitch)
		    else SWBOOLEAN('t',tswitch)
		    else SWINT('i',iswitch,inv,"labelg -i")
		    else SWINT('K',Kswitch,invararg,"labelg -K")
		    else SWRANGE('k',":-",Iswitch,minil,maxil,"labelg -k")
		    else SWRANGE('I',":-",Iswitch,minil,maxil,"labelg -I")
		    else if (sw == 'f')
		    {
			fswitch = TRUE;
			fmt = arg;
			break;
		    }
		    else SWINT('M',Mswitch,secret,"labelg -M")
		    else badargs = TRUE;
		}
	    }
	    else
	    {
		++argnum;
		if      (argnum == 1) infilename = arg;
	        else if (argnum == 2) outfilename = arg;
		else                  badargs = TRUE;
	    }
	}

	if (sswitch && gswitch) 
            gt_abort(">E labelg: -s and -g are incompatible\n");

	if (tswitch && (fswitch || Sswitch))
            gt_abort(">E labelg: -t is incompatible with -S and -f \n");

	if (!Sswitch && iswitch && (inv > NUMINVARS))
	    gt_abort(">E labelg: -i value must be 0..16\n");
	if (Sswitch && iswitch && (inv != 8))
	    gt_abort(">E labelg: only distances (-i8) is available with -S");
	if (tswitch && iswitch)
	    gt_abort(">E labelg: invariants are not available with -t");

	if (iswitch && inv == 0) iswitch = FALSE;

	if (iswitch)
	{
	    if (Iswitch)
	    {
		mininvarlevel = minil;
		maxinvarlevel = maxil;
	    }
	    else
		mininvarlevel = maxinvarlevel = 1;
	    if (!Kswitch) invararg = 3;
	}

	if (!Mswitch) secret = 1;

	if (badargs || argnum > 2)
	{
	    fprintf(stderr,">E Usage: %s\n",USAGE);
	    GETHELP;
	    exit(1);
	}

	if (!qswitch)
	{
	    fprintf(stderr,">A labelg");
	    if (sswitch || gswitch || fswitch || iswitch 
			|| tswitch || Sswitch)
		fprintf(stderr," -");
	    if (sswitch) fprintf(stderr,"s");
	    if (gswitch) fprintf(stderr,"g");
	    if (Sswitch) fprintf(stderr,"S");
	    if (tswitch) fprintf(stderr,"t");
	    if (iswitch)
		fprintf(stderr,"i=%s[%d:%d,%d]",invarproc[inv].name,
		        mininvarlevel,maxinvarlevel,invararg);
	    if (fswitch) fprintf(stderr,"f%s",fmt);
	    if (argnum > 0) fprintf(stderr," %s",infilename);
	    if (argnum > 1) fprintf(stderr," %s",outfilename);
	    fprintf(stderr,"\n");
	    fflush(stderr);
	}

	if (infilename && infilename[0] == '-') infilename = NULL;
	infile = opengraphfile(infilename,&codetype,FALSE,1);
	if (!infile) exit(1);
	if (!infilename) infilename = "stdin";

	if (!outfilename || outfilename[0] == '-')
	{
	    outfilename = "stdout";
	    outfile = stdout;
	}
	else if ((outfile = fopen(outfilename,"w")) == NULL)
	{
	    fprintf(stderr,"Can't open output file %s\n",outfilename);
	    gt_abort(NULL);
	}

	if (uswitch)
	    outcode = 0;
	else if (sswitch || (!gswitch && (codetype&SPARSE6)))
	    outcode = SPARSE6;
	else 
	    outcode = GRAPH6;

	if (!fswitch) fmt = NULL;

	if (!uswitch && (codetype&HAS_HEADER))
	{
	    if (outcode == SPARSE6) writeline(outfile,SPARSE6_HEADER);
	    else    		    writeline(outfile,GRAPH6_HEADER);
	}

	nin = 0;
	orbtotal = 0;
	unorbtotal = 0.0;
	t = CPUTIME;
	if (Sswitch)
	{
	    SG_INIT(sg);
	    SG_INIT(sh);
	    while (TRUE)
	    {
		if (read_sg_loops(infile,&sg,&loops) == NULL) break;
		++nin;
		n = sg.nv;
		m = (n + WORDSIZE - 1) / WORDSIZE;
		SG_ALLOC(sh,n,sg.nde,"labelg");
		for (ii = 0; ii < secret; ++ii)
		    fcanonise_inv_sg(&sg,m,n,&sh,fmt,iswitch?distances_sg:NULL,
		                 mininvarlevel,maxinvarlevel,invararg,loops>0);
		sortlists_sg(&sh);
	        orbtotal += gt_numorbits;
	        unorbtotal += 1.0 / gt_numorbits;
	        if (outcode == SPARSE6) writes6_sg(outfile,&sh);
	        else if (outcode == GRAPH6) writeg6_sg(outfile,&sh);
	    }
	}
	else if (tswitch)
	{
	    SG_INIT(sg);
	    SG_INIT(sh);
	    traces_opts.getcanon = TRUE;
	    traces_opts.writeautoms = FALSE;
	    traces_opts.verbosity = 0;
	    traces_opts.outfile = stdout;
	    
	    while (TRUE)
	    {
		if (read_sg_loops(infile,&sg,&loops) == NULL) break;
		if (loops > 0) gt_abort(">E Traces does not allow loops\n");
		++nin;
		n = sg.nv;
	        DYNALLOC1(int,lab,lab_sz,n,"traces@labelg");
	        DYNALLOC1(int,ptn,ptn_sz,n,"traces@labelg");
	        DYNALLOC1(int,orbits,orbits_sz,n,"traces@labelg");
		SG_ALLOC(sh,n,sg.nde,"labelg");
		for (ii = 0; ii < n; ++ii) { lab[ii] = ii; ptn[ii] = 1; }
		ptn[n-1] = 0;
		for (ii = 0; ii < secret; ++ii)
		    Traces(&sg,lab,ptn,orbits,&traces_opts,&traces_stats,&sh);
		sortlists_sg(&sh);
	        orbtotal += gt_numorbits;
	        unorbtotal += 1.0 / gt_numorbits;
	        if (outcode == SPARSE6) writes6_sg(outfile,&sh);
	        else if (outcode == GRAPH6) writeg6_sg(outfile,&sh);
		SG_FREE(sh);
	    }
	}
	else
	{
	    while (TRUE)
	    {
	        if ((g = readg(infile,NULL,0,&m,&n)) == NULL) break;
	        ++nin;
#if !MAXN
	        DYNALLOC2(graph,h,h_sz,n,m,"labelg");
#endif
		loops = loopcount(g,m,n);
	        for (ii = 0; ii < secret; ++ii)
	            fcanonise_inv(g,m,n,h,fmt,invarproc[inv].entrypoint,
			        mininvarlevel,maxinvarlevel,invararg,loops>0);
	        orbtotal += gt_numorbits;
	        unorbtotal += 1.0 / gt_numorbits;
	        if (outcode == SPARSE6) writes6(outfile,h,m,n);
	        if (outcode == GRAPH6)  writeg6(outfile,h,m,n);
	        FREES(g);
	    }
	}
	t = CPUTIME - t;

#if LONG_LONG_COUNTERS
	if (Oswitch)
	    fprintf(stderr,">C orbit totals = %lld %15.8f\n",
			   orbtotal,unorbtotal);
        if (!qswitch)
            fprintf(stderr,
                    ">Z  %lld graphs labelled from %s to %s in %3.2f sec.\n",
                    nin,infilename,outfilename,t);
#else
	if (Oswitch)
	    fprintf(stderr,">C orbit totals = %ld %15.8f\n",
			   orbtotal,unorbtotal);
        if (!qswitch)
            fprintf(stderr,
                    ">Z  %ld graphs labelled from %s to %s in %3.2f sec.\n",
                    nin,infilename,outfilename,t);
#endif

	exit(0);
}
