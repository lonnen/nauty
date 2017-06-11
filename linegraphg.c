/* linegraphg.c  version 1.1; B D McKay, Sep 2012. */

#define USAGE "linegraphg [-q] [infile [outfile]]"

#define HELPTEXT \
" Take the linegraphs of a file of graphs.\n\
  Input graphs with no edges produce only a warning message.\n\
\n\
    The output file has a header if and only if the input file does.\n\
\n\
    -q  Suppress auxiliary information.\n"

/*************************************************************************/

#include "gtools.h" 

/**************************************************************************/

static void
linegraph(sparsegraph *g, sparsegraph *h)
/* h := linegraph of g */
{
    DYNALLSTAT(size_t,eno,eno_sz);   /* edge number */
    int *ge,*gd,*he,*hd;
    size_t *gv,*hv;
    int gnv,hnv;
    size_t i,j,k,gnde,hnde,xhnde,num;
    size_t hi,lo,mid,v,w;

    sortlists_sg(g);
    SG_VDE(g,gv,gd,ge);
    gnv = g->nv;
    gnde = g->nde;
    DYNALLOC1(size_t,eno,eno_sz,gnde,"linegraphg");

    hnv = gnde/2;
    if (hnv != gnde/2) gt_abort(">E linegraphg: too many input edges\n");

    hnde = 0;
    num = 0;
    for (i = 0; i < gnv; ++i)
    {
        xhnde = hnde;
        hnde += gd[i]*((size_t)gd[i]-1);
        if (hnde < xhnde) gt_abort(">E linegraphg: too many output edges\n");

        for (j = gv[i]; j < gv[i]+gd[i]; ++j)
        {
            if (ge[j] == i)
                gt_abort(">E linegraphg can't handle loops\n");
            else if (ge[j] > i)
                eno[j] = num++;
            else
            {
                lo = gv[ge[j]];
                hi = lo + gd[ge[j]] - 1;
                while (lo <= hi)
                {
                    mid = lo + (hi-lo)/2;
                    if (ge[mid] == i) break;
                    else if  (ge[mid] < i) lo = mid+1;
                    else hi = mid-1;
                }
		if (lo > hi)
		    gt_abort(">E linegraphg : binary search failed\n");
                eno[j] = eno[mid];
            }
        }
    }

    SG_ALLOC(*h,hnv,hnde,"linegraphg");
    h->nv = hnv;
    h->nde = hnde;
    SG_VDE(h,hv,hd,he);

    for (i = 0; i < hnv; ++i) hd[i] = 0;
    for (i = 0; i < gnv; ++i)
    {
        for (j = gv[i]; j < gv[i]+gd[i]; ++j)
            hd[eno[j]] += gd[i]-1;
    }

    hv[0] = 0;
    for (i = 1; i < hnv; ++i) hv[i] = hv[i-1] + hd[i-1];
    for (i = 0; i < hnv; ++i) hd[i] = 0;

    for (i = 0; i < gnv; ++i)
    {
        for (j = gv[i]; j < gv[i]+gd[i]-1; ++j)
        for (k = j+1; k < gv[i]+gd[i]; ++k)
        {
            v = eno[j]; w = eno[k];
            he[hv[v]+(hd[v]++)] = w;
            he[hv[w]+(hd[w]++)] = v;
        }
    }
}

/**************************************************************************/

int
main(int argc, char *argv[])
{
    char *infilename,*outfilename;
    FILE *infile,*outfile;
    boolean badargs,quiet;
    int j,m,n,argnum;
    int codetype,outcode;
    sparsegraph g,h;
    nauty_counter nin,nullgraphs;
    char *arg,sw;
    static graph *gq;
    double t;

    HELP;

    SG_INIT(g);
    SG_INIT(h);

    infilename = outfilename = NULL;
    quiet = FALSE;

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
                SWBOOLEAN('q',quiet)
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

    if (badargs)
    {
        fprintf(stderr,">E Usage: %s\n",USAGE);
        GETHELP;
        exit(1);
    }

    if (!quiet)
    {
        fprintf(stderr,">A linegraphg");
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

    if (codetype&SPARSE6) outcode = SPARSE6;
    else                  outcode = GRAPH6;

    if (codetype&HAS_HEADER)
    {
        if (outcode == SPARSE6) writeline(outfile,SPARSE6_HEADER);
        else                    writeline(outfile,GRAPH6_HEADER);
    }

    nin = 0;
    nullgraphs = 0;
    t = CPUTIME;
    while (read_sg(infile,&g))
    {
        ++nin;

        if (g.nde > 0)
        {
            linegraph(&g,&h);
            if (outcode == SPARSE6) writes6_sg(outfile,&h);
            else                    writeg6_sg(outfile,&h);
        }
        else
            ++nullgraphs;
    }
    t = CPUTIME - t;

    if (!quiet)
    {
        if (nullgraphs > 0)
            fprintf(stderr,">W " COUNTER_FMT " null graphs not written\n",
                    nullgraphs);

        fprintf(stderr,">Z " COUNTER_FMT
                " graphs converted from %s to %s in %3.2f sec.\n",
                nin,infilename,outfilename,t);
    }

    exit(0);
}
