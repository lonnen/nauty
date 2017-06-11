/******************************************************************************
 *                                                                            *
 * This is the main file for traces() version 2.0, which is included into     *
 *   nauty() version 2.5.                                                     *
 *                                                                            *
 *   nauty is Copyright (1984-2014) Brendan McKay. All rights reserved.       *
 *   Subject to the waivers and disclaimers in nauty.h.                       *
 *   Traces is Copyright (2008-2014) Adolfo Piperno. All rights reserved.     *
 *                                                                            *
 *   CHANGE HISTORY                                                           *
 *       28-Dec-12 : final changes for version 2.0                            *
 *       29-Mar-13 : bug correction in automorphism mode                      *
 *       21-May-13 : bug correction (coloured lists)                          *
 *       29-Jun-13 : bug correction (coloured lists and cycles)               *
 *       07-Dec-13 : bug correction in automorphism mode (wrong group size    *
 *                   due to randomness in Schreier-Sims orbit computation)    *
 *                   bug correction (discrete initial partition)              *
 *       15-Feb-14 : CPUDEFS removed (already declared in gtools.h)           *
 *****************************************************************************/

#include "traces.h"
#define SORT_OF_SORT 3
#define SORT_NAME sortindirect
#include "sorttemplates.c"

typedef struct Candidate {
    boolean sortedlab;
    int *invlab;
    int *lab;
    int code;
    int do_it;
    int indnum;
    int name;
    int vertex;
    struct Candidate *next;
    struct searchtrie *stnode;
    unsigned int firstsingcode;
    unsigned int pathsingcode;
    unsigned int singcode;
} Candidate;

typedef struct Partition {
    int *cls;
    int *inv;
    int active;
    int cells;
    int code;
} Partition;

typedef struct trielist {
    struct searchtrie *triearray;
    struct trielist *prev;
    struct trielist *next;
} trielist;

typedef struct TracesVars {
    char digstring[25];
	double autchk;
	double expaths;
	double schreier1;
	double schreier2;
	double schreier3;
    int *currorbit;
    int *orbits;
    int answ;
    int brkstpcount;
    int compstage;
    int cand_level;
    int canlist;
    int digits;
    int expathlength;
    int firstpathlength;
    int fromlevel;
    int indivend;
    int indivstart;
    int linelgth;
    int mark;
    int maxtreelevel;
    int maxspineorblevel;
    int name;
    struct searchtrie *gotonode;
    struct searchtrie *newgotonode;
    struct searchtrie *newst_stage1;
    int newindex;
    int nextlevel;
    int nfix;
    int finalnumcells;
    int samepref;
    int smalldeglevel;
    int specialgens;
    int stackmark;
    int steps;
    int strategy;
    trielist *strielist;
    int strienext;
    int tcell_sz;
    int tcell;
    int tcellevel;
    int tcellexpath_sz;
    int tcellexpath;
    int tcellfromlevel;
    int tolevel_tl;
    int tolevel;
    int treedepth;
    int trienext;
    int triepos;
    TracesOptions *options;
    TracesStats *stats;
    unsigned int singlongcode;
} TracesVars;

typedef struct TracesInfo {
    boolean autofound;
    boolean exitfromref;
    boolean identitygroup;
    boolean minimalinorbits;
    boolean newtc;
    boolean thegraphisparse;
    boolean thegrouphaschanged;
    boolean thereisnextlevel;
    boolean useTempOrbits1;
    boolean useTempOrbits2;
} TracesInfo;

typedef struct TracesSpine {
    boolean thetracexists;
    boolean thexpathexists;
    Candidate *listend;
    Candidate *liststart;
    int ccend;
    int ccstart;
    int idx;
    int listcounter;
    int stpend;
    int stpstart;
    int tgtcell;
    int tgtend;
    int tgtfrom;
    int tgtpos;
    int tgtsize;
    int trcend;
    int trcstart;
    int updates;
    unsigned long keptcounter;
    unsigned long levelcounter;
    Partition *part;
    unsigned int singcode;
} TracesSpine;

typedef struct trie {
	int value;
	struct trie *first_child;
	struct trie *next_sibling;
} trie;

typedef struct searchtrie {
	int index;
	int name;
	int vtx;
    int level;
	struct searchtrie *father;
	struct searchtrie *first_child;
	struct searchtrie *last_child;
	struct searchtrie *next_sibling;
	struct searchtrie *goes_to;
} searchtrie;

typedef struct pair {
	int fst;
	int snd;
} pair;


static int  traces_refine(sparsegraph*, Candidate*, int, int, Partition*,
						  struct TracesVars*, struct TracesInfo*, boolean);
static void traces_refine_notrace(sparsegraph*, Candidate*, int, int, Partition*,
								  struct TracesVars*, struct TracesInfo*);
static void traces_refine_maketrie(sparsegraph*, Candidate*, int, int, Partition*,
								   struct TracesVars*, struct TracesInfo*);
static int  traces_refine_comptrie(sparsegraph*, Candidate*, int, int, Partition*,
								   struct TracesVars*, struct TracesInfo*);
static int  traces_refine_sametrace(sparsegraph*, Candidate*, int, int, Partition*,
									struct TracesVars*, struct TracesInfo*);
static int  traces_refine_refine(sparsegraph*, Candidate*, int, int, Partition*,
								 struct TracesVars*, struct TracesInfo*);
static void quickSort(int*, int);
static struct Candidate* NewCandidate(int, Candidate**, int);
static int FreeList(Candidate*, int);
static int FixBase(int*, struct TracesVars*, Candidate*, int, int);
static void factorial(double*, int*, int);
static void factorial2(double*, int*, int);
static int CheckForAutomorphisms(sparsegraph*, Candidate*, Candidate*,
                                 struct TracesVars*, struct TracesInfo*, int, int, Partition*);
static int CheckForSingAutomorphisms(sparsegraph*, Candidate*, Partition*, Candidate*,
									 struct TracesVars*, struct TracesInfo*,
									 int, int);
static int CheckForMatching(sparsegraph*, Candidate*, Candidate*, Partition*, struct TracesVars*, struct TracesInfo*, int, int);
static void Individualize(Partition*, Candidate*, int, int, int, int);
static boolean TreeFyTwo(sparsegraph*, int, int, Candidate*, Candidate*, Partition*, int, struct TracesVars*, struct TracesInfo*);
static void ExperimentalStep(sparsegraph*, Partition*, Candidate*, TracesVars*, TracesInfo*, int, int);
static boolean TargetCell(sparsegraph*, Candidate*, Partition*, int, TracesVars*, int);
static boolean TargetCellFirstPath(sparsegraph*, Candidate*, Partition*, int, TracesVars*);
static boolean TargetCellExpPath(sparsegraph*, Candidate*, Partition*, int, TracesVars*);
static boolean SelectNextLevel(sparsegraph*, int, TracesVars*);
static void CopyCand(Candidate*, Candidate*, int);
static struct trie* trie_new(int, struct TracesVars*);
static struct trie* trie_make(trie*, int, int, struct TracesVars*);
static struct trie* trie_comp(trie*, int);
static void RemoveFromLevel(int, int, int, boolean);
static void CompStage0(sparsegraph*, Partition*, Partition*, Candidate*, Candidate*,
					   int, int, struct TracesVars*, struct TracesInfo*);
static void CompStage1(sparsegraph*, Partition*, Partition*, Candidate*, Candidate*,
					   int, int,
					   struct TracesVars*, struct TracesInfo*);
static void CompStage2(sparsegraph*, Partition*, Partition*, Candidate*, Candidate*,
					   int, int,
					   struct TracesVars*, struct TracesInfo*);
static void grouporderplus(sparsegraph*, int*, int, Candidate*, schreier*, permnode**,
						   double*, int*, int, TracesVars*, TracesInfo*);

static boolean Prefix(Candidate*, Candidate*, int);
static boolean findperm(permnode*, int*, int, TracesVars*);
static int spinelementorbsize(int*, int*, int, int);
static void Propagate(sparsegraph*, Candidate*, Partition*, int, int, int);
static trielist* searchtrie_new(int, TracesVars*);
static searchtrie* searchtrie_make(Candidate*, Candidate*, int, TracesVars*);
static boolean lookup(searchtrie*);
static int* findcurrorbits(schreier*, int);

static const unsigned int fuzz1[] = {037541, 061532, 005257, 026416};
static const unsigned int fuzz2[] = {006532, 070236, 035523, 062437};

#define FUZZ1(x) ((x) ^ fuzz1[(x)&3])
#define FUZZ2(x) ((x) ^ fuzz2[(x)&3])

#define STATS_INIT stats_arg->grpsize1 = 1; \
stats_arg->grpsize2 = 0; \
stats_arg->numorbits = n; \
stats_arg->treedepth= 0; \
stats_arg->numgenerators = 0; \
stats_arg->numnodes = 1; \
stats_arg->interrupted = 0; \
stats_arg->canupdates = 0; \
stats_arg->peaknodes = 0; \

#define TIME_INIT tv->autchk = 0; \
tv->expaths = 0; \
tv->schreier1 = 0; \
tv->schreier2 = 0; \
tv->schreier3 = 0;

#define MASHCOMM(l, i) ((l) + FUZZ1(i))
#define MASHNONCOMM(l, i) (FUZZ2(l) + (i))
#define MASH(l, i) ((((l) ^ 065435) + (i)) & 077777)
#define MASH1(l, i) ((l + (i*i)) & 077777)
#define CLEANUP(l) ((int)((l) % 0x7FFF))
#define SS(n, sing, plur)  (n), ((n) == 1?(sing):(plur))

#define SETMARK(Arr, Mrk) if (Mrk > (NAUTY_INFINITY-2)) { memset(Arr, 0, n*sizeof(int)); Mrk = 0; } Mrk++;

#define COPYNODE(W, V) \
memcpy(W->lab, V->lab, n*sizeof(int)); \
memcpy(W->invlab, V->invlab, n*sizeof(int)); \
W->code = V->code; \
W->singcode = V->singcode; \
W->do_it = V->do_it;

#define NEXTLINE fprintf(outfile, "\n");

#define PRINTCHAR(c) fprintf(outfile, "%s", c);

#define PRINTCAND(V, Lev) PRINTCHAR(" ") for (tmp=1; tmp<=Lev; tmp++) {fprintf(outfile, tv->digstring, V->lab[Spine[tmp].tgtpos]+labelorg);}

#define PRINTCANDBIG(V, Lev) PRINTCHAR(" ") \
for (tmp=1; tmp<=5; tmp++) {fprintf(outfile, tv->digstring, V->lab[Spine[tmp].tgtpos]+labelorg);} \
fprintf(outfile, "... "); \
for (tmp=Lev-4; tmp<=Lev; tmp++) {fprintf(outfile, tv->digstring, V->lab[Spine[tmp].tgtpos]+labelorg);}

#define LINE(K, c) PRINTCHAR(c) for (tmp=1; tmp<=K; tmp++) {fprintf(outfile, c);}

#define TRACE_CHECK(Tr, Ind, Arg, End) TracePos = Tr+Ind; \
if (newtrace) { \
*TracePos = Arg; \
} \
else { \
if (Ind < *End) { \
if (*TracePos != Arg) { \
if (*TracePos > Arg) { \
return FALSE; \
} \
else { \
*TracePos = Arg; \
newtrace = TRUE; \
} \
} \
} \
else { \
*TracePos = Arg; \
newtrace = TRUE; \
} \
} \
Ind++;

#define SAMETRACE_CHECK(Tr, Ind, Arg, End) TracePos = Tr+Ind; \
if (Ind < *End) { \
if (*TracePos != Arg) { \
return FALSE; \
} \
} \
else { \
return FALSE; \
} \
Ind++;

#define NEWPART(P) P = malloc(sizeof(*(P))); \
if (P == NULL) { \
fprintf(ERRFILE, "\nError, memory not allocated.\n"); \
exit(1); \
} \
P->cls = malloc(n*sizeof(int)); \
if (P->cls == NULL) { \
fprintf(ERRFILE, "\nError, memory not allocated.\n"); \
exit(1); \
} \
P->inv = malloc(n*sizeof(int)); \
if (P->inv == NULL) { \
fprintf(ERRFILE, "\nError, memory not allocated.\n"); \
exit(1); \
} \
P->code = -1; \
P->cells = 0;

#define NEWPARTSPINE(Lev) if (Lev > 3) { \
Spine[Lev].part = malloc(sizeof(*(Spine[Lev].part))); \
if (Spine[Lev].part == NULL) { \
fprintf(ERRFILE, "\nError, memory not allocated.\n"); \
exit(1); \
} \
Spine[Lev].part->cls = Spine[Lev-3].part->cls; \
Spine[Lev].part->inv = Spine[Lev-3].part->inv; \
Spine[Lev-3].part->cls = Spine[Lev-3].part->inv = NULL; \
Spine[Lev].part->code = -1; \
Spine[Lev].part->cells = 0; \
} \
else { \
NEWPART(Spine[Lev].part) \
}

#define FREEPART(Part) if (Part) { \
if (Part->cls) free(Part->cls); \
if (Part->inv) free(Part->inv); \
free(Part); \
}

#define FREECAND(Cand) if (Cand) { \
free(Cand->lab); \
free(Cand->invlab); \
free(Cand); \
}


#define COPYPART(P, Q) memcpy(P->cls, Q->cls, n*sizeof(int)); \
memcpy(P->inv, Q->inv, n*sizeof(int)); \
P->cells = Q->cells; \
P->code = Q->code; \

#define ADDTONEXTLEVEL if (SpineTL->listend) { \
(SpineTL->listend)->next = NewCandidate(n, &GarbList, TRUE); \
if ((tv->compstage < 2) && (SpineTL->listcounter <= (NAUTY_INFINITY-2))) SpineTL->listcounter++; \
SpineTL->listend = (SpineTL->listend)->next; \
CopyCand(SpineTL->listend, NextCand, n); \
} \
else { \
SpineTL->liststart = NewCandidate(n, &GarbList, TRUE); \
if (tv->compstage < 2) SpineTL->listcounter = 1; \
SpineTL->listend = SpineTL->liststart; \
CopyCand(SpineTL->liststart, NextCand, n); \
}

#define ORBITSIZES memset(OrbSize, 0, n*sizeof(int)); \
for (i=0; i<n; i++) { \
OrbSize[tv->orbits[i]]++; \
}

#define CURRORBITSIZES memset(CurrOrbSize, 0, n*sizeof(int)); \
for (i=SpineTL->tgtcell; i<SpineTL->tgtend; i++) { \
CurrOrbSize[tv->currorbit[CurrCand->lab[i]]]++; \
}

#define EXITFROMSTAGE0REFINE PRINT_LINE_PLUS(tv->tolevel) \
if (tv->options->verbosity >= 2) fprintf(outfile, "-=="); \
CurrCand->indnum--; \
RemoveFromLevel(tv->tolevel, tv->treedepth, tv->strategy, FALSE); \
tv->compstage = 1; \
trieroot = trie_new(n, tv); \
trieref = trieroot; \
tv->nextlevel = tv->maxtreelevel = tv->fromlevel; \
ti->thereisnextlevel = TRUE; \
ti->exitfromref = TRUE; \
return;

#define EXITFROMSTAGE0EXPATH2 PRINT_LINE_PLUS(tv->tolevel) \
if (tv->options->verbosity >= 2) fprintf(outfile, "=-="); \
tv->compstage = 1; \
trieroot = trie_new(n, tv); \
trieref = trieroot; \
tv->nextlevel = tv->maxtreelevel = tv->tolevel; \
ti->thereisnextlevel = TRUE; \
ti->exitfromref = FALSE; \
return;

#define EXITFROMSTAGE0EXPATH1 PRINT_RETURN PRINT_LINE_PLUS(tv->tolevel) \
if (tv->options->verbosity >= 2) fprintf(outfile, "==-"); \
if (SpineTL->liststart) { \
AuxCand = SpineTL->liststart; \
SpineTL->liststart = NewCandidate(n, &GarbList, TRUE); \
CopyCand(SpineTL->liststart, NextCand, n); \
SpineTL->liststart->next = AuxCand; \
} \
else { \
SpineTL->liststart = NewCandidate(n, &GarbList, TRUE); \
SpineTL->listend = SpineTL->liststart; \
SpineTL->liststart->next = NULL; \
CopyCand(SpineTL->liststart, NextCand, n); \
} \
tv->compstage = 1; \
trieroot = trie_new(n, tv); \
trieref = trieroot; \
tv->nextlevel = tv->maxtreelevel = tv->tolevel; \
ti->thereisnextlevel = TRUE; \
ti->exitfromref = FALSE; \
return;

#define UPDATE_LINELGTH if (tv->options->verbosity >= 2) { \
if (tv->tolevel < 12) { \
tv->linelgth = (tv->digits+1)*tv->tolevel+16; \
} \
else { \
tv->linelgth = (tv->digits+1)*10+20; \
} \
}

#define PRINT_LINE if ((tv->options->verbosity >= 1) && (tv->strategy == 0)) { \
if (!ti->newtc) \
{ \
if (tv->options->verbosity >= 2) { LINE(tv->linelgth, "-"); \
NEXTLINE} \
} \
ti->newtc = FALSE; \
}


#define PRINT_LINE_PLUS(Lev) if ((tv->options->verbosity >= 1) && (tv->strategy == 0)) { \
if (!ti->newtc) \
{ \
if (tv->options->verbosity >= 2) { LINE(tv->linelgth, "-"); \
fprintf(outfile, " ");} \
fprintf(outfile, "level %d:  %d cell%s; target cell: %d; %d orbit%s; %lu node%s (%lu kept); %d update%s;", \
Lev, SS(Spine[Lev].part->cells, "", "s"), Spine[Lev].tgtsize, SS(tv->stats->numorbits, "", "s"), \
SS(Spine[Lev].levelcounter, "", "s"), Spine[Lev].keptcounter, SS(Spine[Lev].updates, "", "s")); \
NEXTLINE \
} \
ti->newtc = FALSE; \
}

#define PRINT_CANDIDATE(Cand, Lev) { \
for (tmp = Cand->name, cu = 0; tmp > 0; tmp /= 10, ++cu) {} \
for (tmp = Lev, cu1 = 0; tmp > 0; tmp /= 10, ++cu1) {} \
cu = 14-cu-cu1; \
LINE(cu, "-") \
fprintf(outfile, " %d, %d) ", Lev % 10000, Cand->name % 10000000); \
if (Lev < 12) { \
PRINTCAND(Cand, Lev) \
} \
else { \
PRINTCANDBIG(Cand, Lev) \
} \
PRINTCHAR("| ")	\
fprintf(outfile, "{%x, %x} ", Cand->code, Cand->singcode); \
}

#define PRINT_EXPPATHSTEP(Cand, Boo) if (tv->options->verbosity >= 2) { \
if ((tv->tolevel_tl-tv->tolevel < 6) || (NextPart->cells == tv->finalnumcells)) { \
fprintf(outfile, "%d ", Cand->lab[SpineTL_tl->tgtpos]+labelorg); \
if (tv->options->verbosity >= 2) {if (Boo) fprintf(outfile, "{%x} ", Cand->code); else fprintf(outfile, "{interr.(%d)} ", NextPart->cells);} \
else { if (!Boo) fprintf(outfile, "{interr.(%d)} ", NextPart->cells);} \
if ((NextPart->cells == tv->finalnumcells) || (tv->tolevel_tl == tv->smalldeglevel)) { \
fprintf(outfile, "(%d) ", tv->tolevel_tl); \
} \
} \
else { \
if (tv->tolevel_tl-tv->tolevel == 6) { \
fprintf(outfile, "... "); \
} \
} \
}

#define PRINT_RETURN if (tv->options->verbosity >= 2) { \
fprintf(outfile, "\n"); \
}

#define PRINT_FROM_VERB(Verb) if (tv->options->verbosity >= Verb) { \
fprintf(outfile, "FROM: "); \
PRINTCAND(CurrCand, tv->fromlevel) \
fprintf(outfile, " do_it: %d, indnum: %d, stnode->index: %d ", CurrCand->do_it, CurrCand->indnum, CurrCand->stnode->index); \
PRINT_RETURN \
}

#define PRINT_NOTMIN_VERB(Verb) if (tv->options->verbosity >= Verb) { \
fprintf(outfile, " is NOT minimal in orbits (1, %d) [%d]; ", gom_level, CurrCand->lab[Spine[gom_level+1].tgtpos]+labelorg); \
fprintf(outfile, "at lev %d, orb[%d] = %d.\n", gom_level+1, CurrCand->lab[Spine[gom_level+1].tgtpos]+labelorg, tv->currorbit[CurrCand->lab[Spine[gom_level+1].tgtpos]]+labelorg); }

#define PRINT_SKIPPED_VERB(Verb) if (tv->options->verbosity >= Verb) \
fprintf(outfile, " skipped (0) (orbit[%d] = %d)\n", \
NextCand->vertex+labelorg, tv->currorbit[NextCand->vertex]+labelorg);

#define PRINT_REFINE_VERB(Verb) if (tv->options->verbosity >= Verb) \
fprintf(outfile, " REFINE (orbit[%d] = %d)\n", NextCand->vertex+labelorg, tv->currorbit[NextCand->vertex]+labelorg);

#define PRINT_INDIV_VERB(Verb) if (tv->options->verbosity >= Verb) { \
PRINTCAND(CurrCand, tv->fromlevel) \
fprintf(outfile, "| "); \
fprintf(outfile, tv->digstring, NextCand->vertex+labelorg); \
}

#define SPECIALGENERATORS if (tv->options->generators) addpermutation(ring, AUTPERM, n); \
tv->stats->numgenerators++; \
tv->specialgens++; \
if (tv->options->writeautoms) { \
fprintf(outfile, "Gen #%d: ", tv->stats->numgenerators); \
writeperm(outfile, AUTPERM, tv->options->cartesian, tv->options->linelength, n); \
} \
if (tv->options->userautomproc) { \
(*tv->options->userautomproc)(tv->stats->numgenerators, AUTPERM, n); \
}

#define UPDATEMIN(A, B) if (B < A) A = B;

#if !MAXN
DYNALLSTAT(int, AUTPERM, AUTPERM_sz);
DYNALLSTAT(int, BreakSteps, BreakSteps_sz);
DYNALLSTAT(int, CurrOrbSize, CurrOrbSize_sz);
DYNALLSTAT(int, CurrRefCells, CurrRefCells_sz);
DYNALLSTAT(int, fix, fix_sz);
DYNALLSTAT(int, IDENTITY_PERM, IDENTITY_PERM_sz);
DYNALLSTAT(int, Markers, Markers_sz);
DYNALLSTAT(int, MarkHitVtx, MarkHitVtx_sz);
DYNALLSTAT(int, MultRefCells, MultRefCells_sz);
DYNALLSTAT(int, NghCounts, NghCounts_sz);
DYNALLSTAT(int, OrbSize, OrbSize_sz);
DYNALLSTAT(int, RefCells, RefCells_sz);
DYNALLSTAT(searchtrie*, RefPath, RefPath_sz);
DYNALLSTAT(int, SplCls, SplCls_sz);
DYNALLSTAT(int, SplCnt, SplCnt_sz);
DYNALLSTAT(int, SplPos, SplPos_sz);
DYNALLSTAT(int, StackMarkers, StackMarkers_sz);
DYNALLSTAT(int, TheTrace, TheTrace_sz);
DYNALLSTAT(int, TheTraceCC, TheTraceCC_sz);
DYNALLSTAT(int, TheTraceSplNum, TheTraceSplNum_sz);
DYNALLSTAT(int, TheTraceSteps, TheTraceSteps_sz);
DYNALLSTAT(int, TEMPLAB, TEMPLAB_sz);
DYNALLSTAT(int, TEMPINVLAB, TEMPINVLAB_sz);
DYNALLSTAT(int, WorkArray, WorkArray_sz);
DYNALLSTAT(int, WorkArray1, WorkArray1_sz);
DYNALLSTAT(int, WorkArray2, WorkArray2_sz);
DYNALLSTAT(int, WorkArray3, WorkArray3_sz);
DYNALLSTAT(int, WorkArray4, WorkArray4_sz);
DYNALLSTAT(int, WorkArray5, WorkArray5_sz);
DYNALLSTAT(int, TreeStack, TreeStack_sz);
DYNALLSTAT(TracesSpine, Spine, Spine_sz);
DYNALLSTAT(trie*, TrieArray, TrieArray_sz);
#else
static TLS_ATTR int AUTPERM[MAXN];
static TLS_ATTR int BreakSteps[MAXN];
static TLS_ATTR int CurrOrbSize[MAXN];
static TLS_ATTR int CurrRefCells[MAXN];
static TLS_ATTR int fix[MAXN];
static TLS_ATTR int IDENTITY_PERM[MAXN];
static TLS_ATTR int Markers[MAXN];
static TLS_ATTR int MarkHitVtx[MAXN];
static TLS_ATTR int MultRefCells[MAXN];
static TLS_ATTR int NghCounts[MAXN];
static TLS_ATTR int OrbSize[MAXN];
static TLS_ATTR int RefCells[MAXN];
static TLS_ATTR searchtrie* RefPath[MAXN];
static TLS_ATTR int SplCls[MAXN];
static TLS_ATTR int SplCnt[MAXN];
static TLS_ATTR int SplPos[MAXN];
static TLS_ATTR int StackMarkers[MAXN];
static TLS_ATTR int TheTrace[MAXN];
static TLS_ATTR int TheTraceCC[MAXN];
static TLS_ATTR int TheTraceSplNum[MAXN];
static TLS_ATTR int TheTraceSteps[MAXN];
static TLS_ATTR int TEMPLAB[MAXN];
static TLS_ATTR int TEMPINVLAB[MAXN];
static TLS_ATTR int WorkArray[MAXN];
static TLS_ATTR int WorkArray1[MAXN];
static TLS_ATTR int WorkArray2[MAXN];
static TLS_ATTR int WorkArray3[MAXN];
static TLS_ATTR int WorkArray4[MAXN];
static TLS_ATTR int WorkArray5[MAXN];
static TLS_ATTR int TreeStack[MAXN];
static TLS_ATTR TracesSpine Spine[MAXN];
static TLS_ATTR trie* TrieArray[MAXN];
#endif

static TLS_ATTR FILE *outfile;

/* Brendan's SCHREIER */
static TLS_ATTR schreier  *gpB;				/* This will point to the Schreier structure */
static TLS_ATTR permnode  *gensB;			/* This will point to the stored generators */

static TLS_ATTR Candidate *GarbList, *SpOrd, *SpCyc, *SpSwp;
static TLS_ATTR Partition *SpPart1, *SpPart2;
static TLS_ATTR TracesSpine *SpineTL, *SpineFL, *SpineTL_tl;
static TLS_ATTR trie *trieroot, *trieref;
static TLS_ATTR int *TempOrbits = NULL;

static int
given_gens(sparsegraph *g, permnode *gens, int *orbits, boolean digraph)
/* Check if the permutations in the list gens are automorphisms,
 * also set mark and refcount fields and initialise orbits. */
{
	int i, m, n, norbs;
	permnode *pn;
	
	n = g->nv;
	for (i = 0; i < n; ++i) orbits[i] = i;
	memcpy(IDENTITY_PERM, orbits, n*sizeof(int));
	norbs = n;
	
	if (!gens) return norbs;
	
	m = SETWORDSNEEDED(n);
	pn = gens;
	do {
		if (!isautom_sg((graph*)g, pn->p, digraph, m, n)) {
			fprintf(ERRFILE, "Input permutation is not an automorphism\n");
			exit(1);
		}
	    norbs = orbjoin(orbits, pn->p, n);
	    pn->mark = 1;
	    pn->refcount = 0;
	    pn = pn->next;
	} while (pn != gens);
	
	return norbs;
}

void
Traces(sparsegraph *g_arg, int *lab, int *ptn,
	   int *orbits_arg, TracesOptions *options_arg, TracesStats *stats_arg,
	   sparsegraph *canong_arg)
{
	int i, j, cu, cu1;
	int temp, tmp;
    trielist *STStart, *STAux;
    searchtrie *TrieNode;
	
    if (g_arg->nv > (NAUTY_INFINITY-2))
    {
        fprintf(ERRFILE, "Traces: need n <= %d, but n=%d\n\n",
                NAUTY_INFINITY-2, g_arg->nv);
        return;
    }
    
	struct TracesVars *tv = malloc(sizeof(struct TracesVars));
    if (tv == NULL) {
        fprintf(ERRFILE, "\nError, memory not allocated.\n");
        exit(1);
    }
	struct TracesInfo *ti = malloc(sizeof(struct TracesInfo));
    if (ti == NULL) {
        fprintf(ERRFILE, "\nError, memory not allocated.\n");
        exit(1);
    }
	
    Partition *CurrPart, *NextPart;
    Candidate *CurrCand, *NextCand;
    Candidate *BestCand, *AuxCand;
	
	const int n = g_arg->nv;
	const int m = SETWORDSNEEDED(n);
	trieroot = NULL;
	NextCand = GarbList = NULL;
	
    tv->canlist = 0;
    tv->compstage = 0;
    tv->expathlength = n;
    tv->finalnumcells = n;
    tv->firstpathlength = 0;
    tv->linelgth = 0;
    tv->name = 0;
    tv->maxspineorblevel = 0;
    tv->nfix = 0;
    tv->options = options_arg;
    tv->orbits = orbits_arg;
    tv->smalldeglevel = n;
    tv->specialgens = 0;
    tv->stats = stats_arg;
    tv->treedepth = 0;
    tv->gotonode = NULL;
    
	if (tv->options->strategy == 0) {
		tv->steps = n;
		tv->strategy = 0;
	}
	else {
		tv->strategy = 1;
        tv->steps = tv->options->strategy;
        if (tv->steps > n) {
            tv->steps = n;
        }
	}
    
#if !MAXN
    DYNALLOC1(int, AUTPERM, AUTPERM_sz, n, "Traces");
    DYNALLOC1(int, BreakSteps, BreakSteps_sz, n, "Traces");
    DYNALLOC1(int, CurrOrbSize, CurrOrbSize_sz, n, "Traces");
    DYNALLOC1(int, CurrRefCells, CurrRefCells_sz, n, "Traces");
    DYNALLOC1(int, fix, fix_sz, n, "Traces");
    DYNALLOC1(int, IDENTITY_PERM, IDENTITY_PERM_sz, n, "Traces");
	DYNALLOC1(int, Markers, Markers_sz, n, "Traces");
	DYNALLOC1(int, MarkHitVtx, MarkHitVtx_sz, n, "Traces");
    DYNALLOC1(int, MultRefCells, MultRefCells_sz, n, "Traces");
	DYNALLOC1(int, NghCounts, NghCounts_sz, n, "Traces");
    DYNALLOC1(int, OrbSize, OrbSize_sz, n, "Traces");
    DYNALLOC1(int, RefCells, RefCells_sz, n, "Traces");
	DYNALLOC1(int, SplCls, SplCls_sz, n, "Traces");
	DYNALLOC1(int, SplCnt, SplCnt_sz, n, "Traces");
	DYNALLOC1(int, SplPos, SplPos_sz, n, "Traces");
	DYNALLOC1(int, StackMarkers, StackMarkers_sz, n, "Traces");
    DYNALLOC1(int, TheTrace, TheTrace_sz, n+10, "Traces");
    DYNALLOC1(int, TheTraceCC, TheTraceCC_sz, n, "Traces");
    DYNALLOC1(int, TheTraceSplNum, TheTraceSplNum_sz, n, "Traces");
    DYNALLOC1(int, TheTraceSteps, TheTraceSteps_sz, n+10, "Traces");
    DYNALLOC1(int, TEMPLAB, TEMPLAB_sz, n+10, "Traces");
    DYNALLOC1(int, TEMPINVLAB, TEMPINVLAB_sz, n+10, "Traces");
	DYNALLOC1(int, WorkArray, WorkArray_sz, n, "Traces");
    DYNALLOC1(int, WorkArray1, WorkArray1_sz, n, "Traces");
    DYNALLOC1(int, WorkArray2, WorkArray2_sz, n, "Traces");
    DYNALLOC1(int, WorkArray3, WorkArray3_sz, n, "Traces");
    DYNALLOC1(int, WorkArray4, WorkArray4_sz, n, "Traces");
    DYNALLOC1(int, WorkArray5, WorkArray5_sz, n, "Traces");
    DYNALLOC1(int, TreeStack, TreeStack_sz, n, "Traces");
    DYNALLOC1(TracesSpine, Spine, Spine_sz, n, "Traces");
    DYNALLOC1(trie*, TrieArray, TrieArray_sz, n, "Traces");
#endif
    
	outfile = (tv->options->outfile == NULL ? stdout : tv->options->outfile);
    
	SpOrd = SpCyc = SpSwp = NULL;
	SpPart1 = SpPart2 = NULL;
	
    if (tv->options->verbosity >= 2) {
        for (i = n, tv->digits = 0; i > 0; i /= 10, ++tv->digits) {}
        sprintf(tv->digstring, "%s%dd ", "%", tv->digits);
    }
	
#define PERMSTACK WorkArray1
#define CYCLES WorkArray1
#define HitCls WorkArray2
#define CYCOLR WorkArray2
#define HitVtx WorkArray3
#define CYLGTH WorkArray3
#define CStack WorkArray4
#define CYMULT WorkArray4
#define HitCount WorkArray5
#define ElmHitCll WorkArray5
#define CYCPOS WorkArray5
#define CYCHIT AUTPERM
#define LGHATTR RefCells
#define CYCREP MultRefCells
#define TempOrbSize TEMPLAB
#define AutomCount TEMPINVLAB
    
	/* Initialize group and statistics */
    STATS_INIT
    if (tv->options->verbosity >= 2) {
        TIME_INIT
    }
    
	if (tv->options->generators) {
		tv->stats->numorbits = given_gens(g_arg, *tv->options->generators, orbits_arg, tv->options->digraph);
		newgroup(&gpB, NULL, n);
		gensB = *tv->options->generators;
		expandschreier(gpB, &gensB, n);
		ti->thegrouphaschanged = TRUE;
		ti->identitygroup = FALSE;
	}
	else {
		newgroup(&gpB, &gensB, n);
		for (i = 0; i < n; ++i) orbits_arg[i] = i;
		memcpy(IDENTITY_PERM, orbits_arg, n*sizeof(int));
		tv->stats->numorbits = n;
		ti->thegrouphaschanged = FALSE;
		ti->identitygroup = TRUE;
	}
	tv->currorbit = gpB->orbits;
	
	memset(fix, 0, n*sizeof(int));
	memset(TheTraceCC, 0, n*sizeof(int));
	/* ran_init(1234);  any long int as an argument */
    
	/* The graph is sparse? */
	if (g_arg->nde < n || g_arg->nde / n < n / (g_arg->nde / n)) {
		ti->thegraphisparse = TRUE;
	}
	else {
		ti->thegraphisparse = FALSE;
	}
	
	/* Initialize candidate, partition, cells, orbits */
	NEWPART(Spine[0].part);
	Spine[0].part->code = -1;
	
	NEWPART(NextPart);
	CurrPart = Spine[0].part;
	
	CurrCand = NewCandidate(n, &GarbList, TRUE);
	memset(CurrPart->inv, 0, n*sizeof(int));
	
	CurrCand->singcode = 0;
    TempOrbits = NULL;
    STStart = NULL;
    
	if (tv->options->defaultptn) {
		memcpy(CurrCand->lab, IDENTITY_PERM, n*sizeof(int));
		memcpy(CurrCand->invlab, IDENTITY_PERM, n*sizeof(int));
		CurrPart->cells = 1;
		CurrPart->cls[0] = n;
		TheTrace[0] = 0;
	}
	else {
		memcpy(CurrCand->lab, lab, n*sizeof(int));
		CurrPart->cells = 0;
		j = 0;
		for (i = 0; i < n; i++) {
			if (j) {
				CurrPart->inv[i] = j;
			}
			CurrCand->invlab[CurrCand->lab[i]] = i;
			if (!ptn[i]) {
				CurrPart->cls[j] = i-j+1;
				if (CurrPart->cls[j] == 1) {
					CurrCand->singcode = MASHCOMM(CurrCand->singcode, CurrCand->lab[j]);
				}
				TheTrace[CurrPart->cells++] = j;
				j = i+1;
			}
		}
	}
    
	/* Initialization of Spine structure */
	SpineFL = Spine;
	SpineFL->tgtcell = SpineFL->tgtpos = 0;
	SpineFL->tgtend = n;
	SpineFL->trcstart = 0;
	SpineFL->trcend = CurrPart->active = CurrPart->cells;
	SpineFL->ccstart = SpineFL->ccend = 0;
	SpineFL->stpstart = SpineFL->stpend = 0;
	SpineFL->thetracexists = FALSE;
	SpineFL->liststart = SpineFL->listend = CurrCand;
    SpineFL->levelcounter = 1;
    SpineFL->keptcounter = 1;
    SpineFL->updates = 1;
	
	/* Further initializations */
    tv->mark = tv->stackmark = 1073741825;
	tv->maxtreelevel = 0;
	tv->tolevel = 0;
	tv->tcell = 0;
	UPDATE_LINELGTH
    
	/* First refinement */
	traces_refine(g_arg, CurrCand, m, n, CurrPart, tv, ti, FALSE);
    
    CurrCand->name = 0;
    
	if (CurrPart->cells == n) {
		tv->stats->canupdates++;
        /* CANONICAL FORM ? */
		if (tv->options->getcanon) {
			memcpy(lab, CurrCand->lab, n*sizeof(int));
			if (canong_arg) memcpy(CurrPart->inv, CurrCand->invlab, n*sizeof(int));
		}
	}
	else {
        STStart = searchtrie_new(n, tv);
        CurrCand->stnode = tv->strielist->triearray;
		NextCand = NewCandidate(n, &GarbList, TRUE);
		SpineFL->listcounter = 1;
		tv->tcellevel = 0;
		ti->newtc = FALSE;
        ti->exitfromref = FALSE;
        Spine[1].levelcounter = 0;
        Spine[1].updates = 0;
        
		do {
			tv->fromlevel = tv->tolevel;
			SpineFL = Spine+tv->fromlevel;
            
			if (CurrCand) {
				switch (tv->compstage) {
					case 0: CompStage0(g_arg, CurrPart, NextPart, CurrCand, NextCand, m, n, tv, ti);
						break;
					case 1:
                        if (!TempOrbits) {
                            memcpy(WorkArray1, tv->currorbit, n*sizeof(int));
                            TempOrbits = WorkArray1;
                        }
                        memset(TempOrbSize, 0, n*sizeof(int));
                        for (i=SpineTL->tgtcell; i<SpineTL->tgtend; i++) {
                            TempOrbSize[TempOrbits[CurrCand->lab[i]]]++;
                        }
                        CompStage1(g_arg, CurrPart, NextPart, CurrCand, NextCand, m, n, tv, ti);
                        break;
					case 2: CompStage2(g_arg, CurrPart, NextPart, CurrCand, NextCand, m, n, tv, ti);
                        break;
					default:
						break;
				}
			}
			
			/* NEXT CANDIDATE */
			if (ti->thereisnextlevel) {
				if (tv->nextlevel != tv->fromlevel) {
					UPDATE_LINELGTH
				}
				tv->tolevel = tv->nextlevel;
				CurrCand = Spine[tv->nextlevel].liststart;
				CurrPart = Spine[tv->nextlevel].part;
			}
		}
		while (ti->thereisnextlevel);
        
        if (tv->compstage) {
            memset(CurrOrbSize, 0, n*sizeof(int));
            for (i=0; i<n; i++) {
                CurrOrbSize[TempOrbits[i]]++;
            }
        }
        
        if (!tv->options->getcanon) {
			if (tv->compstage) {
                tv->maxtreelevel++;
                TrieNode = Spine[tv->maxtreelevel].liststart->stnode;
                if (TrieNode->father) {
                    TrieNode = TrieNode->father;
                }
                if (!ti->exitfromref) {
                    TrieNode->index = 0;
                    for (i=1; i<AutomCount[0]; i++) {
                        TrieNode->index += CurrOrbSize[AutomCount[i]];
                    }
                }
            }
        }
        
        tv->stats->grpsize1 = 1.0;
        tv->stats->grpsize2 = 0;
        
        if (tv->maxtreelevel) PRINT_LINE_PLUS(tv->maxtreelevel);
        
        AuxCand = Spine[tv->maxtreelevel].liststart;
        while (!AuxCand->do_it) {
            AuxCand = AuxCand->next;
        }
        
        /* GROUP ORDER */
        grouporderplus(g_arg, fix, tv->nfix, AuxCand, gpB, &gensB, &(tv->stats->grpsize1), &(tv->stats->grpsize2),
                       n, tv, ti);
        
        /* CANONICAL FORM ? */
        if (tv->options->getcanon) {
            BestCand = AuxCand;
            AuxCand = AuxCand->next;
            while (AuxCand) {
                if (AuxCand->do_it) {
                    if (comparelab_tr(g_arg, BestCand->lab, BestCand->invlab, AuxCand->lab, AuxCand->invlab) == 1) {
                        BestCand = AuxCand;
                        tv->stats->canupdates++;
                    }
                }
                AuxCand = AuxCand->next;
            }
            if (tv->options->verbosity >= 2) {
                LINE(tv->linelgth+4, "*")
                NEXTLINE
                fprintf(outfile, "Canonical:");
                PRINTCAND(BestCand, tv->maxtreelevel)
                PRINT_RETURN
            }
            memcpy(lab, BestCand->lab, n*sizeof(int));
            if (canong_arg) memcpy(CurrPart->inv, BestCand->invlab, n*sizeof(int));
        }
        
        if (tv->options->verbosity >= 2) {
            if (tv->linelgth < 40) {
                tv->linelgth = 40;
            }
            LINE(tv->linelgth+4, "*")
            NEXTLINE
        }
        
    }
    
    tv->stats->treedepth = tv->treedepth;
    if (Spine[tv->treedepth].part->code == -1) {
        tv->stats->treedepth--;
    }
    
    if (tv->options->verbosity >= 2) {
        fprintf(outfile, "group time: %.2f, %.2f, %.2f; total: %.2f; exp_paths time: %.2f; aut_check time: %.2f\n%lu refinement%s interrupted by trace comparison (%s); special cells: %d\n------", tv->schreier1, tv->schreier2, tv->schreier3, tv->schreier1+tv->schreier2+tv->schreier3,
                tv->expaths, tv->autchk, SS(tv->stats->interrupted, "", "s"), (tv->options->strategy == 0 ? "breadth-first" : "depth-first"), tv->specialgens);
        PRINT_RETURN
    }
    if (tv->options->getcanon && canong_arg) {
        canong_arg->nv  = g_arg->nv;
        canong_arg->nde  = g_arg->nde;
        SG_ALLOC(*canong_arg, g_arg->nv, g_arg->nde, "traces canong");
        updatecan_tr(g_arg, canong_arg, lab, CurrPart->inv, 0);
    }
    
    if (tv->options->generators) {
        deleteunmarked(&gensB);
        *tv->options->generators = gensB;
        freeschreier(&gpB, NULL);
    }
    else {
        freeschreier(&gpB, &gensB);
        schreier_freedyn();
    }
    
    while (STStart) {
        STAux = STStart;
        free(STAux->triearray);
        STStart = STStart->next;
        free(STAux);
    }
    
    tv->canlist = 0;
    for (i=0; i<=tv->treedepth; i++) {
        if (Spine[i].liststart) {
            tv->canlist += FreeList(Spine[i].liststart, TRUE);
            Spine[i].liststart = Spine[i].listend = NULL;
        }
    }
    
    if (GarbList) {
        tv->stats->peaknodes = FreeList(GarbList, FALSE);
    }
    
    FREECAND(NextCand)
    FREECAND(SpOrd)
    FREECAND(SpCyc)
    FREECAND(SpSwp)
    FREEPART(NextPart)
    FREEPART(SpPart1)
    FREEPART(SpPart2)
    
    if (!tv->options->getcanon && trieroot) {
        for (i=0; i<=tv->triepos; i++) {
            free(TrieArray[i]);
        }
    }
    
    for (i=0; i <= tv->treedepth; i++) {
        FREEPART(Spine[i].part)
    }
    
    CurrCand = GarbList = NULL;
    tv->stats->peaknodes += tv->canlist;
    
    free(tv);
    free(ti);
    traces_freedyn();
    return;
}

int traces_refine(sparsegraph *sg,
				  Candidate *Cand,
				  int m,
				  int n,
				  Partition *Part,
				  struct TracesVars* tv,
				  struct TracesInfo *ti,
                  boolean make_code)
{
	int i, j, k, sc, ind0, ind1, ind2, ind3, ind4, jk, tlp1, labi;
	int value, iend, newcell, TcSize;
	int HitClsInd, SplInd, SplCntInd, CStackInd, TraceInd, TraceCCInd, TraceStepsInd;
	size_t j1, iend1;
	unsigned int longcode;
	int newtrace = FALSE;
	int Sparse = TRUE;
	int *lab, *cls, *InvLab, *TracePos, *SplitCell, *LabCell, *TraceEnd, Traceccend, *Tracestpend;
	int BigCell, BigCellPos, BigCellSize;
	boolean TraceCell = FALSE;
    int *nghb;
	
	TcSize = Spine[tv->tolevel].tgtsize;
	
	if (tv->stackmark > (NAUTY_INFINITY-2)) {
		memset(StackMarkers, 0, n*sizeof(int));
		tv->stackmark = 0;
	}
	tv->stackmark++;
	
	SpineTL = Spine+tv->tolevel;
	TraceEnd = &(SpineTL->trcend);
	Traceccend = SpineTL->ccend;
	Tracestpend = &(SpineTL->stpend);
	TraceCCInd = SpineTL->ccstart;
	TraceStepsInd = SpineTL->stpstart;
	
	lab = Cand->lab;
	InvLab = Cand->invlab;
	cls = Part->cls;
	
    UPDATEMIN(Part->active, n-1);
	memcpy(CStack+1, TheTrace+SpineTL->trcstart, (Part->active)*sizeof(int));
	CStackInd = Part->active;
	for (i = 1; i <= CStackInd; i++) {
		StackMarkers[CStack[i]] = tv->stackmark;
	}
	
	longcode = Part->cells;
    TraceInd = SpineTL->trcstart+Part->active;
	
	if (!SpineTL->thetracexists) {
		newtrace = TRUE;
	}
    
	while (CStackInd > 0) {
        
		if (tv->mark > (NAUTY_INFINITY-2)) {
			memset(Markers, 0, n*sizeof(int));
			memset(MarkHitVtx, 0, n*sizeof(int));
			tv->mark = 0;
		}
		tv->mark++;
		
		if (Part->cells == tv->finalnumcells) break;
		
		j = CStackInd;
		k = CStackInd;
		while (--j > 0) {
			if (cls[CStack[j]] < cls[CStack[k]]) {
				k = j;
			}
			if ((cls[CStack[k]] == 1) || (j < CStackInd - 12)) {
				break;
			}
		}
        
		ind0 = CStack[k];
		ind2 = ind0+cls[ind0];
		CStack[k] = CStack[CStackInd--];
		StackMarkers[ind0] = 0;
		
		if (!newtrace) {
			TraceCell = ((ind0 == TheTraceCC[TraceCCInd]) && (TraceCCInd < Traceccend));
        }
		/* Analysis of occurrences of neighbors of the current cell */
		/* The list of cells with neighbors in the current cell is  built */
		if (cls[ind0] == 1) {			/* SINGLETON CURRENT CELL CASE */
			HitClsInd = 0;
			j1 = sg->v[lab[ind0]];
			iend1 = j1+sg->d[lab[ind0]];
			
			for (; j1 < iend1; ++j1) {
				k = sg->e[j1];
				value = Part->inv[InvLab[k]];
				if (cls[value] > 1) {
					if (Markers[value] != tv->mark) {
						HitCls[HitClsInd++] = value;
						Markers[value] = tv->mark;
						ElmHitCll[value] = value;
					}
					HitVtx[ElmHitCll[value]++] = k;
				}
			}
			tv->mark++;
            
			SplInd = 0;
			for (j = 0; j < HitClsInd; j++) {
				ind1 = HitCls[j];
				ElmHitCll[ind1] -= ind1;
				if ((ElmHitCll[ind1] > 0) && (ElmHitCll[ind1] < cls[ind1])) {
					SplCls[SplInd++] = ind1;
				}
			}
			
            /* SINGLETON CC CASE */
			if (SplInd) {
				if (newtrace) {
                    TheTraceCC[TraceCCInd] = ind0;
				}
                if (!TraceCell) {
                    TheTraceCC[TraceCCInd] = ind0;
                    newtrace = TRUE;
                }
                
                TRACE_CHECK(TheTraceSplNum, TraceCCInd, SplInd, &Traceccend)
                
				switch (SplInd) {
					case 0:
					case 1:
						break;
					case 2:
						if (SplCls[0] > SplCls[1]) {
							value = SplCls[0];
							SplCls[0] = SplCls[1];
							SplCls[1] = value;
						}
						break;
					case 3:
					case 4:
					case 5:
					case 6:
					case 7:
					case 8:
						for (k = 1; k < SplInd; ++k) {
							value = SplCls[k];
							i = k - 1;
							while ((i >= 0) && (value < SplCls[i])) {
								SplCls[i + 1] = SplCls[i];
								--i;
							}
							SplCls[i + 1] = value;
						}
						break;
					default:
						quickSort(SplCls, SplInd);
						break;
				}
				
				for (j = 0; j < SplInd; j++) {
					ind1 = SplCls[j];
					i = ind1+cls[ind1]-ElmHitCll[ind1];
                    TRACE_CHECK(TheTrace, TraceInd, i, TraceEnd)
				}
				
				/* REARRANGE THE CELLS */
				for (j = 0; j < SplInd; j++) {
					ind1 = SplCls[j];
					cls[ind1] = cls[ind1]-ElmHitCll[ind1];
					newcell = ind1+cls[ind1];
					cls[newcell] = ElmHitCll[ind1];
					Part->cells++;
					
					if (StackMarkers[ind1] != tv->stackmark) {
						if (cls[newcell] < cls[ind1]) {
							CStack[++CStackInd] = newcell;
							StackMarkers[newcell] = tv->stackmark;
						}
						else {
							CStack[++CStackInd] = ind1;
							StackMarkers[ind1] = tv->stackmark;
						}
					}
					else {
						CStack[++CStackInd] = newcell;
						StackMarkers[newcell] = tv->stackmark;
					}
					
					SplitCell = HitVtx+ind1;
					ind3 = cls[newcell];
					LabCell = lab+newcell;
					
					for (jk = 0; jk < ind3; jk++) {
						k = SplitCell[jk];
						i = LabCell[jk];
						Part->inv[newcell+jk] = newcell;
						lab[InvLab[k]] = i;
						InvLab[i] = InvLab[k];
						LabCell[jk] = k;
						InvLab[k] = newcell+jk;
					}
					if (cls[ind1] == 1) {
						Cand->singcode = MASHCOMM(Cand->singcode, Cand->lab[ind1]);
                    }
					if (cls[newcell] == 1) {
						Cand->singcode = MASHCOMM(Cand->singcode, Cand->lab[newcell]);
                    }
				}
			}
			else {
                if ((!newtrace) && TraceCell) {
                    return 0;
				}
			}
			
		}
		else {
			if (ti->thegraphisparse) {
				if (cls[ind0] == n) {
					memcpy(NghCounts, sg->d, n*sizeof(int));
					HitCls[0] = 0;
					HitClsInd = 1;
					ElmHitCll[0] = 0;
					for (i = 0; i < n; i++) {
						if (sg->d[i]) {
							HitVtx[ElmHitCll[0]++] = i;
						}
					}
				}
				else {
					HitClsInd = 0;
                    for (i = ind0; i < ind2; i++) {
                        labi = lab[i];
                        nghb = sg->e+sg->v[labi];
                        j = sg->d[labi];
                        for (; j--;) {
                            k = nghb[j];
                            if (MarkHitVtx[k] == tv->mark) {
                                NghCounts[k]++;
                            }
                            else {
                                value = Part->inv[InvLab[k]];
                                if (cls[value] > 1) {
                                    MarkHitVtx[k] = tv->mark;
                                    NghCounts[k] = 1;
                                    if (Markers[value] != tv->mark) {
                                        HitCls[HitClsInd++] = value;
                                        Markers[value] = tv->mark;
                                        HitVtx[value] = k;
                                        ElmHitCll[value] = 1;
                                    }
                                    else {
                                        HitVtx[value+ElmHitCll[value]++] = k;
                                    }
                                }
                            }
						}
					}
                    
                }
				
				tv->mark++;
				
                
				SplInd = 0;
				SplCls[0] = n;
				for (j = 0; j < HitClsInd; j++) {
					ind1 = HitCls[j];
					if ((ElmHitCll[ind1] > 0) && (ElmHitCll[ind1] < cls[ind1])) {
						SplCls[SplInd++] = HitCls[j];
					}
					else {
						ind2 = ind1+cls[ind1];
						value = NghCounts[lab[ind1++]];
						for (i = ind1; i < ind2; i++)
						{
							if (NghCounts[lab[i]] != value)
							{
								SplCls[SplInd++] = HitCls[j];
								break;
							}
						}
					}
				}
                
                /* SPARSE CASE */
				if (SplInd) {
					if (newtrace) {
                        TheTraceCC[TraceCCInd] = ind0;
					}
                    if (!TraceCell) {
                        TheTraceCC[TraceCCInd] = ind0;
                        newtrace = TRUE;
                    }
                    TRACE_CHECK(TheTraceSplNum, TraceCCInd, SplInd+n, &Traceccend)
					
					/* Sorting the cells to be split */
					switch (SplInd) {
						case 0:
						case 1:
							break;
						case 2:
							if (SplCls[0] > SplCls[1]) {
								value = SplCls[0];
								SplCls[0] = SplCls[1];
								SplCls[1] = value;
							}
							break;
						case 3:
						case 4:
						case 5:
						case 6:
						case 7:
						case 8:
							for (k = 1; k < SplInd; ++k) {
								value = SplCls[k];
								i = k - 1;
								while ((i >= 0) && (value < SplCls[i])) {
									SplCls[i + 1] = SplCls[i];
									--i;
								}
								SplCls[i + 1] = value;
							}
							break;
						default:
							quickSort(SplCls, SplInd);
							break;
					}
					
					for (sc = 0; sc < SplInd; sc++) {	/* For each cell C to be split */
						ind0 = SplCls[sc];
						ind1 = ind0 + cls[ind0];
						SplCntInd = 0;
						if (ElmHitCll[ind0] < cls[ind0]) {
							SplCnt[SplCntInd++] = 0;
							SplPos[0] = cls[ind0] - ElmHitCll[ind0];
						}
						
						/* According to the numbers of neighbors of C into the current cell */
						/* compute how many vertices in C will be placed into the same new cell */
						iend = ind0 + ElmHitCll[ind0];
						for (i = ind0; i < iend; i++) {
							value = NghCounts[HitVtx[i]];
							if (Markers[value] != tv->mark) {
								Markers[value] = tv->mark;
								SplCnt[SplCntInd++] = value;
								SplPos[value] = 1;
							}
							else {
								SplPos[value]++;
							}
						}
						tv->mark++;
						
						if (SplCntInd) {
                            TRACE_CHECK(TheTraceSteps, TraceStepsInd, SplCntInd+n, Tracestpend)
						}
						
						/* Sort the values deriving from the previous step */
						switch (SplCntInd) {
							case 0:
							case 1:
								break;
							case 2:
								if (SplCnt[0] > SplCnt[1]) {
									value = SplCnt[0];
									SplCnt[0] = SplCnt[1];
									SplCnt[1] = value;
								}
								break;
							case 3:
							case 4:
							case 5:
							case 6:
							case 7:
							case 8:
								for (k = 1; k < SplCntInd; ++k) {
									value = SplCnt[k];
									i = k - 1;
									while ((i >= 0) && (value < SplCnt[i])) {
										SplCnt[i + 1] = SplCnt[i];
										--i;
									}
									SplCnt[i + 1] = value;
								}
								break;
							default:
								quickSort(SplCnt, SplCntInd);
								break;
						}
						
						Part->cells += SplCntInd-1;
						
						/* Split the cell C and update the information for sizes of new cells */
						/* Put the new cells into the stack */
						i = ind0;
						if (StackMarkers[i] != tv->stackmark) {
							BigCellSize = 0;
						}
						for (k = 0; k < SplCntInd; k++) {
							value = SplPos[SplCnt[k]];
							cls[i] = value;
							if ((StackMarkers[ind0] != tv->stackmark) && (value > BigCellSize)) {
								BigCell = i;
								BigCellPos = CStackInd;
								BigCellSize = cls[i];
							}
							SplPos[SplCnt[k]] = i;
							i += value;
							if (i < ind1) {
								CStack[++CStackInd] = i;
								StackMarkers[i] = tv->stackmark;
                                TRACE_CHECK(TheTrace, TraceInd, i, TraceEnd)
							}
						}
						
						if ((StackMarkers[ind0] != tv->stackmark) && (BigCell != ind0)) {
							CStack[BigCellPos] = ind0;
							StackMarkers[BigCell] = 0;
							StackMarkers[ind0] = tv->stackmark;
						}
						/* Permute elements of the cell C */
						iend = ind0 + ElmHitCll[ind0];
						for (i = ind0; i < iend; i++) {
							value = HitVtx[i];
							j = SplPos[NghCounts[value]]++; /* where HitVtx[i] goes */
							k = InvLab[value];				/* where HitVtx[i] is in lab */
							lab[k] = lab[j];
							lab[j] = value;
							InvLab[value] = j;
							InvLab[lab[k]] = k;
							NghCounts[value] = 0;
						}
						
						/* Reconstruct the cell C and update the inverse partition */
						newcell = ind1 - ElmHitCll[ind0];
						i = newcell;
						ind2 = newcell+cls[newcell]-1;
						do {
							Part->inv[i] = newcell;
							if (i == ind2) {
								newcell = i+1;
								if (newcell < n) ind2 = newcell+cls[newcell]-1;
							}
						}
						while (++i < ind1);
						
						for (i = ind0, k = 0; k < SplCntInd; i+=cls[i], k++) {
							if (cls[i] == 1) {
								Cand->singcode = MASHCOMM(Cand->singcode, Cand->lab[i]);
                            }
						}
						
					}
				}
				else {
					if ((!newtrace) && TraceCell) {
                        return 0;
					}
				}
				
			}
			else {
				if (sg->d[lab[ind0]] > n/cls[ind0]) {
					Sparse = FALSE;
				}
				else {
					Sparse = TRUE;
				}
				if (Sparse) {
					/* Counting occurrences of neighbors of the current cell */
					/* The list of cells with neighbors in the current cell is also built */
					if (cls[ind0] == n) {
						memcpy(NghCounts, sg->d, n*sizeof(int));
						HitCls[0] = 0;
						HitClsInd = 1;
					}
					else {
						memset(NghCounts, 0, n*sizeof(int));
						HitClsInd = 0;
						for (i = ind0; i < ind2; i++) {
							j1 = sg->v[lab[i]];
							iend1 = j1+sg->d[lab[i]];
							for (; j1 < iend1; ++j1) {
								k = sg->e[j1];
								(NghCounts[k])++;
								value = Part->inv[InvLab[k]];
								if (Markers[value] != tv->mark) {
									if (cls[value] > 1) HitCls[HitClsInd++] = value;
									Markers[value] = tv->mark;
								}
							}
						}
					}
					
					tv->mark++;
					
                    
					SplInd = 0;
					for (j = 0; j < HitClsInd; j++) {
						ind1 = HitCls[j];
						ind2 = ind1+cls[ind1];
						value = NghCounts[lab[ind1++]];
						for (i = ind1; i < ind2; i++)
						{
							if (NghCounts[lab[i]] != value)
							{
								SplCls[SplInd++] = HitCls[j];
								break;
							}
						}
					}
					
                    /* DENSE-SPARSE CASE */
					if (SplInd) {
						if (newtrace) {
                            TheTraceCC[TraceCCInd] = ind0;
						}
                        if (!TraceCell) {
                            TheTraceCC[TraceCCInd] = ind0;
                            newtrace = TRUE;
                        }
                        TRACE_CHECK(TheTraceSplNum, TraceCCInd, SplInd+2*n, &Traceccend)
						
						/* Sorting the cells to be split */
						switch (SplInd) {
							case 0:
							case 1:
								break;
							case 2:
								if (SplCls[0] > SplCls[1]) {
									value = SplCls[0];
									SplCls[0] = SplCls[1];
									SplCls[1] = value;
								}
								break;
							case 3:
							case 4:
							case 5:
							case 6:
							case 7:
							case 8:
								for (k = 1; k < SplInd; ++k) {
									value = SplCls[k];
									i = k - 1;
									while ((i >= 0) && (value < SplCls[i])) {
										SplCls[i + 1] = SplCls[i];
										--i;
									}
									SplCls[i + 1] = value;
								}
								break;
							default:
								quickSort(SplCls, SplInd);
								break;
						}
						
						for (j = 0; j < SplInd; j++) {	/* For each cell C to be split */
							ind0 = SplCls[j];
							ind1 = ind0+cls[ind0];
							SplCntInd = 0;
							
							/* According to the numbers of neighbors of C into the current cell */
							/* compute how many vertices in C will be placed into the same new cell */
							for (i = ind0; i < ind1; i++) {
								value = NghCounts[lab[i]];
								if (Markers[value] != tv->mark) {
									Markers[value] = tv->mark;
									SplCnt[SplCntInd++] = value;
									SplPos[value] = 1;
								}
								else {
									SplPos[value]++;
								}
							}
							tv->mark++;
							
							if (SplCntInd) {
                                TRACE_CHECK(TheTraceSteps, TraceStepsInd, SplCntInd+2*n, Tracestpend)
							}
							
							/* Sort the values deriving from the previous step */
							switch (SplCntInd) {
								case 0:
								case 1:
									break;
								case 2:
									if (SplCnt[0] > SplCnt[1]) {
										value = SplCnt[0];
										SplCnt[0] = SplCnt[1];
										SplCnt[1] = value;
									}
									break;
								case 3:
								case 4:
								case 5:
								case 6:
								case 7:
								case 8:
									for (k = 1; k < SplCntInd; ++k) {
										value = SplCnt[k];
										i = k - 1;
										while ((i >= 0) && (value < SplCnt[i])) {
											SplCnt[i + 1] = SplCnt[i];
											--i;
										}
										SplCnt[i + 1] = value;
									}
									break;
								default:
									quickSort(SplCnt, SplCntInd);
									break;
							}
							
							Part->cells += SplCntInd-1;
							
							/* Split the cell C and update the information for sizes of new cells */
							/* Put the new cells into the stack */
							i = ind0;
							if (StackMarkers[i] != tv->stackmark) {
								BigCellSize = 0;
							}
							for (k = 0; k < SplCntInd; k++) {
								value = SplPos[SplCnt[k]];
								cls[i] = value;
								
								if ((StackMarkers[ind0] != tv->stackmark) && (value > BigCellSize)) {
									BigCell = i;
									BigCellPos = CStackInd;
									BigCellSize = cls[i];
								}
								SplPos[SplCnt[k]] = i;
								i += value;
								if (i < ind1) {
									CStack[++CStackInd] = i;
									StackMarkers[i] = tv->stackmark;
                                    TRACE_CHECK(TheTrace, TraceInd, i, TraceEnd)
								}
							}
							
							if ((StackMarkers[ind0] != tv->stackmark) && (BigCell != ind0)) {
								CStack[BigCellPos] = ind0;
								StackMarkers[BigCell] = 0;
								StackMarkers[ind0] = tv->stackmark;
							}
							
							/* Permute elements of the cell C */
							i = ind0;
							do {
								SplCnt[SplPos[NghCounts[lab[i]]]++] = lab[i];
							}
							while(++i < ind1);
							
							/* Reconstruct the cell C and update the inverse partition */
							newcell = ind0;
							i = ind0;
							ind2 = newcell+cls[newcell]-1;
							do {
								lab[i] = SplCnt[i];
								InvLab[lab[i]] = i;
								
								Part->inv[i] = newcell;
								if (i == ind2) {
									newcell = i+1;
									if (newcell < n) ind2 = newcell+cls[newcell]-1;
								}
							}
							while (++i < ind1);
							
							for (i = ind0, k = 0; k < SplCntInd; i+=cls[i], k++) {
								if (cls[i] == 1) {
									Cand->singcode = MASHCOMM(Cand->singcode, Cand->lab[i]);
                                }
							}
							
						}
					}
					else {
						if ((!newtrace) && TraceCell) {
                            return 0;
						}
					}
					
				}
				else {
					if (cls[ind0] == n) {
						memcpy(NghCounts, sg->d, n*sizeof(int));
					}
					else {
						memset(NghCounts, 0, n*sizeof(int));
						
						for (i = ind0; i < ind2; i++) {
							j1 = sg->v[lab[i]];
							iend1 = j1+sg->d[lab[i]];
							for (; j1 < iend1; ++j1) {
								k = sg->e[j1];
								(NghCounts[k])++;
							}
						}
					}
					SplInd = 0;
					ind4 = 0;
                    while (ind4 < n) {	/* For each cell C with size(C) > 1 */
						ind1 = ind4+cls[ind4];
						if (cls[ind4] > 1) {
							
                            
							/* Determine whether C must be split */
							SplCntInd = 0;
							value = NghCounts[lab[ind4]];
							for (i = ind4+1; i < ind1; i++) {
								if (NghCounts[lab[i]] != value)
								{
									SplInd++;
									Markers[value] = tv->mark;
									SplCnt[SplCntInd++] = value;
									SplPos[value] = i-ind4;
									do {
										value = NghCounts[lab[i++]];
										if (Markers[value] != tv->mark) {
											Markers[value] = tv->mark;
											SplCnt[SplCntInd++] = value;
											SplPos[value] = 1;
										}
										else {
											SplPos[value]++;
										}
									}
									while(i != ind1);
									break;
								}
							}
							tv->mark++;
							
                            if (SplInd && !TraceCell) newtrace = TRUE;
                            
							if (SplCntInd) {
                                TRACE_CHECK(TheTraceSteps, TraceStepsInd, SplCntInd+3*n, Tracestpend)
								
								/* Sort the values deriving from the previous step */
								switch (SplCntInd) {
									case 0:
									case 1:
										break;
									case 2:
										if (SplCnt[0] > SplCnt[1]) {
											value = SplCnt[0];
											SplCnt[0] = SplCnt[1];
											SplCnt[1] = value;
										}
										break;
									case 3:
									case 4:
									case 5:
									case 6:
									case 7:
									case 8:
										for (k = 1; k < SplCntInd; ++k) {
											value = SplCnt[k];
											i = k - 1;
											while ((i >= 0) && (value < SplCnt[i])) {
												SplCnt[i + 1] = SplCnt[i];
												--i;
											}
											SplCnt[i + 1] = value;
										}
										break;
									default:
										quickSort(SplCnt, SplCntInd);
										break;
								}
								
								Part->cells += SplCntInd-1;
								
								/* Split the cell C and update the information for sizes of new cells */
								/* Put the new cells into the stack */
								i = ind4;
								if (StackMarkers[i] != tv->stackmark) {
									BigCellSize = 0;
								}
								for (k = 0; k < SplCntInd; k++) {
									value = SplPos[SplCnt[k]];
									cls[i] = value;
									
									if ((StackMarkers[ind4] != tv->stackmark) && (value > BigCellSize)) {
										BigCell = i;
										BigCellPos = CStackInd;
										BigCellSize = cls[i];
									}
									SplPos[SplCnt[k]] = i;
									i += value;
									if (i < ind1) {
										CStack[++CStackInd] = i;
										StackMarkers[i] = tv->stackmark;
                                        TRACE_CHECK(TheTrace, TraceInd, i, TraceEnd)
									}
								}
								if ((StackMarkers[ind4] != tv->stackmark) && (BigCell != ind4)) {
									CStack[BigCellPos] = ind4;
									StackMarkers[BigCell] = 0;
									StackMarkers[ind4] = tv->stackmark;
								}
								
								/* Permute elements of the cell C */
								i = ind4;
								do {
									SplCnt[SplPos[NghCounts[lab[i]]]++] = lab[i];
								}
								while(++i < ind1);
								
								/* Reconstruct the cell C and update the inverse partition */
								newcell = ind4;
								i = ind4;
								ind2 = newcell+cls[newcell]-1;
								do {
									lab[i] = SplCnt[i];
									InvLab[lab[i]] = i;
									Part->inv[i] = newcell;
									if (i == ind2) {
										newcell = i+1;
										if (newcell < n) ind2 = newcell+cls[newcell]-1;
									}
								}
								while (++i < ind1);
								
								for (i = ind4; i < ind1; i+=cls[i]) {
									if (cls[i] == 1) {
										Cand->singcode = MASHCOMM(Cand->singcode, Cand->lab[i]);
                                    }
								}
							}
						}
						ind4 = ind1;
					}
                    
                    /* DENSE-DENSE CASE */
					if (SplInd) {
                        if (!TraceCell) {
                            newtrace = TRUE;
                        }
                        TRACE_CHECK(TheTraceSplNum, TraceCCInd, SplInd+3*n, &Traceccend)
						if (newtrace) {
                            TheTraceCC[TraceCCInd-1] = ind0;
						}
					}
					else {
						if ((!newtrace) && TraceCell) {
                            return 0;
						}
					}
					
				}
			}
		}
	}  /* end while (CStackInd > 0) */
    if (make_code) {
        for (i=SpineTL->trcstart; i < TraceInd; i++) {
            ind0 = TheTrace[i];
            longcode = MASHNONCOMM(longcode, Part->inv[ind0]);
            j1 = sg->v[lab[ind0]];
            iend1 = j1+sg->d[lab[ind0]];
            for (; j1 < iend1; ++j1) {
                k = sg->e[j1];
                value = Part->inv[InvLab[k]];
                longcode = MASHCOMM(longcode, value);
            }
        }
    }
    
	Part->code = Cand->code = CLEANUP(longcode);
    tlp1 = tv->tolevel+1;
	if (newtrace) {
        if ((tlp1 < n) && (tlp1 > tv->treedepth)) {
			tv->treedepth = tlp1;
            if (tv->strategy) {
                NEWPART(Spine[tlp1].part)
			}
			else {
                NEWPARTSPINE(tlp1)
			}
			Spine[tlp1].liststart = Spine[tlp1].listend = NULL;
			Spine[tlp1].listcounter = 0;
		}
		*TraceEnd = TraceInd;
		SpineTL->ccend = TraceCCInd;
		*Tracestpend = TraceStepsInd;
		
		SpineTL->thetracexists = TRUE;
		if (tlp1 < n) {
			Spine[tlp1].ccstart = TraceCCInd;
			Spine[tlp1].stpstart = TraceStepsInd;
			Spine[tlp1].thetracexists = FALSE;
			Spine[tlp1].part->code = -1;
		}
		return 2;
	}
	else {
        if (TraceInd < *TraceEnd) {
            return 0;
        }
		if (Cand->code > SpineTL->part->code) {
			return 2;
		}
		else {
			if (Cand->code < SpineTL->part->code) {
				return 0;
			}
		}
		return 1;
	}
}

void traces_refine_notrace(sparsegraph *sg,
						   Candidate *Cand,
						   int m,
						   int n,
						   Partition *Part,
						   struct TracesVars* tv,
						   struct TracesInfo *ti)
{
	int i, i1, j, k, sc, ind0, ind1, ind2, ind3, labi;
	int value, iend, newcell;
	int HitClsInd, SplInd, SplCntInd, CStackInd;
	size_t j1, iend1;
	unsigned int longcode;
	int Split = 0;
	int Sparse = TRUE;
	int *lab, *cls, *InvLab, *SplitCell, *LabCell;
	int BigCell, BigCellPos, BigCellSize;
    int *nghb;
    
 	if (tv->stackmark > (NAUTY_INFINITY-2)) {
		memset(StackMarkers, 0, n*sizeof(int));
		tv->stackmark = 0;
	}
	tv->stackmark++;
	
	lab = Cand->lab;
	InvLab = Cand->invlab;
	cls = Part->cls;
	
	CStackInd = 1;
    CStack[1] = tv->tcellexpath+cls[tv->tcellexpath];
	
	for (i = 1; i <= CStackInd; i++) {
		StackMarkers[CStack[i]] = tv->stackmark;
	}
	
	longcode = Part->cells;
	
	while (CStackInd > 0)
	{
		if (tv->mark > (NAUTY_INFINITY-2)) {
			memset(Markers, 0, n*sizeof(int));
			memset(MarkHitVtx, 0, n*sizeof(int));
			tv->mark = 0;
		}
		tv->mark++;
		
        
		j = CStackInd;
		k = CStackInd;
		while (--j > 0) {
			if (cls[CStack[j]] < cls[CStack[k]]) {
				k = j;
			}
			if ((cls[CStack[k]] == 1) || (j < CStackInd - 12)) {
				break;
			}
		}
		
		ind0 = CStack[k];
		ind2 = ind0+cls[ind0];
		CStack[k] = CStack[CStackInd--];		/* Current Cell */
		longcode = MASHNONCOMM(longcode, ind0);
		StackMarkers[ind0] = 0;
		
		/* Analysis of occurrences of neighbors of the current cell */
		/* The list of cells with neighbors in the current cell is  built */
		if (cls[ind0] == 1) {	/* SINGLETON CURRENT CELL CASE */
			HitClsInd = 0;
			j1 = sg->v[lab[ind0]];
			iend1 = j1+sg->d[lab[ind0]];
			
			for (; j1 < iend1; ++j1) {
				k = sg->e[j1];
				value = Part->inv[InvLab[k]];
				if (cls[value] > 1) {
					if (Markers[value] != tv->mark) {
						HitCls[HitClsInd++] = value;
						Markers[value] = tv->mark;
						ElmHitCll[value] = value;
					}
					HitVtx[ElmHitCll[value]++] = k;
				}
				else {
					longcode = MASHCOMM(longcode, value);
				}
			}
			
			tv->mark++;
			SplInd = 0;
			for (j = 0; j < HitClsInd; j++) {
				ind1 = HitCls[j];
				ElmHitCll[ind1] -= ind1;
				if ((ElmHitCll[ind1] > 0) && (ElmHitCll[ind1] < cls[ind1])) {
					SplCls[SplInd++] = ind1;
				}
			}
			
			switch (SplInd) {
				case 0:
				case 1:
					break;
				case 2:
					if (SplCls[0] > SplCls[1]) {
						value = SplCls[0];
						SplCls[0] = SplCls[1];
						SplCls[1] = value;
					}
					break;
				case 3:
				case 4:
				case 5:
				case 6:
				case 7:
				case 8:
					for (k = 1; k < SplInd; ++k) {
						value = SplCls[k];
						i = k - 1;
						while ((i >= 0) && (value < SplCls[i])) {
							SplCls[i + 1] = SplCls[i];
							--i;
						}
						SplCls[i + 1] = value;
					}
					break;
				default:
					quickSort(SplCls, SplInd);
					break;
			}
			
			/* REARRANGE THE CELL */
			for (j = 0; j < SplInd; j++) {
				ind1 = SplCls[j];
				cls[ind1] = cls[ind1]-ElmHitCll[ind1];
				newcell = ind1+cls[ind1];
				cls[newcell] = ElmHitCll[ind1];
				
				Part->cells++;
				
				if (StackMarkers[ind1] != tv->stackmark) {
					if (cls[newcell] < cls[ind1]) {
						CStack[++CStackInd] = newcell;
						StackMarkers[newcell] = tv->stackmark;
					}
					else {
						CStack[++CStackInd] = ind1;
						StackMarkers[ind1] = tv->stackmark;
					}
				}
				else {
					CStack[++CStackInd] = newcell;
					StackMarkers[newcell] = tv->stackmark;
				}
				
				SplitCell = HitVtx+ind1;
				ind3 = cls[newcell];
				LabCell = lab+newcell;
				
				for (i1 = 0; i1 < ind3; i1++) {
					k = SplitCell[i1];
					i = LabCell[i1];
					Part->inv[newcell+i1] = newcell;
					lab[InvLab[k]] = i;
					InvLab[i] = InvLab[k];
					LabCell[i1] = k;
					InvLab[k] = newcell+i1;
				}
                if (cls[ind1] == 1) {
                    Cand->pathsingcode = MASHCOMM(Cand->pathsingcode, lab[ind1]);
                }
                if (cls[newcell] == 1) {
                    Cand->pathsingcode = MASHCOMM(Cand->pathsingcode, lab[newcell]);
                }
			}
		}
		else {
			if (ti->thegraphisparse) {
				if (cls[ind0] == n) {
					memcpy(NghCounts, sg->d, n*sizeof(int));
					HitCls[0] = 0;
					HitClsInd = 1;
					ElmHitCll[0] = 0;
					for (i = 0; i < n; i++) {
						if (sg->d[i]) {
							HitVtx[ElmHitCll[0]++] = i;
						}
					}
				}
				else {
					HitClsInd = 0;
                    for (i = ind0; i < ind2; i++) {
                        labi = lab[i];
                        nghb = sg->e+sg->v[labi];
                        j = sg->d[labi];
                        for (; j--;) {
                            k = nghb[j];
                            if (MarkHitVtx[k] == tv->mark) {
                                NghCounts[k]++;
                            }
                            else {
                                value = Part->inv[InvLab[k]];
                                if (cls[value] > 1) {
                                    MarkHitVtx[k] = tv->mark;
                                    NghCounts[k] = 1;
                                    if (Markers[value] != tv->mark) {
                                        HitCls[HitClsInd++] = value;
                                        Markers[value] = tv->mark;
                                        HitVtx[value] = k;
                                        ElmHitCll[value] = 1;
                                    }
                                    else {
                                        HitVtx[value+ElmHitCll[value]++] = k;
                                    }
                                }
                                else {
                                    longcode = MASHCOMM(longcode, value);
                                }
                            }
						}
					}
                    
                }
				
				tv->mark++;
				
				SplInd = 0;
				SplCls[0] = n;
				for (j = 0; j < HitClsInd; j++) {
					ind1 = HitCls[j];
					if ((ElmHitCll[ind1] > 0) && (ElmHitCll[ind1] < cls[ind1])) {
						SplCls[SplInd++] = HitCls[j];
					}
					else {
						ind2 = ind1+cls[ind1];
						Split = FALSE;
						value = NghCounts[lab[ind1++]];
						for (i = ind1; i < ind2; i++)
						{
							if (NghCounts[lab[i]] != value)
							{
								SplCls[SplInd++] = HitCls[j];
								Split = TRUE;
								break;
							}
						}
						if (!Split) {
							longcode = MASHCOMM(longcode, ind1);
						}
					}
				}
				/* Sorting the cells to be split */
				switch (SplInd) {
					case 0:
					case 1:
						break;
					case 2:
						if (SplCls[0] > SplCls[1]) {
							value = SplCls[0];
							SplCls[0] = SplCls[1];
							SplCls[1] = value;
						}
						break;
					case 3:
					case 4:
					case 5:
					case 6:
					case 7:
					case 8:
						for (k = 1; k < SplInd; ++k) {
							value = SplCls[k];
							i = k - 1;
							while ((i >= 0) && (value < SplCls[i])) {
								SplCls[i + 1] = SplCls[i];
								--i;
							}
							SplCls[i + 1] = value;
						}
						break;
					default:
						quickSort(SplCls, SplInd);
						break;
				}
				
				for (sc = 0; sc < SplInd; sc++) {	/* For each cell C to be split */
					ind0 = SplCls[sc];
					ind1 = ind0 + cls[ind0];
					SplCntInd = 0;
					if (ElmHitCll[ind0] < cls[ind0]) {
						SplCnt[SplCntInd++] = 0;
						SplPos[0] = cls[ind0] - ElmHitCll[ind0];
					}
					
					/* According to the numbers of neighbors of C into the current cell */
					/* compute how many vertices in C will be placed into the same new cell */
					iend = ind0 + ElmHitCll[ind0];
					for (i = ind0; i < iend; i++) {
						value = NghCounts[HitVtx[i]];
						if (Markers[value] != tv->mark) {
							Markers[value] = tv->mark;
							SplCnt[SplCntInd++] = value;
							SplPos[value] = 1;
						}
						else {
							SplPos[value]++;
						}
					}
					tv->mark++;
					
					/* Sort the values deriving from the previous step */
					switch (SplCntInd) {
						case 0:
						case 1:
							break;
						case 2:
							if (SplCnt[0] > SplCnt[1]) {
								value = SplCnt[0];
								SplCnt[0] = SplCnt[1];
								SplCnt[1] = value;
							}
							break;
						case 3:
						case 4:
						case 5:
						case 6:
						case 7:
						case 8:
							for (k = 1; k < SplCntInd; ++k) {
								value = SplCnt[k];
								i = k - 1;
								while ((i >= 0) && (value < SplCnt[i])) {
									SplCnt[i + 1] = SplCnt[i];
									--i;
								}
								SplCnt[i + 1] = value;
							}
							break;
						default:
							quickSort(SplCnt, SplCntInd);
							break;
					}
					Part->cells += SplCntInd-1;
					
					/* Split the cell C and update the information for sizes of new cells */
					/* Put the new cells into the stack */
					i = ind0;
					if (StackMarkers[i] != tv->stackmark) {
						BigCellSize = 0;
					}
					for (k = 0; k < SplCntInd; k++) {
						value = SplPos[SplCnt[k]];
						cls[i] = value;
						if ((StackMarkers[ind0] != tv->stackmark) && (value > BigCellSize)) {
							BigCell = i;
							BigCellPos = CStackInd;
							BigCellSize = cls[i];
						}
						SplPos[SplCnt[k]] = i;
						i += value;
						if (i < ind1) {
							CStack[++CStackInd] = i;
							StackMarkers[i] = tv->stackmark;
						}
					}
					
					if ((StackMarkers[ind0] != tv->stackmark) && (BigCell != ind0)) {
						CStack[BigCellPos] = ind0;
						StackMarkers[BigCell] = 0;
						StackMarkers[ind0] = tv->stackmark;
					}
					/* Permute elements of the cell C */
					iend = ind0 + ElmHitCll[ind0];
					for (i = ind0; i < iend; i++) {
						value = HitVtx[i];
						j = SplPos[NghCounts[value]]++; /* where HitVtx[i] goes */
						k = InvLab[value];				/* where HitVtx[i] is in lab */
						lab[k] = lab[j];
						lab[j] = value;
						InvLab[value] = j;
						InvLab[lab[k]] = k;
						NghCounts[value] = 0;
					}
					
					/* Reconstruct the cell C and update the inverse partition */
					newcell = ind1 - ElmHitCll[ind0];
					i = newcell;
					ind2 = newcell+cls[newcell]-1;
					do {
						Part->inv[i] = newcell;
						if (i == ind2) {
							newcell = i+1;
							if (newcell < n) ind2 = newcell+cls[newcell]-1;
						}
					}
					while (++i < ind1);
                    
                    for (i = ind0, k = 0; k < SplCntInd; i+=cls[i], k++) {
                        if (cls[i] == 1) {
                            Cand->pathsingcode = MASHCOMM(Cand->pathsingcode, lab[i]);
                        }
                    }
                    
				}
			}
			else {
				if (sg->d[lab[ind0]] > n/cls[ind0]) {
					Sparse = FALSE;
				}
				else {
					Sparse = TRUE;
				}
				if (Sparse) {
					/* Counting occurrences of neighbors of the current cell */
					/* The list of cells with neighbors in the current cell is also built */
					if (cls[ind0] == n) {
						memcpy(NghCounts, sg->d, n*sizeof(int));
						HitCls[0] = 0;
						HitClsInd = 1;
					}
					else {
						memset(NghCounts, 0, n*sizeof(int));
						HitClsInd = 0;
						for (i = ind0; i < ind2; i++) {
							j1 = sg->v[lab[i]];
							iend1 = j1+sg->d[lab[i]];
							for (; j1 < iend1; ++j1) {
								k = sg->e[j1];
								(NghCounts[k])++;
								value = Part->inv[InvLab[k]];
								if (Markers[value] != tv->mark) {
									if (cls[value] > 1) HitCls[HitClsInd++] = value;
									Markers[value] = tv->mark;
								}
							}
						}
					}
					
					tv->mark++;
					
					SplInd = 0;
					for (j = 0; j < HitClsInd; j++) {
						ind1 = HitCls[j];
						ind2 = ind1+cls[ind1];
						value = NghCounts[lab[ind1++]];
						for (i = ind1; i < ind2; i++)
						{
							if (NghCounts[lab[i]] != value)
							{
								SplCls[SplInd++] = HitCls[j];
								break;
							}
						}
					}
					
					/* Sorting the cells to be split */
					switch (SplInd) {
						case 0:
						case 1:
							break;
						case 2:
							if (SplCls[0] > SplCls[1]) {
								value = SplCls[0];
								SplCls[0] = SplCls[1];
								SplCls[1] = value;
							}
							break;
						case 3:
						case 4:
						case 5:
						case 6:
						case 7:
						case 8:
							for (k = 1; k < SplInd; ++k) {
								value = SplCls[k];
								i = k - 1;
								while ((i >= 0) && (value < SplCls[i])) {
									SplCls[i + 1] = SplCls[i];
									--i;
								}
								SplCls[i + 1] = value;
							}
							break;
						default:
							quickSort(SplCls, SplInd);
							break;
					}
					
					for (j = 0; j < SplInd; j++) {	/* For each cell C to be split */
						ind0 = SplCls[j];
						ind1 = ind0+cls[ind0];
						SplCntInd = 0;
						
						/* According to the numbers of neighbors of C into the current cell */
						/* compute how many vertices in C will be placed into the same new cell */
						for (i = ind0; i < ind1; i++) {
							value = NghCounts[lab[i]];
							if (Markers[value] != tv->mark) {
								Markers[value] = tv->mark;
								SplCnt[SplCntInd++] = value;
								SplPos[value] = 1;
							}
							else {
								SplPos[value]++;
							}
						}
						
						tv->mark++;
						
						/* Sort the values deriving from the previous step */
						switch (SplCntInd) {
							case 0:
							case 1:
								break;
							case 2:
								if (SplCnt[0] > SplCnt[1]) {
									value = SplCnt[0];
									SplCnt[0] = SplCnt[1];
									SplCnt[1] = value;
								}
								break;
							case 3:
							case 4:
							case 5:
							case 6:
							case 7:
							case 8:
								for (k = 1; k < SplCntInd; ++k) {
									value = SplCnt[k];
									i = k - 1;
									while ((i >= 0) && (value < SplCnt[i])) {
										SplCnt[i + 1] = SplCnt[i];
										--i;
									}
									SplCnt[i + 1] = value;
								}
								break;
							default:
								quickSort(SplCnt, SplCntInd);
								break;
						}
						
						Part->cells += SplCntInd-1;
						
						/* Split the cell C and update the information for sizes of new cells */
						/* Put the new cells into the stack */
						i = ind0;
						if (StackMarkers[i] != tv->stackmark) {
							BigCellSize = 0;
						}
						for (k = 0; k < SplCntInd; k++) {
							value = SplPos[SplCnt[k]];
							cls[i] = value;
							if ((StackMarkers[ind0] != tv->stackmark) && (value > BigCellSize)) {
								BigCell = i;
								BigCellPos = CStackInd;
								BigCellSize = cls[i];
							}
							SplPos[SplCnt[k]] = i;
							i += value;
							if (i < ind1) {
								CStack[++CStackInd] = i;
								StackMarkers[i] = tv->stackmark;
							}
						}
						if ((StackMarkers[ind0] != tv->stackmark) && (BigCell != ind0)) {
							CStack[BigCellPos] = ind0;
							StackMarkers[BigCell] = 0;
							StackMarkers[ind0] = tv->stackmark;
						}
						
						/* Permute elements of the cell C */
						i = ind0;
						do {
							SplCnt[SplPos[NghCounts[lab[i]]]++] = lab[i];
						}
						while(++i < ind1);
						
						/* Reconstruct the cell C and update the inverse partition */
						newcell = ind0;
						i = ind0;
						ind2 = newcell+cls[newcell]-1;
						do {
							lab[i] = SplCnt[i];
							InvLab[lab[i]] = i;
							Part->inv[i] = newcell;
							
							if (i == ind2) {
								newcell = i+1;
								if (newcell < n) ind2 = newcell+cls[newcell]-1;
							}
						}
						while (++i < ind1);
                        
                        for (i = ind0, k = 0; k < SplCntInd; i+=cls[i], k++) {
                            if (cls[i] == 1) {
                                Cand->pathsingcode = MASHCOMM(Cand->pathsingcode, lab[i]);
                            }
                        }
                        
					}
				}
				else {
					if (cls[ind0] == n) {
						memcpy(NghCounts, sg->d, n*sizeof(int));
					}
					else {
						memset(NghCounts, 0, n*sizeof(int));
						
						for (i = ind0; i < ind2; i++) {
							j1 = sg->v[lab[i]];
							iend1 = j1+sg->d[lab[i]];
							for (; j1 < iend1; ++j1) {
								k = sg->e[j1];
								(NghCounts[k])++;
							}
						}
					}
					
					ind0 = 0;
					while (ind0 < n) {	/* For each cell C with size(C) > 1 */
						ind1 = ind0+cls[ind0];
						if (cls[ind0] > 1) {
							
							/* Determine whether C must be split */
							SplCntInd = 0;
							value = NghCounts[lab[ind0]];
							for (i = ind0+1; i < ind1; i++)
							{
								if (NghCounts[lab[i]] != value)
								{
									Markers[value] = tv->mark;
									SplCnt[SplCntInd++] = value;
									SplPos[value] = i-ind0;
									do {
										value = NghCounts[lab[i++]];
										if (Markers[value] != tv->mark) {
											Markers[value] = tv->mark;
											SplCnt[SplCntInd++] = value;
											SplPos[value] = 1;
										}
										else {
											SplPos[value]++;
										}
									}
									while(i != ind1);
									break;
								}
							}
							
							if (SplCntInd) {
								tv->mark++;
								
								/* Sort the values deriving from the previous step */
								switch (SplCntInd) {
									case 0:
									case 1:
										break;
									case 2:
										if (SplCnt[0] > SplCnt[1]) {
											value = SplCnt[0];
											SplCnt[0] = SplCnt[1];
											SplCnt[1] = value;
										}
										break;
									case 3:
									case 4:
									case 5:
									case 6:
									case 7:
									case 8:
										for (k = 1; k < SplCntInd; ++k) {
											value = SplCnt[k];
											i = k - 1;
											while ((i >= 0) && (value < SplCnt[i])) {
												SplCnt[i + 1] = SplCnt[i];
												--i;
											}
											SplCnt[i + 1] = value;
										}
										break;
									default:
										quickSort(SplCnt, SplCntInd);
										break;
								}
								
								Part->cells += SplCntInd-1;
								
								/* Split the cell C and update the information for sizes of new cells */
								/* Put the new cells into the stack */
								i = ind0;
								if (StackMarkers[i] != tv->stackmark) {
									BigCellSize = 0;
								}
								for (k = 0; k < SplCntInd; k++) {
									value = SplPos[SplCnt[k]];
									cls[i] = value;
									if ((StackMarkers[ind0] != tv->stackmark) && (value > BigCellSize)) {
										BigCell = i;
										BigCellPos = CStackInd;
										BigCellSize = cls[i];
									}
									SplPos[SplCnt[k]] = i;
									i += value;
									if (i < ind1) {
										CStack[++CStackInd] = i;
										StackMarkers[i] = tv->stackmark;
									}
								}
								if ((StackMarkers[ind0] != tv->stackmark) && (BigCell != ind0)) {
									CStack[BigCellPos] = ind0;
									StackMarkers[BigCell] = 0;
									StackMarkers[ind0] = tv->stackmark;
								}
								
								/* Permute elements of the cell C */
								i = ind0;
								do {
									SplCnt[SplPos[NghCounts[lab[i]]]++] = lab[i];
								}
								while(++i < ind1);
								
								/* Reconstruct the cell C and update the inverse partition */
								newcell = ind0;
								i = ind0;
								ind2 = newcell+cls[newcell]-1;
								do {
									lab[i] = SplCnt[i];
									InvLab[lab[i]] = i;
									Part->inv[i] = newcell;
									if (i == ind2) {
										newcell = i+1;
										if (newcell < n) ind2 = newcell+cls[newcell]-1;
									}
								}
								while (++i < ind1);
                                
                                for (i = ind0; i < ind1; i+=cls[i]) {
									if (cls[i] == 1) {
										Cand->pathsingcode = MASHCOMM(Cand->pathsingcode, lab[i]);
                                    }
								}
                                
							}
						}
						ind0 = ind1;
					}
				}
			}
		}
	}  /* end while (CStackInd > 0) */
	
	Cand->code = CLEANUP(longcode);
	return;
}

void traces_refine_maketrie(sparsegraph *sg,
							Candidate *Cand,
							int m,
							int n,
							Partition *Part,
							struct TracesVars* tv,
							struct TracesInfo *ti)
{
	int i, i1, j, k, sc, ind0, ind1, ind2, ind3, labi;
	int value, iend, newcell, TcSize;
	int HitClsInd, SplInd, SplCntInd, CStackInd;
	size_t j1, iend1;
	unsigned int longcode;
	int Split = 0;
	int Sparse = TRUE;
	int *lab, *cls, *InvLab, *SplitCell, *LabCell;
	int BigCell, BigCellPos, BigCellSize;
    int *nghb;
    
	TcSize = Spine[tv->tolevel].tgtsize;
	
	if (tv->stackmark > (NAUTY_INFINITY-2)) {
		memset(StackMarkers, 0, n*sizeof(int));
		tv->stackmark = 0;
	}
	tv->stackmark++;
	
	lab = Cand->lab;
	InvLab = Cand->invlab;
	cls = Part->cls;
	
	CStack[1] = Spine[tv->tolevel_tl].tgtpos;
	CStackInd = 1;
	for (i = 1; i <= CStackInd; i++) {
		StackMarkers[CStack[i]] = tv->stackmark;
	}
	
	longcode = Part->cells;
	
	while (CStackInd > 0)
	{
		if (tv->mark > (NAUTY_INFINITY-2)) {
			memset(Markers, 0, n*sizeof(int));
			memset(MarkHitVtx, 0, n*sizeof(int));
			tv->mark = 0;
		}
		tv->mark++;
		
		if (Part->cells == tv->finalnumcells) break;
		
		j = CStackInd;
		k = CStackInd;
		while (--j > 0) {
			if (cls[CStack[j]] < cls[CStack[k]]) {
				k = j;
			}
			if ((cls[CStack[k]] == 1) || (j < CStackInd - 12)) {
				break;
			}
		}
		
		ind0 = CStack[k];
		ind2 = ind0+cls[ind0];
		CStack[k] = CStack[CStackInd--];
		longcode = MASHNONCOMM(longcode, ind0);
		StackMarkers[ind0] = 0;
		
		/* Analysis of occurrences of neighbors of the current cell */
		/* The list of cells with neighbors in the current cell is  built */
		if (cls[ind0] == 1) {  /* SINGLETON CURRENT CELL CASE */
			HitClsInd = 0;
			j1 = sg->v[lab[ind0]];
			iend1 = j1+sg->d[lab[ind0]];
			
			for (; j1 < iend1; ++j1) {
				k = sg->e[j1];
				value = Part->inv[InvLab[k]];
				if (cls[value] > 1) {
					if (Markers[value] != tv->mark) {
						HitCls[HitClsInd++] = value;
						Markers[value] = tv->mark;
						ElmHitCll[value] = value;
					}
					HitVtx[ElmHitCll[value]++] = k;
				}
				else {
					longcode = MASHCOMM(longcode, value);
				}
			}
			
			tv->mark++;
			
			SplInd = 0;
			for (j = 0; j < HitClsInd; j++) {
				ind1 = HitCls[j];
				ElmHitCll[ind1] -= ind1;
				if ((ElmHitCll[ind1] > 0) && (ElmHitCll[ind1] < cls[ind1])) {
					SplCls[SplInd++] = ind1;
				}
			}
			
			switch (SplInd) {
				case 0:
				case 1:
					break;
				case 2:
					if (SplCls[0] > SplCls[1]) {
						value = SplCls[0];
						SplCls[0] = SplCls[1];
						SplCls[1] = value;
					}
					break;
				case 3:
				case 4:
				case 5:
				case 6:
				case 7:
				case 8:
					for (k = 1; k < SplInd; ++k) {
						value = SplCls[k];
						i = k - 1;
						while ((i >= 0) && (value < SplCls[i])) {
							SplCls[i + 1] = SplCls[i];
							--i;
						}
						SplCls[i + 1] = value;
					}
					break;
				default:
					quickSort(SplCls, SplInd);
					break;
			}
			
			for (j = 0; j < SplInd; j++) {
				ind1 = SplCls[j];
				i = ind1+cls[ind1]-ElmHitCll[ind1];
				trieref = trie_make(trieref, i, n, tv);
			}
			
			/* REARRANGE THE CELL */
			for (j = 0; j < SplInd; j++) {
				ind1 = SplCls[j];
				cls[ind1] = cls[ind1]-ElmHitCll[ind1];
				newcell = ind1+cls[ind1];
				cls[newcell] = ElmHitCll[ind1];
				
				Part->cells++;
				
				if (StackMarkers[ind1] != tv->stackmark) {
					if (cls[newcell] < cls[ind1]) {
						CStack[++CStackInd] = newcell;
						StackMarkers[newcell] = tv->stackmark;
					}
					else {
						CStack[++CStackInd] = ind1;
						StackMarkers[ind1] = tv->stackmark;
					}
				}
				else {
					CStack[++CStackInd] = newcell;
					StackMarkers[newcell] = tv->stackmark;
				}
				
				SplitCell = HitVtx+ind1;
				ind3 = cls[newcell];
				LabCell = lab+newcell;
				
				for (i1 = 0; i1 < ind3; i1++) {
					k = SplitCell[i1];
					i = LabCell[i1];
					Part->inv[newcell+i1] = newcell;
					lab[InvLab[k]] = i;
					InvLab[i] = InvLab[k];
					LabCell[i1] = k;
					InvLab[k] = newcell+i1;
				}
			}
		}
		else {
			if (ti->thegraphisparse) {
				if (cls[ind0] == n) {
					memcpy(NghCounts, sg->d, n*sizeof(int));
					HitCls[0] = 0;
					HitClsInd = 1;
					ElmHitCll[0] = 0;
					for (i = 0; i < n; i++) {
						if (sg->d[i]) {
							HitVtx[ElmHitCll[0]++] = i;
						}
					}
				}
				else {
					HitClsInd = 0;
                    for (i = ind0; i < ind2; i++) {
                        labi = lab[i];
                        nghb = sg->e+sg->v[labi];
                        j = sg->d[labi];
                        for (; j--;) {
                            k = nghb[j];
                            if (MarkHitVtx[k] == tv->mark) {
                                NghCounts[k]++;
                            }
                            else {
                                value = Part->inv[InvLab[k]];
                                if (cls[value] > 1) {
                                    MarkHitVtx[k] = tv->mark;
                                    NghCounts[k] = 1;
                                    if (Markers[value] != tv->mark) {
                                        HitCls[HitClsInd++] = value;
                                        Markers[value] = tv->mark;
                                        HitVtx[value] = k;
                                        ElmHitCll[value] = 1;
                                    }
                                    else {
                                        HitVtx[value+ElmHitCll[value]++] = k;
                                    }
                                }
                                else {
                                    longcode = MASHCOMM(longcode, value);
                                }
                            }
						}
					}
                    
                }
				
				tv->mark++;
				
				SplInd = 0;
				SplCls[0] = n;
				for (j = 0; j < HitClsInd; j++) {
					ind1 = HitCls[j];
					if ((ElmHitCll[ind1] > 0) && (ElmHitCll[ind1] < cls[ind1])) {
						SplCls[SplInd++] = HitCls[j];
					}
					else {
						ind2 = ind1+cls[ind1];
						Split = FALSE;
						value = NghCounts[lab[ind1++]];
						for (i = ind1; i < ind2; i++)
						{
							if (NghCounts[lab[i]] != value)
							{
								SplCls[SplInd++] = HitCls[j];
								Split = TRUE;
								break;
							}
						}
						if (!Split) {
							longcode = MASHCOMM(longcode, ind1);
						}
					}
				}
				/* Sorting the cells to be split */
				switch (SplInd) {
					case 0:
					case 1:
						break;
					case 2:
						if (SplCls[0] > SplCls[1]) {
							value = SplCls[0];
							SplCls[0] = SplCls[1];
							SplCls[1] = value;
						}
						break;
					case 3:
					case 4:
					case 5:
					case 6:
					case 7:
					case 8:
						for (k = 1; k < SplInd; ++k) {
							value = SplCls[k];
							i = k - 1;
							while ((i >= 0) && (value < SplCls[i])) {
								SplCls[i + 1] = SplCls[i];
								--i;
							}
							SplCls[i + 1] = value;
						}
						break;
					default:
						quickSort(SplCls, SplInd);
						break;
				}
				
				for (sc = 0; sc < SplInd; sc++) {	/* For each cell C to be split */
					ind0 = SplCls[sc];
					ind1 = ind0 + cls[ind0];
					SplCntInd = 0;
					if (ElmHitCll[ind0] < cls[ind0]) {
						SplCnt[SplCntInd++] = 0;
						SplPos[0] = cls[ind0] - ElmHitCll[ind0];
					}
					
					/* According to the numbers of neighbors of C into the current cell */
					/* compute how many vertices in C will be placed into the same new cell */
					iend = ind0 + ElmHitCll[ind0];
					for (i = ind0; i < iend; i++) {
						value = NghCounts[HitVtx[i]];
						if (Markers[value] != tv->mark) {
							Markers[value] = tv->mark;
							SplCnt[SplCntInd++] = value;
							SplPos[value] = 1;
						}
						else {
							SplPos[value]++;
						}
					}
					
					tv->mark++;
					
					/* Sort the values deriving from the previous step */
					switch (SplCntInd) {
						case 0:
						case 1:
							break;
						case 2:
							if (SplCnt[0] > SplCnt[1]) {
								value = SplCnt[0];
								SplCnt[0] = SplCnt[1];
								SplCnt[1] = value;
							}
							break;
						case 3:
						case 4:
						case 5:
						case 6:
						case 7:
						case 8:
							for (k = 1; k < SplCntInd; ++k) {
								value = SplCnt[k];
								i = k - 1;
								while ((i >= 0) && (value < SplCnt[i])) {
									SplCnt[i + 1] = SplCnt[i];
									--i;
								}
								SplCnt[i + 1] = value;
							}
							break;
						default:
							quickSort(SplCnt, SplCntInd);
							break;
					}
					Part->cells += SplCntInd-1;
					
					/* Split the cell C and update the information for sizes of new cells */
					/* Put the new cells into the stack */
					i = ind0;
					if (StackMarkers[i] != tv->stackmark) {
						BigCellSize = 0;
					}
					for (k = 0; k < SplCntInd; k++) {
						value = SplPos[SplCnt[k]];
						cls[i] = value;
						if ((StackMarkers[ind0] != tv->stackmark) && (value > BigCellSize)) {
							BigCell = i;
							BigCellPos = CStackInd;
							BigCellSize = cls[i];
						}
						SplPos[SplCnt[k]] = i;
						i += value;
						if (i < ind1) {
							CStack[++CStackInd] = i;
							StackMarkers[i] = tv->stackmark;
							trieref = trie_make(trieref, i, n, tv);
						}
					}
					
					if ((StackMarkers[ind0] != tv->stackmark) && (BigCell != ind0)) {
						CStack[BigCellPos] = ind0;
						StackMarkers[BigCell] = 0;
						StackMarkers[ind0] = tv->stackmark;
					}
					/* Permute elements of the cell C */
					iend = ind0 + ElmHitCll[ind0];
					for (i = ind0; i < iend; i++) {
						value = HitVtx[i];
						j = SplPos[NghCounts[value]]++; /* where HitVtx[i] goes */
						k = InvLab[value];				/* where HitVtx[i] is in lab */
						lab[k] = lab[j];
						lab[j] = value;
						InvLab[value] = j;
						InvLab[lab[k]] = k;
						NghCounts[value] = 0;
					}
					
					/* Reconstruct the cell C and update the inverse partition */
					newcell = ind1 - ElmHitCll[ind0];
					i = newcell;
					ind2 = newcell+cls[newcell]-1;
					do {
						Part->inv[i] = newcell;
						if (i == ind2) {
							newcell = i+1;
							if (newcell < n) ind2 = newcell+cls[newcell]-1;
						}
					}
					while (++i < ind1);
				}
			}
			else {
				if (sg->d[lab[ind0]] > n/cls[ind0]) {
					Sparse = FALSE;
				}
				else {
					Sparse = TRUE;
				}
				if (Sparse) {
					/* Counting occurrences of neighbors of the current cell */
					/* The list of cells with neighbors in the current cell is also built */
					if (cls[ind0] == n) {
						memcpy(NghCounts, sg->d, n*sizeof(int));
						HitCls[0] = 0;
						HitClsInd = 1;
					}
					else {
						memset(NghCounts, 0, n*sizeof(int));
						HitClsInd = 0;
						for (i = ind0; i < ind2; i++) {
							j1 = sg->v[lab[i]];
							iend1 = j1+sg->d[lab[i]];
							for (; j1 < iend1; ++j1) {
								k = sg->e[j1];
								(NghCounts[k])++;
								value = Part->inv[InvLab[k]];
								if (Markers[value] != tv->mark) {
									if (cls[value] > 1) HitCls[HitClsInd++] = value;
									Markers[value] = tv->mark;
								}
							}
						}
					}
					
					tv->mark++;
					
					SplInd = 0;
					for (j = 0; j < HitClsInd; j++) {
						ind1 = HitCls[j];
						ind2 = ind1+cls[ind1];
						value = NghCounts[lab[ind1++]];
						for (i = ind1; i < ind2; i++)
						{
							if (NghCounts[lab[i]] != value)
							{
								SplCls[SplInd++] = HitCls[j];
								break;
							}
						}
					}
					
					/* Sorting the cells to be split */
					switch (SplInd) {
						case 0:
						case 1:
							break;
						case 2:
							if (SplCls[0] > SplCls[1]) {
								value = SplCls[0];
								SplCls[0] = SplCls[1];
								SplCls[1] = value;
							}
							break;
						case 3:
						case 4:
						case 5:
						case 6:
						case 7:
						case 8:
							for (k = 1; k < SplInd; ++k) {
								value = SplCls[k];
								i = k - 1;
								while ((i >= 0) && (value < SplCls[i])) {
									SplCls[i + 1] = SplCls[i];
									--i;
								}
								SplCls[i + 1] = value;
							}
							break;
						default:
							quickSort(SplCls, SplInd);
							break;
					}
					
					for (j = 0; j < SplInd; j++) {	/* For each cell C to be split */
						ind0 = SplCls[j];
						ind1 = ind0+cls[ind0];
						SplCntInd = 0;
						
						/* According to the numbers of neighbors of C into the current cell */
						/* compute how many vertices in C will be placed into the same new cell */
						for (i = ind0; i < ind1; i++) {
							value = NghCounts[lab[i]];
							if (Markers[value] != tv->mark) {
								Markers[value] = tv->mark;
								SplCnt[SplCntInd++] = value;
								SplPos[value] = 1;
							}
							else {
								SplPos[value]++;
							}
						}
						
						tv->mark++;
						
						/* Sort the values deriving from the previous step */
						switch (SplCntInd) {
							case 0:
							case 1:
								break;
							case 2:
								if (SplCnt[0] > SplCnt[1]) {
									value = SplCnt[0];
									SplCnt[0] = SplCnt[1];
									SplCnt[1] = value;
								}
								break;
							case 3:
							case 4:
							case 5:
							case 6:
							case 7:
							case 8:
								for (k = 1; k < SplCntInd; ++k) {
									value = SplCnt[k];
									i = k - 1;
									while ((i >= 0) && (value < SplCnt[i])) {
										SplCnt[i + 1] = SplCnt[i];
										--i;
									}
									SplCnt[i + 1] = value;
								}
								break;
							default:
								quickSort(SplCnt, SplCntInd);
								break;
						}
						
						Part->cells += SplCntInd-1;
						
						/* Split the cell C and update the information for sizes of new cells */
						/* Put the new cells into the stack */
						i = ind0;
						if (StackMarkers[i] != tv->stackmark) {
							BigCellSize = 0;
						}
						for (k = 0; k < SplCntInd; k++) {
							value = SplPos[SplCnt[k]];
							cls[i] = value;
							if ((StackMarkers[ind0] != tv->stackmark) && (value > BigCellSize)) {
								BigCell = i;
								BigCellPos = CStackInd;
								BigCellSize = cls[i];
							}
							SplPos[SplCnt[k]] = i;
							i += value;
							if (i < ind1) {
								CStack[++CStackInd] = i;
								StackMarkers[i] = tv->stackmark;
								trieref = trie_make(trieref, i, n, tv);
							}
						}
						if ((StackMarkers[ind0] != tv->stackmark) && (BigCell != ind0)) {
							CStack[BigCellPos] = ind0;
							StackMarkers[BigCell] = 0;
							StackMarkers[ind0] = tv->stackmark;
						}
						
						/* Permute elements of the cell C */
						i = ind0;
						do {
							SplCnt[SplPos[NghCounts[lab[i]]]++] = lab[i];
						}
						while(++i < ind1);
						
						/* Reconstruct the cell C and update the inverse partition */
						newcell = ind0;
						i = ind0;
						ind2 = newcell+cls[newcell]-1;
						do {
							lab[i] = SplCnt[i];
							InvLab[lab[i]] = i;
							Part->inv[i] = newcell;
							if (i == ind2) {
								newcell = i+1;
								if (newcell < n) ind2 = newcell+cls[newcell]-1;
							}
						}
						while (++i < ind1);
					}
				}
				else {
					if (cls[ind0] == n) {
						memcpy(NghCounts, sg->d, n*sizeof(int));
					}
					else {
						memset(NghCounts, 0, n*sizeof(int));
						
						for (i = ind0; i < ind2; i++) {
							j1 = sg->v[lab[i]];
							iend1 = j1+sg->d[lab[i]];
							for (; j1 < iend1; ++j1) {
								k = sg->e[j1];
								(NghCounts[k])++;
							}
						}
					}
					
					ind0 = 0;
					while (ind0 < n) {	/* For each cell C with size(C) > 1 */
						ind1 = ind0+cls[ind0];
						if (cls[ind0] > 1) {
							
							/* Determine whether C must be split */
							SplCntInd = 0;
							value = NghCounts[lab[ind0]];
							for (i = ind0+1; i < ind1; i++)
							{
								if (NghCounts[lab[i]] != value)
								{
									Markers[value] = tv->mark;
									SplCnt[SplCntInd++] = value;
									SplPos[value] = i-ind0;
									do {
										value = NghCounts[lab[i++]];
										if (Markers[value] != tv->mark) {
											Markers[value] = tv->mark;
											SplCnt[SplCntInd++] = value;
											SplPos[value] = 1;
										}
										else {
											SplPos[value]++;
										}
									}
									while(i != ind1);
									break;
								}
							}
							
							if (SplCntInd) {
								tv->mark++;
								
								/* Sort the values deriving from the previous step */
								switch (SplCntInd) {
									case 0:
									case 1:
										break;
									case 2:
										if (SplCnt[0] > SplCnt[1]) {
											value = SplCnt[0];
											SplCnt[0] = SplCnt[1];
											SplCnt[1] = value;
										}
										break;
									case 3:
									case 4:
									case 5:
									case 6:
									case 7:
									case 8:
										for (k = 1; k < SplCntInd; ++k) {
											value = SplCnt[k];
											i = k - 1;
											while ((i >= 0) && (value < SplCnt[i])) {
												SplCnt[i + 1] = SplCnt[i];
												--i;
											}
											SplCnt[i + 1] = value;
										}
										break;
									default:
										quickSort(SplCnt, SplCntInd);
										break;
								}
								
								Part->cells += SplCntInd-1;
								
								/* Split the cell C and update the information for sizes of new cells */
								/* Put the new cells into the stack */
								i = ind0;
								if (StackMarkers[i] != tv->stackmark) {
									BigCellSize = 0;
								}
								for (k = 0; k < SplCntInd; k++) {
									value = SplPos[SplCnt[k]];
									cls[i] = value;
									if ((StackMarkers[ind0] != tv->stackmark) && (value > BigCellSize)) {
										BigCell = i;
										BigCellPos = CStackInd;
										BigCellSize = cls[i];
									}
									SplPos[SplCnt[k]] = i;
									i += value;
									if (i < ind1) {
										CStack[++CStackInd] = i;
										StackMarkers[i] = tv->stackmark;
									}
								}
								if ((StackMarkers[ind0] != tv->stackmark) && (BigCell != ind0)) {
									CStack[BigCellPos] = ind0;
									StackMarkers[BigCell] = 0;
									StackMarkers[ind0] = tv->stackmark;
									trieref = trie_make(trieref, i, n, tv);
								}
								
								/* Permute elements of the cell C */
								i = ind0;
								do {
									SplCnt[SplPos[NghCounts[lab[i]]]++] = lab[i];
								}
								while(++i < ind1);
								
								/* Reconstruct the cell C and update the inverse partition */
								newcell = ind0;
								i = ind0;
								ind2 = newcell+cls[newcell]-1;
								do {
									lab[i] = SplCnt[i];
									InvLab[lab[i]] = i;
									Part->inv[i] = newcell;
									if (i == ind2) {
										newcell = i+1;
										if (newcell < n) ind2 = newcell+cls[newcell]-1;
									}
								}
								while (++i < ind1);
							}
						}
						ind0 = ind1;
					}
				}
			}
		}
	}
	
	Cand->code = CLEANUP(longcode);
	return;
}

int traces_refine_comptrie(sparsegraph *sg,
						   Candidate *Cand,
						   int m,
						   int n,
						   Partition *Part,
						   struct TracesVars* tv,
						   struct TracesInfo *ti)
{
	int i, i1, j, k, sc, ind0, ind1, ind2, ind3, labi;
	int value, iend, newcell, TcSize;
	int HitClsInd, SplInd, SplCntInd, CStackInd;
	size_t j1, iend1;
	unsigned int longcode;
	int Split = 0;
	int Sparse = TRUE;
	int *lab, *cls, *InvLab, *SplitCell, *LabCell;
	int BigCell, BigCellPos, BigCellSize;
    int *nghb;
    
	TcSize = Spine[tv->tolevel].tgtsize;
	
	if (tv->stackmark > (NAUTY_INFINITY-2)) {
		memset(StackMarkers, 0, n*sizeof(int));
		tv->stackmark = 0;
	}
	tv->stackmark++;
	
	lab = Cand->lab;
	InvLab = Cand->invlab;
	cls = Part->cls;
	
	CStack[1] = Spine[tv->tolevel_tl].tgtpos;
	CStackInd = 1;
	for (i = 1; i <= CStackInd; i++) {
		StackMarkers[CStack[i]] = tv->stackmark;
	}
	
	longcode = Part->cells;
	while (CStackInd > 0)
	{
		if (tv->mark > (NAUTY_INFINITY-2)) {
			memset(Markers, 0, n*sizeof(int));
			memset(MarkHitVtx, 0, n*sizeof(int));
			tv->mark = 0;
		}
		tv->mark++;
		
		if (Part->cells == tv->finalnumcells) break;
		
		j = CStackInd;
		k = CStackInd;
		while (--j > 0) {
			if (cls[CStack[j]] < cls[CStack[k]]) {
				k = j;
			}
			if ((cls[CStack[k]] == 1) || (j < CStackInd - 12)) {
				break;
			}
		}
		
		ind0 = CStack[k];
		ind2 = ind0+cls[ind0];
		CStack[k] = CStack[CStackInd--];		/* Current Cell */
		longcode = MASHNONCOMM(longcode, ind0);
		StackMarkers[ind0] = 0;
		
		/* Analysis of occurrences of neighbors of the current cell */
		/* The list of cells with neighbors in the current cell is  built */
		if (cls[ind0] == 1) {  /* SINGLETON CURRENT CELL CASE */
			HitClsInd = 0;
			j1 = sg->v[lab[ind0]];
			iend1 = j1+sg->d[lab[ind0]];
			
			for (; j1 < iend1; ++j1) {
				k = sg->e[j1];
				value = Part->inv[InvLab[k]];
				if (cls[value] > 1) {
					if (Markers[value] != tv->mark) {
						HitCls[HitClsInd++] = value;
						Markers[value] = tv->mark;
						ElmHitCll[value] = value;
					}
					HitVtx[ElmHitCll[value]++] = k;
				}
				else {
					longcode = MASHCOMM(longcode, value);
				}
			}
			
			tv->mark++;
			
			SplInd = 0;
			for (j = 0; j < HitClsInd; j++) {
				ind1 = HitCls[j];
				ElmHitCll[ind1] -= ind1;
				if ((ElmHitCll[ind1] > 0) && (ElmHitCll[ind1] < cls[ind1])) {
					SplCls[SplInd++] = ind1;
				}
			}
			
			switch (SplInd) {
				case 0:
				case 1:
					break;
				case 2:
					if (SplCls[0] > SplCls[1]) {
						value = SplCls[0];
						SplCls[0] = SplCls[1];
						SplCls[1] = value;
					}
					break;
				case 3:
				case 4:
				case 5:
				case 6:
				case 7:
				case 8:
					for (k = 1; k < SplInd; ++k) {
						value = SplCls[k];
						i = k - 1;
						while ((i >= 0) && (value < SplCls[i])) {
							SplCls[i + 1] = SplCls[i];
							--i;
						}
						SplCls[i + 1] = value;
					}
					break;
				default:
					quickSort(SplCls, SplInd);
					break;
			}
			
			for (j = 0; j < SplInd; j++) {
				ind1 = SplCls[j];
				i = ind1+cls[ind1]-ElmHitCll[ind1];
				trieref = trie_comp(trieref, i);
				if (trieref == NULL) return 0;
			}
			
			/* REARRANGE THE CELL */
			for (j = 0; j < SplInd; j++) {
				ind1 = SplCls[j];
				cls[ind1] = cls[ind1]-ElmHitCll[ind1];
				newcell = ind1+cls[ind1];
				cls[newcell] = ElmHitCll[ind1];
				
				Part->cells++;
				if (StackMarkers[ind1] != tv->stackmark) {
					if (cls[newcell] < cls[ind1]) {
						CStack[++CStackInd] = newcell;
						StackMarkers[newcell] = tv->stackmark;
					}
					else {
						CStack[++CStackInd] = ind1;
						StackMarkers[ind1] = tv->stackmark;
					}
				}
				else {
					CStack[++CStackInd] = newcell;
					StackMarkers[newcell] = tv->stackmark;
				}
				SplitCell = HitVtx+ind1;
				ind3 = cls[newcell];
				LabCell = lab+newcell;
				
				for (i1 = 0; i1 < ind3; i1++) {
					k = SplitCell[i1];
					i = LabCell[i1];
					Part->inv[newcell+i1] = newcell;
					lab[InvLab[k]] = i;
					InvLab[i] = InvLab[k];
					LabCell[i1] = k;
					InvLab[k] = newcell+i1;
				}
			}
		}
		else {
			if (ti->thegraphisparse) {
				if (cls[ind0] == n) {
					memcpy(NghCounts, sg->d, n*sizeof(int));
					HitCls[0] = 0;
					HitClsInd = 1;
					ElmHitCll[0] = 0;
					for (i = 0; i < n; i++) {
						if (sg->d[i]) {
							HitVtx[ElmHitCll[0]++] = i;
						}
					}
				}
				else {
					HitClsInd = 0;
                    for (i = ind0; i < ind2; i++) {
                        labi = lab[i];
                        nghb = sg->e+sg->v[labi];
                        j = sg->d[labi];
                        for (; j--;) {
                            k = nghb[j];
                            if (MarkHitVtx[k] == tv->mark) {
                                NghCounts[k]++;
                            }
                            else {
                                value = Part->inv[InvLab[k]];
                                if (cls[value] > 1) {
                                    MarkHitVtx[k] = tv->mark;
                                    NghCounts[k] = 1;
                                    if (Markers[value] != tv->mark) {
                                        HitCls[HitClsInd++] = value;
                                        Markers[value] = tv->mark;
                                        HitVtx[value] = k;
                                        ElmHitCll[value] = 1;
                                    }
                                    else {
                                        HitVtx[value+ElmHitCll[value]++] = k;
                                    }
                                }
                                else {
                                    longcode = MASHCOMM(longcode, value);
                                }
                            }
						}
					}
                    
                }
				
				tv->mark++;
				
				SplInd = 0;
				SplCls[0] = n;
				for (j = 0; j < HitClsInd; j++) {
					ind1 = HitCls[j];
					if ((ElmHitCll[ind1] > 0) && (ElmHitCll[ind1] < cls[ind1])) {
						SplCls[SplInd++] = HitCls[j];
					}
					else {
						ind2 = ind1+cls[ind1];
						Split = FALSE;
						value = NghCounts[lab[ind1++]];
						for (i = ind1; i < ind2; i++)
						{
							if (NghCounts[lab[i]] != value)
							{
								SplCls[SplInd++] = HitCls[j];
								Split = TRUE;
								break;
							}
						}
						if (!Split) {
							longcode = MASHCOMM(longcode, ind1);
						}
					}
				}
				/* Sorting the cells to be split */
				switch (SplInd) {
					case 0:
					case 1:
						break;
					case 2:
						if (SplCls[0] > SplCls[1]) {
							value = SplCls[0];
							SplCls[0] = SplCls[1];
							SplCls[1] = value;
						}
						break;
					case 3:
					case 4:
					case 5:
					case 6:
					case 7:
					case 8:
						for (k = 1; k < SplInd; ++k) {
							value = SplCls[k];
							i = k - 1;
							while ((i >= 0) && (value < SplCls[i])) {
								SplCls[i + 1] = SplCls[i];
								--i;
							}
							SplCls[i + 1] = value;
						}
						break;
					default:
						quickSort(SplCls, SplInd);
						break;
				}
				
				for (sc = 0; sc < SplInd; sc++) {	/* For each cell C to be split */
					ind0 = SplCls[sc];
					ind1 = ind0 + cls[ind0];
					SplCntInd = 0;
					if (ElmHitCll[ind0] < cls[ind0]) {
						SplCnt[SplCntInd++] = 0;
						SplPos[0] = cls[ind0] - ElmHitCll[ind0];
					}
					
					/* According to the numbers of neighbors of C into the current cell */
					/* compute how many vertices in C will be placed into the same new cell */
					iend = ind0 + ElmHitCll[ind0];
					for (i = ind0; i < iend; i++) {
						value = NghCounts[HitVtx[i]];
						if (Markers[value] != tv->mark) {
							Markers[value] = tv->mark;
							SplCnt[SplCntInd++] = value;
							SplPos[value] = 1;
						}
						else {
							SplPos[value]++;
						}
					}
					
					tv->mark++;
					
					/* Sort the values deriving from the previous step */
					switch (SplCntInd) {
						case 0:
						case 1:
							break;
						case 2:
							if (SplCnt[0] > SplCnt[1]) {
								value = SplCnt[0];
								SplCnt[0] = SplCnt[1];
								SplCnt[1] = value;
							}
							break;
						case 3:
						case 4:
						case 5:
						case 6:
						case 7:
						case 8:
							for (k = 1; k < SplCntInd; ++k) {
								value = SplCnt[k];
								i = k - 1;
								while ((i >= 0) && (value < SplCnt[i])) {
									SplCnt[i + 1] = SplCnt[i];
									--i;
								}
								SplCnt[i + 1] = value;
							}
							break;
						default:
							quickSort(SplCnt, SplCntInd);
							break;
					}
					Part->cells += SplCntInd-1;
					
					/* Split the cell C and update the information for sizes of new cells */
					/* Put the new cells into the stack */
					i = ind0;
					if (StackMarkers[i] != tv->stackmark) {
						BigCellSize = 0;
					}
					for (k = 0; k < SplCntInd; k++) {
						value = SplPos[SplCnt[k]];
						cls[i] = value;
						if ((StackMarkers[ind0] != tv->stackmark) && (value > BigCellSize)) {
							BigCell = i;
							BigCellPos = CStackInd;
							BigCellSize = cls[i];
						}
						SplPos[SplCnt[k]] = i;
						i += value;
						if (i < ind1) {
							CStack[++CStackInd] = i;
							StackMarkers[i] = tv->stackmark;
							trieref = trie_comp(trieref, i);
							if (trieref == NULL) return 0;
						}
					}
					
					if ((StackMarkers[ind0] != tv->stackmark) && (BigCell != ind0)) {
						CStack[BigCellPos] = ind0;
						StackMarkers[BigCell] = 0;
						StackMarkers[ind0] = tv->stackmark;
					}
					/* Permute elements of the cell C */
					iend = ind0 + ElmHitCll[ind0];
					for (i = ind0; i < iend; i++) {
						value = HitVtx[i];
						j = SplPos[NghCounts[value]]++; /* where HitVtx[i] goes */
						k = InvLab[value];				/* where HitVtx[i] is in lab */
						lab[k] = lab[j];
						lab[j] = value;
						InvLab[value] = j;
						InvLab[lab[k]] = k;
						NghCounts[value] = 0;
					}
					
					/* Reconstruct the cell C and update the inverse partition */
					newcell = ind1 - ElmHitCll[ind0];
					i = newcell;
					ind2 = newcell+cls[newcell]-1;
					do {
						Part->inv[i] = newcell;
						if (i == ind2) {
							newcell = i+1;
							if (newcell < n) ind2 = newcell+cls[newcell]-1;
						}
					}
					while (++i < ind1);
				}
			}
			else {
				if (sg->d[lab[ind0]] > n/cls[ind0]) {
					Sparse = FALSE;
				}
				else {
					Sparse = TRUE;
				}
				if (Sparse) {
					/* Counting occurrences of neighbors of the current cell */
					/* The list of cells with neighbors in the current cell is also built */
					if (cls[ind0] == n) {
						memcpy(NghCounts, sg->d, n*sizeof(int));
						HitCls[0] = 0;
						HitClsInd = 1;
					}
					else {
						memset(NghCounts, 0, n*sizeof(int));
						HitClsInd = 0;
						for (i = ind0; i < ind2; i++) {
							j1 = sg->v[lab[i]];
							iend1 = j1+sg->d[lab[i]];
							for (; j1 < iend1; ++j1) {
								k = sg->e[j1];
								(NghCounts[k])++;
								value = Part->inv[InvLab[k]];
								if (Markers[value] != tv->mark) {
									if (cls[value] > 1) HitCls[HitClsInd++] = value;
									Markers[value] = tv->mark;
								}
							}
						}
					}
					
					tv->mark++;
					
					SplInd = 0;
					for (j = 0; j < HitClsInd; j++) {
						ind1 = HitCls[j];
						ind2 = ind1+cls[ind1];
						value = NghCounts[lab[ind1++]];
						for (i = ind1; i < ind2; i++)
						{
							if (NghCounts[lab[i]] != value)
							{
								SplCls[SplInd++] = HitCls[j];
								break;
							}
						}
					}
					
					/* Sorting the cells to be split */
					switch (SplInd) {
						case 0:
						case 1:
							break;
						case 2:
							if (SplCls[0] > SplCls[1]) {
								value = SplCls[0];
								SplCls[0] = SplCls[1];
								SplCls[1] = value;
							}
							break;
						case 3:
						case 4:
						case 5:
						case 6:
						case 7:
						case 8:
							for (k = 1; k < SplInd; ++k) {
								value = SplCls[k];
								i = k - 1;
								while ((i >= 0) && (value < SplCls[i])) {
									SplCls[i + 1] = SplCls[i];
									--i;
								}
								SplCls[i + 1] = value;
							}
							break;
						default:
							quickSort(SplCls, SplInd);
							break;
					}
					
					for (j = 0; j < SplInd; j++) {	/* For each cell C to be split */
						ind0 = SplCls[j];
						ind1 = ind0+cls[ind0];
						SplCntInd = 0;
						
						/* According to the numbers of neighbors of C into the current cell */
						/* compute how many vertices in C will be placed into the same new cell */
						for (i = ind0; i < ind1; i++) {
							value = NghCounts[lab[i]];
							if (Markers[value] != tv->mark) {
								Markers[value] = tv->mark;
								SplCnt[SplCntInd++] = value;
								SplPos[value] = 1;
							}
							else {
								SplPos[value]++;
							}
						}
						
						tv->mark++;
						
						/* Sort the values deriving from the previous step */
						switch (SplCntInd) {
							case 0:
							case 1:
								break;
							case 2:
								if (SplCnt[0] > SplCnt[1]) {
									value = SplCnt[0];
									SplCnt[0] = SplCnt[1];
									SplCnt[1] = value;
								}
								break;
							case 3:
							case 4:
							case 5:
							case 6:
							case 7:
							case 8:
								for (k = 1; k < SplCntInd; ++k) {
									value = SplCnt[k];
									i = k - 1;
									while ((i >= 0) && (value < SplCnt[i])) {
										SplCnt[i + 1] = SplCnt[i];
										--i;
									}
									SplCnt[i + 1] = value;
								}
								break;
							default:
								quickSort(SplCnt, SplCntInd);
								break;
						}
						
						Part->cells += SplCntInd-1;
						
						/* Split the cell C and update the information for sizes of new cells */
						/* Put the new cells into the stack */
						i = ind0;
						if (StackMarkers[i] != tv->stackmark) {
							BigCellSize = 0;
						}
						for (k = 0; k < SplCntInd; k++) {
							value = SplPos[SplCnt[k]];
							cls[i] = value;
							if ((StackMarkers[ind0] != tv->stackmark) && (value > BigCellSize)) {
								BigCell = i;
								BigCellPos = CStackInd;
								BigCellSize = cls[i];
							}
							SplPos[SplCnt[k]] = i;
							i += value;
							if (i < ind1) {
								CStack[++CStackInd] = i;
								StackMarkers[i] = tv->stackmark;
								trieref = trie_comp(trieref, i);
								if (trieref == NULL) return 0;
							}
						}
						if ((StackMarkers[ind0] != tv->stackmark) && (BigCell != ind0)) {
							CStack[BigCellPos] = ind0;
							StackMarkers[BigCell] = 0;
							StackMarkers[ind0] = tv->stackmark;
						}
						
						/* Permute elements of the cell C */
						i = ind0;
						do {
							SplCnt[SplPos[NghCounts[lab[i]]]++] = lab[i];
						}
						while(++i < ind1);
						
						/* Reconstruct the cell C and update the inverse partition */
						newcell = ind0;
						i = ind0;
						ind2 = newcell+cls[newcell]-1;
						do {
							lab[i] = SplCnt[i];
							InvLab[lab[i]] = i;
							Part->inv[i] = newcell;
							
							if (i == ind2) {
								newcell = i+1;
								if (newcell < n) ind2 = newcell+cls[newcell]-1;
							}
						}
						while (++i < ind1);
					}
				}
				else {
					if (cls[ind0] == n) {
						memcpy(NghCounts, sg->d, n*sizeof(int));
					}
					else {
						memset(NghCounts, 0, n*sizeof(int));
						
						for (i = ind0; i < ind2; i++) {
							j1 = sg->v[lab[i]];
							iend1 = j1+sg->d[lab[i]];
							for (; j1 < iend1; ++j1) {
								k = sg->e[j1];
								(NghCounts[k])++;
							}
						}
					}
					
					ind0 = 0;
					while (ind0 < n) {	/* For each cell C with size(C) > 1 */
						ind1 = ind0+cls[ind0];
						if (cls[ind0] > 1) {
							
							/* Determine whether C must be split */
							SplCntInd = 0;
							value = NghCounts[lab[ind0]];
							for (i = ind0+1; i < ind1; i++)
							{
								if (NghCounts[lab[i]] != value)
								{
									Markers[value] = tv->mark;
									SplCnt[SplCntInd++] = value;
									SplPos[value] = i-ind0;
									do {
										value = NghCounts[lab[i++]];
										if (Markers[value] != tv->mark) {
											Markers[value] = tv->mark;
											SplCnt[SplCntInd++] = value;
											SplPos[value] = 1;
										}
										else {
											SplPos[value]++;
										}
									}
									while(i != ind1);
									break;
								}
							}
							
							if (SplCntInd) {
								tv->mark++;
								
								/* Sort the values deriving from the previous step */
								switch (SplCntInd) {
									case 0:
									case 1:
										break;
									case 2:
										if (SplCnt[0] > SplCnt[1]) {
											value = SplCnt[0];
											SplCnt[0] = SplCnt[1];
											SplCnt[1] = value;
										}
										break;
									case 3:
									case 4:
									case 5:
									case 6:
									case 7:
									case 8:
										for (k = 1; k < SplCntInd; ++k) {
											value = SplCnt[k];
											i = k - 1;
											while ((i >= 0) && (value < SplCnt[i])) {
												SplCnt[i + 1] = SplCnt[i];
												--i;
											}
											SplCnt[i + 1] = value;
										}
										break;
									default:
										quickSort(SplCnt, SplCntInd);
										break;
								}
								
								Part->cells += SplCntInd-1;
								
								/* Split the cell C and update the information for sizes of new cells */
								/* Put the new cells into the stack */
								i = ind0;
								if (StackMarkers[i] != tv->stackmark) {
									BigCellSize = 0;
								}
								for (k = 0; k < SplCntInd; k++) {
									value = SplPos[SplCnt[k]];
									cls[i] = value;
									if ((StackMarkers[ind0] != tv->stackmark) && (value > BigCellSize)) {
										BigCell = i;
										BigCellPos = CStackInd;
										BigCellSize = cls[i];
									}
									SplPos[SplCnt[k]] = i;
									i += value;
									if (i < ind1) {
										CStack[++CStackInd] = i;
										StackMarkers[i] = tv->stackmark;
									}
								}
								if ((StackMarkers[ind0] != tv->stackmark) && (BigCell != ind0)) {
									CStack[BigCellPos] = ind0;
									StackMarkers[BigCell] = 0;
									StackMarkers[ind0] = tv->stackmark;
									trieref = trie_comp(trieref, i);
									if (trieref == NULL) return 0;
								}
								
								/* Permute elements of the cell C */
								i = ind0;
								do {
									SplCnt[SplPos[NghCounts[lab[i]]]++] = lab[i];
								}
								while(++i < ind1);
								
								/* Reconstruct the cell C and update the inverse partition */
								newcell = ind0;
								i = ind0;
								ind2 = newcell+cls[newcell]-1;
								do {
									lab[i] = SplCnt[i];
									InvLab[lab[i]] = i;
									Part->inv[i] = newcell;
									if (i == ind2) {
										newcell = i+1;
										if (newcell < n) ind2 = newcell+cls[newcell]-1;
									}
								}
								while (++i < ind1);
							}
						}
						ind0 = ind1;
					}
				}
			}
		}
	}
	
	Cand->code = CLEANUP(longcode);
	return 1;
}

int traces_refine_sametrace(sparsegraph *sg,
							Candidate *Cand,
							int m,
							int n,
							Partition *Part,
							struct TracesVars* tv,
							struct TracesInfo *ti)
{
	int i, j, k, sc, ind0, ind1, ind2, ind3, ind4, jk, labi;
	int value, iend, newcell;
	int HitClsInd, SplInd, SplCntInd, CStackInd, TraceInd, TraceCCInd, TraceStepsInd;
	size_t j1, iend1;
	unsigned int longcode;
	int Sparse = TRUE;
	int *lab, *cls, *InvLab, *TracePos, *SplitCell, *LabCell, *TraceEnd, Traceccend, *Tracestpend;
	int BigCell, BigCellPos, BigCellSize;
	boolean TraceCell = FALSE;
    int *nghb;
    
	if (tv->stackmark > (NAUTY_INFINITY-2)) {
		memset(StackMarkers, 0, n*sizeof(int));
		tv->stackmark = 0;
	}
	tv->stackmark++;
	
	SpineTL = Spine+tv->tolevel;
	TraceEnd = &(SpineTL->trcend);
	Traceccend = SpineTL->ccend;
	Tracestpend = &(SpineTL->stpend);
	TraceCCInd = SpineTL->ccstart;
	TraceStepsInd = SpineTL->stpstart;
	
	lab = Cand->lab;
	InvLab = Cand->invlab;
	cls = Part->cls;
	
	UPDATEMIN(Part->active, n-1);
	memcpy(CStack+1, TheTrace+SpineTL->trcstart, (Part->active)*sizeof(int));
	CStackInd = Part->active;
	for (i = 1; i <= CStackInd; i++) {
		StackMarkers[CStack[i]] = tv->stackmark;
	}
	
	longcode = Part->cells;
	TraceInd = SpineTL->trcstart+Part->active;
	
	while (CStackInd > 0)
	{
        
		if (tv->mark > (NAUTY_INFINITY-2)) {
			memset(Markers, 0, n*sizeof(int));
			memset(MarkHitVtx, 0, n*sizeof(int));
			tv->mark = 0;
		}
		tv->mark++;
		
		if (Part->cells == tv->finalnumcells) break;
		
		j = CStackInd;
		k = CStackInd;
		while (--j > 0) {
			if (cls[CStack[j]] < cls[CStack[k]]) {
				k = j;
			}
			if ((cls[CStack[k]] == 1) || (j < CStackInd - 12)) {
				break;
			}
		}
		
		ind0 = CStack[k];
		ind2 = ind0+cls[ind0];
		CStack[k] = CStack[CStackInd--];
		StackMarkers[ind0] = 0;
        
        
		TraceCell = ((ind0 == TheTraceCC[TraceCCInd]) && (TraceCCInd < Traceccend));
        
		/* Analysis of occurrences of neighbors of the current cell */
		/* The list of cells with neighbors in the current cell is  built */
		if (cls[ind0] == 1) {  /* SINGLETON CURRENT CELL CASE */
			HitClsInd = 0;
			j1 = sg->v[lab[ind0]];
			iend1 = j1+sg->d[lab[ind0]];
			
			for (; j1 < iend1; ++j1) {
				k = sg->e[j1];
				value = Part->inv[InvLab[k]];
				if (cls[value] > 1) {
					if (Markers[value] != tv->mark) {
						HitCls[HitClsInd++] = value;
						Markers[value] = tv->mark;
						ElmHitCll[value] = value;
					}
					HitVtx[ElmHitCll[value]++] = k;
				}
			}
			tv->mark++;
			
			SplInd = 0;
			for (j = 0; j < HitClsInd; j++) {
				ind1 = HitCls[j];
				ElmHitCll[ind1] -= ind1;
				if ((ElmHitCll[ind1] > 0) && (ElmHitCll[ind1] < cls[ind1])) {
					SplCls[SplInd++] = ind1;
				}
			}
            
            /* SINGLETON CC CASE */
			if (SplInd) {
				SAMETRACE_CHECK(TheTraceSplNum, TraceCCInd, SplInd, &Traceccend)
                
				
				switch (SplInd) {
					case 0:
					case 1:
						break;
					case 2:
						if (SplCls[0] > SplCls[1]) {
							value = SplCls[0];
							SplCls[0] = SplCls[1];
							SplCls[1] = value;
						}
						break;
					case 3:
					case 4:
					case 5:
					case 6:
					case 7:
					case 8:
						for (k = 1; k < SplInd; ++k) {
							value = SplCls[k];
							i = k - 1;
							while ((i >= 0) && (value < SplCls[i])) {
								SplCls[i + 1] = SplCls[i];
								--i;
							}
							SplCls[i + 1] = value;
						}
						break;
					default:
						quickSort(SplCls, SplInd);
						break;
				}
				
				for (j = 0; j < SplInd; j++) {
					ind1 = SplCls[j];
					i = ind1+cls[ind1]-ElmHitCll[ind1];
                    SAMETRACE_CHECK(TheTrace, TraceInd, i, TraceEnd)
				}
				
				/* REARRANGE THE CELLS */
				for (j = 0; j < SplInd; j++) {
					ind1 = SplCls[j];
					cls[ind1] = cls[ind1]-ElmHitCll[ind1];
					newcell = ind1+cls[ind1];
					cls[newcell] = ElmHitCll[ind1];
					Part->cells++;
					
					if (StackMarkers[ind1] != tv->stackmark) {
						if (cls[newcell] < cls[ind1]) {
							CStack[++CStackInd] = newcell;
							StackMarkers[newcell] = tv->stackmark;
						}
						else {
							CStack[++CStackInd] = ind1;
							StackMarkers[ind1] = tv->stackmark;
						}
					}
					else {
						CStack[++CStackInd] = newcell;
						StackMarkers[newcell] = tv->stackmark;
					}
					
					SplitCell = HitVtx+ind1;
					ind3 = cls[newcell];
					LabCell = lab+newcell;
					
					for (jk = 0; jk < ind3; jk++) {
						k = SplitCell[jk];
						i = LabCell[jk];
						Part->inv[newcell+jk] = newcell;
						lab[InvLab[k]] = i;
						InvLab[i] = InvLab[k];
						LabCell[jk] = k;
						InvLab[k] = newcell+jk;
					}
					if (cls[ind1] == 1) {
						Cand->singcode = MASHCOMM(Cand->singcode, Cand->lab[ind1]);
                    }
					if (cls[newcell] == 1) {
						Cand->singcode = MASHCOMM(Cand->singcode, Cand->lab[newcell]);
					}
				}
			}
			else {
				if (TraceCell) {
                    return FALSE;
				}
			}
			
		}
		else {
			if (ti->thegraphisparse) {
				if (cls[ind0] == n) {
					memcpy(NghCounts, sg->d, n*sizeof(int));
					HitCls[0] = 0;
					HitClsInd = 1;
					ElmHitCll[0] = 0;
					for (i = 0; i < n; i++) {
						if (sg->d[i]) {
							HitVtx[ElmHitCll[0]++] = i;
						}
					}
				}
				else {
					HitClsInd = 0;
                    for (i = ind0; i < ind2; i++) {
                        labi = lab[i];
                        nghb = sg->e+sg->v[labi];
                        j = sg->d[labi];
                        for (; j--;) {
                            k = nghb[j];
                            if (MarkHitVtx[k] == tv->mark) {
                                NghCounts[k]++;
                            }
                            else {
                                value = Part->inv[InvLab[k]];
                                if (cls[value] > 1) {
                                    MarkHitVtx[k] = tv->mark;
                                    NghCounts[k] = 1;
                                    if (Markers[value] != tv->mark) {
                                        HitCls[HitClsInd++] = value;
                                        Markers[value] = tv->mark;
                                        HitVtx[value] = k;
                                        ElmHitCll[value] = 1;
                                    }
                                    else {
                                        HitVtx[value+ElmHitCll[value]++] = k;
                                    }
                                }
                            }
						}
					}
                    
                }
				
				tv->mark++;
				
				SplInd = 0;
				SplCls[0] = n;
				for (j = 0; j < HitClsInd; j++) {
					ind1 = HitCls[j];
					if ((ElmHitCll[ind1] > 0) && (ElmHitCll[ind1] < cls[ind1])) {
						SplCls[SplInd++] = HitCls[j];
					}
					else {
						ind2 = ind1+cls[ind1];
						value = NghCounts[lab[ind1++]];
						for (i = ind1; i < ind2; i++)
						{
							if (NghCounts[lab[i]] != value)
							{
								SplCls[SplInd++] = HitCls[j];
								break;
							}
						}
					}
				}
				
                /* SPARSE CASE */
				if (SplInd) {
                    SAMETRACE_CHECK(TheTraceSplNum, TraceCCInd, SplInd+n, &Traceccend)
					
					/* Sorting the cells to be split */
					switch (SplInd) {
						case 0:
						case 1:
							break;
						case 2:
							if (SplCls[0] > SplCls[1]) {
								value = SplCls[0];
								SplCls[0] = SplCls[1];
								SplCls[1] = value;
							}
							break;
						case 3:
						case 4:
						case 5:
						case 6:
						case 7:
						case 8:
							for (k = 1; k < SplInd; ++k) {
								value = SplCls[k];
								i = k - 1;
								while ((i >= 0) && (value < SplCls[i])) {
									SplCls[i + 1] = SplCls[i];
									--i;
								}
								SplCls[i + 1] = value;
							}
							break;
						default:
							quickSort(SplCls, SplInd);
							break;
					}
					
					for (sc = 0; sc < SplInd; sc++) {	/* For each cell C to be split */
						ind0 = SplCls[sc];
						ind1 = ind0 + cls[ind0];
						SplCntInd = 0;
						if (ElmHitCll[ind0] < cls[ind0]) {
							SplCnt[SplCntInd++] = 0;
							SplPos[0] = cls[ind0] - ElmHitCll[ind0];
						}
						
						/* According to the numbers of neighbors of C into the current cell */
						/* compute how many vertices in C will be placed into the same new cell */
						iend = ind0 + ElmHitCll[ind0];
						for (i = ind0; i < iend; i++) {
							value = NghCounts[HitVtx[i]];
							if (Markers[value] != tv->mark) {
								Markers[value] = tv->mark;
								SplCnt[SplCntInd++] = value;
								SplPos[value] = 1;
							}
							else {
								SplPos[value]++;
							}
						}
						tv->mark++;
						
						if (SplCntInd) {
                            SAMETRACE_CHECK(TheTraceSteps, TraceStepsInd, SplCntInd+n, Tracestpend)
						}
						
						/* Sort the values deriving from the previous step */
						switch (SplCntInd) {
							case 0:
							case 1:
								break;
							case 2:
								if (SplCnt[0] > SplCnt[1]) {
									value = SplCnt[0];
									SplCnt[0] = SplCnt[1];
									SplCnt[1] = value;
								}
								break;
							case 3:
							case 4:
							case 5:
							case 6:
							case 7:
							case 8:
								for (k = 1; k < SplCntInd; ++k) {
									value = SplCnt[k];
									i = k - 1;
									while ((i >= 0) && (value < SplCnt[i])) {
										SplCnt[i + 1] = SplCnt[i];
										--i;
									}
									SplCnt[i + 1] = value;
								}
								break;
							default:
								quickSort(SplCnt, SplCntInd);
								break;
						}
						
						Part->cells += SplCntInd-1;
						
						/* Split the cell C and update the information for sizes of new cells */
						/* Put the new cells into the stack */
						i = ind0;
						if (StackMarkers[i] != tv->stackmark) {
							BigCellSize = 0;
						}
						for (k = 0; k < SplCntInd; k++) {
							value = SplPos[SplCnt[k]];
							cls[i] = value;
							if ((StackMarkers[ind0] != tv->stackmark) && (value > BigCellSize)) {
								BigCell = i;
								BigCellPos = CStackInd;
								BigCellSize = cls[i];
							}
							SplPos[SplCnt[k]] = i;
							i += value;
							if (i < ind1) {
								CStack[++CStackInd] = i;
								StackMarkers[i] = tv->stackmark;
                                SAMETRACE_CHECK(TheTrace, TraceInd, i, TraceEnd)
							}
						}
						
						if ((StackMarkers[ind0] != tv->stackmark) && (BigCell != ind0)) {
							CStack[BigCellPos] = ind0;
							StackMarkers[BigCell] = 0;
							StackMarkers[ind0] = tv->stackmark;
						}
						/* Permute elements of the cell C */
						iend = ind0 + ElmHitCll[ind0];
						for (i = ind0; i < iend; i++) {
							value = HitVtx[i];
							j = SplPos[NghCounts[value]]++; /* where HitVtx[i] goes */
							k = InvLab[value];				/* where HitVtx[i] is in lab */
							lab[k] = lab[j];
							lab[j] = value;
							InvLab[value] = j;
							InvLab[lab[k]] = k;
							NghCounts[value] = 0;
						}
						
						/* Reconstruct the cell C and update the inverse partition */
						newcell = ind1 - ElmHitCll[ind0];
						i = newcell;
						ind2 = newcell+cls[newcell]-1;
						do {
							Part->inv[i] = newcell;
							if (i == ind2) {
								newcell = i+1;
								if (newcell < n) ind2 = newcell+cls[newcell]-1;
							}
						}
						while (++i < ind1);
						
						for (i = ind0, k = 0; k < SplCntInd; i+=cls[i], k++) {
							if (cls[i] == 1) {
								Cand->singcode = MASHCOMM(Cand->singcode, Cand->lab[i]);
							}
						}
						
					}
				}
				else {
                    if (TraceCell) {
                        return FALSE;
                    }
				}
				
			}
			else {
				if (sg->d[lab[ind0]] > n/cls[ind0]) {
					Sparse = FALSE;
				}
				else {
					Sparse = TRUE;
				}
				if (Sparse) {
					/* Counting occurrences of neighbors of the current cell */
					/* The list of cells with neighbors in the current cell is also built */
					if (cls[ind0] == n) {
						memcpy(NghCounts, sg->d, n*sizeof(int));
						HitCls[0] = 0;
						HitClsInd = 1;
					}
					else {
						memset(NghCounts, 0, n*sizeof(int));
						HitClsInd = 0;
						for (i = ind0; i < ind2; i++) {
							j1 = sg->v[lab[i]];
							iend1 = j1+sg->d[lab[i]];
							for (; j1 < iend1; ++j1) {
								k = sg->e[j1];
								(NghCounts[k])++;
								value = Part->inv[InvLab[k]];
								if (Markers[value] != tv->mark) {
									if (cls[value] > 1) HitCls[HitClsInd++] = value;
									Markers[value] = tv->mark;
								}
							}
						}
					}
					
					tv->mark++;
					
					SplInd = 0;
					for (j = 0; j < HitClsInd; j++) {
						ind1 = HitCls[j];
						ind2 = ind1+cls[ind1];
						value = NghCounts[lab[ind1++]];
						for (i = ind1; i < ind2; i++)
						{
							if (NghCounts[lab[i]] != value)
							{
								SplCls[SplInd++] = HitCls[j];
								break;
							}
						}
					}
					
                    /* DENSE-SPARSE CASE */
					if (SplInd) {
                        SAMETRACE_CHECK(TheTraceSplNum, TraceCCInd, SplInd+2*n, &Traceccend)
						
						/* Sorting the cells to be split */
						switch (SplInd) {
							case 0:
							case 1:
								break;
							case 2:
								if (SplCls[0] > SplCls[1]) {
									value = SplCls[0];
									SplCls[0] = SplCls[1];
									SplCls[1] = value;
								}
								break;
							case 3:
							case 4:
							case 5:
							case 6:
							case 7:
							case 8:
								for (k = 1; k < SplInd; ++k) {
									value = SplCls[k];
									i = k - 1;
									while ((i >= 0) && (value < SplCls[i])) {
										SplCls[i + 1] = SplCls[i];
										--i;
									}
									SplCls[i + 1] = value;
								}
								break;
							default:
								quickSort(SplCls, SplInd);
								break;
						}
						
						for (j = 0; j < SplInd; j++) {	/* For each cell C to be split */
							ind0 = SplCls[j];
							ind1 = ind0+cls[ind0];
							SplCntInd = 0;
							
							/* According to the numbers of neighbors of C into the current cell */
							/* compute how many vertices in C will be placed into the same new cell */
							for (i = ind0; i < ind1; i++) {
								value = NghCounts[lab[i]];
								if (Markers[value] != tv->mark) {
									Markers[value] = tv->mark;
									SplCnt[SplCntInd++] = value;
									SplPos[value] = 1;
								}
								else {
									SplPos[value]++;
								}
							}
							tv->mark++;
							
							if (SplCntInd) {
                                SAMETRACE_CHECK(TheTraceSteps, TraceStepsInd, SplCntInd+2*n, Tracestpend)
							}
							
							/* Sort the values deriving from the previous step */
							switch (SplCntInd) {
								case 0:
								case 1:
									break;
								case 2:
									if (SplCnt[0] > SplCnt[1]) {
										value = SplCnt[0];
										SplCnt[0] = SplCnt[1];
										SplCnt[1] = value;
									}
									break;
								case 3:
								case 4:
								case 5:
								case 6:
								case 7:
								case 8:
									for (k = 1; k < SplCntInd; ++k) {
										value = SplCnt[k];
										i = k - 1;
										while ((i >= 0) && (value < SplCnt[i])) {
											SplCnt[i + 1] = SplCnt[i];
											--i;
										}
										SplCnt[i + 1] = value;
									}
									break;
								default:
									quickSort(SplCnt, SplCntInd);
									break;
							}
							
							Part->cells += SplCntInd-1;
							
							/* Split the cell C and update the information for sizes of new cells */
							/* Put the new cells into the stack */
							i = ind0;
							if (StackMarkers[i] != tv->stackmark) {
								BigCellSize = 0;
							}
							for (k = 0; k < SplCntInd; k++) {
								value = SplPos[SplCnt[k]];
								cls[i] = value;
								if ((StackMarkers[ind0] != tv->stackmark) && (value > BigCellSize)) {
									BigCell = i;
									BigCellPos = CStackInd;
									BigCellSize = cls[i];
								}
								SplPos[SplCnt[k]] = i;
								i += value;
								if (i < ind1) {
									CStack[++CStackInd] = i;
									StackMarkers[i] = tv->stackmark;
                                    SAMETRACE_CHECK(TheTrace, TraceInd, i, TraceEnd)
								}
							}
							
							if ((StackMarkers[ind0] != tv->stackmark) && (BigCell != ind0)) {
								CStack[BigCellPos] = ind0;
								StackMarkers[BigCell] = 0;
								StackMarkers[ind0] = tv->stackmark;
							}
							
							/* Permute elements of the cell C */
							i = ind0;
							do {
								SplCnt[SplPos[NghCounts[lab[i]]]++] = lab[i];
							}
							while(++i < ind1);
							
							/* Reconstruct the cell C and update the inverse partition */
							newcell = ind0;
							i = ind0;
							ind2 = newcell+cls[newcell]-1;
							do {
								lab[i] = SplCnt[i];
								InvLab[lab[i]] = i;
								Part->inv[i] = newcell;
								if (i == ind2) {
									newcell = i+1;
									if (newcell < n) ind2 = newcell+cls[newcell]-1;
								}
							}
							while (++i < ind1);
							
							for (i = ind0, k = 0; k < SplCntInd; i+=cls[i], k++) {
								if (cls[i] == 1) {
									Cand->singcode = MASHCOMM(Cand->singcode, Cand->lab[i]);
								}
							}
							
						}
					}
					else {
                        if (TraceCell) {
                            return FALSE;
                        }
					}
					
				}
				else {
					if (cls[ind0] == n) {
						memcpy(NghCounts, sg->d, n*sizeof(int));
					}
					else {
						memset(NghCounts, 0, n*sizeof(int));
						
						for (i = ind0; i < ind2; i++) {
							j1 = sg->v[lab[i]];
							iend1 = j1+sg->d[lab[i]];
							for (; j1 < iend1; ++j1) {
								k = sg->e[j1];
								(NghCounts[k])++;
							}
						}
					}
					SplInd = 0;
					ind4 = 0;
					while (ind4 < n) {	/* For each cell C with size(C) > 1 */
						ind1 = ind4+cls[ind4];
						if (cls[ind4] > 1) {
							
							/* Determine whether C must be split */
							SplCntInd = 0;
							value = NghCounts[lab[ind4]];
							for (i = ind4+1; i < ind1; i++) {
								if (NghCounts[lab[i]] != value)
								{
									SplInd++;
									Markers[value] = tv->mark;
									SplCnt[SplCntInd++] = value;
									SplPos[value] = i-ind4;
									do {
										value = NghCounts[lab[i++]];
										if (Markers[value] != tv->mark) {
											Markers[value] = tv->mark;
											SplCnt[SplCntInd++] = value;
											SplPos[value] = 1;
										}
										else {
											SplPos[value]++;
										}
									}
									while(i != ind1);
									break;
								}
							}
							tv->mark++;
							
							if (SplCntInd) {
                                SAMETRACE_CHECK(TheTraceSteps, TraceStepsInd, SplCntInd+3*n, Tracestpend)
								
								/* Sort the values deriving from the previous step */
								switch (SplCntInd) {
									case 0:
									case 1:
										break;
									case 2:
										if (SplCnt[0] > SplCnt[1]) {
											value = SplCnt[0];
											SplCnt[0] = SplCnt[1];
											SplCnt[1] = value;
										}
										break;
									case 3:
									case 4:
									case 5:
									case 6:
									case 7:
									case 8:
										for (k = 1; k < SplCntInd; ++k) {
											value = SplCnt[k];
											i = k - 1;
											while ((i >= 0) && (value < SplCnt[i])) {
												SplCnt[i + 1] = SplCnt[i];
												--i;
											}
											SplCnt[i + 1] = value;
										}
										break;
									default:
										quickSort(SplCnt, SplCntInd);
										break;
								}
								
								Part->cells += SplCntInd-1;
								
								/* Split the cell C and update the information for sizes of new cells */
								/* Put the new cells into the stack */
								i = ind4;
								if (StackMarkers[i] != tv->stackmark) {
									BigCellSize = 0;
								}
								for (k = 0; k < SplCntInd; k++) {
									value = SplPos[SplCnt[k]];
									cls[i] = value;
									if ((StackMarkers[ind4] != tv->stackmark) && (value > BigCellSize)) {
										BigCell = i;
										BigCellPos = CStackInd;
										BigCellSize = cls[i];
									}
									SplPos[SplCnt[k]] = i;
									i += value;
									if (i < ind1) {
										CStack[++CStackInd] = i;
										StackMarkers[i] = tv->stackmark;
                                        SAMETRACE_CHECK(TheTrace, TraceInd, i, TraceEnd)
									}
								}
								if ((StackMarkers[ind4] != tv->stackmark) && (BigCell != ind4)) {
									CStack[BigCellPos] = ind4;
									StackMarkers[BigCell] = 0;
									StackMarkers[ind4] = tv->stackmark;
								}
								
								/* Permute elements of the cell C */
								i = ind4;
								do {
									SplCnt[SplPos[NghCounts[lab[i]]]++] = lab[i];
								}
								while(++i < ind1);
								
								/* Reconstruct the cell C and update the inverse partition */
								newcell = ind4;
								i = ind4;
								ind2 = newcell+cls[newcell]-1;
								do {
									lab[i] = SplCnt[i];
									InvLab[lab[i]] = i;
									Part->inv[i] = newcell;
									if (i == ind2) {
										newcell = i+1;
										if (newcell < n) ind2 = newcell+cls[newcell]-1;
									}
								}
								while (++i < ind1);
								
								for (i = ind4, k = 0; k < SplCntInd; i+=cls[i], k++) {
									if (cls[i] == 1) {
										Cand->singcode = MASHCOMM(Cand->singcode, Cand->lab[i]);
									}
								}
								
							}
						}
						ind4 = ind1;
					}
					
                    /* DENSE-DENSE CASE */
					if (SplInd) {
                        SAMETRACE_CHECK(TheTraceSplNum, TraceCCInd, SplInd+3*n, &Traceccend)
                    }
					else {
                        if (TraceCell) {
                            return FALSE;
                        }
					}
                    
				}
			}
		}
	}  /* end while (CStackInd > 0) */
    
    for (i=SpineTL->trcstart; i < TraceInd; i++) {
		ind0 = TheTrace[i];
		longcode = MASHNONCOMM(longcode, Part->inv[ind0]);
        j1 = sg->v[lab[ind0]];
		iend1 = j1+sg->d[lab[ind0]];
		for (; j1 < iend1; ++j1) {
			value = Part->inv[InvLab[sg->e[j1]]];
			longcode = MASHCOMM(longcode, value);
        }
	}
    Part->code = Cand->code = CLEANUP(longcode);
    if ((Cand->code != SpineTL->part->code) || (TraceInd != *TraceEnd)) return FALSE;
	return TRUE;
}

void refine_tr(sparsegraph *sg, int *lab, int *ptn, int *numcells, int *code, TracesOptions *options)
{
	int i, j;
	int TraceEnd;
	
	struct TracesVars tvar, *tv;
	struct TracesInfo tinf, *ti;
	
    Partition *CurrPart;
    Candidate *CurrCand;
	
	const int n = sg->nv;
	const int m = SETWORDSNEEDED(n);
	
    if (n > (NAUTY_INFINITY-2))
    {
        fprintf(ERRFILE, "Traces: need n <= %d, but n=%d\n\n",
                NAUTY_INFINITY-2, n);
        return;
    }
    
#if !MAXN
	DYNALLOC1(int, Markers, Markers_sz, n, "refine_tr");
	DYNALLOC1(int, MarkHitVtx, MarkHitVtx_sz, n, "refine_tr");
	DYNALLOC1(int, NghCounts, NghCounts_sz, n, "refine_tr");
	DYNALLOC1(int, SplPos, SplPos_sz, n, "refine_tr");
	DYNALLOC1(int, SplCls, SplCls_sz, n, "refine_tr");
	DYNALLOC1(int, SplCnt, SplCnt_sz, n, "refine_tr");
	DYNALLOC1(int, StackMarkers, StackMarkers_sz, n, "refine_tr");
	DYNALLOC1(int, TheTrace, TheTrace_sz, n+10, "refine_tr");
	DYNALLOC1(int, TheTraceSteps, TheTraceSteps_sz, n+10, "refine_tr");
	DYNALLOC1(int, TheTraceCC, TheTraceCC_sz, n, "refine_tr");
	DYNALLOC1(int, TheTraceSplNum, TheTraceSplNum_sz, n, "refine_tr");
	DYNALLOC1(int, WorkArray2, WorkArray2_sz, n, "refine_tr");
	DYNALLOC1(int, WorkArray3, WorkArray3_sz, n, "refine_tr");
	DYNALLOC1(int, WorkArray4, WorkArray4_sz, n, "refine_tr");
	DYNALLOC1(int, WorkArray5, WorkArray5_sz, n, "refine_tr");
	DYNALLOC1(TracesSpine, Spine, Spine_sz, n, "refine_tr");
#endif
	
#define HitCls WorkArray2
#define HitVtx WorkArray3
#define CStack WorkArray4
#define ElmHitCll WorkArray5
	
	outfile = (options->outfile == NULL ? stdout : options->outfile);
	
	tv = &tvar;
	ti = &tinf;
	/* The graph is sparse? */
	if (sg->nde < n || sg->nde / n < n / (sg->nde / n)) {
		ti->thegraphisparse = TRUE;
	}
	else {
		ti->thegraphisparse = FALSE;
	}
    
	/* Initialize candidate, partition, cells, orbits */
	CurrPart = malloc(sizeof(Partition));
	if (CurrPart == NULL) {
		fprintf(ERRFILE, "\nError, memory not allocated.\n");
		exit(1);
	}
	CurrPart->cls = malloc(n*sizeof(int));
	if (CurrPart->cls == NULL) {
		fprintf(ERRFILE, "\nError, memory not allocated.\n");
		exit(1);
	}
	CurrPart->inv = malloc(n*sizeof(int));
	if (CurrPart->inv == NULL) {
		fprintf(ERRFILE, "\nError, memory not allocated.\n");
		exit(1);
	}
	CurrPart->code = -1;
	
	CurrCand = malloc(sizeof(Candidate));
	if (CurrCand == NULL) {
		fprintf(ERRFILE, "\nError, memory not allocated.\n");
		exit(1);
	}
	CurrCand->lab = malloc(n*sizeof(*CurrCand->lab));
	if (CurrCand->lab == NULL) {
		fprintf(ERRFILE, "\nError, memory not allocated.\n");
		exit(1);
	}
	CurrCand->invlab = malloc(n*sizeof(*CurrCand->invlab));
	if (CurrCand->invlab == NULL) {
		fprintf(ERRFILE, "\nError, memory not allocated.\n");
		exit(1);
	}
	CurrCand->do_it = TRUE;
	CurrCand->indnum = 0;
	CurrCand->code = 0;
	CurrCand->next = NULL;
	
	memset(CurrPart->cls, 0, n*sizeof(int));
	memset(CurrPart->inv, 0, n*sizeof(int));
	memcpy(CurrCand->lab, lab, n*sizeof(int));
	memcpy(CurrCand->invlab, lab, n*sizeof(int));
	CurrPart->cells = 0;
	
	j = 0;
	for (i = 0; i < n; i++) {
		if (j) {
			CurrPart->inv[i] = j;
		}
		CurrCand->invlab[CurrCand->lab[i]] = i;
		if (!ptn[i]) {
			CurrPart->cls[j] = i-j+1;
			TheTrace[CurrPart->cells++] = j;
			j = i+1;
		}
	}
	
	/* Further initializations */
	tv->mark = tv->stackmark = 1073741825;
	tv->maxtreelevel = 0;
	tv->tolevel = 0;
	
	TraceEnd = traces_refine_refine(sg, CurrCand, m, n, CurrPart, tv, ti);
	
	for (i = CurrPart->active; i < TraceEnd; i++) {
		ptn[TheTrace[i]-1] = 0;
	}
	memcpy(lab, CurrCand->lab, n*sizeof(int));
	*code = CurrCand->code;
	*numcells = CurrPart->cells;
    
	FREECAND(CurrCand)
	FREEPART(CurrPart)
    
#if !MAXN
    DYNFREE(Markers, Markers_sz);
    DYNFREE(MarkHitVtx, MarkHitVtx_sz);
    DYNFREE(NghCounts, NghCounts_sz);
    DYNFREE(SplCls, SplCls_sz);
    DYNFREE(SplCnt, SplCnt_sz);
    DYNFREE(SplPos, SplPos_sz);
    DYNFREE(StackMarkers, StackMarkers_sz);
    DYNFREE(TheTrace, TheTrace_sz);
    DYNFREE(TheTraceCC, TheTraceCC_sz);
    DYNFREE(TheTraceSplNum, TheTraceSplNum_sz);
    DYNFREE(TheTraceSteps, TheTraceSteps_sz);
    DYNFREE(WorkArray2, WorkArray2_sz);
    DYNFREE(WorkArray3, WorkArray3_sz);
    DYNFREE(WorkArray4, WorkArray4_sz);
    DYNFREE(WorkArray5, WorkArray5_sz);
    DYNFREE(Spine, Spine_sz);
#endif
}

int traces_refine_refine(sparsegraph *sg,
						 Candidate *Cand,
						 int m,
						 int n,
						 Partition *Part,
						 struct TracesVars* tv,
						 struct TracesInfo *ti)
{
	int i, j, k, sc, ind0, ind1, ind2, ind3, ind4, jk, labi;
	int value, iend, newcell;
	int HitClsInd, SplInd, SplCntInd, CStackInd, TraceInd, TraceCCInd, TraceStepsInd;
	size_t j1, iend1;
	unsigned int longcode;
	int newtrace = FALSE;
	int Sparse = TRUE;
	int *lab, *cls, *InvLab, *TracePos, *SplitCell, *LabCell, TraceEnd, Traceccend, Tracestpend;
	int BigCell, BigCellPos, BigCellSize;
	boolean TraceCell = FALSE;
    int *nghb;
    
	if (tv->stackmark > (NAUTY_INFINITY-2)) {
		memset(StackMarkers, 0, n*sizeof(int));
		tv->stackmark = 0;
	}
	tv->stackmark++;
	
	TraceEnd = Part->active = Part->cells;
	Traceccend = 0;
	Tracestpend = 0;
	TraceCCInd = 0;
	TraceStepsInd = 0;
	
	lab = Cand->lab;
	InvLab = Cand->invlab;
	cls = Part->cls;
	
	UPDATEMIN(Part->active, n-1);
	memcpy(CStack+1, TheTrace, (Part->active)*sizeof(int));
	CStackInd = Part->active;
	for (i = 1; i <= CStackInd; i++) {
		StackMarkers[CStack[i]] = tv->stackmark;
	}
	
	longcode = Part->cells;
	TraceInd = Part->active;
	
	newtrace = TRUE;
	
	while (CStackInd > 0)
	{
		
		if (tv->mark > (NAUTY_INFINITY-2)) {
			memset(Markers, 0, n*sizeof(int));
			memset(MarkHitVtx, 0, n*sizeof(int));
			tv->mark = 0;
		}
		tv->mark++;
		
		if (Part->cells == n) break;
		
		j = CStackInd;
		k = CStackInd;
		while (--j > 0) {
			if (cls[CStack[j]] < cls[CStack[k]]) {
				k = j;
			}
			if ((cls[CStack[k]] == 1) || (j < CStackInd - 12)) {
				break;
			}
		}
		
		ind0 = CStack[k];
		ind2 = ind0+cls[ind0];
		CStack[k] = CStack[CStackInd--];
		StackMarkers[ind0] = 0;
        if (!newtrace) {
			TraceCell = ((ind0 == TheTraceCC[TraceCCInd]) && (TraceCCInd < Traceccend));
        }
		/* Analysis of occurrences of neighbors of the current cell */
		/* The list of cells with neighbors in the current cell is  built */
		if (cls[ind0] == 1) {  /* SINGLETON CURRENT CELL CASE */
            HitClsInd = 0;
			j1 = sg->v[lab[ind0]];
			iend1 = j1+sg->d[lab[ind0]];
			
			for (; j1 < iend1; ++j1) {
				k = sg->e[j1];
				value = Part->inv[InvLab[k]];
				if (cls[value] > 1) {
					if (Markers[value] != tv->mark) {
						HitCls[HitClsInd++] = value;
						Markers[value] = tv->mark;
						ElmHitCll[value] = value;
					}
					HitVtx[ElmHitCll[value]++] = k;
				}
			}
			tv->mark++;
			
			SplInd = 0;
			for (j = 0; j < HitClsInd; j++) {
				ind1 = HitCls[j];
				ElmHitCll[ind1] -= ind1;
				if ((ElmHitCll[ind1] > 0) && (ElmHitCll[ind1] < cls[ind1])) {
					SplCls[SplInd++] = ind1;
				}
			}
			
			if (SplInd) {
				if (newtrace) {
					TheTraceCC[TraceCCInd] = ind0;
                }
				TRACE_CHECK(TheTraceSplNum, TraceCCInd, SplInd, &Traceccend)
				
				switch (SplInd) {
					case 0:
					case 1:
						break;
					case 2:
						if (SplCls[0] > SplCls[1]) {
							value = SplCls[0];
							SplCls[0] = SplCls[1];
							SplCls[1] = value;
						}
						break;
					case 3:
					case 4:
					case 5:
					case 6:
					case 7:
					case 8:
						for (k = 1; k < SplInd; ++k) {
							value = SplCls[k];
							i = k - 1;
							while ((i >= 0) && (value < SplCls[i])) {
								SplCls[i + 1] = SplCls[i];
								--i;
							}
							SplCls[i + 1] = value;
						}
						break;
					default:
						quickSort(SplCls, SplInd);
						break;
				}
				
				for (j = 0; j < SplInd; j++) {
					ind1 = SplCls[j];
					i = ind1+cls[ind1]-ElmHitCll[ind1];
					TRACE_CHECK(TheTrace, TraceInd, i, &TraceEnd)
				}
				
				/* REARRANGE THE CELLS */
				for (j = 0; j < SplInd; j++) {
					ind1 = SplCls[j];
					cls[ind1] = cls[ind1]-ElmHitCll[ind1];
					newcell = ind1+cls[ind1];
					cls[newcell] = ElmHitCll[ind1];
					Part->cells++;
					
					if (StackMarkers[ind1] != tv->stackmark) {
						if (cls[newcell] < cls[ind1]) {
							CStack[++CStackInd] = newcell;
							StackMarkers[newcell] = tv->stackmark;
						}
						else {
							CStack[++CStackInd] = ind1;
							StackMarkers[ind1] = tv->stackmark;
						}
					}
					else {
						CStack[++CStackInd] = newcell;
						StackMarkers[newcell] = tv->stackmark;
					}
					
					SplitCell = HitVtx+ind1;
					ind3 = cls[newcell];
					LabCell = lab+newcell;
					
					for (jk = 0; jk < ind3; jk++) {
						k = SplitCell[jk];
						i = LabCell[jk];
						Part->inv[newcell+jk] = newcell;
						lab[InvLab[k]] = i;
						InvLab[i] = InvLab[k];
						LabCell[jk] = k;
						InvLab[k] = newcell+jk;
					}
				}
			}
			else {
				if ((!newtrace) && TraceCell) {
					return 0;
				}
			}
			
		}
		else {
			if (ti->thegraphisparse) {
                if (cls[ind0] == n) {
					memcpy(NghCounts, sg->d, n*sizeof(int));
					HitCls[0] = 0;
					HitClsInd = 1;
					ElmHitCll[0] = 0;
					for (i = 0; i < n; i++) {
						if (sg->d[i]) {
							HitVtx[ElmHitCll[0]++] = i;
						}
					}
				}
                else {
					HitClsInd = 0;
                    for (i = ind0; i < ind2; i++) {
                        labi = lab[i];
                        nghb = sg->e+sg->v[labi];
                        j = sg->d[labi];
                        for (; j--;) {
                            k = nghb[j];
                            if (MarkHitVtx[k] == tv->mark) {
                                NghCounts[k]++;
                            }
                            else {
                                value = Part->inv[InvLab[k]];
                                if (cls[value] > 1) {
                                    MarkHitVtx[k] = tv->mark;
                                    NghCounts[k] = 1;
                                    if (Markers[value] != tv->mark) {
                                        HitCls[HitClsInd++] = value;
                                        Markers[value] = tv->mark;
                                        HitVtx[value] = k;
                                        ElmHitCll[value] = 1;
                                    }
                                    else {
                                        HitVtx[value+ElmHitCll[value]++] = k;
                                    }
                                }
                            }
						}
					}
                }
                
				tv->mark++;
				
				SplInd = 0;
				SplCls[0] = n;
				for (j = 0; j < HitClsInd; j++) {
					ind1 = HitCls[j];
					if ((ElmHitCll[ind1] > 0) && (ElmHitCll[ind1] < cls[ind1])) {
						SplCls[SplInd++] = HitCls[j];
					}
					else {
						ind2 = ind1+cls[ind1];
						value = NghCounts[lab[ind1++]];
						for (i = ind1; i < ind2; i++)
						{
							if (NghCounts[lab[i]] != value)
							{
								SplCls[SplInd++] = HitCls[j];
								break;
							}
						}
					}
				}
				
				if (SplInd) {
					if (newtrace) {
                        TheTraceCC[TraceCCInd] = ind0;
                    }
					TRACE_CHECK(TheTraceSplNum, TraceCCInd, SplInd, &Traceccend)
					
					/* Sorting the cells to be split */
					switch (SplInd) {
						case 0:
						case 1:
							break;
						case 2:
							if (SplCls[0] > SplCls[1]) {
								value = SplCls[0];
								SplCls[0] = SplCls[1];
								SplCls[1] = value;
							}
							break;
						case 3:
						case 4:
						case 5:
						case 6:
						case 7:
						case 8:
							for (k = 1; k < SplInd; ++k) {
								value = SplCls[k];
								i = k - 1;
								while ((i >= 0) && (value < SplCls[i])) {
									SplCls[i + 1] = SplCls[i];
									--i;
								}
								SplCls[i + 1] = value;
							}
							break;
						default:
							quickSort(SplCls, SplInd);
							break;
					}
					
					for (sc = 0; sc < SplInd; sc++) {	/* For each cell C to be split */
						ind0 = SplCls[sc];
						ind1 = ind0 + cls[ind0];
						SplCntInd = 0;
						if (ElmHitCll[ind0] < cls[ind0]) {
							SplCnt[SplCntInd++] = 0;
							SplPos[0] = cls[ind0] - ElmHitCll[ind0];
						}
						
						/* According to the numbers of neighbors of C into the current cell */
						/* compute how many vertices in C will be placed into the same new cell */
						iend = ind0 + ElmHitCll[ind0];
						for (i = ind0; i < iend; i++) {
							value = NghCounts[HitVtx[i]];
							if (Markers[value] != tv->mark) {
								Markers[value] = tv->mark;
								SplCnt[SplCntInd++] = value;
								SplPos[value] = 1;
							}
							else {
								SplPos[value]++;
							}
						}
						tv->mark++;
						
						if (SplCntInd) {
							TRACE_CHECK(TheTraceSteps, TraceStepsInd, SplCntInd, &Tracestpend)
						}
						
						/* Sort the values deriving from the previous step */
						switch (SplCntInd) {
							case 0:
							case 1:
								break;
							case 2:
								if (SplCnt[0] > SplCnt[1]) {
									value = SplCnt[0];
									SplCnt[0] = SplCnt[1];
									SplCnt[1] = value;
								}
								break;
							case 3:
							case 4:
							case 5:
							case 6:
							case 7:
							case 8:
								for (k = 1; k < SplCntInd; ++k) {
									value = SplCnt[k];
									i = k - 1;
									while ((i >= 0) && (value < SplCnt[i])) {
										SplCnt[i + 1] = SplCnt[i];
										--i;
									}
									SplCnt[i + 1] = value;
								}
								break;
							default:
								quickSort(SplCnt, SplCntInd);
								break;
						}
						
						Part->cells += SplCntInd-1;
						
						/* Split the cell C and update the information for sizes of new cells */
						/* Put the new cells into the stack */
						i = ind0;
						if (StackMarkers[i] != tv->stackmark) {
							BigCellSize = 0;
						}
						
						for (k = 0; k < SplCntInd; k++) {
							value = SplPos[SplCnt[k]];
							cls[i] = value;
							if ((StackMarkers[ind0] != tv->stackmark) && (value > BigCellSize)) {
								BigCell = i;
								BigCellPos = CStackInd;
								BigCellSize = cls[i];
							}
							SplPos[SplCnt[k]] = i;
							i += value;
							if (i < ind1) {
								CStack[++CStackInd] = i;
								StackMarkers[i] = tv->stackmark;
								TRACE_CHECK(TheTrace, TraceInd, i, &TraceEnd)
							}
						}
						
						if ((StackMarkers[ind0] != tv->stackmark) && (BigCell != ind0)) {
							CStack[BigCellPos] = ind0;
							StackMarkers[BigCell] = 0;
							StackMarkers[ind0] = tv->stackmark;
						}
						/* Permute elements of the cell C */
						iend = ind0 + ElmHitCll[ind0];
						for (i = ind0; i < iend; i++) {
							value = HitVtx[i];
							j = SplPos[NghCounts[value]]++; /* where HitVtx[i] goes */
							k = InvLab[value];				/* where HitVtx[i] is in lab */
							lab[k] = lab[j];
							lab[j] = value;
							InvLab[value] = j;
							InvLab[lab[k]] = k;
							NghCounts[value] = 0;
						}
						
						/* Reconstruct the cell C and update the inverse partition */
						newcell = ind1 - ElmHitCll[ind0];
						i = newcell;
						ind2 = newcell+cls[newcell]-1;
						do {
							Part->inv[i] = newcell;
							if (i == ind2) {
								newcell = i+1;
								if (newcell < n) ind2 = newcell+cls[newcell]-1;
							}
						}
						while (++i < ind1);
					}
				}
				else {
					if ((!newtrace) && TraceCell) {
						return 0;
					}
				}
				
			}
			else {
				if (sg->d[lab[ind0]] > n/cls[ind0]) {
					Sparse = FALSE;
                }
				else {
					Sparse = TRUE;
                }
				if (Sparse) {
					/* Counting occurrences of neighbors of the current cell */
					/* The list of cells with neighbors in the current cell is also built */
					if (cls[ind0] == n) {
						memcpy(NghCounts, sg->d, n*sizeof(int));
						HitCls[0] = 0;
						HitClsInd = 1;
					}
					else {
						memset(NghCounts, 0, n*sizeof(int));
						HitClsInd = 0;
						for (i = ind0; i < ind2; i++) {
							j1 = sg->v[lab[i]];
							iend1 = j1+sg->d[lab[i]];
							for (; j1 < iend1; ++j1) {
								k = sg->e[j1];
								(NghCounts[k])++;
								value = Part->inv[InvLab[k]];
								if (Markers[value] != tv->mark) {
									if (cls[value] > 1) HitCls[HitClsInd++] = value;
									Markers[value] = tv->mark;
								}
							}
						}
					}
					
					tv->mark++;
					
					SplInd = 0;
					for (j = 0; j < HitClsInd; j++) {
						ind1 = HitCls[j];
						ind2 = ind1+cls[ind1];
						value = NghCounts[lab[ind1++]];
						for (i = ind1; i < ind2; i++)
						{
							if (NghCounts[lab[i]] != value)
							{
								SplCls[SplInd++] = HitCls[j];
								break;
							}
						}
					}
					
					if (SplInd) {
						if (newtrace) {
                            TheTraceCC[TraceCCInd] = ind0;
                        }
						TRACE_CHECK(TheTraceSplNum, TraceCCInd, SplInd, &Traceccend)
						
						/* Sorting the cells to be split */
						switch (SplInd) {
							case 0:
							case 1:
								break;
							case 2:
								if (SplCls[0] > SplCls[1]) {
									value = SplCls[0];
									SplCls[0] = SplCls[1];
									SplCls[1] = value;
								}
								break;
							case 3:
							case 4:
							case 5:
							case 6:
							case 7:
							case 8:
								for (k = 1; k < SplInd; ++k) {
									value = SplCls[k];
									i = k - 1;
									while ((i >= 0) && (value < SplCls[i])) {
										SplCls[i + 1] = SplCls[i];
										--i;
									}
									SplCls[i + 1] = value;
								}
								break;
							default:
								quickSort(SplCls, SplInd);
								break;
						}
						
						for (j = 0; j < SplInd; j++) {	/* For each cell C to be split */
							ind0 = SplCls[j];
							ind1 = ind0+cls[ind0];
							SplCntInd = 0;
							
							/* According to the numbers of neighbors of C into the current cell */
							/* compute how many vertices in C will be placed into the same new cell */
							for (i = ind0; i < ind1; i++) {
								value = NghCounts[lab[i]];
								if (Markers[value] != tv->mark) {
									Markers[value] = tv->mark;
									SplCnt[SplCntInd++] = value;
									SplPos[value] = 1;
								}
								else {
									SplPos[value]++;
								}
							}
							tv->mark++;
							
							if (SplCntInd) {
								TRACE_CHECK(TheTraceSteps, TraceStepsInd, SplCntInd, &Tracestpend)
							}
							
							/* Sort the values deriving from the previous step */
							switch (SplCntInd) {
								case 0:
								case 1:
									break;
								case 2:
									if (SplCnt[0] > SplCnt[1]) {
										value = SplCnt[0];
										SplCnt[0] = SplCnt[1];
										SplCnt[1] = value;
									}
									break;
								case 3:
								case 4:
								case 5:
								case 6:
								case 7:
								case 8:
									for (k = 1; k < SplCntInd; ++k) {
										value = SplCnt[k];
										i = k - 1;
										while ((i >= 0) && (value < SplCnt[i])) {
											SplCnt[i + 1] = SplCnt[i];
											--i;
										}
										SplCnt[i + 1] = value;
									}
									break;
								default:
									quickSort(SplCnt, SplCntInd);
									break;
							}
							
							Part->cells += SplCntInd-1;
							
							/* Split the cell C and update the information for sizes of new cells */
							/* Put the new cells into the stack */
							i = ind0;
							if (StackMarkers[i] != tv->stackmark) {
								BigCellSize = 0;
							}
							for (k = 0; k < SplCntInd; k++) {
								value = SplPos[SplCnt[k]];
								cls[i] = value;
								if ((StackMarkers[ind0] != tv->stackmark) && (value > BigCellSize)) {
									BigCell = i;
									BigCellPos = CStackInd;
									BigCellSize = cls[i];
								}
								SplPos[SplCnt[k]] = i;
								i += value;
								if (i < ind1) {
									CStack[++CStackInd] = i;
									StackMarkers[i] = tv->stackmark;
									TRACE_CHECK(TheTrace, TraceInd, i, &TraceEnd)
								}
							}
							
							if ((StackMarkers[ind0] != tv->stackmark) && (BigCell != ind0)) {
								CStack[BigCellPos] = ind0;
								StackMarkers[BigCell] = 0;
								StackMarkers[ind0] = tv->stackmark;
							}
							
							/* Permute elements of the cell C */
							i = ind0;
							do {
								SplCnt[SplPos[NghCounts[lab[i]]]++] = lab[i];
							}
							while(++i < ind1);
							
							/* Reconstruct the cell C and update the inverse partition */
							newcell = ind0;
							i = ind0;
							ind2 = newcell+cls[newcell]-1;
							do {
								lab[i] = SplCnt[i];
								InvLab[lab[i]] = i;
								Part->inv[i] = newcell;
								if (i == ind2) {
									newcell = i+1;
									if (newcell < n) ind2 = newcell+cls[newcell]-1;
								}
							}
							while (++i < ind1);
						}
					}
					else {
						if ((!newtrace) && TraceCell) {
							return 0;
						}
					}
					
				}
				else {
					if (cls[ind0] == n) {
						memcpy(NghCounts, sg->d, n*sizeof(int));
					}
					else {
						memset(NghCounts, 0, n*sizeof(int));
						
						for (i = ind0; i < ind2; i++) {
							j1 = sg->v[lab[i]];
							iend1 = j1+sg->d[lab[i]];
							for (; j1 < iend1; ++j1) {
								k = sg->e[j1];
								(NghCounts[k])++;
							}
						}
					}
					SplInd = 0;
					ind4 = 0;
					while (ind4 < n) {	/* For each cell C with size(C) > 1 */
						ind1 = ind4+cls[ind4];
						if (cls[ind4] > 1) {
							
							/* Determine whether C must be split */
							SplCntInd = 0;
							value = NghCounts[lab[ind4]];
							for (i = ind4+1; i < ind1; i++) {
								if (NghCounts[lab[i]] != value)
								{
									SplInd++;
									Markers[value] = tv->mark;
									SplCnt[SplCntInd++] = value;
									SplPos[value] = i-ind4;
									do {
										value = NghCounts[lab[i++]];
										if (Markers[value] != tv->mark) {
											Markers[value] = tv->mark;
											SplCnt[SplCntInd++] = value;
											SplPos[value] = 1;
										}
										else {
											SplPos[value]++;
										}
									}
									while(i != ind1);
									break;
								}
							}
							tv->mark++;
							
							if (SplCntInd) {
								TRACE_CHECK(TheTraceSteps, TraceStepsInd, SplCntInd, &Tracestpend)
								
								/* Sort the values deriving from the previous step */
								switch (SplCntInd) {
									case 0:
									case 1:
										break;
									case 2:
										if (SplCnt[0] > SplCnt[1]) {
											value = SplCnt[0];
											SplCnt[0] = SplCnt[1];
											SplCnt[1] = value;
										}
										break;
									case 3:
									case 4:
									case 5:
									case 6:
									case 7:
									case 8:
										for (k = 1; k < SplCntInd; ++k) {
											value = SplCnt[k];
											i = k - 1;
											while ((i >= 0) && (value < SplCnt[i])) {
												SplCnt[i + 1] = SplCnt[i];
												--i;
											}
											SplCnt[i + 1] = value;
										}
										break;
									default:
										quickSort(SplCnt, SplCntInd);
										break;
								}
								
								Part->cells += SplCntInd-1;
								
								/* Split the cell C and update the information for sizes of new cells */
								/* Put the new cells into the stack */
								i = ind4;
								if (StackMarkers[i] != tv->stackmark) {
									BigCellSize = 0;
								}
								for (k = 0; k < SplCntInd; k++) {
									value = SplPos[SplCnt[k]];
									cls[i] = value;
									if ((StackMarkers[ind4] != tv->stackmark) && (value > BigCellSize)) {
										BigCell = i;
										BigCellPos = CStackInd;
										BigCellSize = cls[i];
									}
									SplPos[SplCnt[k]] = i;
									i += value;
									if (i < ind1) {
										CStack[++CStackInd] = i;
										StackMarkers[i] = tv->stackmark;
										TRACE_CHECK(TheTrace, TraceInd, i, &TraceEnd)
									}
								}
								if ((StackMarkers[ind4] != tv->stackmark) && (BigCell != ind4)) {
									CStack[BigCellPos] = ind4;
									StackMarkers[BigCell] = 0;
									StackMarkers[ind4] = tv->stackmark;
								}
								
								/* Permute elements of the cell C */
								i = ind4;
								do {
									SplCnt[SplPos[NghCounts[lab[i]]]++] = lab[i];
								}
								while(++i < ind1);
								
								/* Reconstruct the cell C and update the inverse partition */
								newcell = ind4;
								i = ind4;
								ind2 = newcell+cls[newcell]-1;
								do {
									lab[i] = SplCnt[i];
									InvLab[lab[i]] = i;
									Part->inv[i] = newcell;
									if (i == ind2) {
										newcell = i+1;
										if (newcell < n) ind2 = newcell+cls[newcell]-1;
									}
								}
								while (++i < ind1);
							}
						}
						ind4 = ind1;
					}
					
					if (SplInd) {
						if (newtrace) {
                            TheTraceCC[TraceCCInd] = ind0;
                        }
						TheTraceSplNum[TraceCCInd] = SplInd;
						TraceCCInd++;
					}
					else {
						if ((!newtrace) && TraceCell) {
							return 0;
						}
					}
					
				}
			}
		}
	} /* end while (CStackInd > 0) */
	
	for (i=0; i < TraceInd; i++) {
		ind0 = TheTrace[i];
		longcode = MASHNONCOMM(longcode, Part->inv[ind0]);
		j1 = sg->v[lab[ind0]];
		iend1 = j1+sg->d[lab[ind0]];
		for (; j1 < iend1; ++j1) {
			value = Part->inv[InvLab[sg->e[j1]]];
			longcode = MASHCOMM(longcode, value);
		}
	}
	Part->code = Cand->code = CLEANUP(longcode);
	return TraceInd;
}

void CopyVect(int *from, int *to, int u, int n)
{
	int i;
	for (i = u; i<n; i++)
		to[i] = from[i];
	return;
}

void quickSort(int *arr, int elements) {
	
#define MAX_LEVELS 300
	
	int piv, beg[MAX_LEVELS], end[MAX_LEVELS], i = 0, L, R, swap;
	int k, value;
	
	beg[0] = 0;
	end[0] = elements;
	while (i>= 0) {
		L = beg[i];
		R = end[i]-1;
		if (L<R-8) {
			piv = arr[(L+R)/2];
			arr[(L+R)/2] = arr[L];
			arr[L] = piv;
			while (L<R) {
				while (arr[R]>= piv && L<R) R--;
				if (L<R) arr[L++] = arr[R];
				while (arr[L]<= piv && L<R) L++;
				if (L<R) arr[R--] = arr[L];
			}
			arr[L] = piv;
			beg[i+1] = L+1;
			end[i+1] = end[i]; end[i++] = L;
			if (end[i]-beg[i]>end[i-1]-beg[i-1]) {
				swap = beg[i];
				beg[i] = beg[i-1];
				beg[i-1] = swap;
				swap = end[i];
				end[i] = end[i-1];
				end[i-1] = swap;
			}
		}
		else {
			i--;
		}
	}
	for (k = 1; k < elements; ++k) {
		value = arr[k];
		i = k - 1;
		while ((i >= 0) && (value < arr[i])) {
			arr[i + 1] = arr[i];
			--i;
		}
		arr[i + 1] = value;
	}
}

struct Candidate *NewCandidate(int n, Candidate **GarbList, int Mrk)
{
	struct Candidate *Cand;
	
	if (*GarbList) {
		Cand = *GarbList;
		*GarbList = (*GarbList)->next;
	}
	else {
		Cand = malloc(sizeof(*Cand));
		if (Cand == NULL) {
			fprintf(ERRFILE, "\nError, memory not allocated.\n");
			exit(1);
		}
		Cand->lab = malloc(n*sizeof(*Cand->lab));
		if (Cand->lab == NULL) {
			fprintf(ERRFILE, "\nError, memory not allocated.\n");
			exit(1);
		}
		Cand->invlab = malloc(n*sizeof(*Cand->invlab));
		if (Cand->invlab == NULL) {
			fprintf(ERRFILE, "\nError, memory not allocated.\n");
			exit(1);
		}
	}
	Cand->do_it = Mrk;
	Cand->indnum = 0;
	Cand->code = 0;
	Cand->next = NULL;
	Cand->stnode = NULL;
	Cand->sortedlab = FALSE;
	return Cand;
}

int FreeList(Candidate *List, int cond)
{
	Candidate *Temp;
	int conta = 0;
	int conta1 = 0;
	while (List) {
		if (List->do_it == cond) {
			conta1++;
		}
		conta++;
		Temp = List;
		if (List->lab) free(List->lab);
		if (List->invlab) free(List->invlab);
		List = List->next;
		free(Temp);
	}
	
	if (cond) {
		return conta1;
	}
	else {
		return conta;
	}
}

int FixBase(int *fix, struct TracesVars *tv, Candidate *Cand, int from, int to)
{
	int i, j, k, go, nfix;
    nfix = j = 0;
	go = TRUE;
	for (i = from; i < to; i++) {
		k = Cand->lab[Spine[i+1].tgtpos];
		if (go && (nfix < tv->nfix) && (fix[nfix] == k)) {
			j++;
		}
		else {
			fix[nfix] = k;
			if (go) go = FALSE;
		}
		nfix++;
	}
	tv->nfix = nfix;
    return j;
}

void
traces_freedyn(void)
{
	/* Free the static dynamic memory used by Traces */
#if !MAXN
    DYNFREE(AUTPERM, AUTPERM_sz);
    DYNFREE(BreakSteps, BreakSteps_sz);
    DYNFREE(CurrOrbSize, CurrOrbSize_sz);
    DYNFREE(CurrRefCells, CurrRefCells_sz);
    DYNFREE(fix, fix_sz);
    DYNFREE(IDENTITY_PERM, IDENTITY_PERM_sz);
    DYNFREE(Markers, Markers_sz);
    DYNFREE(MarkHitVtx, MarkHitVtx_sz);
    DYNFREE(MultRefCells, MultRefCells_sz);
    DYNFREE(NghCounts, NghCounts_sz);
    DYNFREE(OrbSize, OrbSize_sz);
    DYNFREE(RefCells, RefCells_sz);
    DYNFREE(RefPath, RefPath_sz);
    DYNFREE(Spine, Spine_sz);
    DYNFREE(SplCls, SplCls_sz);
    DYNFREE(SplCnt, SplCnt_sz);
    DYNFREE(SplPos, SplPos_sz);
    DYNFREE(StackMarkers, StackMarkers_sz);
    DYNFREE(TheTrace, TheTrace_sz);
    DYNFREE(TheTraceCC, TheTraceCC_sz);
    DYNFREE(TheTraceSplNum, TheTraceSplNum_sz);
    DYNFREE(TheTraceSteps, TheTraceSteps_sz);
    DYNFREE(TreeStack, TreeStack_sz);
    DYNFREE(TrieArray, TrieArray_sz);
	DYNFREE(TEMPLAB, TEMPLAB_sz);
	DYNFREE(TEMPINVLAB, TEMPINVLAB_sz);
	DYNFREE(WorkArray, WorkArray_sz);
    DYNFREE(WorkArray1, WorkArray1_sz);
    DYNFREE(WorkArray2, WorkArray2_sz);
    DYNFREE(WorkArray3, WorkArray3_sz);
    DYNFREE(WorkArray4, WorkArray4_sz);
    DYNFREE(WorkArray5, WorkArray5_sz);
#endif
}

void factorial(double *size1, int *size2, int k)
{
    int i;
    for(i = k; i; i--) {
        MULTIPLY(*size1, *size2, i);
    }
}

void factorial2(double *size1, int *size2, int k)
{
    int i;
    for(i = k; i; i -= 2) {
        MULTIPLY(*size1, *size2, i);
    }
}

int CheckForAutomorphisms(sparsegraph *g_arg, Candidate *CurrCand, Candidate *NextCand,
                          struct TracesVars* tv, struct TracesInfo* ti,
                          int m, int n, Partition* Part)
{
	Candidate *CheckAutList;
	int i, j, tgt_level;
	int CheckLevel, CheckLevelEnd;
	int temp, tmp, tmp1;
	int start, start1, step, step1, prev, prev1, ngh1, ngh2;
    searchtrie *TrieCandFrom, *TrieCheckFrom;
    
    tv->gotonode = NULL;
    CheckLevel = 0;
    CheckLevelEnd = 0;
    temp = 0;
    
	switch (tv->compstage) {
		case 0:
            if (tv->strategy) {
                CheckLevel = CheckLevelEnd = tv->maxtreelevel;
            }
            else {
                if ((Spine[tv->tolevel].part)->cells == tv->finalnumcells) {
                    CheckLevel = CheckLevelEnd = tv->maxtreelevel;
                }
                else {
                    CheckLevel = 1;
                    if ((Spine[tv->maxtreelevel].part)->cells == tv->finalnumcells) {
                        CheckLevelEnd = tv->maxtreelevel - 1;
                    }
                    else {
                        CheckLevelEnd = tv->maxtreelevel;
                    }
                }
            }
			break;
		case 1:
			CheckLevel = CheckLevelEnd = tv->tolevel;
			break;
		case 2:
			if (m || (tv->tolevel == tv->maxtreelevel+1)) {
				CheckLevel = CheckLevelEnd = tv->maxtreelevel+1;
			}
			else {
				CheckLevel = 1;
				if ((Spine[tv->maxtreelevel].part)->cells == tv->finalnumcells) {
					CheckLevelEnd = tv->maxtreelevel - 1;
				}
				else {
					CheckLevelEnd = tv->maxtreelevel;
				}
			}
			break;
		default:
			break;
	}
    
    while (CheckLevel <= CheckLevelEnd) {
		CheckAutList = Spine[CheckLevel].liststart;
		while (CheckAutList) {
			if (CheckAutList->do_it && lookup(CheckAutList->stnode) && (CheckAutList != NextCand) && (CheckAutList != CurrCand)) {
                if (CheckAutList->code == NextCand->code) {
                    if (Part->cells == n) {
                        for (i = 0; i < n; i++) {
                            AUTPERM[NextCand->lab[i]] = CheckAutList->lab[i];
                        }
                    }
                    else {
                        memcpy(AUTPERM, IDENTITY_PERM, n*sizeof(int));
                        SETMARK(Markers, tv->mark)
                        for (i=0; i<n; i+=Part->cls[i]) {
                            if (Part->cls[i] == 1) {
                                AUTPERM[NextCand->lab[i]] = CheckAutList->lab[i];
                            }
                            else {
                                for (j=i; j<i+Part->cls[i]; j++) {
                                    tmp = CheckAutList->lab[j];
                                    if (g_arg->d[tmp] == 2) {
                                        ngh1 = g_arg->e[g_arg->v[tmp]];
                                        ngh2 = g_arg->e[g_arg->v[tmp]+1];
                                        if ((Markers[tmp] != tv->mark) && ((g_arg->d[ngh1] > 2) || ((g_arg->d[ngh2] > 2)))) {
                                            tmp1 = NextCand->lab[j];
                                            if ((g_arg->d[ngh1] > 2) && ((g_arg->d[ngh2] > 2))) {
                                                if (Part->inv[CheckAutList->invlab[ngh1]] < Part->inv[CheckAutList->invlab[ngh2]]) {
                                                    start = ngh1;
                                                    start1 = NextCand->lab[CheckAutList->invlab[ngh1]];
                                                } else {
                                                    start = ngh2;
                                                    start1 = NextCand->lab[CheckAutList->invlab[ngh2]];
                                                }
                                            }
                                            else {
                                                if (g_arg->d[ngh1] > 2) {
                                                    start = ngh1;
                                                    start1 = NextCand->lab[CheckAutList->invlab[ngh1]];
                                                }
                                                else {
                                                    start = ngh2;
                                                    start1 = NextCand->lab[CheckAutList->invlab[ngh2]];
                                                }
                                            }
                                            step = tmp;
                                            step1 = tmp1;
                                            if (Markers[step] != tv->mark) {
                                                prev = start;
                                                prev1 = start1;
                                                AUTPERM[start1] = start;
                                                do {
                                                    Markers[step] = tv->mark;
                                                    AUTPERM[step1] = step;
                                                    if (g_arg->e[g_arg->v[step]] != prev) {
                                                        prev = step;
                                                        step = g_arg->e[g_arg->v[step]];
                                                    } else {
                                                        prev = step;
                                                        step = g_arg->e[g_arg->v[step]+1];
                                                    }
                                                    if (g_arg->e[g_arg->v[step1]] != prev1) {
                                                        prev1 = step1;
                                                        step1 = g_arg->e[g_arg->v[step1]];
                                                    } else {
                                                        prev1 = step1;
                                                        step1 = g_arg->e[g_arg->v[step1]+1];
                                                    }
                                                } while (g_arg->d[step] == 2);
                                                AUTPERM[step1] = step;
                                            }
                                        }
                                    }
                                    else {
                                        if ((g_arg->d[tmp] == 1) && (g_arg->d[g_arg->e[g_arg->v[tmp]]] > 2)) {
                                            AUTPERM[NextCand->lab[j]] = CheckAutList->lab[j];
                                        }
                                    }
                                }
                            }
                        }
                    }
					if (isautom_sg((graph*)g_arg, AUTPERM, tv->options->digraph, m, n)) {
                        if (!findperm(gensB, AUTPERM, n, tv)) {
                            if (tv->options->verbosity >= 2) tv->schreier3 -= CPUTIME;
                            addgenerator(&gpB, &gensB, AUTPERM, n);
                            if (tv->options->verbosity >= 2) tv->schreier3 += CPUTIME;
                            if (tv->options->verbosity >= 2) {
								fprintf(outfile, "[A (%d, %d)] ", CheckLevel, CheckAutList->name);
							}
							tv->stats->numgenerators++;
							tv->stats->numorbits = orbjoin(tv->orbits, AUTPERM, n);
                            
							ti->thegrouphaschanged = TRUE;
							ti->identitygroup = FALSE;
							if (tv->options->verbosity >= 2 && tv->options->writeautoms) {
								PRINT_RETURN
							}
							if (tv->options->writeautoms) {
								fprintf(outfile, "Gen #%d: ", tv->stats->numgenerators);
								writeperm(outfile, AUTPERM, tv->options->cartesian, tv->options->linelength, n);
							}
							if (tv->options->userautomproc) {
								(*tv->options->userautomproc)(tv->stats->numgenerators, AUTPERM, n);
							}
						}
						else {
							if (tv->options->verbosity >= 2) {
								fprintf(outfile, "[A* (%d, %d)] ", CheckLevel, CheckAutList->name);
							}
						}
                        TrieCandFrom = NULL;
                        TrieCheckFrom = CheckAutList->stnode;
                        if (CurrCand->stnode->level <= 1) {
                            tgt_level = CurrCand->stnode->level + 1;
                            while (TrieCheckFrom->level > tgt_level) {
                                TrieCheckFrom = TrieCheckFrom->father;
                            }
                        }
                        else {
                            if (tv->tolevel <= TrieCheckFrom->level) {
                                tgt_level = tv->tolevel;
                                while (TrieCheckFrom->level != tgt_level) {
                                    TrieCheckFrom = TrieCheckFrom->father;
                                }
                            } else {
                                TrieCandFrom = CurrCand->stnode;
                                tgt_level = TrieCheckFrom->level;
                                while (TrieCandFrom->level != tgt_level) {
                                    TrieCandFrom = TrieCandFrom->father;
                                }
                            }
                        }
                        if (TrieCandFrom) {
                            while (TrieCandFrom->father != TrieCheckFrom->father) {
                                TrieCandFrom = TrieCandFrom->father;
                                TrieCheckFrom = TrieCheckFrom->father;
                            }
                        }
                        else {
                            if ((TrieCheckFrom->level > 1) && (TrieCheckFrom->father != CurrCand->stnode)) {
                                TrieCandFrom = CurrCand->stnode;
                                TrieCheckFrom = TrieCheckFrom->father;
                                while (TrieCandFrom->father != TrieCheckFrom->father) {
                                    TrieCandFrom = TrieCandFrom->father;
                                    TrieCheckFrom = TrieCheckFrom->father;
                                }
                            }
                        }
                        
                        while (TrieCheckFrom->goes_to) {
                            TrieCheckFrom = TrieCheckFrom->goes_to;
                        }
                        
						for (temp=1; temp<=tv->tolevel; temp++) {
							if (CheckAutList->lab[Spine[temp].tgtpos] != NextCand->lab[Spine[temp].tgtpos]) {
								break;
							}
						}
                        if ((temp == tv->tolevel) && TempOrbits) {
                            orbjoin(TempOrbits, AUTPERM, n);
                        }
                        switch (tv->compstage) {
                            case 0:
                                if (tv->strategy && (tv->steps == 1)) {
                                    RemoveFromLevel(temp, tv->maxtreelevel-1, tv->strategy, FALSE);
                                    if (TrieCandFrom) {
                                        TrieCheckFrom->index += TrieCandFrom->index;
                                        TrieCandFrom->goes_to = TrieCheckFrom;
                                    }
                                    else TrieCheckFrom->index++;
                                    NextCand->do_it = FALSE;
                                }
                                else {
                                    if (CheckAutList->lab[Spine[temp].tgtpos] >= NextCand->lab[Spine[temp].tgtpos]) {
                                        CheckAutList->do_it = FALSE;
                                        if (TrieCandFrom) {
                                            TrieCandFrom->index += TrieCheckFrom->index;
                                            tv->newindex = 0;
                                            TrieCheckFrom->goes_to = TrieCandFrom;
                                        }
                                        else {
                                            if (CurrCand->stnode->level > 1) {
                                                tv->newgotonode = TrieCheckFrom;
                                                tv->newindex = TrieCheckFrom->index;
                                            }
                                            else {
                                                tv->newgotonode = NULL;
                                                tv->newindex = 0;
                                            }
                                        }
                                    }
                                    else {
                                        if (TrieCandFrom) {
                                            TrieCheckFrom->index += TrieCandFrom->index;
                                            TrieCandFrom->goes_to = TrieCheckFrom;
                                        } else {
                                            TrieCheckFrom->index++;
                                        }
                                        NextCand->do_it = FALSE;
                                    }
                                }
                                break;
                            case 1:
                                TrieCheckFrom->index ++;
                                tv->gotonode = TrieCheckFrom;
                                break;
                            case 2:
                                if (TrieCandFrom) {
                                    TrieCheckFrom->index += TrieCandFrom->index;
                                    TrieCandFrom->goes_to = TrieCheckFrom;
                                }
                                else {
                                    TrieCheckFrom->index++;
                                }
                                if (temp == tv->maxtreelevel) {
                                    tmp1 = TempOrbits[NextCand->lab[Spine[temp].tgtpos]];
                                    for (i=1; i<AutomCount[0]; i++) if (AutomCount[i] == tmp1) break;
                                    if (i == AutomCount[0]) AutomCount[AutomCount[0]++] = TempOrbits[NextCand->lab[Spine[temp].tgtpos]];
                                }
                                break;
                            default:
                                break;
                        }
						return temp;
					}
				}
			}
			CheckAutList = CheckAutList->next;
		}
		CheckLevel++;
	}
	return FALSE;
}

int CheckForSingAutomorphisms(sparsegraph *g_arg, Candidate *CurrCand, Partition *NextPart, Candidate *NextCand,
							  struct TracesVars* tv, struct TracesInfo* ti,
							  int m, int n)
{
	int i, temp, tmp, tmp1, result, tgt_level;
	TracesSpine *SpineTL;
	Candidate *CheckAutList;
    searchtrie *TrieCandFrom, *TrieCheckFrom;
    
	SpineTL = Spine+tv->tolevel;
	CheckAutList = SpineTL->liststart;
    tv->gotonode = NULL;
    
	result = 0;
	while (CheckAutList != NULL) {
        if (CheckAutList->do_it && (CheckAutList->stnode->father == CurrCand->stnode)) {
            if (CheckAutList->firstsingcode == NextCand->firstsingcode) {
                memcpy(AUTPERM, IDENTITY_PERM, n*sizeof(int));
                if ((tv->tolevel == 1) && (Spine[0].part->cells == 1)) tmp = 2; else tmp = tv->tolevel;
                
                if (TreeFyTwo(g_arg, tmp, tv->tolevel_tl, CheckAutList, NextCand, NextPart, n, tv, ti)) {
                    if (isautom_sg((graph*)g_arg, AUTPERM, tv->options->digraph, m, n)) {
                        if (!findperm(gensB, AUTPERM, n, tv)) {
                            if (tv->options->verbosity >= 2) tv->schreier3 -= CPUTIME;
                            addgenerator(&gpB, &gensB, AUTPERM, n);
                            if (tv->options->verbosity >= 2) tv->schreier3 += CPUTIME;
                            result = CheckAutList->name;
                            if (TempOrbits) orbjoin(TempOrbits, AUTPERM, n);
                            tv->stats->numgenerators++;
                            tv->stats->numorbits = orbjoin(tv->orbits, AUTPERM, n);
                            ti->thegrouphaschanged = TRUE;
                            ti->identitygroup = FALSE;
                            if (tv->options->verbosity >= 2) fprintf(outfile, "[a(%d)] ", CheckAutList->name);
                            if (tv->options->verbosity >= 2 && tv->options->writeautoms) {
                                PRINT_RETURN
                            }
                            if (tv->options->writeautoms) {
                                fprintf(outfile, "Gen #%d: ", tv->stats->numgenerators);
                                writeperm(outfile, AUTPERM, tv->options->cartesian, tv->options->linelength, n);
                            }
                            if (tv->options->userautomproc) {
                                (*tv->options->userautomproc)(tv->stats->numgenerators, AUTPERM, n);
                            }
                        }
                        else {
                            if (tv->options->verbosity >= 2) {
                                fprintf(outfile, "[a*]");
                            }
                            if (TempOrbits) orbjoin(TempOrbits, AUTPERM, n);
                            result = -CheckAutList->name;
                        }
                        
                        
                        
                        TrieCandFrom = NULL;
                        TrieCheckFrom = CheckAutList->stnode;
                        if (CurrCand->stnode->level <= 1) {
                            tgt_level = CurrCand->stnode->level + 1;
                            while (TrieCheckFrom->level > tgt_level) {
                                TrieCheckFrom = TrieCheckFrom->father;
                            }
                        }
                        else {
                            if (tv->tolevel <= TrieCheckFrom->level) {
                                tgt_level = tv->tolevel;
                                while (TrieCheckFrom->level != tgt_level) {
                                    TrieCheckFrom = TrieCheckFrom->father;
                                }
                            } else {
                                TrieCandFrom = CurrCand->stnode;
                                tgt_level = TrieCheckFrom->level;
                                while (TrieCandFrom->level != tgt_level) {
                                    TrieCandFrom = TrieCandFrom->father;
                                }
                            }
                        }
                        if (TrieCandFrom) {
                            while (TrieCandFrom->father != TrieCheckFrom->father) {
                                TrieCandFrom = TrieCandFrom->father;
                                TrieCheckFrom = TrieCheckFrom->father;
                            }
                        }
                        else {
                            if ((TrieCheckFrom->level > 1) && (TrieCheckFrom->father != CurrCand->stnode)) {
                                TrieCandFrom = CurrCand->stnode;
                                TrieCheckFrom = TrieCheckFrom->father;
                                while (TrieCandFrom->father != TrieCheckFrom->father) {
                                    TrieCandFrom = TrieCandFrom->father;
                                    TrieCheckFrom = TrieCheckFrom->father;
                                }
                            }
                        }
                        
                        while (TrieCheckFrom->goes_to) {
                            TrieCheckFrom = TrieCheckFrom->goes_to;
                        }
                        
                        
                        
                        
                        for (temp=1; temp<=tv->tolevel; temp++) {
                            if (CheckAutList->lab[Spine[temp].tgtpos] != NextCand->lab[Spine[temp].tgtpos]) {
                                break;
                            }
                        }
                        
                        switch (tv->compstage) {
                            case 0:
                                if (tv->strategy && (tv->steps == 1)) {
                                    RemoveFromLevel(temp, tv->maxtreelevel-1, tv->strategy, FALSE);
                                    if (TrieCandFrom) {
                                        TrieCheckFrom->index += TrieCandFrom->index;
                                        TrieCandFrom->goes_to = TrieCheckFrom;
                                    }
                                    else TrieCheckFrom->index++;
                                    NextCand->do_it = FALSE;
                                }
                                else {
                                    if (CheckAutList->lab[Spine[temp].tgtpos] >= NextCand->lab[Spine[temp].tgtpos]) {
                                        CheckAutList->do_it = FALSE;
                                        if (TrieCandFrom) {
                                            TrieCandFrom->index += TrieCheckFrom->index;
                                            tv->newindex = 0;
                                            TrieCheckFrom->goes_to = TrieCandFrom;
                                        } else {
                                            if (CurrCand->stnode->level > 1) {
                                                tv->newgotonode = TrieCheckFrom;
                                                tv->newindex = TrieCheckFrom->index;
                                            }
                                            else {
                                                tv->newgotonode = NULL;
                                                tv->newindex = 0;
                                            }
                                        }
                                    }
                                    else {
                                        if (TrieCandFrom) {
                                            TrieCheckFrom->index += TrieCandFrom->index;
                                            TrieCandFrom->goes_to = TrieCheckFrom;
                                        } else {
                                            TrieCheckFrom->index++;
                                        }
                                        NextCand->do_it = FALSE;
                                    }
                                }
                                break;
                            case 1:
                                TrieCheckFrom->index ++;
                                tv->gotonode = TrieCheckFrom;
                                break;
                            case 2:
                                if (TrieCandFrom) {
                                    TrieCheckFrom->index += TrieCandFrom->index;
                                    TrieCandFrom->goes_to = TrieCheckFrom;
                                } else {
                                    TrieCheckFrom->index++;
                                }
                                if (temp == tv->maxtreelevel) {
                                    tmp1 = TempOrbits[NextCand->lab[Spine[temp].tgtpos]];
                                    for (i=1; i<AutomCount[0]; i++) if (AutomCount[i] == tmp1) break;
                                    if (i == AutomCount[0]) AutomCount[AutomCount[0]++] = TempOrbits[NextCand->lab[Spine[temp].tgtpos]];
                                }
                                break;
                            default:
                                break;
                        }
                        
                        return result;
                    }
                }
            }
        }
		CheckAutList = CheckAutList->next;
	}
	return result;
}

int CheckForMatching(sparsegraph *g_arg, Candidate *CurrCand, Candidate *NextCand, Partition *Part, struct TracesVars* tv, struct TracesInfo* ti, int m, int n)
{
	int i, vtx, temp, tmp1, result, tgt_level;
	TracesSpine *SpineTL;
	Candidate *CheckAutList;
    int *cls;
    searchtrie *TrieCandFrom, *TrieCheckFrom;
    
    tv->gotonode = NULL;
	SpineTL = Spine+tv->tolevel;
	CheckAutList = SpineTL->liststart;
	cls = Part->cls;
    
	while (CheckAutList != NULL) {
        if (CheckAutList->do_it && (CheckAutList->singcode == NextCand->singcode)) {
            memcpy(AUTPERM, IDENTITY_PERM, n*sizeof(int));
            for (i=0; i<n; i+=cls[i]) {
                vtx = NextCand->lab[i];
                if ((cls[i] == 1) && (cls[Part->inv[CheckAutList->invlab[vtx]]] == 1)) {
                    AUTPERM[vtx] = CheckAutList->lab[i];
                }
                else {
                    if ((CheckAutList->invlab[vtx]<i) || (CheckAutList->invlab[vtx]>=i+cls[i])) {
                        break;
                    }
                }
            }
            if (i==n && (isautom_sg((graph*)g_arg, AUTPERM, tv->options->digraph, m, n))) {
                
                if (!findperm(gensB, AUTPERM, n, tv)) {
                    if (tv->options->verbosity >= 2) tv->schreier3 -= CPUTIME;
                    addgenerator(&gpB, &gensB, AUTPERM, n);
                    if (tv->options->verbosity >= 2) tv->schreier3 += CPUTIME;
                    result = CheckAutList->name;
                    if (TempOrbits) orbjoin(TempOrbits, AUTPERM, n);
                    tv->stats->numgenerators++;
                    tv->stats->numorbits = orbjoin(tv->orbits, AUTPERM, n);
                    ti->thegrouphaschanged = TRUE;
                    ti->identitygroup = FALSE;
                    if (tv->options->verbosity >= 2) fprintf(outfile, "[M(%d)] ", CheckAutList->name);
                    if (tv->options->verbosity >= 2 && tv->options->writeautoms) {
                        PRINT_RETURN
                    }
                    if (tv->options->writeautoms) {
                        fprintf(outfile, "Gen #%d: ", tv->stats->numgenerators);
                        writeperm(outfile, AUTPERM, tv->options->cartesian, tv->options->linelength, n);
                    }
                    if (tv->options->userautomproc) {
                        (*tv->options->userautomproc)(tv->stats->numgenerators, AUTPERM, n);
                    }
                }
                else {
                    if (tv->options->verbosity >= 2) {
                        fprintf(outfile, "[M*]");
                    }
                    if (TempOrbits) orbjoin(TempOrbits, AUTPERM, n);
                    result = -CheckAutList->name;
                }
                
                TrieCandFrom = NULL;
                TrieCheckFrom = CheckAutList->stnode;
                if (CurrCand->stnode->level <= 1) {
                    tgt_level = CurrCand->stnode->level + 1;
                    while (TrieCheckFrom->level > tgt_level) {
                        TrieCheckFrom = TrieCheckFrom->father;
                    }
                }
                else {
                    if (tv->tolevel <= TrieCheckFrom->level) {
                        tgt_level = tv->tolevel;
                        while (TrieCheckFrom->level != tgt_level) {
                            TrieCheckFrom = TrieCheckFrom->father;
                        }
                    } else {
                        TrieCandFrom = CurrCand->stnode;
                        tgt_level = TrieCheckFrom->level;
                        while (TrieCandFrom->level != tgt_level) {
                            TrieCandFrom = TrieCandFrom->father;
                        }
                    }
                }
                if (TrieCandFrom) {
                    while (TrieCandFrom->father != TrieCheckFrom->father) {
                        TrieCandFrom = TrieCandFrom->father;
                        TrieCheckFrom = TrieCheckFrom->father;
                    }
                }
                else {
                    if ((TrieCheckFrom->level > 1) && (TrieCheckFrom->father != CurrCand->stnode)) {
                        TrieCandFrom = CurrCand->stnode;
                        TrieCheckFrom = TrieCheckFrom->father;
                        while (TrieCandFrom->father != TrieCheckFrom->father) {
                            TrieCandFrom = TrieCandFrom->father;
                            TrieCheckFrom = TrieCheckFrom->father;
                        }
                    }
                }
                
                while (TrieCheckFrom->goes_to) {
                    TrieCheckFrom = TrieCheckFrom->goes_to;
                }
                
                for (temp=1; temp<=tv->tolevel; temp++) {
                    if (CheckAutList->lab[Spine[temp].tgtpos] != NextCand->lab[Spine[temp].tgtpos]) {
                        break;
                    }
                }
                switch (tv->compstage) {
                    case 0:
                        if (tv->strategy && (tv->steps == 1)) {
                            RemoveFromLevel(temp, tv->maxtreelevel-1, tv->strategy, FALSE);
                            if (TrieCandFrom) {
                                TrieCheckFrom->index += TrieCandFrom->index;
                                TrieCandFrom->goes_to = TrieCheckFrom;
                            }
                            else TrieCheckFrom->index++;
                            NextCand->do_it = FALSE;
                        }
                        else {
                            if (CheckAutList->lab[Spine[temp].tgtpos] >= NextCand->lab[Spine[temp].tgtpos]) {
                                CheckAutList->do_it = FALSE;
                                if (TrieCandFrom) {
                                    TrieCandFrom->index += TrieCheckFrom->index;
                                    tv->newindex = 0;
                                    TrieCheckFrom->goes_to = TrieCandFrom;
                                } else {
                                    if (CurrCand->stnode->level > 1) {
                                        tv->newgotonode = TrieCheckFrom;
                                        tv->newindex = TrieCheckFrom->index;
                                    }
                                    else {
                                        tv->newgotonode = NULL;
                                        tv->newindex = 0;
                                    }
                                }
                            }
                            else {
                                if (TrieCandFrom) {
                                    TrieCheckFrom->index += TrieCandFrom->index;
                                    TrieCandFrom->goes_to = TrieCheckFrom;
                                } else {
                                    TrieCheckFrom->index++;
                                }
                                NextCand->do_it = FALSE;
                            }
                        }
                        break;
                    case 1:
                        TrieCheckFrom->index ++;
                        tv->gotonode = TrieCheckFrom;
                        break;
                    case 2:
                        if (TrieCandFrom) {
                            TrieCheckFrom->index += TrieCandFrom->index;
                            TrieCandFrom->goes_to = TrieCheckFrom;
                        } else {
                            TrieCheckFrom->index++;
                        }
                        if (temp == tv->maxtreelevel) {
                            tmp1 = TempOrbits[NextCand->lab[Spine[temp].tgtpos]];
                            for (i=1; i<AutomCount[0]; i++) if (AutomCount[i] == tmp1) break;
                            if (i == AutomCount[0]) AutomCount[AutomCount[0]++] = TempOrbits[NextCand->lab[Spine[temp].tgtpos]];
                        }
                        break;
                    default:
                        break;
                }
                return temp;
            }
        }
        CheckAutList = CheckAutList->next;
	}
    return FALSE;
}

void Individualize(Partition *NextPart, Candidate *NextCand, int K, int Tc, int Cl, int Pos)
{
	int i, j;
	
	NextCand->do_it = TRUE;
	if (NextPart->cls[Tc] > 1) {
		NextPart->cells = Cl+1;
		NextPart->active = 1;
		NextPart->cls[Tc]--;
		NextPart->cls[Pos] = 1;
    }
	NextPart->inv[Pos] = Pos;
	
	j = NextCand->lab[Pos];
	i = NextCand->invlab[K];
	NextCand->lab[Pos] = K;
	NextCand->invlab[K] = Pos;
	NextCand->lab[i] = j;
	NextCand->invlab[j] = i;
	return;
}

boolean TargetCell(sparsegraph *sg, Candidate *TargCand, Partition *Part, int n, struct TracesVars* tv, int Lv) {
	int TCell = -1, TCSize = 1;
	int i;
    
	if (Lv < tv->tcellevel) {
		tv->tcell = Spine[Lv+1].tgtcell;
		return TRUE;
	}
	else {
		while (TCell < 0) {
            for (i = Spine[Lv].tgtcell; i < Spine[Lv].tgtend; i += Part->cls[i]) {
                if ((sg->d[TargCand->lab[i]] > 2) && (sg->d[TargCand->lab[i]] < n-1)) {
                    if (Part->cls[i] > TCSize) {
						TCSize = Part->cls[i];
						TCell = i;
					}
				}
			}
			Lv--;
			if ((Lv < 0) && (TCell < 0)) return FALSE;
		}
        tv->tcell = TCell;
		return TRUE;
	}
}

boolean TargetCellFirstPath(sparsegraph *sg, Candidate *TargCand, Partition *Part, int n, struct TracesVars* tv) {
	int TCell = -1, TCSize = 1;
	int Lv, i, Lev;
    
    Lev = Lv = tv->tolevel_tl++;
    while (TCell < 0) {
        for (i = Spine[Lv].tgtcell; i < Spine[Lv].tgtend; i += Part->cls[i]) {
            if ((sg->d[TargCand->lab[i]] > 2) && (sg->d[TargCand->lab[i]] < n-1)) {
                if (Part->cls[i] > TCSize) {
                    TCSize = Part->cls[i];
                    TCell = i;
                }
            }
        }
        Lv--;
        if ((Lv < 0) && (TCell < 0)) {
            tv->smalldeglevel = Lev;
            tv->finalnumcells = Part->cells;
            return FALSE;
        }
    }
    tv->tcellexpath = TCell;
    Spine[tv->tolevel_tl].tgtcell = tv->tcellexpath;
    Spine[tv->tolevel_tl].tgtsize = TCSize;
    Spine[tv->tolevel_tl].tgtend = Spine[tv->tolevel_tl].tgtcell + TCSize;
    Spine[tv->tolevel_tl].tgtpos = Spine[tv->tolevel_tl].tgtend - 1;
    tv->tcellevel = tv->tolevel_tl;
    
    if (Lv+1 != Lev) {
        BreakSteps[Lev] = ++tv->brkstpcount;
        if (!Spine[tv->tolevel].liststart->firstsingcode) {
            Spine[tv->tolevel].liststart->firstsingcode = Spine[tv->tolevel].liststart->pathsingcode;
        }
    }
    return TRUE;
}

boolean TargetCellExpPath(sparsegraph *sg, Candidate *TargCand, Partition *Part, int n, struct TracesVars* tv) {
	int TCell = -1, TCSize = 1;
	int Lv, i;
    
    if (Part->cells == tv->finalnumcells) {
		return FALSE;
	}
    Lv = tv->tolevel_tl+1;
    SpineTL_tl = Spine+Lv;
	if (tv->tolevel_tl < tv->tcellevel) {
        tv->tcellexpath = Part->inv[SpineTL_tl->tgtcell];
        tv->tcellfromlevel = SpineTL_tl->tgtfrom;
        tv->tolevel_tl++;
        return ((Spine[tv->tolevel_tl].tgtcell >= Spine[tv->tolevel_tl-1].tgtcell) && (Spine[tv->tolevel_tl].tgtend <= Spine[tv->tolevel_tl-1].tgtend));
    }
	else {
        while (TCell < 0) {
            Lv--;
            for (i = Part->inv[Spine[Lv].tgtcell]; i < Spine[Lv].tgtend; i += Part->cls[i]) {
                if ((sg->d[TargCand->lab[i]] > 2) && (sg->d[TargCand->lab[i]] > 2)) {
					if (Part->cls[i] > TCSize) {
						TCSize = Part->cls[i];
						TCell = i;
					}
				}
			}
		}
        tv->tolevel_tl++;
		tv->tcellevel = tv->tolevel_tl;
        Spine[tv->tolevel_tl].tgtcell = tv->tcellexpath = TCell;
        tv->tcellfromlevel = Spine[tv->tolevel_tl].tgtfrom = Lv+1;
        Spine[tv->tolevel_tl].tgtsize = Part->cls[TCell];
		Spine[tv->tolevel_tl].tgtend = TCell+Part->cls[TCell];
		Spine[tv->tolevel_tl].tgtpos = Spine[tv->tolevel_tl].tgtend - 1;
        return ((Spine[tv->tolevel_tl].tgtcell >= Spine[tv->tolevel_tl-1].tgtcell) && (Spine[tv->tolevel_tl].tgtend <= Spine[tv->tolevel_tl-1].tgtend));
	}
}

boolean SelectNextLevel(sparsegraph *g_arg, int n, struct TracesVars *tv) {
    
	switch (tv->compstage) {
		case 2:
			if (tv->smalldeglevel < n) {
                tv->nextlevel = tv->smalldeglevel - 1;
			}
			else {
				tv->nextlevel = tv->maxtreelevel;
			}
			while (tv->nextlevel >=0) {
				if (Spine[tv->nextlevel].liststart) {
					break;
				}
				tv->nextlevel--;
			}
			if (tv->nextlevel < 0) {
				return FALSE;
			}
			break;
		default:
			switch (tv->strategy) {
				case 0:
					tv->nextlevel = tv->fromlevel;
					while (!Spine[tv->nextlevel].liststart) {
						(tv->nextlevel)++;
					}
					if ((Spine[tv->nextlevel].part->cells == tv->finalnumcells) || (tv->nextlevel > tv->maxtreelevel)) {
                        return FALSE;
					}
					break;
				case 1:
					tv->nextlevel = tv->maxtreelevel;
                    if (Spine[tv->nextlevel].part->cells == tv->finalnumcells) {
						(tv->nextlevel)--;
					}
                    while (tv->nextlevel >= 0) {
						if (Spine[tv->nextlevel].liststart) {
                            break;
						}
						tv->nextlevel--;
                    }
					if (tv->nextlevel < 0) {
						return FALSE;
					}
					break;
				default:
					break;
			}
			break;
	}
    return TRUE;
}

boolean TreeFyTwo(sparsegraph *sg, int From, int Num, Candidate *Cand1, Candidate *Cand2, Partition *Part, int n,
                  struct TracesVars* tv, struct TracesInfo *ti)
{
	int k, i, i1, i2;
	size_t j1, iend;
    
	i2=0;
	memcpy(AUTPERM, IDENTITY_PERM, n*sizeof(int));
    i1 = Spine[From].tgtsize;
    memcpy(TreeStack, IDENTITY_PERM+Spine[From].tgtcell, i1*sizeof(int));
    for (i=0; i<i1; i++) {
        AUTPERM[Cand1->lab[TreeStack[i]]] = Cand2->lab[TreeStack[i]];
    }
    while (i2 < i1) {
		k = TreeStack[i2++];
        if (Cand1->lab[k] == Cand2->lab[k]) {
            continue;
        }
        
        j1 = sg->v[Cand1->lab[k]];
		iend = j1+sg->d[Cand1->lab[k]];
        for (; j1 < iend; j1++) {
			k = sg->e[j1];
            if ((k != Cand2->lab[Cand1->invlab[k]]) && (AUTPERM[k] == k)) {
                AUTPERM[k] = Cand2->lab[Cand1->invlab[k]];
                TreeStack[i1++] = Cand1->invlab[k];
            }
		}
    }
    return (i1 != Spine[From].tgtsize);
}

void ExperimentalStep(sparsegraph *g_arg, Partition *NextPart, Candidate *NextCand,
					  TracesVars *tv, TracesInfo *ti, int m, int n) {
	int i, iend, min, tmp;
    
	SpineTL_tl = Spine+tv->tolevel_tl;
	NextPart->active = 1;
	
	/* EXPERIMENTAL PATH INDIVIDUALIZATION AND REFINEMENT */
    if (tv->answ == 2) {
        min = NextCand->lab[tv->tcellexpath];
        tmp = tv->tcellexpath;
        iend = tv->tcellexpath + NextPart->cls[tv->tcellexpath];
        for (i=tv->tcellexpath + 1; i<iend ; i++) {
            if (NextCand->lab[i] < min) {
                min = NextCand->lab[i];
                tmp = i;
            }
        }
    }
    else {
        tmp = tv->tcellexpath+KRAN(NextPart->cls[tv->tcellexpath]);
    }
    if (NextPart->cls[tv->tcellexpath] == 2) {
        NextCand->pathsingcode = MASHCOMM(NextCand->pathsingcode, NextCand->lab[tv->tcellexpath]);
        NextCand->pathsingcode = MASHCOMM(NextCand->pathsingcode, NextCand->lab[tv->tcellexpath+1]);
    }
    else {
        NextCand->pathsingcode = MASHCOMM(NextCand->pathsingcode, NextCand->lab[tmp]);
    }
    
    Individualize(NextPart, NextCand, NextCand->lab[tmp], tv->tcellexpath, NextPart->cells, tv->tcellexpath + NextPart->cls[tv->tcellexpath]-1);
    
	tv->stats->numnodes++;
	if (tv->compstage == 0) {
		traces_refine_notrace(g_arg,
							  NextCand,
							  m,
							  n,
							  NextPart, tv, ti);
	}
	else {
		if (tv->tolevel_tl == tv->maxtreelevel+1) {
			trieref = trieroot;
			tv->answ = traces_refine_comptrie(g_arg,
											  NextCand,
											  m,
											  n,
											  NextPart, tv, ti);
			if (tv->answ == 0 ) {
				tv->stats->interrupted++;
			}
		}
		else {
			traces_refine_notrace(g_arg,
								  NextCand,
								  m,
								  n,
								  NextPart, tv, ti);
		}
	}
}

trie* trie_new(int n, struct TracesVars* tv)
{
	TrieArray[0] = malloc(n*sizeof(trie));
	if (TrieArray[0] == NULL) {
		fprintf(ERRFILE, "\nError, memory not allocated.\n");
		exit(1);
	}
	TrieArray[0][0].first_child = TrieArray[0][0].next_sibling = NULL;
	tv->triepos = 0;
	tv->trienext = 1;
	return TrieArray[0];
}

struct trie *trie_make(trie *t, int value, int n, struct TracesVars* tv)
{
	trie *t1;
	t1 = t;
	if (tv->trienext == n) {
		tv->trienext = 0;
		tv->triepos++;
		TrieArray[tv->triepos] = malloc(n*sizeof(trie));
		if (TrieArray[tv->triepos] == NULL) {
			fprintf(ERRFILE, "\nError, memory not allocated.\n");
			exit(1);
		}
	}
	if (t->first_child) {
		t = t->first_child;
		if (value < t->value) {
			t1->first_child = &TrieArray[tv->triepos][tv->trienext++];
			t1->first_child->next_sibling = t;
			t1->first_child->first_child = NULL;
			t = t1->first_child;
			t->value = value;
			return t;
		}
		while (value > t->value) {
			t1 = t;
			if (t->next_sibling) {
				t = t->next_sibling;
			}
			else break;
		}
		if (value == t->value) {
			return t;
		}
		t1->next_sibling = &TrieArray[tv->triepos][tv->trienext++];
		t1->next_sibling->first_child = t1->next_sibling->next_sibling = NULL;
		if (t != t1) {
			t1->next_sibling->next_sibling = t;
		}
		t = t1->next_sibling;
	}
	else {
		t->first_child = &TrieArray[tv->triepos][tv->trienext++];
		t = t->first_child;
		t->first_child = t->next_sibling = NULL;
	}
	t->value = value;
	return t;
}

struct trie *trie_comp(trie *t, int value)
{
	if (t->first_child) {
		t = t->first_child;
		while (t) {
			if  (value != t->value) {
				t = t->next_sibling;
			}
			else {
				break;
			}
		}
		return t;
	}
	else {
		return NULL;
	}
}


void CopyCand(Candidate *W, Candidate *V, int n) {
	memcpy(W->lab, V->lab, n*sizeof(int));
	memcpy(W->invlab, V->invlab, n*sizeof(int));
	W->name = V->name;
    W->vertex = V->vertex;
	W->code = V->code;
	W->singcode = V->singcode;
	W->firstsingcode = V->firstsingcode;
	W->do_it = V->do_it;
    W->sortedlab = FALSE;
}

void RemoveFromLevel(int from, int to, int strategy, boolean reinit) {
	int i;
    
	for (i=from; i<=to; i++) {
		if (Spine[i].listend) {
			(Spine[i].listend)->next = GarbList;
			GarbList = Spine[i].liststart;
			Spine[i].liststart = Spine[i].listend = NULL;
		}
        if (strategy == 0 || reinit) {
            Spine[i].listcounter = 0;
            if (i>from) {
                Spine[i].thetracexists = FALSE;
                Spine[i].part->code = -1;
            }
        }
	}
}

void CompStage0(sparsegraph *g_arg, Partition *CurrPart, Partition *NextPart, Candidate *CurrCand, Candidate *NextCand,
				int m, int n, struct TracesVars* tv, struct TracesInfo *ti)
{
	int i, j, i1, j2, k, cu, cu1;
	int temp, tmp, auxcode, search_vtx, gom_level;
	boolean closeloop, firstsing;
	Candidate *SpTLliststart, *AuxCand;
    searchtrie *TreeNode, *TreeNode1, *TreeNode2;
	
    PRINT_FROM_VERB(3)
    
	if (TargetCell(g_arg, CurrCand, CurrPart, n, tv, tv->tolevel)) {
        ++tv->tolevel;
		SpineTL = Spine+tv->tolevel;
		SpineTL->tgtcell = tv->tcell;
        SpineTL->tgtsize = CurrPart->cls[tv->tcell];
		SpineTL->tgtend = tv->tcell+SpineTL->tgtsize;
		SpineTL->tgtpos = SpineTL->tgtend - 1;
	}
	else {
		tv->smalldeglevel = tv->tolevel;
        tv->finalnumcells = CurrPart->cells;
		ti->thereisnextlevel = SelectNextLevel(g_arg, n, tv);
		return;
	}
    
    tv->newgotonode = NULL;
	
	/*  CANDIDATE */
	temp = CurrCand->lab[Spine[1].tgtpos];
	k = SpineTL->tgtend;
	
    TreeNode = CurrCand->stnode;
    while (TreeNode) {
        if (TreeNode->goes_to) {
            CurrCand->do_it = FALSE;
            break;
        }
        TreeNode = TreeNode->father;
    }
    
	if (CurrCand->do_it) {
        if ((tv->orbits[temp] == temp) || tv->tolevel == 1) {
            ti->minimalinorbits = TRUE;
            if (((Spine[tv->fromlevel].liststart != Spine[tv->fromlevel].listend) && (CurrPart->cls[tv->tcell] > 6))
                || tv->strategy
                || (tv->expathlength <=10)) {
                TempOrbits = NULL;
                tv->samepref = FixBase(fix, tv, CurrCand, 0, tv->fromlevel);
                if (tv->samepref != tv->nfix || ti->thegrouphaschanged) {
                    if (tv->options->verbosity >= 2) tv->schreier1 -= CPUTIME;
                    gom_level = getorbitsmin(fix, tv->nfix, gpB, &gensB, &tv->currorbit,
                                             CurrCand->lab+tv->tcell, CurrPart->cls[tv->tcell], n, TRUE);
                    if (tv->options->verbosity >= 2) tv->schreier1 += CPUTIME;
                    ti->thegrouphaschanged = FALSE;
                    
                    if (gom_level < tv->nfix) {
                        PRINT_NOTMIN_VERB(3)
                        
                        TreeNode = CurrCand->stnode;
                        j2 = CurrCand->lab[Spine[gom_level+1].tgtpos];
                        i1 = tv->currorbit[j2];
                        for (j=0; j < tv->nfix - gom_level; j++) {
                            TreeNode = TreeNode->father;
                        }
                        TreeNode1 = TreeNode->first_child;
                        while (TreeNode1) {
                            if (TreeNode1->vtx == i1) {
                                break;
                            }
                            TreeNode1 = TreeNode1->next_sibling;
                        }
                        if (TreeNode1) {
                            while (TreeNode1->goes_to) {
                                TreeNode1 = TreeNode1->goes_to;
                            }
                            TreeNode2 = TreeNode1->next_sibling;
                            while (TreeNode2->vtx != j2) {
                                TreeNode2 = TreeNode2->next_sibling;
                            }
                            TreeNode1->index += TreeNode2->index;
                            TreeNode2->goes_to = TreeNode1;
                            
                            ti->minimalinorbits = FALSE;
                        }
                        else {
                            tv->currorbit = getorbits(fix, tv->nfix, gpB, &gensB, n);
                        }
                    }
                }
                else {
                    tv->currorbit = findcurrorbits(gpB, tv->nfix);
                }
            }
            else {
                TempOrbits = WorkArray1;
                memcpy(TempOrbits, IDENTITY_PERM, n*sizeof(int));
                tv->currorbit = TempOrbits;
            }
            
            if (ti->minimalinorbits) {
                memcpy(NextCand->lab, CurrCand->lab, n*sizeof(int));
                memcpy(NextCand->invlab, CurrCand->invlab, n*sizeof(int));
                auxcode = CurrCand->code;
                SpineTL->trcstart = CurrPart->cells;
                TheTrace[SpineTL->trcstart] = SpineTL->tgtpos;
                if (!CurrCand->sortedlab) {
                    quickSort(CurrCand->lab+tv->tcell, CurrPart->cls[tv->tcell]);
                    for (i=tv->tcell; i<tv->tcell+CurrPart->cls[tv->tcell]; i++) {
                        CurrCand->invlab[CurrCand->lab[i]] = i;
                    }
                    CurrCand->sortedlab = TRUE;
                }
                
                tv->indivstart = tv->tcell+CurrCand->indnum;
                tv->indivend = tv->indivstart+tv->steps;
                if (tv->indivend > SpineTL->tgtend) {
                    tv->indivend = SpineTL->tgtend;
                }
                
                temp = CurrCand->lab[tv->indivstart];
                
                for (k = tv->indivstart; k < tv->indivend; k++) {
                    CurrCand->indnum++;
                    NextCand->singcode = CurrCand->singcode;
                    NextCand->vertex = CurrCand->lab[k];
                    NextCand->name = ++tv->name;
                    if (NextCand->name == (NAUTY_INFINITY-2)) {
                        NextCand->name = tv->name = 1;
                    }
                    
                    PRINT_INDIV_VERB(3)
                    if (tv->currorbit[NextCand->vertex] != NextCand->vertex) {
                        PRINT_SKIPPED_VERB(3)
                        
                        search_vtx = tv->currorbit[NextCand->vertex];
                        TreeNode = CurrCand->stnode;
                        if (TreeNode->first_child) {
                            TreeNode = TreeNode->first_child;
                            while (TreeNode) {
                                if (TreeNode->vtx == search_vtx) {
                                    break;
                                }
                                TreeNode = TreeNode->next_sibling;
                            }
                            if (TreeNode) {
                                while (TreeNode->goes_to) {
                                    TreeNode = TreeNode->goes_to;
                                }
                                TreeNode->index++;
                                continue;
                            }
                        }
                    }
                    PRINT_REFINE_VERB(3)
                    
                    memcpy(NextPart->cls, CurrPart->cls, n*sizeof(int));
                    memcpy(NextPart->inv, CurrPart->inv, n*sizeof(int));
                    
                    if (NextPart->cls[tv->tcell] == 2) {
                        NextCand->singcode = MASHCOMM(NextCand->singcode, CurrCand->lab[tv->tcell]);
                        NextCand->singcode = MASHCOMM(NextCand->singcode, CurrCand->lab[tv->tcell+1]);
                    }
                    else {
                        NextCand->singcode = MASHCOMM(NextCand->singcode, NextCand->vertex);
                    }
                    
                    Individualize(NextPart, NextCand, NextCand->vertex, tv->tcell, CurrPart->cells, SpineTL->tgtpos);
                    tv->stats->numnodes++;
                    
                    tv->answ = traces_refine(g_arg,
                                             NextCand,
                                             m,
                                             n,
                                             NextPart, tv, ti, TRUE);
                    
                    if (!tv->strategy && !tv->options->getcanon && (NextPart->cells == tv->finalnumcells) && (tv->tolevel > 1)) {
                        temp = 0;
                        for (i=0; i<tv->tolevel; i++) {
                            temp += Spine[i].listcounter;
                        }
                        if (temp > 5) {
                            EXITFROMSTAGE0REFINE
                        }
                    }
                    switch (tv->answ) {
                        case 0:				/* Interrupted refinement: do not add to the list */
                            tv->stats->interrupted++;
                            SpineTL->levelcounter++;
                            break;
                        case 1 :			/* The same trace has been found once more : add to the list */
                            SpineTL->levelcounter++;
                            
                            NextCand->do_it = TRUE;
                            if (tv->options->verbosity >= 2) PRINT_CANDIDATE(NextCand, tv->tolevel);
                            
                            tv->tolevel_tl = tv->tolevel;
                            NextCand->pathsingcode = NextCand->singcode;
                            NextCand->firstsingcode = 0;
                            
                            if (tv->steps > 1) {
                                if (tv->tolevel < tv->tcellevel) {
                                    closeloop = CheckForMatching(g_arg, CurrCand, NextCand, NextPart, tv, ti, m, n);
                                    if (NextCand->do_it) {
                                        firstsing = TRUE;
                                        if (tv->options->verbosity >= 2) tv->expaths -= CPUTIME;
                                        
                                        /* EXPERIMENTAL PATH */
                                        while (NextPart->cells < tv->finalnumcells) {
                                            if (firstsing && BreakSteps[tv->tolevel]) {
                                                firstsing = FALSE;
                                                NextCand->firstsingcode = NextCand->pathsingcode;
                                                if (CheckForSingAutomorphisms(g_arg, CurrCand, NextPart, NextCand, tv, ti, m, n))
                                                    if (!NextCand->do_it) {
                                                        break;
                                                    }
                                            }
                                            if (!TargetCellExpPath(g_arg, NextCand, NextPart, n, tv)) {
                                                NextCand->firstsingcode = NextCand->pathsingcode;
                                            }
                                            ExperimentalStep(g_arg, NextPart, NextCand, tv, ti, m, n);
                                            PRINT_EXPPATHSTEP(NextCand, TRUE)
                                        }
                                        
                                        if (NextPart->cells == tv->finalnumcells) {
                                            UPDATEMIN(tv->expathlength, tv->tolevel_tl);
                                        }
                                        if (tv->options->verbosity >= 2) tv->expaths += CPUTIME;
                                    }
                                    else {
                                        if (closeloop < tv->tolevel) k = SpineTL->tgtend;
                                        PRINT_RETURN
                                        break;
                                    }
                                }
                                
                                if (!tv->strategy && !tv->options->getcanon && (NextPart->cells == tv->finalnumcells) && (tv->tolevel_tl == tv->tolevel + 1)) {
                                    tv->maxtreelevel = tv->tolevel_tl;
                                    if (tv->tolevel == 1) {
                                        tv->newst_stage1 = searchtrie_make(CurrCand, NextCand, n, tv);
                                        EXITFROMSTAGE0EXPATH1
                                    }
                                    else {
                                        temp = 0;
                                        for (i=0; i<tv->tolevel; i++) {
                                            temp += Spine[i].listcounter;
                                        }
                                        if (temp > 5) {
                                            tv->newst_stage1 = searchtrie_make(CurrCand, NextCand, n, tv);
                                            EXITFROMSTAGE0EXPATH1
                                        }
                                    }
                                }
                                
                                /* ANY AUTOMORPHISM? */
                                if (tv->options->verbosity >= 2) tv->autchk -= CPUTIME;
                                tv->newindex = 0;
                                if (NextCand->do_it) {
                                    closeloop = CheckForAutomorphisms(g_arg, CurrCand, NextCand, tv, ti, m, n, NextPart);
                                    if (!NextCand->do_it && closeloop < tv->tolevel) k = SpineTL->tgtend;
                                }
                                if (tv->options->verbosity >= 2) tv->autchk += CPUTIME;
                                
                                if (NextCand->do_it) {
                                    ADDTONEXTLEVEL;
                                    SpineTL->keptcounter++;
                                    searchtrie_make(CurrCand, SpineTL->listend, n, tv);
                                }
                            }
                            else {
                                if (BreakSteps[tv->tolevel]) {
                                    NextCand->firstsingcode = NextCand->pathsingcode;
                                    if (CheckForSingAutomorphisms(g_arg, CurrCand, NextPart, NextCand, tv, ti, m, n))
                                        if (!NextCand->do_it) {
                                            PRINT_RETURN
                                            break;
                                        }
                                }
                                
                                /* ANY AUTOMORPHISM? */
                                if (tv->options->verbosity >= 2) tv->autchk -= CPUTIME;
                                tv->newindex = 0;
                                if (NextCand->do_it) {
                                    closeloop = CheckForAutomorphisms(g_arg, CurrCand, NextCand, tv, ti, m, n, NextPart);
                                    if (!NextCand->do_it && closeloop < tv->tolevel) k = SpineTL->tgtend;
                                }
                                if (tv->options->verbosity >= 2) tv->autchk += CPUTIME;
                                
                                if (NextCand->do_it) {
                                    ADDTONEXTLEVEL;
                                    SpineTL->keptcounter++;
                                    searchtrie_make(CurrCand, SpineTL->listend, n, tv);
                                }
                            }
                            PRINT_RETURN
                            break;
                        case 2 :	/* Delete the old list and start a new one: a better trace has been found */
                            if (NextPart->cells == tv->finalnumcells) tv->stats->canupdates++;
                            
                            if (tv->tolevel > tv->treedepth) {
                                tv->treedepth = tv->tolevel;
                                if (tv->strategy) {
                                    NEWPART(SpineTL->part);
                                }
                                else {
                                    NEWPARTSPINE(tv->tolevel);
                                }
                            }
                            
                            if (!tv->strategy && (tv->tolevel > 1) && !SpineTL->liststart) {
                                memcpy(NextCand->lab, TEMPLAB, n*sizeof(int));
                                memcpy(NextCand->invlab, TEMPINVLAB, n*sizeof(int));
                                tv->maxtreelevel = tv->tolevel;
                                
                                SpineTL->liststart = NewCandidate(n, &GarbList, TRUE);
                                SpineTL->listend = SpineTL->liststart;
                                
                                if (NextPart->cells != tv->finalnumcells) NextCand->code = auxcode;
                                CopyCand(SpineTL->liststart, NextCand, n);
                                COPYPART(SpineTL->part, NextPart);
                                tv->newindex = 0;
                                tv->newst_stage1 = searchtrie_make(CurrCand, SpineTL->listend, n, tv);
                                
                                SpineTL->listcounter = 1;
                                SpTLliststart = SpineTL->liststart;
                                
                                i = tv->tolevel;
                                if (tv->brkstpcount) {
                                    while ((i<n) && !BreakSteps[i]) {
                                        i++;
                                    }
                                    if (i<n) SpineTL->liststart->firstsingcode = Spine[i].singcode;
                                }
                                
                                SpineTL->updates = 1;
                                SpineTL->levelcounter = 1;
                                SpineTL->keptcounter = 1;
                                PRINT_LINE_PLUS(tv->fromlevel)
                                
                                if (tv->options->verbosity >= 2) PRINT_CANDIDATE(SpineTL->liststart, tv->tolevel);
                                PRINT_RETURN;
                                
                                if (!tv->strategy && !tv->options->getcanon && (tv->tolevel+1 == tv->firstpathlength)) {
                                    if (tv->tolevel == 1) {
                                        EXITFROMSTAGE0EXPATH2;
                                    }
                                    else {
                                        temp = 0;
                                        for (i=0; i<tv->tolevel; i++) {
                                            temp += Spine[i].listcounter;
                                        }
                                        if (temp > 5) {
                                            EXITFROMSTAGE0EXPATH2;
                                        }
                                    }
                                }
                            }
                            else {
                                tv->tcellevel = tv->maxtreelevel = tv->tolevel;
                                SpineTL->levelcounter++;
                                SpineTL->updates++;
                                SpineTL->keptcounter = 1;
                                
                                RemoveFromLevel(tv->tolevel, tv->treedepth, tv->strategy, TRUE);
                                SpineTL->liststart = NewCandidate(n, &GarbList, TRUE);
                                SpineTL->listend = SpineTL->liststart;
                                
                                CopyCand(SpineTL->liststart, NextCand, n);
                                COPYPART(SpineTL->part, NextPart);
                                tv->newindex = 0;
                                
                                tv->newst_stage1 = searchtrie_make(CurrCand, SpineTL->listend, n, tv);
                                
                                SpineTL->listcounter = 1;
                                SpTLliststart = SpineTL->liststart;
                                
                                SpTLliststart->pathsingcode = SpineTL->singcode = SpTLliststart->singcode;
                                SpTLliststart->firstsingcode = 0;
                                
                                PRINT_LINE
                                if (tv->options->verbosity >= 2) PRINT_CANDIDATE(SpTLliststart, tv->tolevel)
                                    
                                    tv->tolevel_tl = tv->tolevel;
                                
                                memset(BreakSteps, 0, n*sizeof(int));
                                tv->brkstpcount = 0;
                                
                                if (tv->steps > 1) {
                                    if (tv->options->verbosity >= 2) tv->expaths -= CPUTIME;
                                    
                                    /* EXPERIMENTAL PATH */
                                    while (NextPart->cells < tv->finalnumcells) {
                                        if (!TargetCellFirstPath(g_arg, SpTLliststart, NextPart, n, tv)) {
                                            tv->tolevel_tl--;
                                            if (tv->tolevel_tl > 6) {
                                                PRINT_EXPPATHSTEP(SpTLliststart, TRUE)
                                            }
                                            break;
                                        }
                                        ExperimentalStep(g_arg, NextPart, SpTLliststart, tv, ti, m, n);
                                        PRINT_EXPPATHSTEP(SpTLliststart, TRUE)
                                        
                                        Spine[tv->tolevel_tl].singcode = SpTLliststart->pathsingcode;
                                    }
                                    if (NextPart->cells == tv->finalnumcells) {
                                        UPDATEMIN(tv->expathlength, tv->tolevel_tl);
                                    }
                                    
                                    if (tv->options->verbosity >= 2) tv->expaths += CPUTIME;
                                    
                                    if (tv->tolevel_tl > tv->smalldeglevel) {
                                        tv->smalldeglevel = tv->tolevel_tl;
                                    }
                                    if ((tv->finalnumcells < n) && (tv->tolevel_tl < tv->smalldeglevel)) {
                                        tv->smalldeglevel = tv->tolevel_tl;
                                    }
                                    
                                    tv->firstpathlength = tv->tolevel_tl;
                                    PRINT_RETURN
                                    if (!tv->strategy && !tv->options->getcanon && (NextPart->cells == tv->finalnumcells) && (tv->tolevel_tl == tv->tolevel + 1)) {
                                        tv->maxtreelevel = tv->tolevel_tl;
                                        if (tv->tolevel == 1) {
                                            EXITFROMSTAGE0EXPATH2
                                        }
                                        else {
                                            temp = 0;
                                            for (i=0; i<tv->tolevel; i++) {
                                                temp += Spine[i].listcounter;
                                            }
                                            if (temp > 5) {
                                                EXITFROMSTAGE0EXPATH2
                                            }
                                        }
                                    }
                                    memcpy(TEMPLAB, SpTLliststart->lab, n*sizeof(int));
                                    memcpy(TEMPINVLAB, SpTLliststart->invlab, n*sizeof(int));
                                }
                                else {
                                    PRINT_RETURN
                                }
                            }
                            break;
                        default:
                            break;
                    }
                } /* end for */
            }
        }
    }
    
	/* REMOVE CURRENT CANDIDATE */
	if (SpineFL->liststart && (k >= SpineTL->tgtend)) {
		SpineFL->liststart = CurrCand->next;
		if (CurrCand->next == NULL) {
			SpineFL->listend = NULL;
		}
		SpineFL->listcounter--;
		CurrCand->next = GarbList;
		GarbList = CurrCand;
	}
	ti->thereisnextlevel = SelectNextLevel(g_arg, n, tv);
	return;
}

void CompStage1(sparsegraph *g_arg, Partition *CurrPart, Partition *NextPart, Candidate *CurrCand, Candidate *NextCand,
				int m, int n,
				struct TracesVars* tv, struct TracesInfo *ti)
{
	int i, k, cu, cu1, tmp, gom_level, search_vtx;
    searchtrie *TreeNode;
    
    CurrCand->stnode = tv->newst_stage1;
    
    tv->tolevel++;
	SpineTL = Spine+tv->tolevel;
	tv->tcell = SpineTL->tgtcell;
    SpineTL->levelcounter = 0;
    SpineTL->keptcounter = 0;
    SpineTL->updates = 1;
    
	
	if (tv->options->verbosity >= 2) {
		LINE(tv->linelgth-3, "=");
		NEXTLINE
	}
	
	memset(RefCells, 0, n*sizeof(int));
	memset(MultRefCells, 0, n*sizeof(int));
	ti->thegrouphaschanged = TRUE;
	
	/*  CANDIDATE */
	memcpy(NextCand->lab, CurrCand->lab, n*sizeof(int));
	memcpy(NextCand->invlab, CurrCand->invlab, n*sizeof(int));
	NextCand->do_it = TRUE;
	SpineTL->trcstart = CurrPart->cells;
	
	tv->indivstart = tv->tcell;
	tv->indivend = SpineTL->tgtend;
	if (g_arg->d[CurrCand->lab[tv->indivstart]] == 1) {
		tv->indivstart = SpineTL->tgtend-1;
	}
    
    FixBase(fix, tv, NextCand, 0, tv->fromlevel);
    
	if (!ti->identitygroup) {
        if (tv->options->verbosity >= 2) tv->schreier2 -= CPUTIME;
        tv->currorbit = getorbits(fix, tv->nfix, gpB, &gensB, n);
        if (tv->options->verbosity >= 2) tv->schreier2 += CPUTIME;
    }
    else {
        if (n / CurrPart->cls[tv->tcell] < 256) {
            memcpy(tv->currorbit, IDENTITY_PERM, n*sizeof(int));
        }
        else {
            for (k = tv->indivstart; k < tv->indivend; k++) {
                tv->currorbit[CurrCand->lab[k]] = CurrCand->lab[k];
            }
        }
    }

    if (!CurrCand->sortedlab) {
        quickSort(CurrCand->lab+tv->tcell, CurrPart->cls[tv->tcell]);
        for (i=tv->tcell; i<tv->tcell+CurrPart->cls[tv->tcell]; i++) {
            CurrCand->invlab[CurrCand->lab[i]] = i;
        }
        CurrCand->sortedlab = TRUE;
    }
    
	for (k = tv->indivstart; k < tv->indivend; k++) {
        NextCand->vertex = CurrCand->lab[k];
        NextCand->name = ++tv->name;
        if (NextCand->name == (NAUTY_INFINITY-2)) {
            NextCand->name = tv->name = 1;
        }
        if (tv->currorbit[CurrCand->lab[k]] != CurrCand->lab[k]) {
            search_vtx = tv->currorbit[NextCand->vertex];
            TreeNode = CurrCand->stnode;
            if (TreeNode->first_child) {
                TreeNode = TreeNode->first_child;
                while (TreeNode) {
                    if (TreeNode->vtx == search_vtx) {
                        break;
                    }
                    TreeNode = TreeNode->next_sibling;
                }
                if (TreeNode) {
                    while (TreeNode->goes_to) {
                        TreeNode = TreeNode->goes_to;
                    }
                    TreeNode->index++;
                }
            }
			continue;
		}
        PRINT_REFINE_VERB(3)
		memcpy(NextPart->cls, CurrPart->cls, n*sizeof(int));
		memcpy(NextPart->inv, CurrPart->inv, n*sizeof(int));
		
		Individualize(NextPart, NextCand, CurrCand->lab[k], tv->tcell, CurrPart->cells, SpineTL->tgtpos);
		
		tv->stats->numnodes++;
        SpineTL->levelcounter++;
        tv->tolevel_tl = tv->tolevel;
		trieref = trieroot;
        SpineTL->levelcounter++;
        
		traces_refine_maketrie(g_arg,
							   NextCand,
							   m,
							   n,
							   NextPart, tv, ti);
		
		RefCells[CurrCand->lab[k]] = NextPart->cells;
		if (NextPart->cells == tv->finalnumcells) {
			tv->answ=1;
			if (tv->options->verbosity >= 2) PRINT_CANDIDATE(NextCand, tv->tolevel)
                
                if (NextPart->cells < tv->finalnumcells) {
                    if (!Spine[tv->tolevel].part) {
                        NEWPART(Spine[tv->tolevel].part)
                    }
                }
            
			/* ANY AUTOMORPHISM? */
			if (tv->options->verbosity >= 2) tv->autchk -= CPUTIME;
            if (NextPart->cells == tv->finalnumcells) {
                CheckForAutomorphisms(g_arg, CurrCand, NextCand, tv, ti, m, n, NextPart);
            }
            else {
                CheckForSingAutomorphisms(g_arg, CurrCand, NextPart, NextCand, tv, ti, m, n);
            }
            
			if (tv->options->verbosity >= 2) tv->autchk += CPUTIME;
			
			PRINT_RETURN
			
			/* ADD TO NEXT LEVEL */
			if (tv->answ == 1) {
                SpineTL->keptcounter++;
                if (!Spine[tv->tolevel].listend) COPYPART(Spine[tv->tolevel].part, NextPart);
				ADDTONEXTLEVEL;
                searchtrie_make(CurrCand, SpineTL->listend, n, tv);
			}
		}
	} /* end for */
	
	for (k = tv->indivstart; k < tv->indivend; k++) {
		MultRefCells[RefCells[tv->currorbit[CurrCand->lab[k]]] % tv->finalnumcells]++;
	}
	if (tv->options->verbosity >= 2) {
        for (k=1; k<n; k++) {
			if (MultRefCells[k]) {
				fprintf(outfile, tv->digstring, k);
				fprintf(outfile, "cells: %d\n", MultRefCells[k]);
			}
		}
	}
	
#if !MAXN
    DYNALLOC1(searchtrie*, RefPath, RefPath_sz, tv->tolevel, "Traces-CS1");
#endif
    
    TreeNode = CurrCand->stnode;
    while (TreeNode) {
        RefPath[TreeNode->level] = TreeNode;
        TreeNode = TreeNode->father;
    }
    
	/* REMOVE CURRENT CANDIDATE */
	SpineFL->liststart = CurrCand->next;
	if (CurrCand->next == NULL) {
		SpineFL->listend = NULL;
		SpineFL->listcounter = 1;
	}
	SpineFL->listcounter--;
	CurrCand->next = GarbList;
	GarbList = CurrCand;
	
	if (tv->options->verbosity >= 2) {
		LINE(tv->linelgth, "=");
		NEXTLINE
	}
	tv->compstage = 2;
	tv->steps = n;
    
	if (tv->options->verbosity >= 2) tv->schreier1 -= CPUTIME;
	gom_level = getorbitsmin(fix, tv->nfix, gpB, &gensB, &tv->currorbit,
							 CurrCand->lab+tv->tcell, CurrPart->cls[tv->tcell], n, TRUE);
	if (tv->options->verbosity >= 2) tv->schreier1 += CPUTIME;
	ORBITSIZES
	ti->thereisnextlevel = SelectNextLevel(g_arg, n, tv);
	SpineTL->part->cells = tv->finalnumcells;
    
    AutomCount[0] = 2;
    AutomCount[1] = CurrCand->vertex;
    
    return;
}

void CompStage2(sparsegraph *g_arg, Partition *CurrPart, Partition *NextPart, Candidate *CurrCand, Candidate *NextCand,
				int m, int n,
				struct TracesVars* tv, struct TracesInfo *ti)
{
	int i, j, i1, j2, k, cu, cu1, vertex, search_vtx, gom_level;
	int temp, tmp, autom;
	Candidate *AuxCand;
    searchtrie *TreeNode, *TreeNode1, *TreeNode2;
    int *CuOrb;
    boolean schreierwrong;
	
	autom = 0;
    schreierwrong = FALSE;
	
    TreeNode = CurrCand->stnode;
    tv->cand_level = 0;
    
    while (TreeNode) {
        if (TreeNode->goes_to) {
            CurrCand->do_it = FALSE;
        }
        if (!tv->cand_level && TreeNode == RefPath[TreeNode->level]) {
            tv->cand_level = TreeNode->level;
        }
        TreeNode = TreeNode->father;
    }
    if (tv->cand_level+1 == tv->maxtreelevel) {
        ti->useTempOrbits1 = TRUE;
    }
    else {
        ti->useTempOrbits1 = FALSE;
    }
    if (tv->cand_level == tv->fromlevel) {
        ti->useTempOrbits2 = TRUE;
    }
    else {
        ti->useTempOrbits2 = FALSE;
    }
    
    PRINT_FROM_VERB(3)
    
    if (CurrCand->do_it) {
        if (tv->tolevel == 0) {
            tv->fromlevel = tv->tolevel;
            SpineFL = Spine+tv->fromlevel;
            vertex = Spine[tv->maxtreelevel+1].liststart->lab[Spine[1].tgtpos];
            k = n;
            
            if (TargetCell(g_arg, CurrCand, CurrPart, n, tv, tv->tolevel)) {
                ++tv->tolevel;
                SpineTL = Spine+tv->tolevel;
                SpineTL->tgtcell = tv->tcell;
                SpineTL->tgtsize = CurrPart->cls[tv->tcell];
                SpineTL->tgtend = tv->tcell+SpineTL->tgtsize;
                SpineTL->tgtpos = SpineTL->tgtend - 1;
            }
            else {
                tv->smalldeglevel = tv->tolevel;
                tv->finalnumcells = CurrPart->cells;
                return;
            }
            
            memcpy(NextCand->lab, CurrCand->lab, n*sizeof(int));
            memcpy(NextCand->invlab, CurrCand->invlab, n*sizeof(int));
            SpineTL->trcstart = CurrPart->cells;
            TheTrace[SpineTL->trcstart] = SpineTL->tgtpos;
            
            tv->indivstart = tv->tcell+CurrCand->indnum;
            tv->indivend = tv->indivstart+tv->steps;
            if (tv->indivend > SpineTL->tgtend) {
                tv->indivend = SpineTL->tgtend;
            }
            memset(CurrRefCells, 0, n*sizeof(int));
            ti->thegrouphaschanged = TRUE;
            
            if (!CurrCand->sortedlab) {
                quickSort(CurrCand->lab+tv->tcell, CurrPart->cls[tv->tcell]);
                for (i=tv->tcell; i<tv->tcell+CurrPart->cls[tv->tcell]; i++) {
                    CurrCand->invlab[CurrCand->lab[i]] = i;
                }
                CurrCand->sortedlab = TRUE;
            }
            
            for (k = tv->indivstart; k < tv->indivend; k++) {
                if ((tv->orbits[CurrCand->lab[k]] == CurrCand->lab[k]) && (OrbSize[tv->orbits[CurrCand->lab[k]]] >= OrbSize[tv->orbits[vertex]])) {
                    
                    CurrCand->indnum++;
                    NextCand->singcode = CurrCand->singcode;
                    NextCand->vertex = CurrCand->lab[k];
                    NextCand->name = ++tv->name;
                    if (NextCand->name == (NAUTY_INFINITY-2)) {
                        NextCand->name = tv->name = 1;
                    }
                    
                    if (ti->thegrouphaschanged) {
                        if (tv->fromlevel == tv->maxtreelevel) {
                            CURRORBITSIZES
                        }
                        ti->thegrouphaschanged = FALSE;
                    }
                    
                    if (tv->currorbit[CurrCand->lab[k]] != CurrCand->lab[k]) {
                        continue;
                    }
                    
                    memcpy(NextPart->cls, CurrPart->cls, n*sizeof(int));
                    memcpy(NextPart->inv, CurrPart->inv, n*sizeof(int));
                    if (NextPart->cls[tv->tcell] == 2) {
                        NextCand->singcode = MASHCOMM(NextCand->singcode, CurrCand->lab[tv->tcell]);
                        NextCand->singcode = MASHCOMM(NextCand->singcode, CurrCand->lab[tv->tcell+1]);
                    }
                    else {
                        NextCand->singcode = MASHCOMM(NextCand->singcode, CurrCand->lab[k]+labelorg);
                    }
                    
                    Individualize(NextPart, NextCand, CurrCand->lab[k], tv->tcell, CurrPart->cells, SpineTL->tgtpos);
                    
                    tv->stats->numnodes++;
                    Spine[tv->tolevel+1].levelcounter++;
                    if (tv->fromlevel == tv->maxtreelevel) {
                        tv->tolevel_tl = tv->tolevel;
                        trieref = trieroot;
                        
                        tv->answ = traces_refine_comptrie(g_arg,
                                                          NextCand,
                                                          m,
                                                          n,
                                                          NextPart, tv, ti);
                        
                        if (tv->answ) {
                            if (NextPart->cells != n) {
                                CurrRefCells[NextPart->cells] += CurrOrbSize[CurrCand->lab[k]];
                                if (CurrRefCells[NextPart->cells] > MultRefCells[NextPart->cells]) {
                                    k = n;
                                    break;
                                }
                                continue;
                            }
                            if (tv->options->verbosity >= 2) PRINT_CANDIDATE(NextCand, tv->tolevel)
                                }
                    }
                    else {
                        tv->answ = traces_refine_sametrace(g_arg,
                                                           NextCand,
                                                           m,
                                                           n,
                                                           NextPart, tv, ti);
                        if (tv->answ) {
                            if (tv->options->verbosity >= 2) PRINT_CANDIDATE(NextCand, tv->tolevel)
                                if (tv->tolevel == tv->maxtreelevel) {
                                    tv->tolevel_tl = tv->tolevel;
                                    if (tv->options->verbosity >= 2) tv->expaths -= CPUTIME;
                                    TargetCellExpPath(g_arg, NextCand, NextPart, n, tv);
                                    ExperimentalStep(g_arg, NextPart, NextCand, tv, ti, m, n);
                                    PRINT_EXPPATHSTEP(NextCand, tv->answ)
                                    if (NextPart->cells == tv->finalnumcells) {
                                        UPDATEMIN(tv->expathlength, tv->tolevel_tl);
                                    }
                                    
                                    if (tv->options->verbosity >= 2) tv->expaths += CPUTIME;
                                    if (!tv->answ) {
                                        PRINT_RETURN
                                    }
                                }
                        }
                    }
                    if (tv->answ) {
                        if (NextPart->cells == tv->finalnumcells) {
                            if (tv->options->verbosity >= 2) tv->autchk -= CPUTIME;
                            temp = (tv->tolevel_tl == tv->tolevel+1);
                            autom = CheckForAutomorphisms(g_arg, CurrCand, NextCand,
                                                          tv, ti, temp, n, NextPart);
                            if (tv->options->verbosity >= 2) tv->autchk += CPUTIME;
                            
                            if (ti->thegrouphaschanged) {
                                ORBITSIZES
                            }
                        }
                        PRINT_RETURN
                        
                        /* ADD TO NEXT LEVEL */
                        if ((NextPart->cells < tv->finalnumcells) || (tv->tolevel != tv->maxtreelevel) || (tv->tolevel_tl != tv->tolevel+1)) {
                            ADDTONEXTLEVEL;
                            searchtrie_make(CurrCand, SpineTL->listend, n, tv);
                        }
                    }
                    else {
                        tv->stats->interrupted++;
                    }
                    if (tv->fromlevel == tv->maxtreelevel) {
                        k = n;
                        break;
                    }
                }
            } /* end for */
        }
        else {
            temp = CurrCand->lab[Spine[1].tgtpos];
            vertex = Spine[tv->maxtreelevel+1].liststart->lab[Spine[1].tgtpos];
            k = n;
            
            if (tv->cand_level || ((tv->orbits[temp] == temp) && (OrbSize[tv->orbits[CurrCand->lab[Spine[1].tgtpos]]] >= OrbSize[tv->orbits[vertex]])))
            {
                tv->fromlevel = tv->tolevel;
                SpineFL = Spine+tv->fromlevel;
                
                if (TargetCell(g_arg, CurrCand, CurrPart, n, tv, tv->tolevel)) {
                    tv->tcellevel = ++tv->tolevel;
                    SpineTL = Spine+tv->tolevel;
                    SpineTL->tgtcell = tv->tcell;
                    SpineTL->tgtsize = CurrPart->cls[tv->tcell];
                    SpineTL->tgtend = tv->tcell+SpineTL->tgtsize;
                    SpineTL->tgtpos = SpineTL->tgtend - 1;
                }
                else {
                    tv->smalldeglevel = tv->tolevel;
                    tv->finalnumcells = CurrPart->cells;
                    return;
                }
                ti->minimalinorbits = TRUE;
                
                if (!ti->identitygroup) {
                    
                    if (ti->useTempOrbits1 && ti->useTempOrbits2) {
                        CuOrb = TempOrbits;
                    }
                    else {
                        FixBase(fix, tv, CurrCand, 0, tv->fromlevel);
                        if (ti->useTempOrbits1 && tv->fromlevel == tv->maxtreelevel) {
                            tv->currorbit = getorbits(fix, tv->nfix, gpB, &gensB, n);
                            CuOrb = tv->currorbit;
                        } else {
                            if (tv->options->verbosity >= 2) tv->schreier1 -= CPUTIME;
                            gom_level = getorbitsmin(fix, tv->nfix, gpB, &gensB, &tv->currorbit,
                                                     CurrCand->lab+tv->tcell, CurrPart->cls[tv->tcell], n, TRUE);
                            if (tv->options->verbosity >= 2) tv->schreier1 += CPUTIME;
                            CuOrb = tv->currorbit;
                            if (gom_level < tv->nfix) {
                                PRINT_NOTMIN_VERB(3)
                                if (ti->useTempOrbits1) {
                                    for (i=1; i<AutomCount[0]; i++) if (AutomCount[i] == CuOrb[tv->currorbit[CurrCand->vertex]]) break;
                                    if (i < AutomCount[0]) {
                                        AutomCount[AutomCount[0]++] = CurrCand->vertex;
                                    }
                                    ti->minimalinorbits = FALSE;
                                }
                                else {
                                    TreeNode = CurrCand->stnode;
                                    j2 = CurrCand->lab[Spine[gom_level+1].tgtpos];
                                    i1 = tv->currorbit[j2];
                                    for (j=0; j < tv->nfix - gom_level; j++) {
                                        TreeNode = TreeNode->father;
                                    }
                                    TreeNode1 = TreeNode->first_child;
                                    while (TreeNode1) {
                                        if (TreeNode1->vtx == i1) {
                                            break;
                                        }
                                        TreeNode1 = TreeNode1->next_sibling;
                                    }
                                    schreierwrong = FALSE;
                                    if (TreeNode1) {
                                        while (TreeNode1->goes_to) {
                                            TreeNode1 = TreeNode1->goes_to;
                                        }
                                        TreeNode2 = TreeNode->first_child;
                                        while (TreeNode2->vtx != j2) {
                                            TreeNode2 = TreeNode2->next_sibling;
                                        }
                                        
                                        TreeNode1->index += TreeNode2->index;
                                        TreeNode2->goes_to = TreeNode1;
                                        
                                        ti->minimalinorbits = FALSE;
                                    }
                                    else {
                                        tv->currorbit = getorbits(fix, tv->nfix, gpB, &gensB, n);
                                        schreierwrong = TRUE;
                                    }
                                }
                            }
                        }
                    }
                    ti->thegrouphaschanged = FALSE;
                }
                else {
                    CuOrb = IDENTITY_PERM;
                }
                
                if (ti->minimalinorbits) {
                    memcpy(NextCand->lab, CurrCand->lab, n*sizeof(int));
                    memcpy(NextCand->invlab, CurrCand->invlab, n*sizeof(int));
                    SpineTL->trcstart = CurrPart->cells;
                    TheTrace[SpineTL->trcstart] = SpineTL->tgtpos;
                    
                    tv->indivstart = tv->tcell+CurrCand->indnum;
                    tv->indivend = tv->indivstart+tv->steps;
                    if (tv->indivend > SpineTL->tgtend) {
                        tv->indivend = SpineTL->tgtend;
                    }
                    memset(CurrRefCells, 0, n*sizeof(int));
                    if (!ti->identitygroup) ti->thegrouphaschanged = TRUE;
                    
                    if (!CurrCand->sortedlab) {
                        quickSort(CurrCand->lab+tv->tcell, CurrPart->cls[tv->tcell]);
                        for (i=tv->tcell; i<tv->tcell+CurrPart->cls[tv->tcell]; i++) {
                            CurrCand->invlab[CurrCand->lab[i]] = i;
                        }
                        CurrCand->sortedlab = TRUE;
                    }
                    
                    for (k = tv->indivstart; k < tv->indivend; k++) {
                        CurrCand->indnum++;
                        NextCand->singcode = CurrCand->singcode;
                        NextCand->vertex = CurrCand->lab[k];
                        NextCand->name = ++tv->name;
                        if (NextCand->name == (NAUTY_INFINITY-2)) {
                            NextCand->name = tv->name = 1;
                        }
                        
                        if (ti->thegrouphaschanged) {
                            if (tv->fromlevel == tv->maxtreelevel) {
                                CURRORBITSIZES
                            }
                            ti->thegrouphaschanged = FALSE;
                        }
                        
                        if (!schreierwrong) {
                            if (CuOrb[CurrCand->lab[k]] != CurrCand->lab[k]) {
                                continue;
                            }
                        }
                        
                        memcpy(NextPart->cls, CurrPart->cls, n*sizeof(int));
                        memcpy(NextPart->inv, CurrPart->inv, n*sizeof(int));
                        if (NextPart->cls[tv->tcell] == 2) {
                            NextCand->singcode = MASHCOMM(NextCand->singcode, CurrCand->lab[tv->tcell]);
                            NextCand->singcode = MASHCOMM(NextCand->singcode, CurrCand->lab[tv->tcell+1]);
                        }
                        else {
                            NextCand->singcode = MASHCOMM(NextCand->singcode, CurrCand->lab[k]);
                        }
                        
                        Individualize(NextPart, NextCand, CurrCand->lab[k], tv->tcell, CurrPart->cells, SpineTL->tgtpos);
                        
                        tv->stats->numnodes++;
                        Spine[tv->tolevel+1].levelcounter++;
                        if (tv->fromlevel == tv->maxtreelevel) {
                            tv->tolevel_tl = tv->tolevel;
                            trieref = trieroot;
                            
                            tv->answ = traces_refine_comptrie(g_arg,
                                                              NextCand,
                                                              m,
                                                              n,
                                                              NextPart, tv, ti);
                            if (tv->answ) {
                                if (NextPart->cells != tv->finalnumcells) {
                                    CurrRefCells[NextPart->cells] += CurrOrbSize[CurrCand->lab[k]];
                                    if (CurrRefCells[NextPart->cells] > MultRefCells[NextPart->cells]) {
                                        k = n;
                                        break;
                                    }
                                    continue;
                                }
                                if (tv->options->verbosity >= 2) PRINT_CANDIDATE(NextCand, tv->tolevel);
                            }
                        }
                        else {
                            tv->answ = traces_refine_sametrace(g_arg,
                                                               NextCand,
                                                               m,
                                                               n,
                                                               NextPart, tv, ti);
                            if (tv->answ) {
                                if (tv->options->verbosity >= 2) PRINT_CANDIDATE(NextCand, tv->tolevel)
                                    if (tv->tolevel == tv->maxtreelevel) {
                                        tv->tolevel_tl = tv->tolevel;
                                        if (tv->options->verbosity >= 2) tv->expaths -= CPUTIME;
                                        if (TargetCellExpPath(g_arg, NextCand, NextPart, n, tv)) {
                                            ExperimentalStep(g_arg, NextPart, NextCand, tv, ti, m, n);
                                            PRINT_EXPPATHSTEP(NextCand, tv->answ)
                                            if (NextPart->cells == tv->finalnumcells) {
                                                UPDATEMIN(tv->expathlength, tv->tolevel_tl);
                                            }
                                            
                                        }
                                        if (tv->options->verbosity >= 2) tv->expaths += CPUTIME;
                                        if (!tv->answ) {
                                            PRINT_RETURN
                                        }
                                    }
                            }
                        }
                        if (tv->answ) {
                            if (NextPart->cells == tv->finalnumcells) {
                                if (tv->options->verbosity >= 2) tv->autchk -= CPUTIME;
                                temp = (tv->tolevel_tl == tv->tolevel+1);
                                autom = CheckForAutomorphisms(g_arg, CurrCand, NextCand,
                                                              tv, ti, temp, n, NextPart);
                                if (tv->options->verbosity >= 2) tv->autchk += CPUTIME;
                                if (autom) {
                                    for (i=autom; i<=tv->maxtreelevel; i++) {
                                        AuxCand = Spine[i].liststart;
                                        while (AuxCand && Prefix(AuxCand, NextCand, autom)) {
                                            AuxCand->do_it = FALSE;
                                            AuxCand = AuxCand->next;
                                        }
                                    }
                                    if (autom == tv->tolevel) {
                                        autom = 0;
                                    }
                                }
                                
                                if (ti->thegrouphaschanged) {
                                    ORBITSIZES
                                }
                            }
                            PRINT_RETURN
                            
                            /* ADD TO NEXT LEVEL */
                            if ((NextPart->cells < tv->finalnumcells) || (tv->tolevel != tv->maxtreelevel) || (tv->tolevel_tl != tv->tolevel+1)) {
                                ADDTONEXTLEVEL;
                                searchtrie_make(CurrCand, SpineTL->listend, n, tv);
                            }
                        }
                        else {
                            tv->stats->interrupted++;
                        }
                        if (autom) {
                            k = n;
                            autom = 0;
                            break;
                        }
                        if (tv->fromlevel == tv->maxtreelevel) {
                            k = n;
                            break;
                        }
                    } /* end for */
                    TreeNode = RefPath[tv->maxtreelevel];
                    search_vtx = TreeNode->vtx;
                }
            }
            else SpineTL = &Spine[tv->tolevel+1];
        }
        
    }
    
	/* REMOVE CURRENT CANDIDATE */
	if (!CurrCand->do_it || k >= SpineTL->tgtend) {
		SpineFL->liststart = CurrCand->next;
		if (CurrCand->next == NULL) {
			SpineFL->listend = NULL;
		}
		CurrCand->next = GarbList;
		GarbList = CurrCand;
	}
    ti->thereisnextlevel = SelectNextLevel(g_arg, n, tv);
    return;
}

void
grouporderplus(sparsegraph *sg, int *fix, int nfix, Candidate *Cand,
               schreier *gp, permnode **ring,
			   double *grpsize1, int *grpsize2, int n, TracesVars *tv, TracesInfo *ti) {
    
	int i, i1, j, j0, j2, k, k1, k2, w, w1, w2, c, c1, c2, c3, c4, n1, n2;
    int prev, step, start, counts, StInd, CyInd, cycnum, orbrep, cellord, LghAttrInd;
	int tmp, temp, halfsize, vtx, cell1;
    int *cylab, *cycol;
	TracesSpine SpineSDegLev;
    searchtrie *TrieNode;
    
    *grpsize1 = 1.0; *grpsize2 = 0;
    
    TrieNode = Spine[tv->maxtreelevel].liststart->stnode;
    if (tv->options->verbosity >= 2) {
        fprintf(outfile, "-->> ");
        while (TrieNode->father) {
            fprintf(outfile, "%d (%d), ", TrieNode->name, TrieNode->index);
            if (TrieNode->father->name) {
                MULTIPLY(tv->stats->grpsize1, tv->stats->grpsize2, TrieNode->index);
            } else {
                temp = spinelementorbsize(tv->orbits, Spine[tv->maxtreelevel].liststart->lab+Spine[1].tgtcell, Spine[1].tgtsize, TrieNode->vtx);
                MULTIPLY(*grpsize1, *grpsize2, temp);
                fprintf(outfile, "orbcount vtx %d from cell %d (%d): %d\n",
                        TrieNode->vtx+labelorg, Spine[1].tgtcell, Spine[1].tgtsize, temp);
            }
            TrieNode = TrieNode->father;
        }
    }
    else {
        while (TrieNode->father) {
            if (TrieNode->father->name) {
                MULTIPLY(tv->stats->grpsize1, tv->stats->grpsize2, TrieNode->index);
            } else {
                temp = spinelementorbsize(tv->orbits, Spine[tv->maxtreelevel].liststart->lab+Spine[1].tgtcell, Spine[1].tgtsize, TrieNode->vtx);
                MULTIPLY(*grpsize1, *grpsize2, temp);
            }
            TrieNode = TrieNode->father;
        }
    }
    
    if (tv->smalldeglevel < n) {
        SpineSDegLev = Spine[tv->smalldeglevel];
        
        /* Deg 2 and at least one nghb with deg > 2 */
        SETMARK(StackMarkers, tv->stackmark)
        for (i=0; i<n; i += SpineSDegLev.part->cls[i]) {
            SETMARK(Markers, tv->mark)
            if (SpineSDegLev.part->cls[i] > 1) {
                tmp = SpineSDegLev.liststart->lab[i];
                /* tmp has deg 2 and at least one nghb with deg > 2 */
                if ((sg->d[tmp] == 2) && ((sg->d[sg->e[sg->v[tmp]]] > 2) || ((sg->d[sg->e[sg->v[tmp]+1]] > 2)))) {
                    n1 = sg->e[sg->v[tmp]];
                    n2 = sg->e[sg->v[tmp]+1];
                    if (sg->d[n1] > 2) {
                        if (sg->d[n2] > 2) {
                            if (Cand->invlab[n1] < Cand->invlab[n2]) {
                                start = n1;
                            }
                            else {
                                start = n2;
                            }
                        }
                        else {
                            start = n1;
                        }
                    }
                    else {
                        start = n2;
                    }
                    counts = 0;
                    StInd = 0;
                    for (j=i; j<i+SpineSDegLev.part->cls[i]; j++) {
                        step = SpineSDegLev.liststart->lab[j];
                        if (Markers[step] != tv->mark) {
                            prev = start;
                            counts++;
                            do {
                                Markers[step] = tv->mark;
                                PERMSTACK[StInd++] = step;
                                if (sg->e[sg->v[step]] != prev) {
                                    prev = step;
                                    step = sg->e[sg->v[step]];
                                } else {
                                    prev = step;
                                    step = sg->e[sg->v[step]+1];
                                }
                            } while (sg->d[step] == 2);
                            if (sg->d[step] == 1) {
                                PERMSTACK[StInd++] = step;
                            }
                        }
                    }
                    if (counts == SpineSDegLev.part->cls[i]) {
                        factorial(grpsize1, grpsize2, SpineSDegLev.part->cls[i]);
                        if (tv->options->generators || tv->options->writeautoms || tv->options->userautomproc) {
                            memcpy(AUTPERM, IDENTITY_PERM, n*sizeof(int));
                            for (k=0; k<StInd/counts; k++) {
                                i1 = PERMSTACK[k];
                                for (j0=0; j0<counts-1; j0++) {
                                    AUTPERM[PERMSTACK[j0*(StInd/counts)+k]] = PERMSTACK[(j0+1)*(StInd/counts)+k];
                                }
                                AUTPERM[PERMSTACK[j0*(StInd/counts)+k]] =i1;
                            }
                        }
                        SPECIALGENERATORS
                        if (counts > 2) {
                            if (tv->options->generators || tv->options->writeautoms || tv->options->userautomproc) {
                                memcpy(AUTPERM, IDENTITY_PERM, n*sizeof(int));
                                for (k=0; k<StInd/counts; k++) {
                                    i1 = PERMSTACK[k];
                                    for (j0=0; j0<1; j0++) {
                                        AUTPERM[PERMSTACK[j0*(StInd/counts)+k]] = PERMSTACK[(j0+1)*(StInd/counts)+k];
                                    }
                                    AUTPERM[PERMSTACK[j0*(StInd/counts)+k]] =i1;
                                }
                            }
                            SPECIALGENERATORS
                        }
                        /* the next loop for orbits and canonical label */
                        for (k=0; k<StInd/counts; k++) {
                            orbrep = n;
                            for (j0=0; j0<counts; j0++) {
                                tmp = PERMSTACK[j0*(StInd/counts)+k];
                                if (tv->orbits[tmp] < orbrep) {
                                    orbrep = tv->orbits[tmp];
                                }
                            }
                            cellord = SpineSDegLev.part->inv[SpineSDegLev.liststart->invlab[PERMSTACK[k]]];
                            StackMarkers[orbrep] = tv->stackmark;
                            for (j0=0; j0<counts; j0++) {
                                tmp = PERMSTACK[j0*(StInd/counts)+k];
                                SpineSDegLev.liststart->lab[cellord+j0] = tmp;
                                SpineSDegLev.liststart->invlab[tmp] = cellord+j0;
                                if (StackMarkers[tv->orbits[tmp]] != tv->stackmark) {
                                    tv->stats->numorbits--;
                                    StackMarkers[tv->orbits[tmp]] = tv->stackmark;
                                }
                                tv->orbits[tmp] = orbrep;
                            }
                        }
                    }
                    else {
                        factorial2(grpsize1, grpsize2, SpineSDegLev.part->cls[i]);
                        for (j=0; j<counts; j++) {
                            j0 = j*(StInd/counts);
                            k1 = (j+1)*(StInd/counts);
                            if (tv->options->generators || tv->options->writeautoms || tv->options->userautomproc) {
                                memcpy(AUTPERM, IDENTITY_PERM, n*sizeof(int));
                                for (k=j0, i1=k1-1; k<k1; k++, i1--) {
                                    AUTPERM[PERMSTACK[k]] = PERMSTACK[i1];
                                }
                            }
                            SPECIALGENERATORS
                        }
                        if (counts > 1) {
                            if (tv->options->generators || tv->options->writeautoms || tv->options->userautomproc) {
                                memcpy(AUTPERM, IDENTITY_PERM, n*sizeof(int));
                                for (k=0; k<StInd/counts; k++) {
                                    i1 = PERMSTACK[k];
                                    for (j0=0; j0<counts-1; j0++) {
                                        AUTPERM[PERMSTACK[j0*(StInd/counts)+k]] = PERMSTACK[(j0+1)*(StInd/counts)+k];
                                    }
                                    AUTPERM[PERMSTACK[j0*(StInd/counts)+k]] =i1;
                                }
                            }
                            SPECIALGENERATORS
                        }
                        if (counts > 2) {
                            if (tv->options->generators || tv->options->writeautoms || tv->options->userautomproc) {
                                memcpy(AUTPERM, IDENTITY_PERM, n*sizeof(int));
                                for (k=0; k<StInd/counts; k++) {
                                    i1 = PERMSTACK[k];
                                    for (j0=0; j0<1; j0++) {
                                        AUTPERM[PERMSTACK[j0*(StInd/counts)+k]] = PERMSTACK[(j0+1)*(StInd/counts)+k];
                                    }
                                    AUTPERM[PERMSTACK[j0*(StInd/counts)+k]] =i1;
                                }
                            }
                            SPECIALGENERATORS
                        }
                        /* the next loop for orbits and canonical label */
                        i1 = StInd/counts;
                        for (j=0; j<i1/2; j++) {
                            orbrep = n;
                            for (k=0, k1=i1-1; k<StInd; k+=i1, k1+=i1) {
                                tmp = PERMSTACK[k+j];
                                temp = PERMSTACK[k1-j];
                                if (tv->orbits[tmp] < orbrep) {
                                    orbrep = tv->orbits[tmp];
                                }
                                if (tv->orbits[temp] < orbrep) {
                                    orbrep = tv->orbits[temp];
                                }
                            }
                            cellord = SpineSDegLev.part->inv[SpineSDegLev.liststart->invlab[PERMSTACK[j]]];
                            StackMarkers[orbrep] = tv->stackmark;
                            for (k=0, k1=i1-1; k<StInd; k+=i1, k1+=i1) {
                                tmp = PERMSTACK[k+j];
                                SpineSDegLev.liststart->lab[cellord] = tmp;
                                SpineSDegLev.liststart->invlab[tmp] = cellord++;
                                if (StackMarkers[tv->orbits[tmp]] != tv->stackmark) {
                                    tv->stats->numorbits--;
                                    StackMarkers[tv->orbits[tmp]] = tv->stackmark;
                                }
                                tv->orbits[tmp] = orbrep;
                                
                                temp = PERMSTACK[k1-j];
                                SpineSDegLev.liststart->lab[cellord] = temp;
                                SpineSDegLev.liststart->invlab[temp] = cellord++;
                                if (StackMarkers[tv->orbits[temp]] != tv->stackmark) {
                                    tv->stats->numorbits--;
                                    StackMarkers[tv->orbits[temp]] = tv->stackmark;
                                }
                                tv->orbits[temp] = orbrep;
                            }
                        }
                        if (i1 % 2) {
                            orbrep = n;
                            for (k=j; k<StInd; k+=i1) {
                                tmp = PERMSTACK[k];
                                if (tv->orbits[tmp] < orbrep) {
                                    orbrep = tv->orbits[tmp];
                                }
                            }
                            cellord = SpineSDegLev.part->inv[SpineSDegLev.liststart->invlab[PERMSTACK[j]]];
                            StackMarkers[orbrep] = tv->stackmark;
                            j0 = 0;
                            for (k=j; k<StInd; k+=i1) {
                                tmp = PERMSTACK[k];
                                SpineSDegLev.liststart->lab[cellord+j0] = tmp;
                                SpineSDegLev.liststart->invlab[tmp] = cellord+j0++;
                                if (StackMarkers[tv->orbits[tmp]] != tv->stackmark) {
                                    tv->stats->numorbits--;
                                    StackMarkers[tv->orbits[tmp]] = tv->stackmark;
                                }
                                tv->orbits[tmp] = orbrep;
                            }
                        }
                    }
                    Propagate(sg, SpineSDegLev.liststart, SpineSDegLev.part, i, i+SpineSDegLev.part->cls[i], n);
                }
            }
        }
        
        /* Deg 2 and at least one nghb with deg = 1 */
        SETMARK(Markers, tv->mark)
        for (i=0; i<n; i += SpineSDegLev.part->cls[i]) {
            if (SpineSDegLev.part->cls[i] > 1) {
                tmp = SpineSDegLev.liststart->lab[i];
                /* tmp has deg 2 and at least one nghb with == 1 */
                if ((sg->d[tmp] == 2) && ((sg->d[sg->e[sg->v[tmp]]] == 1) || ((sg->d[sg->e[sg->v[tmp]+1]] == 1)))) {
                    counts = 0;
                    StInd = 0;
                    for (j=i; j<i+SpineSDegLev.part->cls[i]; j++) {
                        step = SpineSDegLev.liststart->lab[j];
                        if (Markers[step] != tv->mark) {
                            n1 = sg->e[sg->v[step]];
                            n2 = sg->e[sg->v[step]+1];
                            if (sg->d[n1] == 1) {
                                if (sg->d[n2] == 1) {
                                    if (Cand->invlab[n1] < Cand->invlab[n2]) {
                                        start = n1;
                                    }
                                    else {
                                        start = n2;
                                    }
                                }
                                else {
                                    start = n1;
                                }
                            }
                            else {
                                start = n2;
                            }
                            PERMSTACK[StInd++] = start;
                            prev = start;
                            counts++;
                            do {
                                Markers[step] = tv->mark;
                                PERMSTACK[StInd++] = step;
                                if (sg->e[sg->v[step]] != prev) {
                                    prev = step;
                                    step = sg->e[sg->v[step]];
                                } else {
                                    prev = step;
                                    step = sg->e[sg->v[step]+1];
                                }
                            } while (sg->d[step] == 2);
                            PERMSTACK[StInd++] = step;
                        }
                    }
                    if (counts == SpineSDegLev.part->cls[i]) {
                        if (SpineSDegLev.part->inv[Cand->invlab[PERMSTACK[0]]] != SpineSDegLev.part->inv[Cand->invlab[PERMSTACK[StInd/counts-1]]]) {
                            factorial(grpsize1, grpsize2, SpineSDegLev.part->cls[i]);
                        }
                        else {
                            factorial2(grpsize1, grpsize2, 2*SpineSDegLev.part->cls[i]);
                            for (j=0; j<counts; j++) {
                                if (tv->options->generators || tv->options->writeautoms || tv->options->userautomproc) {
                                    j0 = j*(StInd/counts);
                                    k1 = (j+1)*(StInd/counts);
                                    memcpy(AUTPERM, IDENTITY_PERM, n*sizeof(int));
                                    for (k=j0, i1=k1-1; k<k1; k++, i1--) {
                                        AUTPERM[PERMSTACK[k]] = PERMSTACK[i1];
                                    }
                                }
                                SPECIALGENERATORS
                            }
                        }
                        if (tv->options->generators || tv->options->writeautoms || tv->options->userautomproc) {
                            memcpy(AUTPERM, IDENTITY_PERM, n*sizeof(int));
                            for (k=0; k<StInd/counts; k++) {
                                i1 = PERMSTACK[k];
                                for (j0=0; j0<counts-1; j0++) {
                                    AUTPERM[PERMSTACK[j0*(StInd/counts)+k]] = PERMSTACK[(j0+1)*(StInd/counts)+k];
                                }
                                AUTPERM[PERMSTACK[j0*(StInd/counts)+k]] =i1;
                            }
                        }
                        SPECIALGENERATORS
                        if (counts > 2) {
                            if (tv->options->generators || tv->options->writeautoms || tv->options->userautomproc) {
                                memcpy(AUTPERM, IDENTITY_PERM, n*sizeof(int));
                                for (k=0; k<StInd/counts; k++) {
                                    i1 = PERMSTACK[k];
                                    for (j0=0; j0<1; j0++) {
                                        AUTPERM[PERMSTACK[j0*(StInd/counts)+k]] = PERMSTACK[(j0+1)*(StInd/counts)+k];
                                    }
                                    AUTPERM[PERMSTACK[j0*(StInd/counts)+k]] =i1;
                                }
                            }
                            SPECIALGENERATORS
                        }
                    }
                    else {
                        factorial2(grpsize1, grpsize2, SpineSDegLev.part->cls[i]);
                        for (j=0; j<counts; j++) {
                            if (tv->options->generators || tv->options->writeautoms || tv->options->userautomproc) {
                                j0 = j*(StInd/counts);
                                k1 = (j+1)*(StInd/counts);
                                memcpy(AUTPERM, IDENTITY_PERM, n*sizeof(int));
                                for (k=j0, i1=k1-1; k<k1; k++, i1--) {
                                    AUTPERM[PERMSTACK[k]] = PERMSTACK[i1];
                                }
                            }
                            SPECIALGENERATORS
                        }
                        if (counts > 1) {
                            if (tv->options->generators || tv->options->writeautoms || tv->options->userautomproc) {
                                memcpy(AUTPERM, IDENTITY_PERM, n*sizeof(int));
                                for (k=0; k<StInd/counts; k++) {
                                    i1 = PERMSTACK[k];
                                    for (j0=0; j0<counts-1; j0++) {
                                        AUTPERM[PERMSTACK[j0*(StInd/counts)+k]] = PERMSTACK[(j0+1)*(StInd/counts)+k];
                                    }
                                    AUTPERM[PERMSTACK[j0*(StInd/counts)+k]] =i1;
                                }
                            }
                            SPECIALGENERATORS
                        }
                        if (counts > 2) {
                            if (tv->options->generators || tv->options->writeautoms || tv->options->userautomproc) {
                                memcpy(AUTPERM, IDENTITY_PERM, n*sizeof(int));
                                for (k=0; k<StInd/counts; k++) {
                                    i1 = PERMSTACK[k];
                                    for (j0=0; j0<1; j0++) {
                                        AUTPERM[PERMSTACK[j0*(StInd/counts)+k]] = PERMSTACK[(j0+1)*(StInd/counts)+k];
                                    }
                                    AUTPERM[PERMSTACK[j0*(StInd/counts)+k]] =i1;
                                }
                            }
                            SPECIALGENERATORS
                        }
                    }
                    i1 = StInd/counts;
                    if (SpineSDegLev.part->inv[Cand->invlab[PERMSTACK[0]]] != SpineSDegLev.part->inv[Cand->invlab[PERMSTACK[StInd/counts-1]]]) {
                        for (j=0; j<i1; j++) {
                            orbrep = n;
                            for (k=0; k<StInd; k+=i1) {
                                tmp = PERMSTACK[k+j];
                                if (tv->orbits[tmp] < orbrep) {
                                    orbrep = tv->orbits[tmp];
                                }
                            }
                            cellord = SpineSDegLev.part->inv[SpineSDegLev.liststart->invlab[PERMSTACK[j]]];
                            StackMarkers[orbrep] = tv->stackmark;
                            for (k=0; k<StInd; k+=i1) {
                                tmp = PERMSTACK[k+j];
                                SpineSDegLev.liststart->lab[cellord] = tmp;
                                SpineSDegLev.liststart->invlab[tmp] = cellord++;
                                if (StackMarkers[tv->orbits[tmp]] != tv->stackmark) {
                                    tv->stats->numorbits--;
                                    StackMarkers[tv->orbits[tmp]] = tv->stackmark;
                                }
                                tv->orbits[tmp] = orbrep;
                            }
                        }
                    }
                    else {
                        for (j=0; j<i1/2; j++) {
                            orbrep = n;
                            for (k=0, k1=i1-1; k<StInd; k+=i1, k1+=i1) {
                                tmp = PERMSTACK[k+j];
                                temp = PERMSTACK[k1-j];
                                if (tv->orbits[tmp] < orbrep) {
                                    orbrep = tv->orbits[tmp];
                                }
                                if (tv->orbits[temp] < orbrep) {
                                    orbrep = tv->orbits[temp];
                                }
                            }
                            cellord = SpineSDegLev.part->inv[SpineSDegLev.liststart->invlab[PERMSTACK[j]]];
                            StackMarkers[orbrep] = tv->stackmark;
                            for (k=0, k1=i1-1; k<StInd; k+=i1, k1+=i1) {
                                tmp = PERMSTACK[k+j];
                                SpineSDegLev.liststart->lab[cellord] = tmp;
                                SpineSDegLev.liststart->invlab[tmp] = cellord++;
                                if (StackMarkers[tv->orbits[tmp]] != tv->stackmark) {
                                    tv->stats->numorbits--;
                                    StackMarkers[tv->orbits[tmp]] = tv->stackmark;
                                }
                                tv->orbits[tmp] = orbrep;
                                
                                temp = PERMSTACK[k1-j];
                                SpineSDegLev.liststart->lab[cellord] = temp;
                                SpineSDegLev.liststart->invlab[temp] = cellord++;
                                if (StackMarkers[tv->orbits[temp]] != tv->stackmark) {
                                    tv->stats->numorbits--;
                                    StackMarkers[tv->orbits[temp]] = tv->stackmark;
                                }
                                tv->orbits[temp] = orbrep;
                            }
                        }
                        if (i1 % 2) {
                            orbrep = n;
                            for (k=j; k<StInd; k+=i1) {
                                tmp = PERMSTACK[k];
                                if (tv->orbits[tmp] < orbrep) {
                                    orbrep = tv->orbits[tmp];
                                }
                            }
                            cellord = SpineSDegLev.part->inv[SpineSDegLev.liststart->invlab[PERMSTACK[j]]];
                            StackMarkers[orbrep] = tv->stackmark;
                            j0 = 0;
                            for (k=j; k<StInd; k+=i1) {
                                tmp = PERMSTACK[k];
                                SpineSDegLev.liststart->lab[cellord+j0] = tmp;
                                SpineSDegLev.liststart->invlab[tmp] = cellord+j0++;
                                if (StackMarkers[tv->orbits[tmp]] != tv->stackmark) {
                                    tv->stats->numorbits--;
                                    StackMarkers[tv->orbits[tmp]] = tv->stackmark;
                                }
                                tv->orbits[tmp] = orbrep;
                            }
                        }
                    }
                    Propagate(sg, SpineSDegLev.liststart, SpineSDegLev.part, i, i+SpineSDegLev.part->cls[i], n);
                }
            }
        }
        
        /* Cycles */
        for (i=0; i<n; i += SpineSDegLev.part->cls[i]) {
            if (SpineSDegLev.part->cls[i] > 1) {
                tmp = SpineSDegLev.liststart->lab[i];
                if (sg->d[tmp] == 2) {
                    SETMARK(Markers, tv->mark)
                    CyInd = StInd = cycnum = 0;
                    for (j=i; j<i+SpineSDegLev.part->cls[i]; j++) {
                        start = SpineSDegLev.liststart->lab[j];
                        if (Markers[start] != tv->mark) {
                            counts = 1;
                            CYCLES[StInd] = start;
                            CYCOLR[StInd++] = SpineSDegLev.part->inv[SpineSDegLev.liststart->invlab[start]];
                            Markers[start] = tv->mark;
                            k = SpineSDegLev.liststart->invlab[sg->e[sg->v[start]]];
                            k1 = SpineSDegLev.liststart->invlab[sg->e[sg->v[start]+1]];
                            if (SpineSDegLev.part->inv[k] < SpineSDegLev.part->inv[k1]) {
                                step = sg->e[sg->v[start]];
                            } else {
                                step = sg->e[sg->v[start]+1];
                            }
                            prev = start;
                            do {
                                counts++;
                                Markers[step] = tv->mark;
                                CYCLES[StInd] = step;
                                CYCOLR[StInd++] = SpineSDegLev.part->inv[SpineSDegLev.liststart->invlab[step]];
                                if (sg->e[sg->v[step]] != prev) {
                                    prev = step;
                                    step = sg->e[sg->v[step]];
                                } else {
                                    prev = step;
                                    step = sg->e[sg->v[step]+1];
                                }
                            } while (step != start);
                            CYLGTH[CyInd++] = counts;
                            cycnum++;
                        }
                    }
                    
                    CYCPOS[0] = 0;
                    WorkArray[0] = CYLGTH[0];
                    for (j=1; j<CyInd; j++) {
                        CYCPOS[j] = CYCPOS[j-1]+CYLGTH[j-1];
                        WorkArray[CYCPOS[j]] = CYLGTH[j];
                    }
                    sortindirect(CYCPOS, WorkArray, CyInd);
                    
                    /* the next loop for orbits and canonical label */
                    SETMARK(Markers, tv->mark)
                    
                    c1 = WorkArray[CYCPOS[0]];
                    LghAttrInd = 0;
                    for (c=0; c<CyInd; c++) {
                        if (WorkArray[CYCPOS[c]] != c1) {
                            c1 = WorkArray[CYCPOS[c]];
                        }
                        cylab = CYCLES+CYCPOS[c];
                        cycol = CYCOLR+CYCPOS[c];
                        for (c2=0; c2<WorkArray[CYCPOS[c]]; c2++) {
                            LGHATTR[CYCPOS[c]+c2] = WorkArray[CYCPOS[c]];
                            if (Markers[cycol[c2]] != tv->mark) {
                                Markers[cycol[c2]] = tv->mark;
                                CYCHIT[cycol[c2]] = cycol[c2];
                                SpineSDegLev.liststart->lab[CYCHIT[cycol[c2]]] = cylab[c2];
                                SpineSDegLev.liststart->invlab[cylab[c2]] = CYCHIT[cycol[c2]];
                            }
                            else {
                                SpineSDegLev.liststart->lab[++CYCHIT[cycol[c2]]] = cylab[c2];
                                SpineSDegLev.liststart->invlab[cylab[c2]] = CYCHIT[cycol[c2]];
                            }
                        }
                    }
                    
                    SETMARK(Markers, tv->mark)
                    c2 = LGHATTR[CYCPOS[0]];
                    c = k = 0;
                    while (c<CyInd) {
                        if (LGHATTR[CYCPOS[c]] != c2) {
                            c2 = LGHATTR[CYCPOS[c]];
                            for (c3 = k; c3<c; c3++) {
                                for (c4=CYCPOS[c3]; c4<CYCPOS[c3]+LGHATTR[CYCPOS[c3]]; c4++) {
                                    tv->orbits[CYCLES[c4]] = CYCREP[CYCOLR[c4]];
                                }
                            }
                            k = c;
                            SETMARK(Markers, tv->mark)
                        }
                        for (c1=CYCPOS[c]; c1<CYCPOS[c]+LGHATTR[CYCPOS[c]]; c1++) {
                            if (Markers[CYCOLR[c1]] != tv->mark) {
                                Markers[CYCOLR[c1]] = tv->mark;
                                CYCREP[CYCOLR[c1]] = CYCLES[c1];
                            }
                            else {
                                if (CYCLES[c1] < CYCREP[CYCOLR[c1]]) {
                                    CYCREP[CYCOLR[c1]] = CYCLES[c1];
                                }
                            }
                        }
                        c++;
                    }
                    for (c3 = k; c3<c; c3++) {
                        for (c4=CYCPOS[c3]; c4<CYCPOS[c3]+LGHATTR[CYCPOS[c3]]; c4++) {
                            tv->orbits[CYCLES[c4]] = CYCREP[CYCOLR[c4]];
                        }
                    }
                    
                    k = 0;
                    for (i1=0; i1<CyInd; i1++) {
                        k1 = CYCOLR[k];
                        k2 = CYCOLR[k+1];
                        for (j=1; j<=CYLGTH[i1]/2; j++) {
                            w1 = CYCOLR[j+k];
                            w2 = CYCOLR[j+1+k];
                            if ((w1 == k1) && (w2 == k2)) {
                                if (tv->options->generators || tv->options->writeautoms || tv->options->userautomproc) {
                                    memcpy(AUTPERM, IDENTITY_PERM, n*sizeof(int));
                                }
                                for (w=0; w<CYLGTH[i1]; w++) {
                                    if (CYCOLR[w+k] == CYCOLR[((w+j) % CYLGTH[i1]) + k]) {
                                        if (tv->options->generators || tv->options->writeautoms || tv->options->userautomproc) {
                                            AUTPERM[CYCLES[w+k]] = CYCLES[((w+j) % CYLGTH[i1]) + k];
                                        }
                                    } else {
                                        break;
                                    }
                                }
                                SPECIALGENERATORS
                                if (w == CYLGTH[i1]) {
                                    MULTIPLY(*grpsize1, *grpsize2, CYLGTH[i1]/j);
                                    break;
                                }
                            }
                        }
                        if (SpineSDegLev.part->cls[k1] >= SpineSDegLev.part->cls[k2]) {
                            for (j=CYLGTH[i1]-1; j>0; j--) {
                                w1 = CYCOLR[j % CYLGTH[i1] + k];
                                w2 = CYCOLR[(j-1) % CYLGTH[i1] + k];
                                if ((w1 == k1) && (w2 == k2)) {
                                    if (tv->options->generators || tv->options->writeautoms || tv->options->userautomproc) {
                                        memcpy(AUTPERM, IDENTITY_PERM, n*sizeof(int));
                                        for (w=0; w<CYLGTH[i1]; w++) {
                                            AUTPERM[CYCLES[w+k]] = CYCLES[((j-w+(w>j)*CYLGTH[i1]) % CYLGTH[i1]) + k];
                                        }
                                    }
                                    SPECIALGENERATORS
                                    MULTIPLY(*grpsize1, *grpsize2, 2);
                                    break;
                                }
                            }
                        }
                        else {
                            j=CYLGTH[i1]-1;
                            w2 = CYCOLR[j % CYLGTH[i1] + k];
                            if (w2 == k2) {
                                if (tv->options->generators || tv->options->writeautoms || tv->options->userautomproc) {
                                    memcpy(AUTPERM, IDENTITY_PERM, n*sizeof(int));
                                    for (w=1; w<CYLGTH[i1]; w++) {
                                        AUTPERM[CYCLES[w+k]] = CYCLES[CYLGTH[i1]-w+k];
                                    }
                                }
                                SPECIALGENERATORS
                                MULTIPLY(*grpsize1, *grpsize2, 2);
                            }
                        }
                        k += CYLGTH[i1];
                    }
                    k = 0;
                    for (i1=0; i1<CyInd; i1++) {
                        if (CYLGTH[i1] > 0) {
                            CYMULT[0] = k;
                            k1 = k;
                            counts = 1;
                            for (j0=i1+1; j0<CyInd; j0++) {
                                k1 += abs(CYLGTH[j0]);
                                if (CYLGTH[j0] == CYLGTH[i1]) {
                                    CYMULT[counts++] = k1;
                                    CYLGTH[j0] = -CYLGTH[j0];
                                }
                            }
                            if (counts > 1) {
                                if (tv->options->generators || tv->options->writeautoms || tv->options->userautomproc) {
                                    memcpy(AUTPERM, IDENTITY_PERM, n*sizeof(int));
                                    for (j0=0; j0<CYLGTH[i1]; j0++) {
                                        for (j2 = 0; j2<counts-1; j2++) {
                                            AUTPERM[CYCLES[CYMULT[j2]+j0]] = CYCLES[CYMULT[j2+1]+j0];
                                        }
                                        AUTPERM[CYCLES[CYMULT[j2]+j0]] = CYCLES[CYMULT[0]+j0];
                                    }
                                }
                                SPECIALGENERATORS
                                if (counts > 2) {
                                    if (tv->options->generators || tv->options->writeautoms || tv->options->userautomproc) {
                                        memcpy(AUTPERM, IDENTITY_PERM, n*sizeof(int));
                                        for (j0=0; j0<CYLGTH[i1]; j0++) {
                                            AUTPERM[CYCLES[CYMULT[0]+j0]] = CYCLES[CYMULT[1]+j0];
                                            AUTPERM[CYCLES[CYMULT[1]+j0]] = CYCLES[CYMULT[0]+j0];
                                        }
                                    }
                                    SPECIALGENERATORS
                                }
                                factorial(grpsize1, grpsize2, counts);
                            }
                        }
                        k += abs(CYLGTH[i1]);
                        CYLGTH[i1] = -CYLGTH[i1];
                    }
                    Propagate(sg, SpineSDegLev.liststart, SpineSDegLev.part, i, i+SpineSDegLev.part->cls[i], n);
                }
            }
        }
        
        /* Deg 1 */
        SETMARK(Markers, tv->mark)
        for (i=0; i<n; i += SpineSDegLev.part->cls[i]) {
			if (SpineSDegLev.part->cls[i] > 1) {
				tmp = SpineSDegLev.liststart->lab[i];
                /* tmp had degree 1, and its nghb too */
                if ((sg->d[tmp] == 1) && (sg->d[sg->e[sg->v[tmp]]] == 1)) {
                    cell1 = SpineSDegLev.part->inv[Cand->invlab[sg->e[sg->v[tmp]]]];
                    if (i < cell1) {
                        if (tv->options->generators || tv->options->writeautoms || tv->options->userautomproc) {
                            memcpy(AUTPERM, IDENTITY_PERM, n*sizeof(int));
                            for (j=i; j<i+SpineSDegLev.part->cls[i]-1; j++) {
                                AUTPERM[SpineSDegLev.liststart->lab[j]] = SpineSDegLev.liststart->lab[j+1];
                                AUTPERM[sg->e[sg->v[SpineSDegLev.liststart->lab[j]]]] = sg->e[sg->v[SpineSDegLev.liststart->lab[j+1]]];
                            }
                            AUTPERM[SpineSDegLev.liststart->lab[j]] = tmp;
                            AUTPERM[sg->e[sg->v[SpineSDegLev.liststart->lab[j]]]] = sg->e[sg->v[tmp]];
                        }
                        SPECIALGENERATORS
                        if (SpineSDegLev.part->cls[i] > 2) {
                            if (tv->options->generators || tv->options->writeautoms || tv->options->userautomproc) {
                                memcpy(AUTPERM, IDENTITY_PERM, n*sizeof(int));
                                AUTPERM[tmp] = SpineSDegLev.liststart->lab[i+1];
                                AUTPERM[SpineSDegLev.liststart->lab[i+1]] = tmp;
                                AUTPERM[sg->e[sg->v[tmp]]] = sg->e[sg->v[SpineSDegLev.liststart->lab[i+1]]];
                                AUTPERM[sg->e[sg->v[SpineSDegLev.liststart->lab[i+1]]]] = sg->e[sg->v[tmp]];
                            }
                            SPECIALGENERATORS
                        }
                        factorial(grpsize1, grpsize2, SpineSDegLev.part->cls[i]);
                        c1 = cell1;
                        for (j=i; j<i+SpineSDegLev.part->cls[i]; j++) {
                            vtx = SpineSDegLev.liststart->lab[j];
                            SpineSDegLev.liststart->lab[c1] = sg->e[sg->v[vtx]];
                            SpineSDegLev.liststart->invlab[sg->e[sg->v[vtx]]] = c1;
                            c1++;
                        }
                        orbrep = n;
                        for (j=i; j<i+SpineSDegLev.part->cls[i]; j++) {
                            if (tv->orbits[SpineSDegLev.liststart->lab[j]] < orbrep) {
                                orbrep = tv->orbits[SpineSDegLev.liststart->lab[j]];
                            }
                        }
                        for (j=i; j<i+SpineSDegLev.part->cls[i]; j++) {
                            tv->orbits[SpineSDegLev.liststart->lab[j]] = orbrep;
                        }
                        tv->stats->numorbits -= SpineSDegLev.part->cls[i] - 1;
                    }
                    else {
                        if (i > cell1) {
                            orbrep = n;
                            for (j=i; j<i+SpineSDegLev.part->cls[i]; j++) {
                                if (tv->orbits[SpineSDegLev.liststart->lab[j]] < orbrep) {
                                    orbrep = tv->orbits[SpineSDegLev.liststart->lab[j]];
                                }
                            }
                            for (j=i; j<i+SpineSDegLev.part->cls[i]; j++) {
                                tv->orbits[SpineSDegLev.liststart->lab[j]] = orbrep;
                            }
                            tv->stats->numorbits -= SpineSDegLev.part->cls[i] - 1;
                        }
                        if (i == cell1) {
                            factorial2(grpsize1, grpsize2, SpineSDegLev.part->cls[i]);
                            /* the cell has size two */
                            if (SpineSDegLev.part->cls[i] == 2) {
                                if (tv->options->generators || tv->options->writeautoms || tv->options->userautomproc) {
                                    memcpy(AUTPERM, IDENTITY_PERM, n*sizeof(int));
                                    AUTPERM[SpineSDegLev.liststart->lab[i]] = SpineSDegLev.liststart->lab[i+1];
                                    AUTPERM[SpineSDegLev.liststart->lab[i+1]] = tmp;
                                }
                                SPECIALGENERATORS
                            }
                            else {
                                /* the cell has size greater than two */
                                if (tv->options->generators || tv->options->writeautoms || tv->options->userautomproc) {
                                    memcpy(AUTPERM, IDENTITY_PERM, n*sizeof(int));
                                    SETMARK(Markers, tv->mark)
                                    halfsize = SpineSDegLev.part->cls[i]/2;
                                    i1 = 0;
                                    for (j=i; j<i+SpineSDegLev.part->cls[i]; j++) {
                                        if (Markers[SpineSDegLev.liststart->lab[j]] != tv->mark) {
                                            Markers[sg->e[sg->v[SpineSDegLev.liststart->lab[j]]]] = tv->mark;
                                            PERMSTACK[i1] = SpineSDegLev.liststart->lab[j];
                                            PERMSTACK[i1+halfsize] = sg->e[sg->v[SpineSDegLev.liststart->lab[j]]];
                                            i1++;
                                        }
                                    }
                                    temp = PERMSTACK[0];
                                    for (j=0; j<SpineSDegLev.part->cls[i]-1; j++) {
                                        AUTPERM[PERMSTACK[j]] = PERMSTACK[j+1];
                                    }
                                    AUTPERM[PERMSTACK[j]] = temp;
                                }
                                
                                SPECIALGENERATORS
                                if (tv->options->generators || tv->options->writeautoms || tv->options->userautomproc) {
                                    memcpy(AUTPERM, IDENTITY_PERM, n*sizeof(int));
                                    memcpy(PERMSTACK+halfsize, PERMSTACK+halfsize+1, (halfsize-1)*sizeof(int));
                                    temp = PERMSTACK[1];
                                    for (j=1; j<SpineSDegLev.part->cls[i]-2; j++) {
                                        AUTPERM[PERMSTACK[j]] = PERMSTACK[j+1];
                                    }
                                    AUTPERM[PERMSTACK[j]] = temp;
                                }
                                SPECIALGENERATORS
                            }
                            SETMARK(Markers, tv->mark)
                            for (j=i; j<i+SpineSDegLev.part->cls[i]; j++) {
                                temp = SpineSDegLev.liststart->lab[j];
                                if (Markers[temp] != tv->mark) {
                                    tmp = SpineSDegLev.liststart->lab[j+1];
                                    Markers[sg->e[sg->v[temp]]] = tv->mark;
                                    i1 = SpineSDegLev.liststart->invlab[sg->e[sg->v[temp]]];
                                    SpineSDegLev.liststart->lab[j+1] = sg->e[sg->v[temp]];
                                    SpineSDegLev.liststart->invlab[sg->e[sg->v[temp]]] = j+1;
                                    SpineSDegLev.liststart->lab[i1] = tmp;
                                    SpineSDegLev.liststart->invlab[tmp] = i1;
                                }
                            }
                            orbrep = n;
                            for (j=i; j<i+SpineSDegLev.part->cls[i]; j++) {
                                if (tv->orbits[SpineSDegLev.liststart->lab[j]] < orbrep) {
                                    orbrep = tv->orbits[SpineSDegLev.liststart->lab[j]];
                                }
                            }
                            for (j=i; j<i+SpineSDegLev.part->cls[i]; j++) {
                                tv->orbits[SpineSDegLev.liststart->lab[j]] = orbrep;
                            }
                            tv->stats->numorbits -= SpineSDegLev.part->cls[i] - 1;
                        }
                    }
				}
                /* tmp had degree 1, its nghb does not have degreee 1 */
				else if (sg->d[tmp] == 1) {
					factorial(grpsize1, grpsize2, SpineSDegLev.part->cls[i]);
                    if (tv->options->generators || tv->options->writeautoms || tv->options->userautomproc) {
                        memcpy(AUTPERM, IDENTITY_PERM, n*sizeof(int));
                        for (j=i; j<i+SpineSDegLev.part->cls[i]-1; j++) {
                            AUTPERM[SpineSDegLev.liststart->lab[j]] = SpineSDegLev.liststart->lab[j+1];
                        }
                        AUTPERM[SpineSDegLev.liststart->lab[j]] = tmp;
                    }
                    SPECIALGENERATORS
					if (SpineSDegLev.part->cls[i] > 2) {
                        if (tv->options->generators || tv->options->writeautoms || tv->options->userautomproc) {
                            memcpy(AUTPERM, IDENTITY_PERM, n*sizeof(int));
                            AUTPERM[SpineSDegLev.liststart->lab[i]] = SpineSDegLev.liststart->lab[i+1];
                            AUTPERM[SpineSDegLev.liststart->lab[i+1]] = tmp;
                        }
                        SPECIALGENERATORS
					}
                    orbrep = n;
                    for (j=i; j<i+SpineSDegLev.part->cls[i]; j++) {
                        if (tv->orbits[SpineSDegLev.liststart->lab[j]] < orbrep) {
                            orbrep = tv->orbits[SpineSDegLev.liststart->lab[j]];
                            if (tv->orbits[SpineSDegLev.liststart->lab[j]] != SpineSDegLev.liststart->lab[j]) {
                            }
                        }
                    }
                    for (j=i; j<i+SpineSDegLev.part->cls[i]; j++) {
                        tv->orbits[SpineSDegLev.liststart->lab[j]] = orbrep;
                    }
                    
                    tv->stats->numorbits -= SpineSDegLev.part->cls[i] - 1;
                }
			}
		}
        
        /* Deg 0 or n-1 */
        for (i=0; i<n; i += SpineSDegLev.part->cls[i]) {
            if (SpineSDegLev.part->cls[i] > 1) {
                tmp = SpineSDegLev.liststart->lab[i];
                /* tmp has deg 0 or n-1 */
                if ((sg->d[tmp] == 0) || (sg->d[tmp] == n-1)) {
                    if (tv->options->generators || tv->options->writeautoms || tv->options->userautomproc) {
                        memcpy(AUTPERM, IDENTITY_PERM, n*sizeof(int));
                        for (j=i; j<i+SpineSDegLev.part->cls[i]-1; j++) {
                            AUTPERM[SpineSDegLev.liststart->lab[j]] = SpineSDegLev.liststart->lab[j+1];
                        }
                        AUTPERM[SpineSDegLev.liststart->lab[j]] = tmp;
                    }
                    SPECIALGENERATORS
                    if (SpineSDegLev.part->cls[i] > 2) {
                        if (tv->options->generators || tv->options->writeautoms || tv->options->userautomproc) {
                            memcpy(AUTPERM, IDENTITY_PERM, n*sizeof(int));
                            AUTPERM[tmp] = SpineSDegLev.liststart->lab[i+1];
                            AUTPERM[SpineSDegLev.liststart->lab[i+1]] = tmp;
                        }
                        SPECIALGENERATORS
                    }
                    factorial(grpsize1, grpsize2, SpineSDegLev.part->cls[i]);
                    orbrep = n;
                    for (j=i; j<i+SpineSDegLev.part->cls[i]; j++) {
                        if (tv->orbits[SpineSDegLev.liststart->lab[j]] < orbrep) {
                            orbrep = tv->orbits[SpineSDegLev.liststart->lab[j]];
                        }
                    }
                    for (j=i; j<i+SpineSDegLev.part->cls[i]; j++) {
                        tv->orbits[SpineSDegLev.liststart->lab[j]] = orbrep;
                    }
                    
                    tv->stats->numorbits -= SpineSDegLev.part->cls[i] - 1;
                }
            }
        }
	}

    SETMARK(Markers, tv->mark)
    i1=0;
    for (c1=0; c1<n; c1++) {
        if (Markers[tv->orbits[c1]] != tv->mark) {
            i1++;
            Markers[tv->orbits[c1]] = tv->mark;
        }
    }
    tv->stats->numorbits = i1;
    return;
}

boolean Prefix(Candidate *Cand1, Candidate *Cand2, int k)
{
	int i;
	for (i=1; i<=k; i++) {
		if (Cand1->lab[Spine[k].tgtpos] != Cand2->lab[Spine[k].tgtpos]) {
			break;
		}
	}
	return (i>k);
}


boolean findperm(permnode *pn, int *p, int n, TracesVars *tv)
{
    permnode *rn;
    if (!pn) {
        return FALSE;
    }
    rn = pn;
    do {
        if (!memcmp(rn->p, p, n*sizeof(int))) {
            return TRUE;
        }
        rn = rn->next;
    } while (rn != pn);
    return FALSE;
}


int spinelementorbsize(int *orbits, int *lab, int size, int elem)
{
    int i, j, val;
    j = 0;
    val = orbits[elem];
    for (i = 0; i < size; ++i) {
        if (orbits[lab[i]] == val) ++j;
    }
    return j;
}

void Propagate(sparsegraph *sg, Candidate *Cand, Partition *Part, int start, int end, int n) {
    int i, vtx, n1, n2;
    
    for (i=start; i<end; i++) {
        Part->cls[i] = 1;
        Part->inv[i] = i;
    }
    Part->cells = Part->cells + end - start - 1;
    vtx = Cand->lab[start];
    if (sg->d[vtx] == 2) {
        n1 = sg->e[sg->v[vtx]];
        n2 = sg->e[sg->v[vtx]+1];
        if ((Part->cls[Part->inv[Cand->invlab[n1]]] > 1) && (sg->d[n1] < n-1))
            Propagate(sg, Cand, Part, Part->inv[Cand->invlab[n1]], Part->inv[Cand->invlab[n1]]+Part->cls[Part->inv[Cand->invlab[n1]]], n);
        if ((Part->cls[Part->inv[Cand->invlab[n2]]] > 1) && (sg->d[n2] < n-1))
            Propagate(sg, Cand, Part, Part->inv[Cand->invlab[n2]], Part->inv[Cand->invlab[n2]]+Part->cls[Part->inv[Cand->invlab[n2]]], n);
        return;
    }
}

trielist *searchtrie_new(int n, struct TracesVars *tv)
{
	tv->strielist = malloc(sizeof(struct trielist));
	if (tv->strielist == NULL) {
		fprintf(ERRFILE, "\nError, memory not allocated.\n");
		exit(1);
	}
    tv->strielist->prev = tv->strielist->next = NULL;
    tv->strielist->triearray = malloc(n*sizeof(searchtrie));
	if (tv->strielist->triearray == NULL) {
		fprintf(ERRFILE, "\nError, memory not allocated.\n");
		exit(1);
	}
	tv->strielist->triearray[0].father = tv->strielist->triearray[0].first_child = NULL;
	tv->strielist->triearray[0].next_sibling = tv->strielist->triearray[0].last_child = NULL;
    tv->strielist->triearray[0].goes_to = NULL;
    tv->strielist->triearray[0].index = 1;
    tv->strielist->triearray[0].name = tv->strielist->triearray[0].level = 0;
    tv->strielist->triearray[0].vtx = n;
    
	tv->strienext = 1;
    return tv->strielist;
}

searchtrie *searchtrie_make(Candidate *CurrCand, Candidate *NextCand, int n, struct TracesVars *tv)
{
    searchtrie *st;
	if (tv->strienext == n) {
		tv->strienext = 0;
		tv->strielist->next = malloc(sizeof(struct trielist));
		if (tv->strielist->next == NULL) {
			fprintf(ERRFILE, "\nError, memory not allocated.\n");
			exit(1);
		}
        tv->strielist->next->prev = tv->strielist;
        tv->strielist = tv->strielist->next;
        tv->strielist->next = NULL;
        tv->strielist->triearray = malloc(n*sizeof(searchtrie));
        if (tv->strielist->triearray == NULL) {
            fprintf(ERRFILE, "\nError, memory not allocated.\n");
            exit(1);
        }
    }
    st = &(tv->strielist->triearray[tv->strienext]);
	st->father = CurrCand->stnode;
	st->name = NextCand->name;
	st->index = tv->newindex+1;
	st->vtx = NextCand->vertex;
	st->level = tv->tolevel;
    st->first_child = st->next_sibling = st->last_child = st->goes_to = NULL;
    if (st->father) {
        if (st->father->first_child) {
            st->father->last_child->next_sibling = st;
            st->father->last_child = st;
        } else {
            st->father->first_child = st->father->last_child = st;
        }
    }
    NextCand->stnode = st;
    if (tv->newgotonode) {
        tv->newgotonode->goes_to = st;
    }
    if (tv->gotonode) {
        st->goes_to = tv->gotonode;
        tv->gotonode = NULL;
    }
    tv->strienext++;
    return st;
}

boolean lookup(searchtrie *t) {
    searchtrie *TreeNode;
    TreeNode = t;
    while (TreeNode->level >= 1) {
        if (TreeNode->goes_to) {
            return FALSE;
        }
        TreeNode = TreeNode->father;
    }
    return TRUE;
}

int *findcurrorbits(schreier *gp, int k)
{
    int i;
    schreier *sh;
    
    sh = gp;
    for (i = 0; i < k; i++) {
        sh = sh->next;
    }
    return sh->orbits;
}
