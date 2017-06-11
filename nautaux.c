/*****************************************************************************
*                                                                            *
* Auxiliary procedures for use with nauty 1.9.                               *
* None of these procedures are needed by nauty or by dreadnaut.              *
*                                                                            *
*   Copyright (1984-2000) Brendan McKay.  All rights reserved.               *
*   Subject to waivers and disclaimers in nauty.h.                           *
*                                                                            *
*   CHANGE HISTORY                                                           *
*       26-Apr-89 : initial creation for Version 1.5.                        *
*       14-Oct-90 : renamed to version 1.6 (no changes to this file)         *
*        5-Jun-93 : renamed to version 1.7+ (no changes to this file)        *
*       18-Aug-93 : renamed to version 1.8 (no changes to this file)         *
*       17-Sep-93 : renamed to version 1.9 (no changes to this file)         *
*       24-Jan-00 : renamed to version 2.0 (no changes to this file)
*                                                                            *
*****************************************************************************/

#define  EXTDEFS 1
#include "naututil.h"          /* which includes "nauty.h" and <stdio.h> */
#include "nautaux.h"

#if  MAXM==1
#define M 1
#else
#define M m
#endif

static permutation workperm[MAXN+2];   /* used for scratch purposes */
static set workset[MAXM];              /* used for scratch purposes */

/*****************************************************************************
*                                                                            *
*  ptncode(g,lab,ptn,level,m,n) returns a long integer invariant which       *
*  depends on the (assumed equitable) partition at the stated level, and     *
*  the number of edges betwen the various cells.                             *
*  Neither nauty nor dreadnaut use this.                                     *
*                                                                            *
*  GLOBALS ACCESSED: bit<r>,setinter()                                       *
*                                                                            *
*****************************************************************************/

long
ptncode(g,lab,ptn,level,m,n)
graph *g;
register nvector *lab,*ptn;
int level,m,n;
{
        register int i;
        register long code;
        int v1,v2,nc,inter,cellend;

        /* find all cells: put starts in workperm[0..n] */

        i = nc = 0;
        code = 0;

        while (i < n)
        {
            workperm[nc++] = i;
            code = ((code << 13) ^ (code >> 19)) + i;
            while (ptn[i] > level)
                ++i;
            ++i;
        }
        workperm[nc] = n;

        for (v2 = 0; v2 < nc; ++v2)
        {
            EMPTYSET(workset,m);
            for (i = workperm[v2], cellend = workperm[v2+1] - 1;
                                                          i <= cellend; ++i)
                ADDELEMENT(workset,lab[i]);
            for (v1 = 0; v1 < nc; ++v1)
            {
                i = workperm[v1];
                cellend = workperm[v1+1] - 1;
                inter = setinter(workset,GRAPHROW(g,lab[i],M),M);
                code = ((code << 13) ^ (code >> 19)) + inter;
            }
        }

        return(code);
}

/*****************************************************************************
*                                                                            *
*  equitable(g,lab,ptn,level,m,n) checks that the partition at the given     *
*  level is equitable.  Neither nauty nor dreadnaut use this.                *
*                                                                            *
*  GLOBALS ACCESSED: bit<r>,setinter()                                       *
*                                                                            *
*****************************************************************************/

boolean
equitable(g,lab,ptn,level,m,n)
graph *g;
register nvector *lab,*ptn;
int level,m,n;
{
        register int i;
        int v1,v2,nc,inter,cellend;
        boolean ok;

        /* find all cells: put starts in workperm[0..n] */

        i = nc = 0;

        while (i < n)
        {
            workperm[nc++] = i;
            while (ptn[i] > level)
                ++i;
            ++i;
        }
        workperm[nc] = n;

        ok = TRUE;
        for (v2 = 0; v2 < nc; ++v2)
        {
            EMPTYSET(workset,m);
            for (i = workperm[v2], cellend = workperm[v2+1] - 1;
                                                          i <= cellend; ++i)
                ADDELEMENT(workset,lab[i]);
            for (v1 = 0; v1 < nc; ++v1)
            {
                i = workperm[v1];
                cellend = workperm[v1+1] - 1;
                if (i == cellend)
                    continue;
                inter = setinter(workset,GRAPHROW(g,lab[i],M),M);
                while (++i <= cellend)
                     if (setinter(workset,GRAPHROW(g,lab[i],M),M) != inter)
                         ok = FALSE;
            }
        }

        return(ok);
}

/*****************************************************************************
*                                                                            *
*  component(g,v,c,m,n) determines the set of all vertices that can be       *
*  reached along directed paths starting at vertex v, including v itself.    *
*  This set is returned as c, unless c is null.  The size of the set is      *
*  returned as the function value.                                           *
*                                                                            *
*  GLOBALS ACCESSED: bit<r>,setinter(),nextelement()                         *
*                                                                            *
*****************************************************************************/

int
component(g,v,cmpt,m,n)
graph *g;
int v,m,n;
set *cmpt;
{
        register int i,z;
        set newverts[MAXM],*gx;
        int head,tail,x;

        EMPTYSET(workset,m);
        ADDELEMENT(workset,v);
        head = 0;
        tail = 1;
        workperm[head] = v;

        while (head < tail && tail < n)
        {
            x = workperm[head++];
            gx = GRAPHROW(g,x,m);
            for (i = m; --i >= 0;)
            {
                newverts[i] = gx[i] & ~workset[i];
                workset[i] |= gx[i];
            }
            for (z = nextelement(newverts,m,-1); z >= 0;
                                                 z = nextelement(newverts,m,z))
                workperm[tail++] = z;
        }

        if (cmpt != NILSET)
            for (i = m; --i >= 0;)
                cmpt[i] = workset[i];

        return(tail);
}
