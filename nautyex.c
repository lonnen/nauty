/* This program prints generators for the automorphism group of an
   n-vertex polygon, where n is a number supplied by the user.
   It needs to be linked with nauty.c and nautil.c.
*/

#include <stdio.h>

#define MAXN 100
#include "nauty.h"

main()
{
        graph g[MAXN*MAXM];
        nvector lab[MAXN],ptn[MAXN],orbits[MAXN];
        static DEFAULTOPTIONS(options);
        statsblk(stats);
        setword workspace[50*MAXM];

        int n,m,v;
        set *gv;

        options.writemarkers = FALSE;

        printf("\nenter n : ");
        if (scanf("%d",&n) == 1)
        {
            if (n < 1 || n > MAXN)
            {
                printf("n must be in the range 1..%d\n",MAXN);
                exit(1);
            }

            m = (n + WORDSIZE - 1) / WORDSIZE;

	    nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);

            for (v = 0; v < n; ++v)
            {
                gv = GRAPHROW(g,v,m);

                EMPTYSET(gv,m);
                ADDELEMENT(gv,(v+n-1)%n);
                ADDELEMENT(gv,(v+1)%n);
            }

            printf("Generators for Aut(C[%d]):\n",n);
            nauty(g,lab,ptn,NILSET,orbits,&options,&stats,
                            workspace,50*MAXM,m,n,NILGRAPH);
        }
}
