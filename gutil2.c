/* gutil2.c: Some more graph utilities. */

#include "gtools.h"
#include "gutils.h"

/**************************************************************************/

long
digoncount(graph *g, int m, int n)
/* Number of digons (cycles of length 2). Useful for digraphs. */
{
    int i,j;
    set *gi;
    setword w;
    long ans;

    ans = 0;

    if (m == 1)
    {
        for (i = 0; i < n; ++i)
        {
            w = g[i] & BITMASK(i);
            while (w)
            {
                TAKEBIT(j,w);
                if ((g[j] & bit[i])) ++ans;
            }
        }
    }
    else
    {
        for (i = 0, gi = g; i < n; ++i, gi += m)
        {
            for (j = i; (j = nextelement(gi,m,j)) > 0; )
                if (ISELEMENT(g+j*m,i)) ++ans;
        }
    }

    return ans;
}

/**************************************************************************/

int
loopcount(graph *g, int m, int n)
/* Number of loops */
{
    set *gi;
    int i,nl;

    nl = 0;
    for (i = 0, gi = g; i < n; ++i, gi += m)
        if (ISELEMENT(gi,i)) ++nl;

    return nl;
}

/**************************************************************************/

long
pathcount1(graph *g, int start, setword body, setword last)
/* Number of paths in g starting at start, lying within body and
   ending in last.  {start} and last should be disjoint subsets of body. */
{
    long count;
    setword gs,w;
    int i;

    gs = g[start];
    w = gs & last;
    count = POPCOUNT(w);

    body &= ~bit[start];
    w = gs & body;
    while (w)
    {
        TAKEBIT(i,w);
        count += pathcount1(g,i,body,last&~bit[i]);
    }

    return count;
}

/**************************************************************************/

long
cyclecount1(graph *g, int n)
/* The total number of cycles in g (assumed no loops), m=1 only */
{
    setword body,nbhd;
    long total;
    int i,j;

    body = ALLMASK(n);
    total = 0;

    for (i = 0; i < n-2; ++i)
    {
        body ^= bit[i];
        nbhd = g[i] & body;
        while (nbhd)
        {
            TAKEBIT(j,nbhd);
            total += pathcount1(g,j,body,nbhd);
        }
    }

    return total;
}

/**************************************************************************/

long
cyclecount(graph *g, int m, int n)
/* The total number of cycles in g (assumed no loops) */
{
    if (n == 0) return 0;
    if (m == 1) return cyclecount1(g,n);

    gt_abort(">E cycle counting is only implemented for n <= WORDSIZE\n");
    return 0;
}

/**************************************************************************/

long
indpathcount1(graph *g, int start, setword body, setword last)
/* Number of induced paths in g starting at start, extravertices within
 * body and ending in last.
 * {start}, body and last should be disjoint. */
{
    long count;
    setword gs,w;
    int i;

    gs = g[start];
    w = gs & last;
    count = POPCOUNT(w);

    w = gs & body;
    while (w)
    {
        TAKEBIT(i,w);
        count += indpathcount1(g,i,body&~gs,last&~bit[i]&~gs);
    }

    return count;
}

/**************************************************************************/

long
indcyclecount1(graph *g, int n)
/* The total number of induced cycles in g (assumed no loops), m=1 only */
{
    setword body,last,cni;
    long total;
    int i,j;

    body = ALLMASK(n);
    total = 0;

    for (i = 0; i < n-2; ++i)
    {
        body ^= bit[i];
        last = g[i] & body;
        cni = g[i] | bit[i];
        while (last)
        {
            TAKEBIT(j,last);
            total += indpathcount1(g,j,body&~cni,last);
        }
    }

    return total;
}

/**************************************************************************/

long
indcyclecount(graph *g, int m, int n)
/* The total number of induced cycles in g (assumed no loops) */
{
    if (n == 0) return 0;
    if (m == 1) return indcyclecount1(g,n);

    gt_abort(
      ">E induced cycle counting is only implemented for n <= WORDSIZE\n");
    return 0;
}
   
/**************************************************************************/

long
numind3sets1(graph *g, int n)
/* The number of triangles in g; undirected only */
{
    int i,j;
    setword gi,w;
    long total;

    total = 0;
    for (i = 2; i < n; ++i)
    {
        gi = ALLMASK(i) & ~g[i];
        while (gi)
        {
            TAKEBIT(j,gi);
            w = gi & ~g[j];
            total += POPCOUNT(w);
        }
    }

    return total;
}

/**************************************************************************/

long
numind3sets(graph *g, int m, int n)
/* The number of independent 3-sets in g; undirected only */
{
    if (m == 1) return numind3sets1(g,n);

    gt_abort(">E numind3sets is only implemented for n <= WORDSIZE\n");
    return 0;
}

/**************************************************************************/

long
numtriangles1(graph *g, int n)
/* The number of triangles in g; undirected only */
{
    int i,j;
    setword gi,w;
    long total;

    total = 0;
    for (i = 0; i < n-2; ++i)
    {
        gi = g[i] & BITMASK(i);
        while (gi)
        {
            TAKEBIT(j,gi);
            w = g[j] & gi;
            total += POPCOUNT(w);
        }
    }

    return total;
}

/**************************************************************************/

long
numtriangles(graph *g, int m, int n)
/* The number of triangles in g; undirected only */
{
    int i,j,k,kw;
    setword *gi,*gj,w;
    long total;

    if (m == 1) return numtriangles1(g,n);

    total = 0;
    for (i = 0, gi = g; i < n-2; ++i, gi += m)
        for (j = i; (j = nextelement(gi,m,j)) > 0; )
        {
            gj = GRAPHROW(g,j,m);
            kw = SETWD(j);
            w = gi[kw] & gj[kw] & BITMASK(SETBT(j));
            if (w) total += POPCOUNT(w);
            for (k = kw+1; k < m; ++k)
            {
                w = gi[k] & gj[k];
                total += POPCOUNT(w);
            }
        }

    return total;
}    

/**************************************************************************/

long
numdirtriangles1(graph *g, int n)
/* Number of directed triangles; m=1 only. */
{
    long total;
    int i,j,k;
    setword biti,bm,wi,wj;

    total = 0;
    for (i = 0; i < n; ++i)
    {
        biti = bit[i];
        bm = BITMASK(i);
        wi = g[i] & bm;
        while (wi)
        {
            TAKEBIT(j,wi);
            wj = g[j] & bm;
            while (wj)
            {
                TAKEBIT(k,wj);
                if ((g[k] & biti)) ++total;
            }
        }
    }

    return total;
}

/**************************************************************************/

long
numdirtriangles(graph *g, int m, int n)
/* The number of directed triangles in g */
{
    long total;
    int i,j,k;
    set *gi,*gj;

    if (m == 1) return numdirtriangles1(g,n);

    total = 0;
    for (i = 0, gi = g; i < n-2; ++i, gi += m)
        for (j = i; (j = nextelement(gi,m,j)) >= 0;)
        {
            gj = GRAPHROW(g,j,m);
            for (k = i; (k = nextelement(gj,m,k)) >= 0; )
                if (k != j && ISELEMENT(GRAPHROW(g,k,m),i)) ++total;
        }

    return total;
} 

/**************************************************************************/

void
commonnbrs(graph *g, int *minadj, int *maxadj, int *minnon, int *maxnon,
       int m, int n)
/* Count the common neighbours of pairs of vertices, and give the minimum
   and maximum for adjacent and non-adjacent vertices.  Undirected only.
   Null minimums are n+1 and null maximums are -1.
*/
{
    int j,k;
    int mina,maxa,minn,maxn;
    int cn;
    set *gi,*gj;
    setword w;

    if (n == 0)
    {
        *minadj = *maxadj = *minnon = *maxnon = 0;
        return;
    }

    mina = minn = n+1;
    maxa = maxn = -1;

    for (j = 0, gj = g; j < n; ++j, gj += m)
    for (gi = g; gi != gj; gi += m)
    {
        cn = 0;
        for (k = 0; k < m; ++k)
        {
            w = gi[k] & gj[k];
            if (w) cn += POPCOUNT(w);
        }

        if (ISELEMENT(gi,j))
        {
            if (cn < mina) mina = cn;
            if (cn > maxa) maxa = cn;
        }
        else
        {
            if (cn < minn) minn = cn;
            if (cn > maxn) maxn = cn;
        }
    }

    *minadj = mina;
    *maxadj = maxa;
    *minnon = minn;
    *maxnon = maxn;
}

/**************************************************************************/

void
delete1(graph *g, graph *h, int v, int n)
/* Delete vertex v from g, result in h */
{
    setword mask1,mask2,gi;
    int i;

    mask1 = ALLMASK(v);
    mask2 = BITMASK(v);

    for (i = 0; i < v; ++i)
    {
        gi = g[i];
        h[i] = (gi & mask1) | ((gi & mask2) << 1);
    }
    for (i = v; i < n-1; ++i)
    {
        gi = g[i+1];
        h[i] = (gi & mask1) | ((gi & mask2) << 1);
    }
}

/**************************************************************************/

void
contract1(graph *g, graph *h, int v, int w, int n)  
/* Contract distinct vertices v and w (not necessarily adjacent)
   with result in h.  No loops are created. */
{
    int x,y;
    setword bitx,bity,mask1,mask2;
    int i;

    if (w < v)
    {
        x = w; 
        y = v;
    }
    else
    {
        x = v; 
        y = w;
    }

    bitx = bit[x];
    bity = bit[y];
    mask1 = ALLMASK(y);
    mask2 = BITMASK(y);

    for (i = 0; i < n; ++i)
        if (g[i] & bity)
            h[i] = (g[i] & mask1) | bitx | ((g[i] & mask2) << 1);
        else
            h[i] = (g[i] & mask1) | ((g[i] & mask2) << 1);

    h[x] |= h[y];
    for (i = y+1; i < n; ++i) h[i-1] = h[i];
    h[x] &= ~bitx;
}

/**************************************************************************/

static TLS_ATTR int knm[18][16];  /* knm[n,m] = conncontent(K_n - m*K_2) */
static TLS_ATTR boolean knm_computed = FALSE;

int
conncontent(graph *g, int m, int n)
/* number of connected spanning subgraphs with an even number
   of edges minus the number with an odd number of edges */
{
    graph h[WORDSIZE];
    setword gj;
    int i,j,v1,v2,x,y;
    int minv,mindeg,deg,goodv;
    long ne;
    
    if (m > 1) NAUTY_ABORT("conncontent only implemented for m=1");

 /* First handle tiny graphs */

    if (n <= 3)
    {
        if (n == 1) return 1;
        if (n == 2) return (g[0] ? -1 : 0);
        if (!g[0] || !g[1] || !g[2]) return 0;    /* disconnected */
        if (g[0]^g[1]^g[2]) return 1;             /* path */
        return 2;                                 /* triangle */
    }

 /* Now compute
      ne = number of edges
      mindeg = minimum degree
      minv = a vertex of minimum degree
      goodv = a vertex with a clique neighbourhood (-1 if none)
 */

    mindeg = n;
    ne = 0;
    goodv = -1;
    for (j = 0; j < n; ++j)
    {
        gj = g[j];
        deg = POPCOUNT(gj);
        ne += deg;
        if (deg < mindeg)
        {
            mindeg = deg;
            minv = j;
            if (deg == 1) goodv = j;
        }
        if (deg >= 3 && deg <= 4 && goodv < 0)
        {
            while (gj)
            {
                TAKEBIT(i,gj);
                if (gj & ~g[i]) break;
            }
            if (!gj) goodv = j;
        }
    }
    ne /= 2;

/* Cases of isolated vertex or tree */

    if (mindeg == 0) return 0;

#if 0
    if (mindeg == 1 && ne == n-1)
    {
        if (isconnected1(g,n)) return ((n&1) ? 1 : -1);
        else                   return 0;
    }
#endif

/* Cases of clique and near-clique */

    if (mindeg == n-1)
    {
        j = -1;
        for (i = 2; i < n; ++i) j *= -i;
        return j;
    }

    if (mindeg == n-2 && n < 16)
    {
        if (!knm_computed)
        {
            knm_computed = TRUE;
            knm[1][0] = 1;
            for (i = 2; i < 16; ++i)
            {
                knm[i][0] = -knm[i-1][0] * (i-1);
                for (j = 1; j+j <= i; ++j)
                    knm[i][j] = knm[i][j-1] + knm[i-1][j-1];
            }
        }
        return knm[n][(n*n-n)/2-ne];
    }

/*  Case of vertex with clique neighbourhood */

    if (goodv >= 0)
    {
        delete1(g,h,goodv,n);
        return -POPCOUNT(g[goodv]) * conncontent(h,m,n-1);
    }

/* Case of minimum degree 2 */

    if (mindeg == 2)
    {
        x = FIRSTBITNZ(g[minv]);
        y = FIRSTBITNZ(g[minv]^bit[x]);
        if (x > minv) --x;
        if (y > minv) --y;
        delete1(g,h,minv,n);
        v1 = conncontent(h,m,n-1);
        if (h[x] & bit[y]) return -2*v1;   /* adjacent neighbours */

        h[x] |= bit[y];
        h[y] |= bit[x];
        v2 = conncontent(h,m,n-1);
        return -v1 - v2;
    }

/* Case of more than 2/3 dense but not complete */

    if (3*ne > n*n-n) 
    {
        j = FIRSTBITNZ(g[minv] ^ bit[minv] ^ ALLMASK(n));   /* non-neighbour */

        g[minv] ^= bit[j];
        g[j] ^= bit[minv];
        v1 = conncontent(g,m,n);
        g[minv] ^= bit[j];
        g[j] ^= bit[minv];

        contract1(g,h,minv,j,n);
        v2 = conncontent(h,m,n-1);
    
        return v1 + v2;
    }
 
/* All remaining cases */

    j = FIRSTBITNZ(g[minv]);     /* neighbour */

    g[minv] ^= bit[j];
    g[j] ^= bit[minv];
    v1 = conncontent(g,m,n);
    g[minv] ^= bit[j];
    g[j] ^= bit[minv];
    
    contract1(g,h,minv,j,n);
    v2 = conncontent(h,m,n-1);

    return v1 - v2;
}

boolean
stronglyconnected(graph *g, int m, int n)
/* test if digraph g is strongly connected */
{
    int sp,v,vc;
    int numvis;
    set *gv;
#if MAXN
    int num[MAXN],lowlink[MAXN],stack[MAXN];
#else
    DYNALLSTAT(int,num,num_sz);
    DYNALLSTAT(int,lowlink,lowlink_sz);
    DYNALLSTAT(int,stack,stack_sz);
#endif

#if !MAXN
    DYNALLOC1(int,num,num_sz,n,"stronglyconnected");
    DYNALLOC1(int,lowlink,lowlink_sz,n,"stronglyconnected");
    DYNALLOC1(int,stack,stack_sz,n,"stronglyconnected");
#endif

    if (n == 0) return FALSE;

    num[0] = 0;
    for (v = 1; v < n; ++v) num[v] = -1;
    lowlink[0] = 0;
    numvis = 1;
    sp = 0;
    stack[0] = 0;
    v = 0;
    vc = -1;
    gv = (set*)g;

    for (;;)
    {
        vc = nextelement(gv,m,vc);
        if (vc < 0)
        {
            if (sp == 0) break;
            if (lowlink[v] == num[v])  return FALSE;
            vc = v;
            v = stack[--sp];
            gv = GRAPHROW(g,v,m);
            if (lowlink[vc] < lowlink[v]) lowlink[v] = lowlink[vc];
        }
        else if (num[vc] < 0)
        {
            stack[++sp] = vc;
            v = vc;
            gv = GRAPHROW(g,v,m);
            vc = -1;
            lowlink[v] = num[v] = numvis++;
        }
        else if (vc != v)
        {
            if (num[vc] < lowlink[v])  lowlink[v] = num[vc];
        }
    }

    return numvis == n;
}

/**************************************************************************/

long
numsquares(graph *g, int m, int n)
/* Number of 4-cycles. Undirected graphs only, loops ok. */
{
    setword w,bitij;
    int i,j,k;
    graph *gi,*gj;
    unsigned long total,t;
    boolean iloop,jloop;

    total = 0;

    if (m == 1)
    {
        for (j = 1; j < n; ++j)
        for (i = 0; i < j; ++i)
        {
            w = (g[i] & g[j]) & ~(bit[i] | bit[j]);
            t = POPCOUNT(w);
            total += t*(t-1)/2;
        }
        return (long)(total / 2);
    }
    else
    {
        for (j = 1, gj = g + m; j < n; ++j, gj += m)
        {
            jloop = ISELEMENT(gj,j);
            if (jloop) DELELEMENT(gj,j);
            for (i = 0, gi = g; i < j; ++i, gi += m)
            {
                iloop = ISELEMENT(gi,i);
                if (iloop) DELELEMENT(gi,i);
                t = 0;
                for (k = 0; k < m; ++k)
                    t += POPCOUNT(gi[k]&gj[k]);
                total += t*(t-1)/2;
                if (iloop) ADDELEMENT(gi,i);
            }
            if (jloop) ADDELEMENT(gj,j);
        }
        return (long)(total / 2);
    }
}

/**************************************************************************/

long
numdiamonds(graph *g, int m, int n)
/* Number of diamonds (squares with diagonal). Undirected only, no loops. */
{
    int i,j,k;
    setword w;
    long ans,c;
    set *gi,*gj;

    ans = 0;

    if (m == 1)
    {
        for (i = 0; i < n; ++i)
        {
            w = g[i] & BITMASK(i);
            while (w)
            {
                TAKEBIT(j,w);
                c = POPCOUNT(g[i]&g[j]);
                ans += c*(c-1)/2;
            }
        }
    }
    else
    {
        for (i = 0, gi = g; i < n; ++i, gi += m)
        {
            for (j = i; (j = nextelement(gi,m,j)) >= 0; )
            {
                gj = g + m*j;
                c = 0;
                for (k = 0; k < m; ++k) c += POPCOUNT(gi[k]&gj[k]);
                ans += c*(c-1)/2;
            }
        }
    }

    return ans;
}

long
numpentagons(graph *g, int m, int n)
/* Number of pentagons (undirected only, no loops) */
{
    int i,j,k,s;
    set *gi,*gj,*gs;
    setword wgi;
    unsigned long count,wsi,wsj,wsij;

    count = 0;

    if (m == 1)
    {
        for (i = 0; i < n; ++i)
        {
            wgi = g[i] & BITMASK(i);
            while (wgi)
            {
                TAKEBIT(j,wgi);
                for (s = 0; s < n; ++s)
                if (s != i && s != j)
                {
                    wsi = POPCOUNT(g[s]&g[i]&~bit[j]);
                    wsj = POPCOUNT(g[s]&g[j]&~bit[i]);
                    count += wsi*wsj - POPCOUNT(g[s]&g[i]&g[j]);
                }
            }
        }
    }
    else
    {
        for (i = 0, gi = g; i < n-1; ++i, gi += m)
        for (j = i; (j = nextelement(gi,m,j)) >= 0; )
        {
            gj = g + j*(size_t)m;

            for (s = 0, gs = g; s < n; ++s, gs += m)
            if (s != i && s != j)
            {
                wsi = wsj = wsij = 0;
                for (k = 0; k < m; ++k)
                {
                    wsi += POPCOUNT(gi[k]&gs[k]);
                    wsj += POPCOUNT(gj[k]&gs[k]);
                    wsij += POPCOUNT(gi[k]&gj[k]&gs[k]);
                }
                if (ISELEMENT(gs,j)) --wsi;
                if (ISELEMENT(gs,i)) --wsj;
                count += wsi*wsj - wsij;
            }
        }
    }

    return (long)(count/5);
}

/**************************************************************************/

int
ktreeness1(graph *g, int n)
/* Version of ktreeness() for m=1. */
{
    int i,v,k,deg[WORDSIZE];
    setword left,degk,w;

    k = n+1;
    for (i = 0; i < n; ++i)
    {
        deg[i] = POPCOUNT(g[i]);
        if (deg[i] < k)
        {
            k = deg[i];
            degk = bit[i];
        }
        else if (deg[i] == k)
            degk |= bit[i];
    }

    if (k == n-1) return n;
    if (k == 0) return 0;

    left = ALLMASK(n);

    while (degk != left && degk != 0)
    {
        TAKEBIT(v,degk);
        if ((g[v] & degk)) return 0;
        left &= ~bit[v];
        w = g[v] & left;
        while (w)
        {
            TAKEBIT(i,w);
            if ((g[i] & w) != w) return 0;
            --deg[i];
            if (deg[i] == k) degk |= bit[i];
        }
    }

    if (degk == 0)                  return 0;
    else if (POPCOUNT(left) == k+1) return k;
    else                            return 0;
}

int
ktreeness(graph *g, int m, int n)
/* Largest k for which g is a k-tree. The complete graph K_n
   is both an (n-1)-tree and an n-tree. This function returns n. */
{
    int i,j,v,k,ndegk,nleft;
    graph *gi;
    DYNALLSTAT(int,deg,deg_sz);
    DYNALLSTAT(set,degk,degk_sz);
    DYNALLSTAT(set,left,left_sz);
    DYNALLSTAT(set,w,w_sz);

    if (m == 1) return ktreeness1(g,n);

    DYNALLOC1(int,deg,deg_sz,n,"ktreeness");
    DYNALLOC1(set,degk,degk_sz,m,"ktreeness");
    DYNALLOC1(set,left,left_sz,m,"ktreeness");
    DYNALLOC1(set,w,w_sz,m,"ktreeness");

    k = n+1;
    for (i = 0, gi = g; i < n; ++i, gi += m)
    {
        SETSIZE(deg[i],gi,m);
        if (deg[i] < k)
        {
            k = deg[i];
            EMPTYSET(degk,m);        
            ADDELEMENT(degk,i);
            ndegk = 1; 
        }
        else if (deg[i] == k)
        {
            ADDELEMENT(degk,i);
            ++ndegk;
        }
    }

    if (k == n-1) return n;
    if (k == 0) return 0;

    FILLSET(left,m,n);
    nleft = n;

    while (ndegk != nleft && ndegk > 0)
    {
        v = nextelement(degk,m,-1);
        DELELEMENT(degk,v);
        --ndegk;
        gi = g + m*v;
        for (i = 0; i < m; ++i)
           if ((gi[i] & degk[i])) return 0;
        DELELEMENT(left,v);
        --nleft;
        for (i = 0; i < m; ++i) w[i] = gi[i] & left[i];

        for (i = -1; (i = nextelement(w,m,i)) >= 0; )
        {
            DELELEMENT(w,i);
            gi = g + m*i;
            for (j = 0; j < m; ++j)
                if ((gi[j] & w[j]) != w[j]) return 0;
            --deg[i];
            if (deg[i] == k) { ADDELEMENT(degk,i); ++ndegk; }
        }
    }

    if (ndegk == 0)        return 0;
    else if (nleft == k+1) return k;
    else                   return 0;
}
