/*****************************************************************************
*    This is the header file for Version 2.0 of nauty().                     *
*****************************************************************************/

#ifndef  NAUTYH_READ    /* only process this file once */

#include <stdio.h>

/* Exactly one of the symbols with names starting with "SYS_" should
   have the value 1.  All the others should have the value 0. */

#define SYS_VAXBSD     0   /* older BSD unix 4.2 on a VAX */
#define SYS_CRAY       0   /* Cray UNIX, portable or standard C */
#define SYS_APOLLO     0   /* DOMAIN C on Apollo */
#define SYS_BSDUNIX    0   /* other BSD unix or Linux */
#define SYS_UNIX       1   /* miscellaneous non-BSD unix, including A/UX */
#define SYS_ALPHA      0   /* Alpha UNIX */
#define SYS_ALPHA32    0   /* Alpha UNIX using 32-bit setwords */
#define SYS_VAXVMS     0   /* VAX11C or GCC on a VAX under VMS */
#define SYS_MACLSC     0   /* Lightspeed C on an Apple Macintosh (Version 1) */
#define SYS_MACTHINK   0   /* THINK C (v4-5) or CodeWarrior (v6+) on a Mac */
#define SYS_MACAZT     0   /* Aztec C on an Apple Macintosh */
#define SYS_MACMPW     0   /* Apple Macintosh Workshop C */
#define SYS_AMIGALC    0   /* Lattice C on a Commodore Amiga */
#define SYS_AMIGAAZT   0   /* Aztec C on a Commodore Amiga */
#define SYS_PCMS4      0   /* Microsoft C 4.0 on IBM PC */
#define SYS_PCMS5      0   /* Microsoft C 5.1 or Quick C 2.5 on IBM PC */
#define SYS_PCTURBO    0   /* Turbo C on IBM PC */
#define SYS_IBMC       0   /* IBM C Set/2 under OS/2 */
#define SYS_MISC       0   /* anything else */

/*****************************************************************************
*                                                                            *
*    AUTHOR: Brendan D. McKay                                                *
*            Computer Science Department                                     *
*            Australian National University                                  *
*            Canberra, Australia 0200                                        *
*            phone:  +61 2 6249 3845    fax:  +61 2 6249 0010                *
*            email:  bdm@cs.anu.edu.au                                       *
*                                                                            *
*   Copyright (1984-2000) Brendan McKay.  All rights reserved.  Permission   *
*   is hereby given for use and/or distribution with the exception of        *
*   sale for profit or application with nontrivial military significance.    *
*   You must not remove this copyright notice, and you must document any     *
*   changes that you make to this program.                                   *
*   This software is subject to this copyright only, irrespective of         *
*   any copyright attached to any package of which this is a part.           *
*                                                                            *
*   This program is only provided "as is".  No responsibility will be taken  *
*   by the author, his employer or his pet rabbit* for any misfortune which  *
*   befalls you because of its use.  I don't think it will delete all your   *
*   files, burn down your computer room or turn your children against you,   *
*   but if it does: stiff cheddar.  On the other hand, I very much welcome   *
*   bug reports, or at least I would if there were any bugs.                 *
*                                                       * RIP, 1989          *
*                                                                            *
*   If you wish to acknowledge use of this program in published articles,    *
*   please do so by citing the User's Guide:                                 *
*                                                                            *
*     B. D. McKay, nauty User's Guide (Version 1.5), Technical Report        *
*         TR-CS-90-02, Australian National University, Department of         *
*         Computer Science, 1990.                                            *
*                                                                            *
*   CHANGE HISTORY                                                           *
*       10-Nov-87 : final changes for version 1.2                            *
*        5-Dec-87 : renamed to version 1.3 (no changes to this file)         *
*       28-Sep-88 : added PC Turbo C support, making version 1.4             *
*       23-Mar-89 : changes for version 1.5 :                                *
*                   - reworked M==1 code                                     *
*                   - defined NAUTYVERSION string                            *
*                   - made NAUTYH_READ to allow this file to be read twice   *
*                   - added optional ANSI function prototypes                *
*                   - added validity check for WORDSIZE                      *
*                   - added new fields to optionblk structure                *
*                   - updated DEFAULTOPTIONS to add invariants fields        *
*                   - added (set*) cast to definition of GRAPHROW            *
*                   - added definition of ALLOCS and FREES                   *
*       25-Mar-89 : - added declaration of new function doref()              *
*                   - added UNION macro                                      *
*       29-Mar-89 : - reduced the default MAXN for small machines            *
*                   - removed OUTOFSPACE (no longer used)                    *
*                   - added SETDIFF and XOR macros                           *
*        2-Apr-89 : - extended statsblk structure                            *
*        4-Apr-89 : - added IS_* macros                                      *
*                   - added ERRFILE definition                               *
*                   - replaced statsblk.outofspace by statsblk.errstatus     *
*        5-Apr-89 : - deleted definition of np2vector (no longer used)       *
*                   - introduced EMPTYSET macro                              *
*       12-Apr-89 : - eliminated MARK, UNMARK and ISMARKED (no longer used)  *
*       18-Apr-89 : - added MTOOBIG and CANONGNIL                            *
*       12-May-89 : - made ISELEM1 and ISELEMENT return 0 or 1               *
*        2-Mar-90 : - added EXTPROC macro and used it                        *
*       12-Mar-90 : - added SYS_CRAY, with help from N. Sloane and A. Grosky *
*                   - added dummy groupopts field to optionblk               *
*                   - select some ANSI things if __STDC__ exists             *
*       20-Mar-90 : - changed default MAXN for Macintosh versions            *
*                   - created SYS_MACTHINK for Macintosh THINK compiler      *
*       27-Mar-90 : - split SYS_MSDOS into SYS_PCMS4 and SYS_PCMS5           *
*       13-Oct-90 : changes for version 1.6:                                 *
*                   - fix definition of setword for WORDSIZE==64             *
*       14-Oct-90 : - added SYS_APOLLO version to avoid compiler bug         *
*       15-Oct-90 : - improve detection of ANSI conformance                  *
*       17-Oct-90 : - changed temp name in EMPTYSET to avoid A/UX bug        *
*       16-Apr-91 : changes for version 1.7:                                 *
*                   - made version SYS_PCTURBO use free(), not cfree()       *
*        2-Sep-91 : - noted that SYS_PCMS5 also works for Quick C            *
*                   - moved MULTIPLY to here from nauty.c                    *
*       12-Jun-92 : - changed the top part of this comment                   *
*       27-Aug-92 : - added version SYS_IBMC, thanks to Ivo Duentsch         *
*        5-Jun-93 : - renamed to version 1.7+, only change in naututil.h     *
*       29-Jul-93 : changes for version 1.8:                                 *
*                   - fixed error in default 64-bit version of FIRSTBIT      *
*                     (not used in any version before ALPHA)                 *
*                   - installed ALPHA version (thanks to Gordon Royle)       *
*                   - defined ALLOCS,FREES for SYS_IBMC                      *
*        3-Sep-93 : - make calloc void* in ALPHA version                     *
*       17-Sep-93 : - renamed to version 1.9,                                *
*                        changed only dreadnaut.c and nautinv.c              *
*       24-Feb-94 : changes for version 1.10:                                *
*                   - added version SYS_AMIGAAZT, thanks to Carsten Saager   *
*                     (making 1.9+)                                          *
*       19-Apr-95 : - added prototype wrapper for C++,                       *
*                     thanks to Daniel Huson                                 *
*        5-Mar-96 : - added SYS_ALPHA32 version (32-bit setwords on Alpha)   *
*       13-Jul-96 : changes for version 2.0:                                 *
*                   - added dynamic allocation                               *
*                   - ERRFILE must be defined                                *
*                   - added FLIPELEM1 and FLIPELEMENT macros                 *
*       13-Aug-96 : - added SWCHUNK? macros                                  *
*                   - added TAKEBIT macro                                    *
*       28-Nov-96 : - include sys/types.h if not ANSI (tentative!)           *
*       24-Jan-97 : - and stdlib.h if ANSI                                   *
*                   - removed use of cfree() from UNIX variants              *
*       25-Jan-97 : - changed options.getcanon from boolean to int           *
*                     Backwards compatibility is ok, as boolean and int      *
*                     are the same.  Now getcanon=2 means to get the label   *
*                     and not care about the group.  Sometimes faster.       *
*        6-Feb-97 : - Put in #undef for FALSE and TRUE to cope with          *
*                     compilers that illegally predefine them.               *
*                   - declared nauty_null and nautil_null                    *
*        2-Jul-98 : - declared ALLBITS                                       *
*       21-Oct-98 : - allow WORDSIZE==64 using unsigned long long            *
*                   - added BIGNAUTY option for really big graphs            *
*       11-Dec-99 : - made bit, leftbit and bytecount static in each file    *
*        9-Jan-00 : - declared nauty_check() and nautil_check()              *
*       12-Feb-00 : - Used #error for compile-time checks                    *
*                   - Added DYNREALLOC                                       *
*        4-Mar-00 : - declared ALLMASK(n)                                    *
*       27-May-00 : - declared CONDYNFREE                                    *
*       28-May-99 : - declared nautil_freedyn()                              *
*                                                                            *
*****************************************************************************/

#define IS_UNIX  (SYS_VAXBSD | SYS_BSDUNIX | SYS_CRAY | SYS_UNIX \
   | SYS_APOLLO | SYS_ALPHA | SYS_ALPHA32)
#define IS_MAC   (SYS_MACMPW | SYS_MACAZT | SYS_MACLSC | SYS_MACTHINK)
#define IS_PC    (SYS_PCTURBO | SYS_PCMS4 | SYS_PCMS5)
#define IS_AMIGA (SYS_AMIGALC | SYS_AMIGAAZT)

/*****************************************************************************
*                                                                            *
*   16-bit, 32-bit and 64-bit versions can be selected by defining WORDSIZE. *
*   The largest graph that can be handled has MAXN vertices.                 *
*   Both WORDSIZE and MAXN can be defined on the command line.               *
*   WORDSIZE must be 16, 32 or 64; MAXN must be <= INFINITY - 2;             *
*                                                                            *
*   With a very slight loss of efficiency (depending on platform), nauty     *
*   can be compiled to dynamically allocate arrays.  Predefine MAXN=0 to     *
*   achieve this effect, which is default behaviour from version 2.0.        *
*   In that case, graphs of size up to INFINITY-2 can be handled if the      *
*   the memory is available.                                                 *
*                                                                            *
*   If only very small graphs need to be processed, use MAXN<=WORDSIZE       *
*   since this causes substantial code optimizations.                        *
*   The symbol EXTDEFS should be defined (with arbitrary value) before this  *
*   file is included, except for one of each group of files linked together. *
*   The file containing the main program is suggested.  This ensures that    *
*   arrays bit, bytecount and leftbit, and parameter labelorg, are only      *
*   initialised once, which some linkers think is necessary.                 *
*                                                                            *
*   Conventions and Assumptions:                                             *
*                                                                            *
*    A 'setword' is the chunk of memory that is occupied by one part of      *
*    a set.  This is assumed to be >= WORDSIZE bits in size.                 *
*                                                                            *
*    The rightmost (loworder) WORDSIZE bits of setwords are numbered         *
*    0..WORDSIZE-1, left to right.  It is necessary that the 2^WORDSIZE      *
*    setwords with the other bits zero are totally ordered under <,=,>.      *
*    This needs care on a 1's-complement machine.                            *
*                                                                            *
*    The int variables m and n have consistent meanings throughout.          *
*    Graphs have n vertices always, and sets have m setwords always.         *
*                                                                            *
*    A 'set' consists of m contiguous setwords, whose bits are numbered      *
*    0,1,2,... from left (high-order) to right (low-order), using only       *
*    the rightmost WORDSIZE bits of each setword.  It is used to             *
*    represent a subset of {0,1,...,n-1} in the usual way - bit number x     *
*    is 1 iff x is in the subset.  Bits numbered n or greater, and           *
*    unnumbered bits, are assumed permanently zero.                          *
*                                                                            *
*    A 'graph' consists of n contiguous sets.  The i-th set represents       *
*    the vertices adjacent to vertex i, for i = 0,1,...,n-1.                 *
*                                                                            *
*    A 'permutation' is an array of n short ints (ints in the case that      *
*    BIGNAUTY is defined) , repesenting a permutation of the set             *
*    {0,1,...,n-1}.  The value of the i-th entry is the number to which      *
*    i is mapped.                                                            *
*                                                                            *
*    If g is a graph and p is a permutation, then g^p is the graph in        *
*    which vertex i is adjacent to vertex j iff vertex p[i] is adjacent      *
*    to vertex p[j] in g.                                                    *
*                                                                            *
*    An 'nvector' is any array of n ints.  Actually, nvector==int.           *
*                                                                            *
*    A partition nest is represented by a pair (lab,ptn), where lab and ptn  *
*    are nvectors.  The "partition at level x" is the partition whose cells  *
*    are {lab[i],lab[i+1],...,lab[j]}, where [i,j] is a maximal subinterval  *
*    of [0,n-1] such that ptn[k] > x for i <= k < j and ptn[j] <= x.         *
*    The partition at level 0 is given to nauty by the user.  This is        *
*    refined for the root of the tree, which has level 1.                    *
*                                                                            *
*****************************************************************************/

#if  SYS_VAXBSD
#define NAUTYVERSION "2.0beta (VAXBSD)"
#endif
#if  SYS_BSDUNIX
#define NAUTYVERSION "2.0beta (BSDUNIX)"
#endif
#if  SYS_UNIX
#define NAUTYVERSION "2.0beta (UNIX)"
#endif
#if  SYS_APOLLO
#define NAUTYVERSION "2.0beta (APOLLO)"
#endif
#if  SYS_VAXVMS
#define NAUTYVERSION "2.0beta (VAXVMS)"
#endif
#if  SYS_MACLSC
#define NAUTYVERSION "2.0beta (MACLSC)"
#endif
#if  SYS_MACTHINK
#define NAUTYVERSION "2.0beta (MACTHINK)"
#endif
#if  SYS_MACAZT
#define NAUTYVERSION "2.0beta (MACAZT)"
#endif
#if  SYS_MACMPW
#define NAUTYVERSION "2.0beta (MACMPW)"
#endif
#if  SYS_AMIGALC
#define NAUTYVERSION "2.0beta (AMIGALC)"
#endif
#if  SYS_AMIGAAZT
#define NAUTYVERSION "2.0beta (AMIGAAZT)"
#endif
#if  SYS_PCMS4
#define NAUTYVERSION "2.0beta (PCMS4)"
#endif
#if  SYS_PCMS5
#define NAUTYVERSION "2.0beta (PCMS5)"
#endif
#if  SYS_PCTURBO
#define NAUTYVERSION "2.0beta (PCTURBO)"
#endif
#if  SYS_CRAY
#define NAUTYVERSION "2.0beta (CRAY)"
#endif
#if  SYS_IBMC
#define NAUTYVERSION "2.0beta (IBM C Set/2 for OS/2)"
#endif
#if  SYS_ALPHA
#define NAUTYVERSION "2.0beta (Alpha)"
#endif
#if  SYS_ALPHA32
#define NAUTYVERSION "2.0beta (Alpha 32 bits)"
#endif
#if  SYS_MISC
#define NAUTYVERSION "2.0beta (MISC)"
#endif

#ifdef BIGNAUTY
#define NAUTYVERSIONID 2101      /* 1000 times the version number */
#else                            /* + 1 for BIGNAUTY */
#define NAUTYVERSIONID 2100
#endif

/* WORDSIZE is the number of set elements per setword (16, 32 or 64) */

#ifdef  WORDSIZE

#if  ((WORDSIZE != 16) & (WORDSIZE != 32) & (WORDSIZE != 64))
#error "WORDSIZE must be 16, 32 or 64"
#endif

#ifdef BIGNAUTY
#if  WORDSIZE == 16
#error "BIGNAUTY requires WORDSIZE 32 or 64"
#endif
#endif

#else  /* WORDSIZE undefined */

#if  IS_PC
#ifdef BIGNAUTY
#define WORDSIZE 32
#else
#define WORDSIZE 16
#endif
#else
#if  (SYS_CRAY | SYS_ALPHA)
#define WORDSIZE 64
#else
#define WORDSIZE 32
#endif
#endif

#endif  /* WORDSIZE */

#ifndef  MAXN                   /* maximum allowed n value */
			        /* Use 0 for dynamic sizing. */
#define MAXN 0
#define MAXM 0
#else
#define MAXM ((MAXN+WORDSIZE-1)/WORDSIZE)  /* max setwords in a set */
#endif  /* MAXN */

    /* set operations (setadd is its address, pos is the bit number): */
#if  MAXM==1

#define SETWD(pos) 0
#define SETBT(pos) (pos)

#else   /* MAXM > 1 */

#if  WORDSIZE==64
#define SETWD(pos) ((pos)>>6)    /* number of setword containing bit pos */
#define SETBT(pos) ((pos)&077)   /* position within setword of bit pos */
#else
#if  WORDSIZE==32
#define SETWD(pos) ((pos)>>5)
#define SETBT(pos) ((pos)&037)
#else   /* WORDSIZE==16 */
#define SETWD(pos) ((pos)>>4)
#define SETBT(pos) ((pos)&017)
#endif
#endif

#endif  /* MAXM */

#if  WORDSIZE==64
#define TIMESWORDSIZE(w) ((w)<<6)    /* w*WORDSIZE */
#else
#if  WORDSIZE==32
#define TIMESWORDSIZE(w) ((w)<<5)
#else  /* WORDSIZE==16 */
#define TIMESWORDSIZE(w) ((w)<<4)
#endif
#endif

#define ADDELEM1(setadd,pos)  (*(setadd) |= bit[pos])
#define DELELEM1(setadd,pos)  (*(setadd) &= ~bit[pos])
#define FLIPELEM1(setadd,pos) (*(setadd) ^= bit[pos])
#define ISELEM1(setadd,pos)   ((*(setadd) & bit[pos]) != 0)

#define ADDELEMENT(setadd,pos)  ((setadd)[SETWD(pos)] |= bit[SETBT(pos)])
#define DELELEMENT(setadd,pos)  ((setadd)[SETWD(pos)] &= ~bit[SETBT(pos)])
#define FLIPELEMENT(setadd,pos) ((setadd)[SETWD(pos)] ^= bit[SETBT(pos)])
#define ISELEMENT(setadd,pos)   (((setadd)[SETWD(pos)] & bit[SETBT(pos)]) != 0)

#if  MAXM==1            /* initialise set to empty */
#define EMPTYSET(setadd,m) *(setadd) = 0;
#else
#define EMPTYSET(setadd,m) \
    {register setword *es; \
    for (es = (setword*)(setadd)+(m); --es >= (setword*)(setadd);) *es=0;}
#endif

#if  MAXM==1            /* obsolete version of EMPTYSET */
#define MAKEEMPTY(setadd,m,i) *(setadd) = 0
#else
#define MAKEEMPTY(setadd,m,i) for(i=m;--i>=0;)setadd[i]=0
#endif

#define NOTSUBSET(word1,word2) ((word1) & ~(word2))  /* test if the 1-bits
                    in setword word1 do not form a subset of those in word2  */
#define INTERSECT(word1,word2) ((word1) &= (word2))  /* AND word2 into word1 */
#define UNION(word1,word2)     ((word1) |= (word2))  /* OR word2 into word1 */
#define SETDIFF(word1,word2)   ((word1) &= ~(word2)) /* - word2 into word1 */
#define XOR(word1,word2)       ((word1) ^= (word2))  /* XOR word2 into word1 */
#define ZAPBIT(word,x) ((word) &= ~bit[x])  /* delete bit x in setword */
#define TAKEBIT(iw,w) {iw = FIRSTBIT(w); w ^= bit[iw];}

#if  WORDSIZE==64
#if  SYS_CRAY
#define POPCOUNT(x) _popcnt(x)
#define FIRSTBIT(x) _leadz(x)
#define BITMASK(x)  _mask(65+(x))
#else
#define POPCOUNT(x) (bytecount[(x)>>56 & 0377] + bytecount[(x)>>48 & 0377] \
                   + bytecount[(x)>>40 & 0377] + bytecount[(x)>>32 & 0377] \
                   + bytecount[(x)>>24 & 0377] + bytecount[(x)>>16 & 0377] \
                   + bytecount[(x)>>8 & 0377]  + bytecount[(x) & 0377])
        /* number of 1-bits in a setword */
#define FIRSTBIT(x) ((x) & 01777777777740000000000 ? \
                       (x) & 01777770000000000000000 ? \
                         (x) & 01774000000000000000000 ? \
                         0+leftbit[((x)>>56) & 0377] : \
                         8+leftbit[(x)>>48] \
                       : (x) & 07760000000000000 ? \
                         16+leftbit[(x)>>40] : \
                         24+leftbit[(x)>>32] \
                     : (x) & 037777600000 ? \
                         (x) & 037700000000 ? \
                         32+leftbit[(x)>>24] : \
                         40+leftbit[(x)>>16] \
                       : (x) & 0177400 ? \
                         48+leftbit[(x)>>8] : \
                         56+leftbit[x])
        /* get number of first 1-bit in non-zero setword (0..WORDSIZE-1) */
#define BITMASK(x)  (0777777777777777777777 >> (x)) /* setword whose rightmost
  WORDSIZE-x-1 (numbered) bits are 1 and the rest 0 (0 <= x < WORDSIZE) */
#endif /* CRAY */
#if  (SYS_CRAY | SYS_ALPHA)
#define SWCHUNK0(w) ((long)((w)>>48)&0177777L)
#define SWCHUNK1(w) ((long)((w)>>32)&0177777L)
#define SWCHUNK2(w) ((long)((w)>>16)&0177777L)
#define SWCHUNK3(w) ((long)(w)&0177777L)
#else
#define SWCHUNK0(w) ((long)((w)>>48)&0177777LL)
#define SWCHUNK1(w) ((long)((w)>>32)&0177777LL)
#define SWCHUNK2(w) ((long)((w)>>16)&0177777LL)
#define SWCHUNK3(w) ((long)(w)&0177777LL)
#endif
#else
#if  WORDSIZE==32
#define POPCOUNT(x) (bytecount[(x)>>24 & 0377] + bytecount[(x)>>16 & 0377] \
                        + bytecount[(x)>>8 & 0377] + bytecount[(x) & 0377])
#define FIRSTBIT(x) ((x) & 037777600000 ? ((x) & 037700000000 ? \
                     leftbit[((x)>>24) & 0377] : 8+leftbit[(x)>>16]) \
                    : ((x) & 0177400 ? 16+leftbit[(x)>>8] : 24+leftbit[x]))
#define BITMASK(x)  (017777777777 >> (x)) /* setword whose rightmost
  WORDSIZE-x-1 (numbered) bits are 1 and the rest 0 (0 <= x < WORDSIZE) */
#define SWCHUNK0(w) ((long)((w)>>16)&0177777L)
#define SWCHUNK1(w) ((long)(w)&0177777L)
#else   /* WORDSIZE==16 */
#define POPCOUNT(x) (bytecount[(x)>>8 & 0377] + bytecount[(x) & 0377])
#define FIRSTBIT(x) ((x) & 0177400 ? leftbit[((x)>>8) & 0377] : 8+leftbit[x])
#define BITMASK(x)  (077777 >> (x))
#define SWCHUNK0(w) ((long)(w)&0177777L)
#endif
#endif

#define ALLBITS (~(setword)0)   /* All numbered bits. */
#define ALLMASK(n) (~(setword)BITMASK((n)-1))   /* First n>0 bits */

#if  MAXM==1
#define GRAPHROW(g,v,m) ((set*)(g) + (v))
#else
#define GRAPHROW(g,v,m) ((set*)(g) + (long)(v) * (long)(m))
    /* address of row v of graph g.  v and m are ints.  Beware of v*m being
       too large for an int. */
#endif

    /* various constants: */
#ifdef FALSE
#undef FALSE
#undef TRUE
#endif
#define FALSE    0
#define TRUE     1
#ifdef BIGNAUTY
#define INFINITY 07777777       /* positive int greater than MAXN+2 */
typedef int shortish;
#else
#define INFINITY 077777         /* positive short int greater than MAXN+2 */
typedef short shortish;
#endif

    /* typedefs for sets, graphs, permutations, etc.: */
typedef int boolean;    /* boolean MUST be the same as int */

#if WORDSIZE == 64
#if (SYS_ALPHA | SYS_CRAY)
typedef unsigned long setword;
#else
typedef unsigned long long setword;
#endif

#else  /* WORDSIZE != 64 */

#if  SYS_ALPHA32
typedef unsigned int setword;
#else
#if  WORDSIZE==32
typedef unsigned long setword;
#else
typedef unsigned short setword;
#endif
#endif

#endif  /* WORDSIZE == 64 */

#if  SYS_VAXBSD
#define UPROC int       /* avoids compiler bug in BSD 4.2 on VAX */
#else
#define UPROC void      /* type of user-defined procedures */
#endif

typedef setword set,graph;
typedef int nvector,np2vector;
typedef shortish permutation;

typedef struct
{
    int group_xyz;            /* dummy field to make it non-empty */
} groupblk;

typedef struct
{
    int getcanon;             /* make canong and canonlab? */
#define LABELONLY 2   /* new value allowed at version 2.0 */
    boolean digraph;          /* multiple edges or loops? */
    boolean writeautoms;      /* write automorphisms? */
    boolean writemarkers;     /* write stats on pts fixed, etc.? */
    boolean defaultptn;       /* set lab,ptn,active for single cell? */
    boolean cartesian;        /* use cartesian rep for writing automs? */
    int linelength;           /* max chars/line (excl. '\n') for output */
    FILE *outfile;            /* file for output, if any */
    UPROC (*userrefproc)();   /* replacement for usual refine procedure */
    UPROC (*userautomproc)(); /* procedure called for each automorphism */
    UPROC (*userlevelproc)(); /* procedure called for each level */
    UPROC (*usernodeproc)();  /* procedure called for each node */
    UPROC (*usertcellproc)(); /* replacement for targetcell procedure */
    UPROC (*invarproc)();     /* procedure to compute vertex-invariant */
    int tc_level;             /* max level for smart target cell choosing */
    int mininvarlevel;        /* min level for invariant computation */
    int maxinvarlevel;        /* max level for invariant computation */
    int invararg;             /* value passed to (*invarproc)() */
    groupblk *groupopts;      /* placeholder for future group options */
} optionblk;

#ifndef CONSOLWIDTH
#if (SYS_MACAZT | SYS_MACMPW)
#define CONSOLWIDTH 72
#else
#if  IS_AMIGA
#define CONSOLWIDTH 75
#else
#define CONSOLWIDTH 78
#endif
#endif
#endif

#define NILFUNCTION ((UPROC(*)())NULL)     /* nil pointer to user-function */
#define NILSET      ((set*)NULL)           /* nil pointer to set */
#define NILGRAPH    ((graph*)NULL)         /* nil pointer to graph */

#define DEFAULTOPTIONS(options) optionblk options = {0,FALSE,TRUE,TRUE,\
    TRUE,FALSE,CONSOLWIDTH,(FILE*)NULL,NILFUNCTION,NILFUNCTION,NILFUNCTION,\
    NILFUNCTION,NILFUNCTION,NILFUNCTION,0,0,0,0,(groupblk*)NULL}

#if  IS_AMIGA
#define PUTC(c,f) fputc(c,f)    /* best way to write character c to file f */
#else
#define PUTC(c,f) putc(c,f)
#endif

typedef struct
{
    double grpsize1;        /* size of group is */
    int grpsize2;           /*    grpsize1 * 10^grpsize2 */
#if  !SYS_MACAZT
#define groupsize1 grpsize1     /* for backwards compatibility */
#define groupsize2 grpsize2
#endif
    int numorbits;          /* number of orbits in group */
    int numgenerators;      /* number of generators found */
    int errstatus;          /* if non-zero : an error code */
#define outofspace errstatus;   /* for backwards compatibility */
    long numnodes;          /* total number of nodes */
    long numbadleaves;      /* number of leaves of no use */
    int maxlevel;           /* maximum depth of search */
    long tctotal;           /* total size of all target cells */
    long canupdates;        /* number of updates of best label */
    long invapplics;        /* number of applications of invarproc */
    long invsuccesses;      /* number of successful applics of invarproc() */
    int invarsuclevel;      /* least level where invarproc worked */
} statsblk;

/* codes for errstatus field: */
#define NTOOBIG      1      /* n > MAXN or n > WORDSIZE*m */
#define MTOOBIG      2      /* m > MAXM */
#define CANONGNIL    3      /* canong = NILGRAPH, but getcanon = TRUE */

/* manipulation of real approximation to group size */
#define MULTIPLY(s1,s2,i) if ((s1 *= i) >= 1e10) {s1 /= 1e10; s2 += 10;}

#define ANSI_STDC 0
#ifdef __STDC__
#if    __STDC__ == 1
#undef  ANSI_STDC
#define ANSI_STDC 1
#endif
#endif

#if ANSI_STDC
#include <stddef.h>
#include <stdlib.h>
#endif

/* ALLOCS(x,y) should return a pointer (any pointer type) to x*y units of new
   storage, not necessarily initialised.  A "unit" of storage is defined by
   the sizeof operator.   x and y are integer values of type int or larger, 
   but x*y may well be too large for an int.  The macro should cast to the
   correct type for the call.  On failure, ALLOCS(x,y) should return a NULL 
   pointer.  FREES(p) should free storage previously allocated by ALLOCS, 
   where p is the value that ALLOCS returned.  */

#if  SYS_MACAZT
extern char *NewPtr();
extern void DisposPtr();
#define ALLOCS(x,y)  NewPtr((long)(x) * (long)(y))
#define FREES(p)     DisposPtr(p)
#endif

#if  (IS_UNIX | SYS_PCTURBO | SYS_PCMS4 | SYS_VAXVMS \
              | SYS_MACLSC | SYS_MISC | SYS_IBMC)
#if  ANSI_STDC
#if  SYS_APOLLO
#include <sys/types.h>          /* Domain C has no stddef.h */
#endif  /* SYS_APOLLO */
#define ALLOCS(x,y)  calloc(x,y)
#define FREES(p)  free(p)
#else  /* not ANSI_STDC */
#include <sys/types.h>
#if  (SYS_ALPHA | SYS_ALPHA32)
extern void *calloc();
#else
extern char *calloc();
#endif  /* SYS_ALPHA */
extern void free();
#define ALLOCS(x,y)  calloc(x,y)
#define FREES(p)  free(p)
#endif  /* ANSI_STDC */
#endif  /* IS_UNIX ... */

#if SYS_PCMS5
#include <malloc.h>
#define ALLOCS(x,y)  calloc(x,y)
#define FREES(p)  free(p)
#endif

#if  SYS_MACTHINK
#include <stddef.h>
extern void *calloc(size_t,size_t);
extern void free(void*);
#define ALLOCS(x,y)  calloc(x,y)
#define FREES(p)  free(p)
#endif

#if  (SYS_MACMPW | IS_AMIGA)
extern char *calloc();
extern void free();
#define ALLOCS(x,y)  calloc(x,y)
#define FREES(p)  free(p)
#endif


/* The following macros are used by nauty if MAXN=0.  They dynamically
   allocate arrays of size dependent on m or n.  For each array there
   should be two static variables:
     type *name;
     size_t name_sz;
   "name" will hold a pointer to an allocated array.  "name_sz" will hold
   the size of the allocated array in units of sizeof(type).  DYNALLSTAT
   declares both variables and initialises name_sz=0.  DYNALLOC1 and
   DYNALLOC2 test if there is enough space allocated, and if not free
   the existing space and allocate a bigger space.  The allocated space
   is not initialised.
   
   In the case of DYNALLOC1, the space is allocated using
       ALLOCS(sz,sizeof(type)).
   In the case of DYNALLOC2, the space is allocated using
       ALLOCS(sz1,sz2*sizeof(type)).

   DYNREALLOC is like DYNALLOC1 except that the old contents are copied
   into the new space.  realloc() is assumed.  This is not currently
   used by nauty or dreadnaut.

   DYNFREE frees any allocated array and sets name_sz back to 0.
   CONDYNFREE does the same, but only if name_sz exceeds some limit.
*/

#define DYNALLSTAT(type,name,name_sz) static type *name; static size_t name_sz=0
#define DYNALLOC1(type,name,name_sz,sz,msg) \
 if ((size_t)(sz) > name_sz) \
 { if (name_sz) FREES(name); name_sz = (sz); \
 if ((name=(type*)ALLOCS(sz,sizeof(type))) == NULL) alloc_error(msg);}
#define DYNALLOC2(type,name,name_sz,sz1,sz2,msg) \
 if ((size_t)(sz1)*(size_t)(sz2) > name_sz) \
 { if (name_sz) FREES(name); name_sz = (size_t)(sz1)*(size_t)(sz2); \
 if ((name=(type*)ALLOCS((sz1),(sz2)*sizeof(type))) == NULL) alloc_error(msg);}
#define DYNREALLOC(type,name,name_sz,sz,msg) \
 {if ((size_t)(sz) > name_sz) \
 { if ((name = (type*)realloc(name,(sz)*sizeof(type))) == NULL) \
      alloc_error(msg); else name_sz = (sz);}}
#define DYNFREE(name,name_sz) if (name_sz) {FREES(name); name_sz = 0;}
#define CONDYNFREE(name,name_sz,minsz) \
 if (name_sz > (size_t)(szlim)) {FREES(name); name_sz = 0;}

/* File to write error messages to (used as first argument to fprintf()). */
#define ERRFILE stderr

#ifdef OLDEXTDEFS
#define EXTDEF_CLASS
#ifdef EXTDEFS
#define EXTDEF_TYPE 1
#else
#define EXTDEF_TYPE 2
#endif
#else
#define EXTDEF_CLASS static
#define EXTDEF_TYPE 2
#endif

#ifdef EXTDEFS
extern int labelorg;
#else
int labelorg = 0;       /* number of least-numbered vertex */
#endif

#if  EXTDEF_TYPE==1
extern setword bit[];
extern int bytecount[];
extern int leftbit[];

#else
    /* array giving setwords with single 1-bit */
#if  WORDSIZE==64
#if  !(SYS_ALPHA | SYS_CRAY)
EXTDEF_CLASS
setword bit[] = {01000000000000000000000LL,0400000000000000000000LL,
                 0200000000000000000000LL,0100000000000000000000LL,
                 040000000000000000000LL,020000000000000000000LL,
                 010000000000000000000LL,04000000000000000000LL,
                 02000000000000000000LL,01000000000000000000LL,
                 0400000000000000000LL,0200000000000000000LL,
                 0100000000000000000LL,040000000000000000LL,
                 020000000000000000LL,010000000000000000LL,
                 04000000000000000LL,02000000000000000LL,
                 01000000000000000LL,0400000000000000LL,0200000000000000LL,
                 0100000000000000LL,040000000000000LL,020000000000000LL,
                 010000000000000LL,04000000000000LL,02000000000000LL,
                 01000000000000LL,0400000000000LL,0200000000000LL,
		 0100000000000LL,040000000000LL,020000000000LL,010000000000LL,
		 04000000000LL,02000000000LL,01000000000LL,0400000000LL,
		 0200000000LL,0100000000LL,040000000LL,020000000LL,
		 010000000LL,04000000LL,02000000LL,01000000LL,0400000LL,
		 0200000LL,0100000LL,040000LL,020000LL,010000LL,04000LL,
                 02000LL,01000LL,0400LL,0200LL,0100LL,040LL,020LL,010LL,
		 04LL,02LL,01LL};
#else
EXTDEF_CLASS
setword bit[] = {01000000000000000000000,0400000000000000000000,
                 0200000000000000000000,0100000000000000000000,
                 040000000000000000000,020000000000000000000,
                 010000000000000000000,04000000000000000000,
                 02000000000000000000,01000000000000000000,
                 0400000000000000000,0200000000000000000,
                 0100000000000000000,040000000000000000,020000000000000000,
                 010000000000000000,04000000000000000,02000000000000000,
                 01000000000000000,0400000000000000,0200000000000000,
                 0100000000000000,040000000000000,020000000000000,
                 010000000000000,04000000000000,02000000000000,
                 01000000000000,0400000000000,0200000000000,0100000000000,
                 040000000000,020000000000,010000000000,04000000000,
                 02000000000,01000000000,0400000000,0200000000,0100000000,
                 040000000,020000000,010000000,04000000,02000000,01000000,
                 0400000,0200000,0100000,040000,020000,010000,04000,
                 02000,01000,0400,0200,0100,040,020,010,04,02,01};
#endif
#else
#if  WORDSIZE==32
EXTDEF_CLASS
setword bit[] = {020000000000,010000000000,04000000000,02000000000,
                 01000000000,0400000000,0200000000,0100000000,040000000,
                 020000000,010000000,04000000,02000000,01000000,0400000,
                 0200000,0100000,040000,020000,010000,04000,02000,01000,
                 0400,0200,0100,040,020,010,04,02,01};
#else   /* WORDSIZE==16 */
EXTDEF_CLASS
setword bit[] = {0100000,040000,020000,010000,04000,02000,01000,0400,0200,
                 0100,040,020,010,04,02,01};
#endif
#endif

    /*  array giving number of 1-bits in bytes valued 0..255: */
EXTDEF_CLASS
int bytecount[] = {0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,
                   1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
                   1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
                   2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
                   1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
                   2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
                   2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
                   3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
                   1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
                   2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
                   2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
                   3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
                   2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
                   3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
                   3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
                   4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8};

    /* array giving position (1..7) of high-order 1-bit in byte: */
EXTDEF_CLASS
int leftbit[] =   {8,7,6,6,5,5,5,5,4,4,4,4,4,4,4,4,
                   3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,
                   2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
                   2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
                   1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                   1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                   1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                   1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
#endif  /* EXTDEFS */

#if  ANSI_STDC
#ifndef  ANSIPROT
#define ANSIPROT 1
#endif
#endif

#ifdef  ANSIPROT
#define EXTPROC(func,args) extern func args;
#else
#define EXTPROC(func,args) extern func();
#endif

/* The following is for C++ programs that read nauty.h.  Compile nauty
   itself using C, not C++.  */

#ifdef __cplusplus
extern "C" {
#endif

EXTPROC(void alloc_error,(char*))
EXTPROC(int bestcell,(graph*,nvector*,nvector*,int,int,int,int))
EXTPROC(void breakout,(nvector*,nvector*,int,int,int,set*,int))
EXTPROC(boolean cheapautom,(nvector*,int,boolean,int))
EXTPROC(void doref,(graph*,nvector*,nvector*,int,int*,int*,permutation*,set*,
                  int*,UPROC(*)(),UPROC(*)(),int,int,int,boolean,int,int))
EXTPROC(boolean isautom,(graph*,permutation*,boolean,int,int))
EXTPROC(int itos,(int,char*))
EXTPROC(void fmperm,(permutation*,set*,set*,int,int))
EXTPROC(void fmptn,(nvector*,nvector*,int,set*,set*,int,int))
EXTPROC(void longprune,(set*,set*,set*,set*,int))
EXTPROC(void nauty,(graph*,nvector*,nvector*,set*,nvector*,optionblk*,
                  statsblk*,set*,int,int,int,graph*))
EXTPROC(void nauty_null,(void))
EXTPROC(void nautil_null,(void))
EXTPROC(int nextelement,(set*,int,int))
EXTPROC(int orbjoin,(nvector*,permutation*,int))
EXTPROC(void permset,(set*,set*,int,permutation*))
EXTPROC(void putstring,(FILE*,char*))
EXTPROC(UPROC refine,(graph*,nvector*,nvector*,int,int*,
                    permutation*,set*,int*,int,int))
EXTPROC(UPROC refine1,(graph*,nvector*,nvector*,int,int*,
                    permutation*,set*,int*,int,int))
EXTPROC(void shortprune,(set*,set*,int))
EXTPROC(UPROC targetcell,(graph*,nvector*,nvector*,int,int,set*,int*,
                        int*,int,int,int,int))
EXTPROC(int testcanlab,(graph*,graph*,nvector*,int*,int,int))
EXTPROC(void updatecan,(graph*,graph*,permutation*,int,int,int))
EXTPROC(void writeperm,(FILE*,permutation*,boolean,int,int))
EXTPROC(void nauty_check,(int,int,int,int))
EXTPROC(void nautil_check,(int,int,int,int))
EXTPROC(void nautil_freedyn,(void))

#ifdef __cplusplus
}
#endif

#define  NAUTYH_READ 1
#endif  /* NAUTYH_READ */
