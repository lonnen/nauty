/*****************************************************************************
* This is the main header file for gtools.  nauty version 2.8.7
* Subject to the copyright notice in nauty.h.                                *
* gtools.h.  Generated from gtools-h.in by configure.
*****************************************************************************/

/* The parts between the ==== lines are modified by configure when
creating gtools.h out of gtools-h.in.  If configure is not being
used, it is necessary to check they are correct.
====================================================================*/

#ifndef  _GTOOLS_H_    /* only process this file once */
#define  _GTOOLS_H_

#define HAVE_ERRNO_H  1      /* <errno.h> exists */
#define HAVE_PERROR  1          /* perror() exists */
#define HAVE_PIPE  1          /* pipe() exists */
#define HAVE_WAIT  1          /* wait() exists */
#define HAVE_WAIT_H  1     /* <sys/wait.h> exists */
#define HAVE_POPEN  1          /* popen() and pclose() exist */
#define POPEN_DEC  1         /* popen() is declared */
#define FTELL_DEC  1         /* ftell() is declared */
#define FDOPEN_DEC  1        /* fdopen() is declared */
#define SORTPROG  "sort"         /* name of sort program */
#define SORT_NEWKEY 1  /* if -k is supported */
#define SORT_SIZE 1         /* if -S is supported */
#define HAVE_PID_T 1    /* pid_t is defined */
#define PUTENV_DEC 1   /* putenv() is declared */
#define SETENV_DEC 1   /* setenv() is declared */
#define HAVE_PUTENV 1   /* putenv() exists */
#define HAVE_SETENV 1   /* setenv() exists */
#define HAVE_FORK 1   /* fork() exists */
#define HAVE_SIGNAL_H  1      /* <signal.h> exists */
#define HAVE_FSEEKO 1  /* fseeko() and ftello() exist */
#define HAVE_SIGACTION 1  /* sigaction() exists */
#define HAVE_SIGPROCMASK 1  /* sigprocmask() exists */
#define ALLOW_INTERRUPT 1 /* no --disable-interrupt */
#define HAVE_GUNZIP 1 /* gunzip program is available */

/* ++++++ This file is automatically generated, don't edit it by hand! ++++++ */

/*==================================================================*/

#ifndef MAXN 
#define MAXN  0
#endif

#define SIZELEN(n) ((n)<=SMALLN?1:((n)<=SMALLISHN?4:8))
        /* length of size code in bytes */
#define G6BODYLEN(n) \
   (((size_t)(n)/12)*((size_t)(n)-1) + (((size_t)(n)%12)*((size_t)(n)-1)+11)/12)
#define G6LEN(n) (SIZELEN(n) + G6BODYLEN(n))
  /* exact graph6 string length excluding \n\0 
     This twisted expression works up to n=227023 in 32-bit arithmetic
     and for larger n if size_t has 64 bits.  */
#define D6BODYLEN(n) \
   ((n)*(size_t)((n)/6) + (((n)*(size_t)((n)%6)+5)/6))
#define D6LEN(n) (1 + SIZELEN(n) + D6BODYLEN(n))
  /* exact digraph6 string length excluding \n\0 
     This twisted expression works up to n=160529 in 32-bit arithmetic
     and for larger n if size_t has 64 bits.  */

#include "naututil.h"      /* which includes stdio.h */
#include "nausparse.h"

#if HAVE_ERRNO_H
#include <errno.h>
#else
extern int errno;
#endif

#if HAVE_WAIT_H && !defined(AVOID_SYS_WAIT_H)
#include <sys/wait.h>
#endif

#if HAVE_SIGNAL_H
#include <signal.h>
#endif

#if HAVE_PERROR
#define NAUTY_ABORT(msg) do {if (errno != 0) perror(msg); exit(1);} while(0)
#else
#define NAUTY_ABORT(msg) do {exit(1);} while(0)
#endif

#define gt_abort_1(fmt,x) \
  { char msg_[257]; snprintf(msg_,256,fmt,x); gt_abort(msg_); }
#define gt_abort_2(fmt,x,y) \
  { char msg_[257]; snprintf(msg_,256,fmt,x,y); gt_abort(msg_); }
#define gt_abort_3(fmt,x,y,z) \
  { char msg_[257]; snprintf(msg_,256,fmt,x,y,z); gt_abort(msg_); }

/* Here we set environment variables that determine the sorting order
   for the shortg program.  Older docs for sort say that it uses
   LC_COLLATE, but the POSIX description of locales says that the
   LC_ALL variable takes precedence over LC_COLLATE.  To be safe,
   we will define both.  Also, define this to be nothing if the
   variable KEEP_SORT_LOCALE is defined. */
#ifdef KEEP_SORT_LOCALE
#define SET_C_COLLATION
#else
#if PUTENV_DEC && HAVE_PUTENV
#define SET_C_COLLATION putenv("LC_ALL=C"); putenv("LC_COLLATE=C")
#elif SETENV_DEC && HAVE_SETENV
#define SET_C_COLLATION setenv("LC_ALL","C",1); setenv("LC_COLLATE","C",1)
#elif HAVE_PUTENV
int putenv(char*);
#define SET_C_COLLATION putenv("LC_ALL=C"); putenv("LC_COLLATE=C")
#elif HAVE_SETENV
int setenv(const char*,const char*,int);
#define SET_C_COLLATION setenv("LC_ALL","C",1); setenv("LC_COLLATE","C",1)
#else
#define SET_C_COLLATION
#endif
#endif

#if HAVE_FLOCKFILE
#define FLOCKFILE(f) flockfile(f)
#define FUNLOCKFILE(f) funlockfile(f)
#else
#define FLOCKFILE(f)
#define FUNLOCKFILE(f)
#endif

#if HAS_STDIO_UNLOCK && !defined(NAUTY_IN_MAGMA) && !defined(IS_JAVA)
#define GETC(f) getc_unlocked(f)
#undef PUTC
#define PUTC(c,f) putc_unlocked(c,f)
#else
#define GETC(f) getc(f)
#undef PUTC
#define PUTC(c,f) putc(c,f)
#endif

#define BIAS6 63
#define MAXBYTE 126
#define SMALLN 62
#define SMALLISHN 258047
#define TOPBIT6 32
#define C6MASK 63

#define GRAPH6_HEADER ">>graph6<<"
#define SPARSE6_HEADER ">>sparse6<<"
#define DIGRAPH6_HEADER ">>digraph6<<"
#define PLANARCODE_HEADER ">>planar_code<<"
#define PLANARCODELE_HEADER ">>planar_code le<<"
#define PLANARCODEBE_HEADER ">>planar_code be<<"
#define EDGECODE_HEADER ">>edge_code<<"

#define GRAPH6         1
#define SPARSE6        2
#define PLANARCODE     4
#define PLANARCODELE   8
#define PLANARCODEBE  16
#define EDGECODE      32
#define INCSPARSE6    64
#define PLANARCODEANY (PLANARCODE|PLANARCODELE|PLANARCODEBE)
#define DIGRAPH6     128
#define UNKNOWN_TYPE 256
#define HAS_HEADER   512

#define NODIGRAPHSYET(code) if (((code)&DIGRAPH6)) \
  gt_abort(">E Sorry, this program doesn't support digraphs yet.\n")

#define ARG_OK 0
#define ARG_MISSING 1
#define ARG_TOOBIG 2
#define ARG_ILLEGAL 3

#define MAXARG 2140000000L
#define NOLIMIT (MAXARG+31L)
#if SIZEOF_LONG == 8
#define MAXLONGARG 9220000000000000000L
#else
#define MAXLONGARG MAXARG
#endif

#define SWBOOLEAN(c,boool) if (sw==c) boool=TRUE;
#define SWCOUNT(c,count) if (sw==c) ++count;
#define SWINT(c,boool,val,id) if (sw==c) \
        {boool=TRUE;arg_int(&arg,&val,id);}
#define SWLONG(c,boool,val,id) if (sw==c) \
        {boool=TRUE;arg_long(&arg,&val,id);}
#define SWULL(c,boool,val,id) if (sw==c) \
        {boool=TRUE;arg_ull(&arg,&val,id);}
#define SWRANGE(c,sep,boool,val1,val2,id) if (sw==c) \
        {boool=TRUE;arg_range(&arg,sep,&val1,&val2,id);}
#define SWREAL(c,boool,val,id) if (sw==c) \
        {boool=TRUE;arg_double(&arg,&val,id);}
#define SWREALRANGE(c,sep,boool,val1,val2,id) if (sw==c) \
        {boool=TRUE;arg_doublerange(&arg,sep,&val1,&val2,id);}
#define SWSEQUENCE(c,sep,boool,val,maxvals,numvals,id) if (sw==c) \
        {boool=TRUE;arg_sequence(&arg,sep,val,maxvals,&numvals,id);}
#define SWSEQUENCEMIN(c,sep,boool,val,minvals,maxvals,numvals,id) if (sw==c) \
        {boool=TRUE;arg_sequence_min(&arg,sep,val,minvals,maxvals,&numvals,id);}

#ifdef HELPUSECMD
#ifdef HELPTEXT2 
#define PUTHELPTEXT printf("\nUsage: %s %s\n\n%s",argv[0],USAGE,HELPTEXT1);\
                    printf("%s",HELPTEXT2);
#else
#define PUTHELPTEXT printf("\nUsage: %s %s\n\n%s",argv[0],USAGE,HELPTEXT)
#endif
#else
#ifdef HELPTEXT2 
#define PUTHELPTEXT printf("\nUsage: %s\n\n%s",USAGE,HELPTEXT1);\
                    printf("%s",HELPTEXT2);
#else
#define PUTHELPTEXT printf("\nUsage: %s\n\n%s",USAGE,HELPTEXT)
#endif
#endif

#define HELP if (argc > 1 && (strcmp(argv[1],"-help")==0 \
                           || strcmp(argv[1],"/?")==0 \
                           || strcmp(argv[1],"--help")==0)) \
       { PUTHELPTEXT; return 0; }

#define PUTVERSION if (argc > 1 && (strcmp(argv[1],"-version")==0 \
                           || strcmp(argv[1],"--version")==0)) \
       { printf("Nauty&Traces version %.4f (%d bits)\n",\
       NAUTYVERSIONID/10000.0,WORDSIZE); return 0; }

#define GETHELP \
fprintf(stderr,"   Use %s -help to see more detailed instructions.\n",argv[0])

#define alloc_error gt_abort

#define CATMSG0(fmt) sprintf(msg+strlen(msg),fmt)
#define CATMSG1(fmt,x1) sprintf(msg+strlen(msg),fmt,x1)
#define CATMSG2(fmt,x1,x2) sprintf(msg+strlen(msg),fmt,x1,x2)
#define CATMSG3(fmt,x1,x2,x3) sprintf(msg+strlen(msg),fmt,x1,x2,x3)
#define CATMSG4(fmt,x1,x2,x3,x4) sprintf(msg+strlen(msg),fmt,x1,x2,x3,x4)
#define CATMSG5(fmt,x1,x2,x3,x4,x5) sprintf(msg+strlen(msg),fmt,x1,x2,x3,x4,x5)
#define CATMSG6(fmt,x1,x2,x3,x4,x5,x6) \
                sprintf(msg+strlen(msg),fmt,x1,x2,x3,x4,x5,x6)
#define CATMSG7(fmt,x1,x2,x3,x4,x5,x6,x7) \
                sprintf(msg+strlen(msg),fmt,x1,x2,x3,x4,x5,x6,x7)
#define CATMSG8(fmt,x1,x2,x3,x4,x5,x6,x7,x8) \
                sprintf(msg+strlen(msg),fmt,x1,x2,x3,x4,x5,x6,x7,x8)
#define CATMSG9(fmt,x1,x2,x3,x4,x5,x6,x7,x8,x9) \
                sprintf(msg+strlen(msg),fmt,x1,x2,x3,x4,x5,x6,x7,x8,x9)

/************************************************************************/

/* ++++++ This file is automatically generated, don't edit it by hand! ++++++ */

#ifdef __cplusplus
extern "C" {
#endif
 
extern void gtools_check(int,int,int,int);
extern FILE *opengraphfile(char*,int*,boolean,long);
extern void writeline(FILE*,char*);
extern char *gtools_getline(FILE*);     /* formerly getline() */
extern int graphsize(char*);
extern void encodegraphsize(int,char**);
extern void stringcounts(char*,int*,size_t*);
extern void stringtograph(char*,graph*,int);
extern void stringtograph_inc(char*,graph*,int,graph*,int);
extern size_t edgecount(char*);
extern int checkgline(char*);
extern graph *readgg(FILE*,graph*,int,int*,int*,boolean*);
extern graph *readg(FILE*,graph*,int,int*,int*);
extern graph *readgg_inc(FILE*,graph*,int,int*,int*,graph*,int,int,boolean*);
extern graph *readg_inc(FILE*,graph*,int,int*,int*,graph*,int,int);
extern char *ntog6(graph*,int,int);
extern char *ntos6(graph*,int,int);
extern char *ntod6(graph*,int,int);
extern char *sgtos6(sparsegraph*);
extern char *sgtog6(sparsegraph*);
extern char *sgtod6(sparsegraph*);
extern void writeg6(FILE*,graph*,int,int);
extern void writed6(FILE*,graph*,int,int);
extern void writes6(FILE*,graph*,int,int);
extern void writeg6_sg(FILE*,sparsegraph*);
extern void writes6_sg(FILE*,sparsegraph*);
extern void writed6_sg(FILE*,sparsegraph*);
extern char *ntois6(graph*,graph*,int,int);
extern void writeis6(FILE*,graph*,graph*,int,int);
extern void writepc_sg(FILE*,sparsegraph*);
extern void stringtosparsegraph(char*,sparsegraph*,int*);
extern sparsegraph *read_sg(FILE*,sparsegraph*);
extern sparsegraph *read_sg_loops(FILE*,sparsegraph*,int*);
extern sparsegraph *read_sgg_loops(FILE*,sparsegraph*,int*,boolean*);
extern sparsegraph *readpc_sg(FILE*,sparsegraph*);
extern sparsegraph *readpcle_sg(FILE*,sparsegraph*);
extern char *getecline(FILE*);
extern void writelast(FILE*);
extern int longvalue(char**,long*);
extern int ullvalue(char**,unsigned long long*);
extern void arg_int(char**,int*,char*);
extern void arg_long(char**,long*,char*);
extern void arg_ull(char**,unsigned long long*,char*);
extern void arg_range(char**,char*,long*,long*,char*);
extern int doublevalue(char**,double*);
extern void arg_double(char**,double*,char*);
extern void arg_doublerange(char**,char*,double*,double*,char*);
extern void arg_sequence(char**,char*,long*,int,int*,char*);
extern void arg_sequence_min(char**,char*,long*,int,int,int*,char*);

extern void writerange(FILE*,int,long,long);
extern void gt_abort(const char*);
extern char *stringcopy(char*);
extern boolean strhaschar(char*,int);

extern void fcanonise(graph*,int,int,graph*,char*,boolean);
extern void fcanonise_inv
             (graph*,int,int,graph*,char*,void(*)(graph*,int*,int*,int,
               int,int,int*,int,boolean,int,int),int,int,int,boolean);
extern void fcanonise_inv_sg
           (sparsegraph*,int,int,sparsegraph*,char*,void(*)(graph*,int*,int*,
             int,int,int,int*,int,boolean,int,int),int,int,int,boolean);
extern void setlabptn(int*,int*,int*,int);
extern int breakcellwt(int*,int*,int*,int,int);

extern void fgroup(graph*,int,int,char*,int*,int*);
extern void fgroup_inv
             (graph*,int,int,char*,int*,int*,void(*)(graph*,int*,int*,int,
                int,int,int*,int,boolean,int,int),int,int,int);
extern int istransitive(graph*,int,int,graph*);
extern void tg_canonise(graph*,graph*,int,int);

extern TLS_ATTR int readg_code;
extern TLS_ATTR char *readg_line;
extern TLS_ATTR size_t ogf_linelen;
extern TLS_ATTR boolean is_pipe;

#ifdef __cplusplus
}
#endif

#ifdef CPUDEFS
CPUDEFS
#endif

/* ++++++ This file is automatically generated, don't edit it by hand! ++++++ */

#endif /* _GTOOLS_H_  */
