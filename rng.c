/* rng.c
   This file contains the code for a high-quality random number
   generator written by Don Knuth.  The auxilliary routine
   ran_arr_cycle() has been modified, and ran_init() is new.

   To use it:

      0.  #include "rng.h" (or "naututil.h" if you are using nauty)

      1.  Call ran_init(seed), where seed is any long integer.
          This step is optional, but if you don't use it you
          will always get the same sequence of random numbers.

      2.  For each random number, use the NEXTRAN macro.  It will
          give a random value in the range 0..2^30-1.  Alternatively,
          KRAN(k) will have a random value in the range 0..k-1.
          KRAN(k) actually gives you NEXTRAN mod k, so it is not
          totally uniform if k is very large.  In that case, you
          can use the slightly slower GETKRAN(k,var) to set the
          variable var to a better random number from 0..k-1.

    Brendan McKay, December 1997.

    The parameter QUALITY should be increased if the best quality is
    essential even at the expence of increased time.  Don suggests
    the value 1009.
*/

/*    This program by D E Knuth is in the public domain and freely copyable
 *    AS LONG AS YOU MAKE ABSOLUTELY NO CHANGES!
 *
 * From BDM:  Whoops, too late.  Somebadlyneededwhitespacewas
 * addedtocorrectfortheKnuthagrophobia.  Sorry Don!
 *
 *    It is explained in Seminumerical Algorithms, 3rd edition, Section 3.6
 *    (or in the errata to the 2nd edition --- see
 *        http://www-cs-faculty.stanford.edu/~knuth/taocp.html
 *    in the changes to pages 171 and following of Volume 2).              */

/*    If you find any bugs, please report them immediately to
 *                 taocp@cs.stanford.edu
 *    (and you will be rewarded if the bug is genuine). Thanks!            */

/************ see the book for explanations and caveats! *******************/

/* the old C calling conventions are used here, for reasons of portability */

#define KK 100                      /* the long lag */
#define LL  37                      /* the short lag */
#define MM (1L<<30)                 /* the modulus */
#define mod_diff(x,y) (((x)-(y))&(MM-1)) /* subtraction mod MM */

static long ran_x[KK];              /* the generator state */

/* void ran_array(long aa[],int n) */
void ran_array(aa,n)    /* put n new random numbers in aa */
	long *aa;   /* destination */
	int n;      /* array length (must be at least KK) */
{
	register int i,j;

	for (j = 0; j < KK; j++) aa[j] = ran_x[j];
	for (;j < n; j++) aa[j] = mod_diff(aa[j-KK],aa[j-LL]);
	for (i = 0; i < LL; i++,j++) ran_x[i] = mod_diff(aa[j-KK],aa[j-LL]);
	for (;i < KK; i++,j++) ran_x[i] = mod_diff(aa[j-KK],ran_x[i-LL]);
}

#define TT  70   /* guaranteed separation between streams */
#define is_odd(x)  ((x)&1)          /* units bit of x */
#define evenize(x) ((x)&(MM-2))     /* make x even */

/* void ran_start(long seed) */
void ran_start(seed)    /* do this before using ran_array */
long seed;              /* selector for different streams */
{
	register int t,j;
	long x[KK+KK-1];              /* the preparation buffer */
	register long ss = evenize(seed+2);

	for (j = 0;j < KK;j++) {
	    x[j] = ss;                      /* bootstrap the buffer */
	    ss <<= 1;
	    if (ss >= MM) ss -= MM-2; /* cyclic shift 29 bits */
	}
	for (;j < KK+KK-1; j++) x[j] = 0;
	x[1]++;              /* make x[1] (and only x[1]) odd */
	ss = seed&(MM-1);
	t = TT-1;

	while (t) {
	    for (j = KK-1; j > 0; j--) x[j+j] = x[j];  /* "square" */
	    for (j = KK+KK-2; j > KK-LL; j -= 2) x[KK+KK-1-j] = evenize(x[j]);
	    for (j = KK+KK-2; j >= KK; j--) if (is_odd(x[j])) {
	        x[j-(KK-LL)] = mod_diff(x[j-(KK-LL)],x[j]);
	        x[j-KK] = mod_diff(x[j-KK],x[j]);
	    }
	    if (is_odd(ss)) {              /* "multiply by z" */
	        for (j = KK; j > 0; j--)  x[j] = x[j-1];
	        x[0] = x[KK];            /* shift the buffer cyclically */
	        if (is_odd(x[KK])) x[LL] = mod_diff(x[LL],x[KK]);
	    }
	    if (ss) ss >>= 1; else t--;
	}
	for (j = 0; j < LL; j++) ran_x[j+KK-LL] = x[j];
	for (;j < KK; j++) ran_x[j-LL] = x[j];
}

/* the following routines are from exercise 3.6--15 */
/* after calling ran_start, get new randoms by, e.g., "x=ran_arr_next()" */

#define QUALITY 139 /* at least KK */
#if QUALITY < KK
#error  "QUALITY >= KK is essential!"
#endif
static long ran_arr_buf[QUALITY+1];
static long ran_arr_sentinel = -1;
long *ran_arr_ptr = &ran_arr_sentinel; /* the next random number, or -1 */

void
ran_init(seed)    /* Added by BDM: use instead of ran_start. */
long seed;
{
	int i;

	seed = (unsigned long)seed % (MM-2);

 	ran_start(seed);
 	for (i = 0; i < 10; ++i)
	     ran_array(ran_arr_buf,QUALITY);
	ran_arr_ptr = ran_arr_buf;
	ran_arr_buf[KK] = -1;
}

long ran_arr_cycle()    /* Modified by BDM to automatically initialise */
                        /* if no explicit initialisation has been done */
{
	if (ran_arr_ptr == &ran_arr_sentinel)
	    ran_init(12345L);
	else
	    ran_array(ran_arr_buf,QUALITY);

	ran_arr_buf[KK] = -1;
	ran_arr_ptr = ran_arr_buf+1;
	return ran_arr_buf[0];
}
