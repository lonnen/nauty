# makefile for nauty 2.0

CC=gcc
CFLAGS=-O4 
LDFLAGS=

SMALL=-DMAXN=WORDSIZE
DYNAMIC=-DMAXN=0
BIG=-DMAXN=0 -DBIGNAUTY
L1=-DMAXN=WORDSIZE -DWORDSIZE=64
L=-DMAXN=0 -DWORDSIZE=64
FOURK=-DMAXN=4096

nauty: dreadnaut dreadnaut1 dreadnautL1 dreadnautB
	echo ""

rng.o: rng.c
	${CC} -c ${CFLAGS} rng.c

dreadnaut: dreadnaut.c naututil.o nauty.o nautil.o nautinv.o rng.o
	${CC} -o dreadnaut ${CFLAGS} ${DYNAMIC} dreadnaut.c \
	    naututil.o nauty.o nautil.o nautinv.o rng.o ${LDFLAGS}
naututil.o: nauty.h naututil.h naututil.c
	${CC} -c ${CFLAGS} ${DYNAMIC} naututil.c
nautil.o: nauty.h nautil.c
	${CC} -c ${CFLAGS} ${DYNAMIC} nautil.c
nauty.o: nauty.h nauty.c
	${CC} -c ${CFLAGS} ${DYNAMIC} nauty.c
nautinv.o: nauty.h naututil.h nautinv.c
	${CC} -c ${CFLAGS} ${DYNAMIC} nautinv.c
nautaux.o: nautaux.h nauty.h naututil.h nautaux.c
	${CC} -c ${CFLAGS} ${DYNAMIC} nautaux.c

dreadnaut1: dreadnaut.c naututil1.o nauty1.o nautil1.o nautinv1.o rng.o
	${CC} -o dreadnaut1 ${CFLAGS} ${SMALL} dreadnaut.c \
	    naututil1.o nauty1.o nautil1.o nautinv1.o rng.o ${LDFLAGS}
naututil1.o: nauty.h naututil.h naututil.c
	${CC} -c ${CFLAGS} ${SMALL} -o naututil1.o naututil.c
nautil1.o: nauty.h nautil.c
	${CC} -c ${CFLAGS} ${SMALL} -o nautil1.o nautil.c
nauty1.o: nauty.h nauty.c
	${CC} -c ${CFLAGS} ${SMALL} -o nauty1.o nauty.c
nautinv1.o: nauty.h naututil.h nautinv.c
	${CC} -c ${CFLAGS} ${SMALL} -o nautinv1.o nautinv.c
nautaux1.o: nautaux.h nauty.h naututil.h nautaux.c
	${CC} -c ${CFLAGS} ${SMALL} -o nautaux1.o nautaux.c

dreadnautB: dreadnaut.c naututilB.o nautyB.o nautilB.o nautinvB.o rng.o
	${CC} -o dreadnautB ${CFLAGS} ${BIG} dreadnaut.c \
	    naututilB.o nautyB.o nautilB.o nautinvB.o rng.o ${LDFLAGS}
naututilB.o: nauty.h naututil.h naututil.c
	${CC} -c ${CFLAGS} ${BIG} -o naututilB.o naututil.c
nautilB.o: nauty.h nautil.c
	${CC} -c ${CFLAGS} ${BIG} -o nautilB.o nautil.c
nautyB.o: nauty.h nauty.c
	${CC} -c ${CFLAGS} ${BIG} -o nautyB.o nauty.c
nautinvB.o: nauty.h naututil.h nautinv.c
	${CC} -c ${CFLAGS} ${BIG} -o nautinvB.o nautinv.c
nautauxB.o: nautaux.h nauty.h naututil.h nautaux.c
	${CC} -c ${CFLAGS} ${BIG} -o nautauxB.o nautaux.c

dreadnaut4K: dreadnaut.c naututil4K.o nauty4K.o nautil4K.o nautinv4K.o rng.o
	${CC} -o dreadnaut4K ${CFLAGS} ${FOURK} dreadnaut.c \
	    naututil4K.o nauty4K.o nautil4K.o nautinv4K.o rng.o ${LDFLAGS}
naututil4K.o: nauty.h naututil.h naututil.c
	${CC} -c ${CFLAGS} ${FOURK} -o naututil4K.o naututil.c
nautil4K.o: nauty.h nautil.c
	${CC} -c ${CFLAGS} ${FOURK} -o nautil4K.o nautil.c
nauty4K.o: nauty.h nauty.c
	${CC} -c ${CFLAGS} ${FOURK} -o nauty4K.o nauty.c
nautinv4K.o: nauty.h naututil.h nautinv.c
	${CC} -c ${CFLAGS} ${FOURK} -o nautinv4K.o nautinv.c
nautaux4K.o: nautaux.h nauty.h naututil.h nautaux.c
	${CC} -c ${CFLAGS} ${FOURK} -o nautaux4K.o nautaux.c

nautyex: nauty.h nautyex.c nauty.o nautil.o
	${CC} -o nautyex ${CFLAGS} nautyex.c nauty.o nautil.o ${LDFLAGS}

dreadnautL1: dreadnaut.c naututilL1.o nautyL1.o nautilL1.o nautinvL1.o rng.o
	${CC} -o dreadnautL1 ${CFLAGS} ${L1} dreadnaut.c \
	    naututilL1.o nautyL1.o nautilL1.o nautinvL1.o rng.o ${LDFLAGS}
naututilL1.o: nauty.h naututil.h naututil.c
	${CC} -c ${CFLAGS} ${L1} -o naututilL1.o naututil.c
nautilL1.o: nauty.h nautil.c
	${CC} -c ${CFLAGS} ${L1} -o nautilL1.o nautil.c
nautyL1.o: nauty.h nauty.c
	${CC} -c ${CFLAGS} ${L1} -o nautyL1.o nauty.c
nautinvL1.o: nauty.h naututil.h nautinv.c
	${CC} -c ${CFLAGS} ${L1} -o nautinvL1.o nautinv.c
nautauxL1.o: nautaux.h nauty.h naututil.h nautaux.c
	${CC} -c ${CFLAGS} ${L1} -o nautauxL1.o nautaux.c

dreadnautL: dreadnaut.c naututilL.o nautyL.o nautilL.o nautinvL.o rng.o
	${CC} -o dreadnautL ${CFLAGS} ${L} dreadnaut.c \
	    naututilL.o nautyL.o nautilL.o nautinvL.o rng.o ${LDFLAGS}
naututilL.o: nauty.h naututil.h naututil.c
	${CC} -c ${CFLAGS} ${L} -o naututilL.o naututil.c
nautilL.o: nauty.h nautil.c
	${CC} -c ${CFLAGS} ${L} -o nautilL.o nautil.c
nautyL.o: nauty.h nauty.c
	${CC} -c ${CFLAGS} ${L} -o nautyL.o nauty.c
nautinvL.o: nauty.h naututil.h nautinv.c
	${CC} -c ${CFLAGS} ${L} -o nautinvL.o nautinv.c
nautauxL.o: nautaux.h nauty.h naututil.h nautaux.c
	${CC} -c ${CFLAGS} ${L} -o nautauxL.o nautaux.c

nauty-clean:
	rm -f nauty*.o nautil*.o nautaux*.o nautinv*.o naututil*.o

nauty-dist: 
	mkdir nauty20
	cp dreadnaut.c nautaux.h nauty.h nautil.c nautyex.c \
	   nautinv.c naututil.c nautyex.c read.me naututil.h nautaux.c \
	   nauty.c rng.c rng.h oldmanual.ps read.me  nauty20
	sed -e 's///' <makefile >nauty20/makefile
	tar cvf nauty20.tar nauty20
	gzip nauty20.tar
	rm -r nauty20

# Now follows make scripts for gtools, which is distributed with nauty
# but not in the same package.

gtools: copyg listg labelg dretog amtog geng complg shortg readg NRswitchg \
     deledgeg countg pickg
	echo ""

gtools.h : nauty.h naututil.h
	touch gtools.h

gtools.o : gtools.h gtools.c
	${CC} -c ${CFLAGS} gtools.c

gtnauty.o : gtools.h gtnauty.c
	${CC} -c ${CFLAGS} gtnauty.c

gutil1.o : gtools.h gutils.h gutil1.c
	${CC} -c ${CFLAGS} gutil1.c

copyg : gtools.h copyg.c gtools.o
	${CC} -o copyg ${CFLAGS} copyg.c gtools.o ${LDFLAGS}

listg : gtools.h listg.c gtools.o nautil.o
	${CC} -o listg ${CFLAGS} listg.c gtools.o nautil.o ${LDFLAGS}

labelg : gtools.h labelg.c gtools.o gtnauty.o nauty.o nautil.o
	${CC} -o labelg ${CFLAGS} labelg.c \
		gtools.o gtnauty.o nauty.o nautil.o ${LDFLAGS}

shortg : gtools.h shortg.c gtools.o gtnauty.o nauty.o nautil.o
	${CC} -o shortg ${CFLAGS} shortg.c \
		gtools.o gtnauty.o nauty.o nautil.o ${LDFLAGS}

dretog : gtools.h dretog.c gtools.o naututil.o nautil.o rng.o
	${CC} -o dretog ${CFLAGS} dretog.c \
		gtools.o naututil.o nautil.o rng.o ${LDFLAGS}

amtog : gtools.h amtog.c gtools.o
	${CC} -o amtog ${CFLAGS} amtog.c gtools.o ${LDFLAGS}

geng : gtools.h geng.c gtools.o nauty1.o nautil1.o
	${CC} -o geng ${CFLAGS} geng.c gtools.o nauty1.o nautil1.o ${LDFLAGS}

complg : gtools.h complg.c gtools.o gtnauty.o nauty.o nautil.o
	${CC} -o complg ${CFLAGS} complg.c \
		gtools.o gtnauty.o nauty.o nautil.o ${LDFLAGS}

NRswitchg : gtools.h NRswitchg.c gtools.o gtnauty.o nauty.o nautil.o
	${CC} -o NRswitchg ${CFLAGS} NRswitchg.c gtools.o gtnauty.o \
		 nauty.o nautil.o ${LDFLAGS}

deledgeg : gtools.h deledgeg.c gtools.o gtnauty.o nauty.o nautil.o
	${CC} -o deledgeg ${CFLAGS} deledgeg.c gtools.o gtnauty.o \
		 nauty.o nautil.o ${LDFLAGS}

pickg : gtools.h testg.c splay.c gtools.o gtnauty.o nauty.o nautil.o gutil1.o
	${CC} -o pickg ${CFLAGS} testg.c gtools.o gtnauty.o gutil1.o \
		 nauty.o nautil.o ${LDFLAGS}

countg : gtools.h testg.c splay.c gtools.o gtnauty.o nauty.o nautil.o gutil1.o
	${CC} -o countg ${CFLAGS} testg.c gtools.o gtnauty.o gutil1.o \
		 nauty.o nautil.o ${LDFLAGS}

readg: readg.c
	${CC} -o readg ${CFLAGS} readg.c ${LDFLAGS}

sumlines: sumlines.c
	${CC} -o sumlines ${CFLAGS} sumlines.c

gtools-dist:
	mkdir gtools10
	cp gtools.h gtools.c geng.c gtnauty.c complg.c labelg.c shortg.c \
	   listg.c readg.c dretog.c amtog.c copyg.c formats.txt \
           NRswitchg.c deledgeg.c testg.c splay.c gutil1.c sumlines.c \
	   gtools10
	sed -e 's///' <makefile >gtools10/makefile
	tar cvf gtools10.tar gtools10
	gzip gtools10.tar
	rm -r gtools10

all:  nauty gtools
	echo ""
