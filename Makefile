CC=gcc
CFLAGS=-O3 -Wall
# don't have time to chase all those up
# CFLAGS=-O3
DBGCFLAGS=-g -Wall
TDBGCFLAGS=-g -Wall -DDBG # True debug flags!

EXES=chipmonk chipmonk_t chipmonk_d gfmatchup gfmatchup_t gfmatchup_d chipmonk2 chipmonk2_t chipmonk2_d hblop hblop_d gffsimp astarp gffsidn

# production binary
chipmonk: chipmonk.c
	${CC} ${CFLAGS} -o $@ $^
# testing mode binary
chipmonk_t: chipmonk.c
	${CC} ${DBGCFLAGS} -o $@ $^
# testing mode binary
chipmonk_d: chipmonk.c
	${CC} ${TDBGCFLAGS} -o $@ $^

# chipmonk was getting crowded. In any case, a good deal of it was about matching gff files against
# each other, in the hope of establishing clean set operations. So, really, a dedicated program should exist for this
# and hence, gfmatchup.c
# production binary
gfmatchup: gfmatchup.c
	${CC} ${CFLAGS} -o $@ $^
# testing mode binary
gfmatchup_t: gfmatchup.c
	${CC} ${DBGCFLAGS} -o $@ $^
# testing mode binary
gfmatchup_d: gfmatchup.c
	${CC} ${TDBGCFLAGS} -o $@ $^

# OK yeast has a small ggf, but when dealing with big gffs, you want a simple way to grep it properly
# sometimes you just want to be able to retrieve the sequence for a certain feature. gffsimp is an attempt to that.
gffsimp: gffsimp.c
	${CC} ${DBGCFLAGS} -o $@ $^
# the -u option is for a single column of names ... gffsidn changes this to dual names, helps cross referencing. So we have two columns in the name file. Actually three. sorry won;t change name becaue of that.
# # this just affects the words_t struct
gffsidn: gffsidn.c
	${CC} ${DBGCFLAGS} -o $@ $^

# derived from gffsimp, reading in array star mirna reports
astarp: astarp.c
	${CC} ${DBGCFLAGS} -o $@ $^

# I then abandoned gfmatchup, ane went back to chipmonk,
# chipmonk2 is an attempt to integrate
# # primarily this is just about getting the blastoutput file format in .. -b
chipmonk2: chipmonk2.c
	${CC} ${CFLAGS} -o $@ $^
# testing mode binary
chipmonk2_t: chipmonk2.c
	${CC} ${DBGCFLAGS} -o $@ $^
# testing mode binary
chipmonk2_d: chipmonk2.c
	${CC} ${TDBGCFLAGS} -o $@ $^

hblop: hblop.c
	${CC} ${CFLAGS} -o $@ $^
# testing mode binary
hblop_t: hblop.c
	${CC} ${DBGCFLAGS} -o $@ $^
# testing mode binary
hblop_d: hblop.c
	${CC} ${TDBGCFLAGS} -o $@ $^

.PHONY: clean

clean:
	rm -f ${EXES}
