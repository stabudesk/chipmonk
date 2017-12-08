CC=gcc
CFLAGS=-O3 -Wall
# don't have time to chase all those up
# CFLAGS=-O3
DBGCFLAGS=-g -Wall
TDBGCFLAGS=-g -Wall -DDBG # True debug flags!

EXES=chipmonk chipmonk_t chipmonk_d gfmatchup gfmatchup_t gfmatchup_d chipmonk2 chipmonk2_t chipmonk2_d

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

.PHONY: clean

clean:
	rm -f ${EXES}
