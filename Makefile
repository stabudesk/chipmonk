CC=gcc
CFLAGS=-O3
DBGCFLAGS=-g -Wall
TDBGCFLAGS=-g -Wall -DDBG # True debug flags!

EXES=chipmonk

# production binary
chipmonk: chipmonk.c
	${CC} ${CFLAGS} -o $@ $^

# testing mode binary
chipmonk_t: chipmonk.c
	${CC} ${DBGCFLAGS} -o $@ $^

# testing mode binary
chipmonk_d: chipmonk.c
	${CC} ${TDBGCFLAGS} -o $@ $^

.PHONY: clean

clean:
	rm -f ${EXES}
