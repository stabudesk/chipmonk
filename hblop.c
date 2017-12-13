/* hblop.c derived from chipmonk.c: just matching blastoutput and a program to "tackle" BED (genomic features file) files.
   Copyright (C) 2017 Ramon Fallon, University of St Andrews.

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<unistd.h> // required for optopt, opterr and optarg.
#include <locale.h>

#ifdef DBG
#define GBUF 4
#define WBUF 4
#else
#define GBUF 32
#define WBUF 32
#endif

#define SSZ 2 /* CG count, first, AT count second, third are the anomalous characters */
#define NUMBUCKETS 20
#define GF22IDCNUM 9 /* fr the ID field 9th col: GF22 ID Col Number */
#define CUTOFFPCT 50.0

// the following is the way we cut out later columns that we have chosen to ignore
#define MXCOL2VIEW 4

#define CONDREALLOC(x, b, c, a, t); \
	if((x)>=((b)-1)) { \
		(b) += (c); \
		(a)=realloc((a), (b)*sizeof(t)); \
		memset(((a)+(b)-(c)), 0, (c)*sizeof(t)); \
	}

#define CONDREALLOC2(x, b, c, a1, a2, t); \
	if((x)==((b)-1)) { \
		(b) += (c); \
		(a1)=realloc((a1), (b)*sizeof(t)); \
		memset((a1)+(b)-(c), '\0', (c)*sizeof(t)); \
		(a2)=realloc((a2), (b)*sizeof(t)); \
		memset((a2)+(b)-(c), '\0', (c)*sizeof(t)); \
	}

typedef unsigned char boole;
typedef enum { /* feature category could really be called ycat because it's for the -y (gf22) file */
	ABS=1, /* an absent feature */
	TMR=2, /* a telomere */
	FTH=3, /* fthing */
	YTH=4, /* ything */
	NCA=5, /* no categry assigned */
	ALL=6 /* the fc member won't get this */
} fcat;
#define FCQUAN 6
char *fcnames[FCQUAN]={"ABS", "TMR", "FTH", "YTH", "NCA", "ALL"};
fcat getfc(char *cnam) /* this function has special permission to be up here */
{
	int i;
	fcat ret;
	for(i=0;i<FCQUAN;++i)
		if(!strcmp(fcnames[i], cnam)) {
			ret=i+1;
			break;
		}
	return ret;
}

typedef struct /* opt_t, a struct for the options */
{
	boole dflg; /* details / information only */
	char *bstr; /* blast output */
	char *ystr; /* the gf22_t type */
} opt_t;

typedef struct /* wseq_t */
{
	size_t *wln;
	size_t wsbuf;
	size_t quan;
	size_t lbuf; /* a buffer for the number of lines */
	size_t numl; /* number of lines, i.e. rows */
	size_t *wpla; /* words per line array: the number of words on each line */
} wseq_t;

typedef struct /* gf22_t: the y option, based on rmf_t */
{
	char *n;
	size_t nsz; /* size of the name r ID field */
	long c[2]; /* coords: 1) start 2) cols 4 and 5 */
	char *t; /* the typestring ... 3rd column */
	char *i; /* the ID string ... last column */
	char sd; /* strand + or - */
	size_t isz; /* size of iD string */
	size_t tsz; /* size of type stringstring */
	fcat fc;
	unsigned char blofound; // relates to specific -b -y operation ... did this get a hit? if so skip othe hits; zero if not.
} gf22_t;

typedef struct /* blop_t: the b option, bast output format 7 */
{
	char *n;
	size_t nsz; /* size of the name r ID field */
	char *tc;
	size_t tcsz; /* target string ..yep */
	char *qfs; /* query Feature substring */
	long q[2]; /* query feature beginning and end */
	char *fs; /* target feature substring */
	size_t fssz; /* size of fs */
	size_t qfssz; /* size of fs */
	int al; /* alignment length */
	float eval;
	float pcti;
	long c[2]; /* coords on the target : 1) start 2) cols 8 and 9 */
	unsigned qfnoc; /* query feature name number of occurrences */
} blop_t;

struct strchainode /* gf22snod struct */
{
    gf22_t *gf22; /* ptr to a single element, not an array, TODO void* it. */
    struct strchainode *n;
};
typedef struct strchainode gf22snod; /* yes, leave this alone, it's the way a struct can have a ptr ot its own type! */

struct strchainodeblo /* blosnod struct */
{
    blop_t *blop; /* ptr to a single element, not an array, TODO void* it. */
    struct strchainodeblo *n;
	unsigned char ocn; /* occurrence number, For we dont want the first occurence ... beware though, if there is only one occurrence it'll be the firt one */
};
typedef struct strchainodeblo blosnod; /* yes, leave this alone, it's the way a struct can have a ptr ot its own type! */

int catchopts(opt_t *opts, int oargc, char **oargv)
{
	int c;
	opterr = 0;

	/* Note: see the data structure for docs */
	while ((c = getopt (oargc, oargv, "dsni:f:u:p:g:r:q:y:a:h:l:z:b:x:")) != -1)
		switch (c) {
			case 'd':
				opts->dflg = 1;
				break;
			case 'b':
				opts->bstr = optarg;
				break;
			case 'y':
				opts->ystr = optarg;
				break;
			case '?':
				fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
				return 1;
			default:
				fprintf (stderr, "Wrong arguments. Please launch without arguments to see help file.\n");
				exit(EXIT_FAILURE);
		}
	return 0;
}

unsigned hashit(char *str, unsigned tsz) /* Dan Bernstein's one */
{
    unsigned long hash = 5381;
    int c;

    char *tstr=str;
    while ((c = *tstr++))
        hash = ((hash << 5) + hash) + c; /*  hash * 33 + c */

    return hash % tsz;
}

gf22snod **gf22tochainharr(gf22_t *gf22, unsigned numsq, unsigned tsz)
{
    unsigned i;

    gf22snod **stab=malloc(tsz*sizeof(gf22snod *));
    for(i=0;i<tsz;++i) 
        stab[i]=NULL; /* _is_ a valid ptr, but it's unallocated. Initialization is possible though. */
    gf22snod *tsnod0, *tsnod2;

    unsigned tint;
    for(i=0; i<numsq; ++i) {
        tint=hashit(gf22[i].i, tsz);
        if( (stab[tint] == NULL) ) {
            stab[tint]=malloc(sizeof(gf22snod));
            stab[tint]->gf22=gf22+i;
            stab[tint]->n=NULL;
            continue;
        }
        tsnod2=stab[tint];
        while( (tsnod2 != NULL) ){
            if(!strcmp(tsnod2->gf22->i, gf22[i].i)) {
                goto nxt;
            }
            tsnod0=tsnod2;
            tsnod2=tsnod2->n;
        }
        tsnod0->n=malloc(sizeof(gf22snod));
        tsnod0->n->gf22=gf22+i;
        tsnod0->n->n=NULL;
nxt:        continue;
    }
    return stab;
}

blosnod **blotochainharr(blop_t *blop, unsigned numsq, unsigned tsz, unsigned char docn /* desired occurrence number */)
{
    unsigned i;

    blosnod **stab=malloc(tsz*sizeof(blosnod *));
    for(i=0;i<tsz;++i)
        stab[i]=NULL; /* _is_ a valid ptr, but it's unallocated. Initialization is possible though. */
    blosnod *tsnod0, *tsnod2;

    unsigned tint;
    for(i=0; i<numsq; ++i) {
        if( blop[i].fssz<3) // pretty much a "no label" blast hit
			continue; 
        tint=hashit(blop[i].qfs, tsz);
        if( (stab[tint] == NULL) ) {
            stab[tint]=malloc(sizeof(blosnod));
            stab[tint]->blop=blop+i;
            stab[tint]->ocn=1;
            stab[tint]->n=NULL;
            continue;
        }
        tsnod2=stab[tint];
        while( (tsnod2 != NULL) ){
            if(!strcmp(tsnod2->blop->qfs, blop[i].qfs)) {
				if(stab[tint]->ocn<docn) { /* less than the desired occurrence */
					stab[tint]->blop=blop+i;
					stab[tint]->ocn++;
				}
				goto nxt;
			}
            tsnod0=tsnod2;
            tsnod2=tsnod2->n;
        }
        tsnod0->n=malloc(sizeof(blosnod));
        tsnod0->n->blop=blop+i;
        tsnod0->n->n=NULL;
nxt:        continue;
    }
    return stab;
}

wseq_t *create_wseq_t(size_t initsz)
{
	wseq_t *words=malloc(sizeof(wseq_t));
	words->wsbuf = initsz;
	words->quan = initsz;
	words->wln=calloc(words->wsbuf, sizeof(size_t));
	words->lbuf=WBUF;
	words->numl=0;
	words->wpla=calloc(words->lbuf, sizeof(size_t));
	return words;
}

void free_wseq(wseq_t *wa)
{
	free(wa->wln);
	free(wa->wpla);
	free(wa);
}

void prtgf22chaharr(gf22snod **stab, unsigned tsz)
{
    unsigned i;
    gf22snod *tsnod2;
    for(i=0;i<tsz;++i) {
        printf("Tablepos %i: ", i); 
        tsnod2=stab[i];
        while(tsnod2) {
            printf("\"%s\" ", tsnod2->gf22->i); 
            tsnod2=tsnod2->n;
        }
        printf("\n"); 
    }
    return;
}

void prtblochaharr(blosnod **stab, unsigned tsz)
{
    unsigned i;
    blosnod *tsnod2;
    for(i=0;i<tsz;++i) {
        printf("Tablepos %i: ", i); 
        tsnod2=stab[i];
        while(tsnod2) {
            printf("%u|%s ", tsnod2->blop->qfnoc, tsnod2->blop->qfs); 
            tsnod2=tsnod2->n;
        }
		putchar('\n');
    }
    printf("\nINFO: table format entry was numoccurrences|queryfeature\n"); 
    return;
}

void freegf22chainharr(gf22snod **stab, size_t tsz)
{
    int i;
    gf22snod *tsnod0, *tsnod2;
    for(i=0; i<tsz; ++i) {
        if( (stab[i] != NULL) ) {
            while( (stab[i]->n != NULL) ) {
                tsnod0=stab[i];
                tsnod2=stab[i]->n;
                while((tsnod2->n != NULL) ){
                    tsnod0=tsnod2;
                    tsnod2=tsnod2->n;
                }
                free(tsnod0->n);
                tsnod0->n=NULL;
            }
            free(stab[i]);
        }
    }
    free(stab);
    return;
}

void freeblochainharr(blosnod **stab, size_t tsz)
{
    int i;
    blosnod *tsnod0, *tsnod2;
    for(i=0; i<tsz; ++i) {
        if( (stab[i] != NULL) ) {
            while( (stab[i]->n != NULL) ) {
                tsnod0=stab[i];
                tsnod2=stab[i]->n;
                while((tsnod2->n != NULL) ){
                    tsnod0=tsnod2;
                    tsnod2=tsnod2->n;
                }
                free(tsnod0->n);
                tsnod0->n=NULL;
            }
            free(stab[i]);
        }
    }
    free(stab);
    return;
}

void prtblop(char *fname, blop_t *blop, int mb)
{
	int i;
	for(i=0;i<mb;++i) // note how we cut out the spurious parts of the motif string to leave it pure and raw (slightly weird why two-char deletion is necessary.
		printf("%i\t%s\t%i\t%s\tqfs:%s\ttfs:%s\t%li\t%li\t%4.2f\t%4.2f\n", blop[i].qfnoc, blop[i].n, blop[i].al, blop[i].tc, blop[i].qfs, blop[i].fs, blop[i].c[0], blop[i].c[1], blop[i].pcti, blop[i].eval);

	printf("You have just seen the %i entries of blast output file called \"%s\".\n", mb, fname); 
	return;
}

void prtgf22(char *fname, gf22_t *gf22, int m7)
{
	int i;
	for(i=0;i<m7;++i) // note how we cut out the spurious parts of the motif string to leave it pure and raw (slightly weird why two-char deletion is necessary.
		printf("%s\tYGD\t%s\t%li\t%li\t.\t%c\t.\tID=%s;NAME=%s\t%i\n", gf22[i].n, gf22[i].t, gf22[i].c[0], gf22[i].c[1], gf22[i].sd, gf22[i].i, gf22[i].i, gf22[i].fc);

#ifdef DBG
	printf("You have just seen the %i entries of basic gff2 file called \"%s\".\n", m7, fname); 
#endif
	return;
}

void prtgf22foz(char *fname, gf22_t *gf22, int m, int n, fcat fc, char *label) /* print the absent ones */
{
	int i, quan=0;
	long rangerep=0; /* representing rang ... */
	printf("%s file called %s is %i rows by %i columns and has following features:\n", label, fname, m, n); 
	printf("You can direct these name into a file wheereupon, you might like to edit it.\n");
	printf("Then re-present these names to ths program with the -u option, whereupon only those names will be looked at\n");
	printf("Note: only feature names without dot and underscore will be printed\n");
	for(i=0;i<m;++i)
		if(gf22[i].fc == fc){
			printf("%s\t%li\t%li\t%c\t%s\t%s\n", gf22[i].n, gf22[i].c[0], gf22[i].c[1], gf22[i].sd, gf22[i].t, gf22[i].i);
			rangerep += gf22[i].c[1] - gf22[i].c[0];
			quan++;
		}
	printf("Total number of annots in %s feature category= %i, representing range of %li basepairs.\n", fcnames[fc-1], quan, rangerep); 

	return;
}

blop_t *processblop(char *fname, int *m, int *n) /* blast output formatting */
{
	/* declarations */
	FILE *fp=fopen(fname,"r");
	int i;
	size_t couc /*count chars per line */, couw=0 /* count words */, oldcouw = 0, nethersz, qfpresz /*query feature pre size */;
	int c;
	int idanomals=0; /* this is for testing col 8 ... ID and string often the same, lie to keep it that way too. This will count when not */
	boole inword=0;
	wseq_t *wa=create_wseq_t(GBUF);
	size_t bwbuf=WBUF;
	char *bufword=calloc(bwbuf, sizeof(char)); /* this is the string we'll keep overwriting. */
	int cncols /* canonical numcols ... we will use the first line */;
	char *tmpc=NULL, *tmpc2=NULL, *tmpc3=NULL, *tmpc4=NULL, *tmpc5=NULL;
	char *oldqfs=calloc(32, sizeof(char));
	unsigned oldqfnoc=0;

	blop_t *blop=calloc(GBUF, sizeof(blop_t));

	while( (c=fgetc(fp)) != EOF) { /* grab a char */
		if( (c== '\n') | (c == ' ') | (c == '\t') | (c=='#')) { /* word closing events */
			if( inword==1) { /* first word closing event */
				wa->wln[couw]=couc;
				bufword[couc++]='\0';
				bufword = realloc(bufword, couc*sizeof(char)); /* normalize */
				/* for the struct, we want to know if it's the first word in a line, like so: */
				if(couw==oldcouw) {
					tmpc=strchr(bufword, ':');
					tmpc2=strchr(tmpc, '|'); // hacky ... avoiding strtok by grabbing r after 2nd |
					tmpc3=strchr(tmpc2, '_'); // token before sart point
					tmpc4=strchr(tmpc3, ':'); // token before end point
					tmpc5=strchr(tmpc4, '|'); // token before end point
					if( (tmpc == NULL) | (tmpc2==NULL) | (tmpc2 <= tmpc) | (tmpc3 == NULL) | (tmpc4==NULL) | (tmpc5 ==NULL) )
						printf("query name did not have a : or | or | appeared before : ... you'll get a segfault.\n");
					blop[wa->numl].n=malloc(couc*sizeof(char));
					blop[wa->numl].nsz=couc;
					strcpy(blop[wa->numl].n, bufword);
					// OK. get the ID string
					// qfpresz=(size_t)(tmpc2-tmpc)-1UL;
					qfpresz=(size_t)(tmpc2-tmpc);
					blop[wa->numl].qfssz=qfpresz;
					blop[wa->numl].qfs=calloc(qfpresz, sizeof(char));
					strncpy(blop[wa->numl].qfs, tmpc+1, qfpresz-1UL);
					if(!strcmp(oldqfs, blop[wa->numl].qfs)) { // number of occurrences
						oldqfnoc++;
						blop[wa->numl].qfnoc = oldqfnoc;
					} else { 
						strcpy(oldqfs, blop[wa->numl].qfs);
						oldqfnoc=0;
					}
					// OK get start and end pos. Notice we will overwrite the bugword here, but it is not being utilised after this.
					// printf("qb=%.*s/qe=%.*s\n", (int)(tmpc4-tmpc3-1), tmpc3+1, (int)(tmpc5-tmpc4-1), tmpc4+1);
					tmpc3[(int)(tmpc4-tmpc3)]='\0';
					tmpc4[(int)(tmpc5-tmpc4)]='\0';
					blop[wa->numl].q[0] = atol(tmpc3+1);
					blop[wa->numl].q[1] = atol(tmpc4+1);
#ifdef DBG2
					printf("qb=%li/qe=%li\n", blop[wa->numl].q[0] , blop[wa->numl].q[1]);
#endif

				} else if((couw-oldcouw)==8) { /* fourth col */
					blop[wa->numl].c[0]=atol(bufword)-1L; // change to zero indexing
				} else if((couw-oldcouw)==9) { /* it's not the first word, and it's 1st and second col */
					blop[wa->numl].c[1]=atol(bufword); // no 0 indexing change required here.
				} else if((couw-oldcouw)==2) { /* 3rd col is pct id */
					blop[wa->numl].pcti=atof(bufword); // no 0 indexing change required here.
				} else if((couw-oldcouw)==10) { /* 3rd col is pct id */
					blop[wa->numl].eval=atof(bufword); // no 0 indexing change required here.
				} else if((couw-oldcouw)==3) { /* 4th col is alignment length */
					blop[wa->numl].al=atoi(bufword); // no 0 indexing change required here.
				} else if( (couw-oldcouw)==1) { // the type string
					tmpc=strrchr(bufword, '|');
					nethersz=(size_t)(tmpc-bufword);
					blop[wa->numl].tcsz=nethersz+1UL;
					blop[wa->numl].fssz=couc - nethersz;
					blop[wa->numl].tc=calloc(blop[wa->numl].tcsz, sizeof(char)); // uninit vgerrs clear by changing from malloc to calloc on these.
					blop[wa->numl].fs=calloc(blop[wa->numl].fssz, sizeof(char));
					strncpy(blop[wa->numl].tc, bufword, nethersz);
					strcpy(blop[wa->numl].fs, tmpc+1);
				}
				couc=0;
				couw++;
			}
			if(c=='#') { /* comment case */
				while( (c=fgetc(fp)) != '\n') ;
				continue;
			} else if(c=='\n') { /* end of a line */
				if(wa->numl == wa->lbuf-1) { /* enought space in our current array? */
					wa->lbuf += WBUF;
					wa->wpla=realloc(wa->wpla, wa->lbuf*sizeof(size_t));
					blop=realloc(blop, wa->lbuf*sizeof(blop_t));
					memset(wa->wpla+(wa->lbuf-WBUF), 0, WBUF*sizeof(size_t));
					memset(blop+(wa->lbuf-WBUF), 0, WBUF*sizeof(blop_t));
				}
				wa->wpla[wa->numl] = couw-oldcouw; /* number of words in current line */
				/* extra bit to catch canonical colnums */
				if(wa->numl==0)
					cncols=couw-oldcouw;
				oldcouw=couw; /* restart words per line count */
				wa->numl++; /* brand new line coming up */
			}
			inword=0;
		} else if(inword==0) { /* deal with first character of new word, + and - also allowed */
			if(couw == wa->wsbuf-1) {
				wa->wsbuf += GBUF;
				wa->wln=realloc(wa->wln, wa->wsbuf*sizeof(size_t));
				blop=realloc(blop, wa->wsbuf*sizeof(blop_t));
				for(i=wa->wsbuf-GBUF;i<wa->wsbuf;++i) {
					wa->wln[i]=0;
					blop[i].qfnoc=0;
				}
			}
			couc=0;
			bwbuf=WBUF;
			bufword=realloc(bufword, bwbuf*sizeof(char)); /* don't bother with memset, it's not necessary */
			bufword[couc++]=c; /* no need to check here, it's the first character */
			inword=1;
		} else {
			if(couc == bwbuf-1) { /* the -1 so that we can always add and extra (say 0) when we want */
				bwbuf += WBUF;
				bufword = realloc(bufword, bwbuf*sizeof(char));
			}
			bufword[couc++]=c;
		}


	} /* end of big for statement */
	fclose(fp);
	free(bufword);

	/* normalization stage */
	wa->quan=couw;
	wa->wln = realloc(wa->wln, wa->quan*sizeof(size_t)); /* normalize */
	blop = realloc(blop, wa->quan*sizeof(blop_t)); /* normalize */
	wa->wpla= realloc(wa->wpla, wa->numl*sizeof(size_t));

	*m= wa->numl;
#ifdef DBG2
	for(i=1;i<wa->numl;++i)
		// if(cncols != wa->wpla[i])
		if(GF22IDCNUM != wa->wpla[i])
			printf("Warning: Numcols is not uniform at %i words per line on all lines. This file has one with %zu.\n", cncols, wa->wpla[i]); 
#endif
	*n= cncols; 
	free_wseq(wa);
	free(oldqfs);

	if(idanomals)
		printf("Warning, idanomals was not zero, it was %i, so some IDs and are not exactly the same as NAMEs.\n", idanomals); 

	return blop;
}

gf22_t *processgf22(char *fname, int *m, int *n) /* this is dummy nmae .. it's for the special gff format for ;yze annottion */
{
	/* In order to make no assumptions, the file is treated as lines containing the same amount of words each,
	 * except for lines starting with #, which are ignored (i.e. comments). These words are checked to make sure they contain only floating number-type
	 * characters [0123456789+-.] only, one string variable is continually written over and copied into a growing floating point array each time */

	/* declarations */
	FILE *fp=fopen(fname,"r");
	int i, j;
	size_t couc /*count chars per line */, couw=0 /* count words */, oldcouw = 0;
	int c;
	int idanomals=0; /* this is for testing col 8 ... ID and string often the same, lie to keep it that way too. This will count when not */
	boole inword=0;
	wseq_t *wa=create_wseq_t(GBUF);
	size_t bwbuf=WBUF;
	char *bufword=calloc(bwbuf, sizeof(char)); /* this is the string we'll keep overwriting. */
	char *sm=NULL /*first semicolon marker*/, *fem=NULL /* first equals marker */, *lem=NULL /* last equals marker */;
	int cncols /* canonical numcols ... we will use the first line */;
	char tch0[16]={0}; /* temporary char holder for attaching numbers to absentidcols */
	int numnoc9=0;

	gf22_t *gf22=malloc(GBUF*sizeof(gf22_t));

	while( (c=fgetc(fp)) != EOF) { /* grab a char */
		if( (c== '\n') | (c == ' ') | (c == '\t') | (c=='#')) { /* word closing events */
			if( inword==1) { /* first word closing event */
				wa->wln[couw]=couc;
				bufword[couc++]='\0';
				bufword = realloc(bufword, couc*sizeof(char)); /* normalize */
				/* for the struct, we want to know if it's the first word in a line, like so: */
				if(couw==oldcouw) {
					gf22[wa->numl].n=malloc(couc*sizeof(char));
					gf22[wa->numl].nsz=couc;
					strcpy(gf22[wa->numl].n, bufword);
					gf22[wa->numl].blofound=0; //seemingly irrelevant I know .. only an initialisation anyhow, useful with -b -y
				} else if((couw-oldcouw)==3) { /* fourth col */
					gf22[wa->numl].c[couw-oldcouw-3]=atol(bufword)-1L; // change to zero indexing
				} else if((couw-oldcouw)==4) { /* it's not the first word, and it's 1st and second col */
					gf22[wa->numl].c[couw-oldcouw-3]=atol(bufword); // no 0 indexing change required here.
				} else if((couw-oldcouw)==6 )  { /* the strand: simply grab first character, it will be + or - */
					gf22[wa->numl].sd=bufword[0];
				} else if( (couw-oldcouw)==2) { // the type string
					gf22[wa->numl].t=malloc(couc*sizeof(char));
					gf22[wa->numl].tsz=couc;
					strcpy(gf22[wa->numl].t, bufword);
				} else if( (couw-oldcouw)==8) { // the ID string at the end
					fem=strchr(bufword, '=');
					lem=strrchr(bufword, '=');
					sm=strchr(bufword, ';');
					gf22[wa->numl].isz=(size_t)(sm-fem);
					gf22[wa->numl].i=malloc(gf22[wa->numl].isz*sizeof(char));
					memcpy(gf22[wa->numl].i, fem+1, (gf22[wa->numl].isz-1)*sizeof(char)); // strncpy writes an extra bit for \0
					gf22[wa->numl].i[gf22[wa->numl].isz-1]='\0'; // null terminate
					/* grab the fcat now, just use first letter */
					for(j=0; j<FCQUAN; j++) {
						if(gf22[wa->numl].i[0] == 'T') {
							gf22[wa->numl].fc=TMR; /* Telomere with any luck */
							break;
						} else if(gf22[wa->numl].i[0] == 'F') {
							gf22[wa->numl].fc=FTH; /* f thingie */
							break;
						} else if(gf22[wa->numl].i[0] == 'Y') {
							gf22[wa->numl].fc=YTH; /* y thingie */
							break;
						} else 
							gf22[wa->numl].fc=NCA; /* no particular category assigned */
					}
					if( strcmp(lem+1, gf22[wa->numl].i))
						idanomals++;
				}
				couc=0;
				couw++;
			}
			if(c=='#') { /* comment case */
				while( (c=fgetc(fp)) != '\n') ;
				continue;
			} else if(c=='\n') { /* end of a line */
				if(wa->numl == wa->lbuf-1) { /* enought space in our current array? */
					wa->lbuf += WBUF;
					wa->wpla=realloc(wa->wpla, wa->lbuf*sizeof(size_t));
					gf22=realloc(gf22, wa->lbuf*sizeof(gf22_t));
					memset(wa->wpla+(wa->lbuf-WBUF), 0, WBUF*sizeof(size_t));
				}
				wa->wpla[wa->numl] = couw-oldcouw; /* number of words in current line */
				/* extra bit to catch canonical colnums */
				if(wa->numl==0)
					cncols=couw-oldcouw;
				else if (couw-oldcouw== GF22IDCNUM-1) {
					// printf("couw-oldcouw=%i\n", couw-oldcouw); 
					gf22[wa->numl].isz=17L;
					gf22[wa->numl].i=malloc(gf22[wa->numl].isz*sizeof(char));
					strcpy(gf22[wa->numl].i, "ABSENTIDCOL"); //* absent ID column */
					gf22[wa->numl].fc=ABS;
					sprintf(tch0, "%05i", numnoc9);
					strcat(gf22[wa->numl].i, tch0);
					numnoc9++;
				}
				oldcouw=couw; /* restart words per line count */
				wa->numl++; /* brand new line coming up */
			}
			inword=0;
		} else if(inword==0) { /* deal with first character of new word, + and - also allowed */
			if(couw == wa->wsbuf-1) {
				wa->wsbuf += GBUF;
				wa->wln=realloc(wa->wln, wa->wsbuf*sizeof(size_t));
				gf22=realloc(gf22, wa->wsbuf*sizeof(gf22_t));
				for(i=wa->wsbuf-GBUF;i<wa->wsbuf;++i)
					wa->wln[i]=0;
			}
			couc=0;
			bwbuf=WBUF;
			bufword=realloc(bufword, bwbuf*sizeof(char)); /* don't bother with memset, it's not necessary */
			bufword[couc++]=c; /* no need to check here, it's the first character */
			inword=1;
		} else {
			if(couc == bwbuf-1) { /* the -1 so that we can always add and extra (say 0) when we want */
				bwbuf += WBUF;
				bufword = realloc(bufword, bwbuf*sizeof(char));
			}
			bufword[couc++]=c;
		}

	} /* end of big for statement */
	fclose(fp);
	free(bufword);

	/* normalization stage */
	wa->quan=couw;
	wa->wln = realloc(wa->wln, wa->quan*sizeof(size_t)); /* normalize */
	gf22 = realloc(gf22, wa->quan*sizeof(gf22_t)); /* normalize */
	wa->wpla= realloc(wa->wpla, wa->numl*sizeof(size_t));

	*m= wa->numl;
#ifdef DBG2
	for(i=1;i<wa->numl;++i)
		// if(cncols != wa->wpla[i])
		if(GF22IDCNUM != wa->wpla[i])
			printf("Warning: Numcols is not uniform at %i words per line on all lines. This file has one with %zu.\n", cncols, wa->wpla[i]); 
#endif
	*n= cncols; 
	free_wseq(wa);

	if(idanomals)
		printf("Warning, idanomals was not zero, it was %i, so some IDs and are not exactly the same as NAMEs.\n", idanomals); 

	return gf22;
}

void mgf222blop_d(char *fname0, char *fname2, gf22_t *gf22, int m7, blosnod **stab, unsigned tsz, blop_t *blop, int mb) /* match blastoutput to the gf22 ... but only is ABSENTCOLS see alternate function for full gf22 search */
{
	unsigned char found;
	int numfound2=0;
	int i, k=0, kk=0, ki=0;
	char *tmpc=NULL;
    blosnod *tsnodd0=NULL, *tsnodd2=NULL;
    unsigned tint;
	unsigned gbuf=GBUF, gbuf2=GBUF;
	int *nfiarr=malloc(gbuf*sizeof(int)); /* the not found index array */
	int *inadarr=malloc(gbuf2*sizeof(int)); /* inadequate finds: due to pct identity being too low or something */
	float pcthresh=85.0, evthresh=0.00001;
	printf("%s gf22 file with %i records wll not be matched against blast outp file %s with %i records at %2.1f identity:\n", fname0, m7, fname2, mb, pcthresh);
	for(i=0;i<m7;++i) {
		if(gf22[i].fc != ABS) { // in this version of thefunction we're not interested in non-ABS
			printf("%s\tYGD\t%s\t%li\t%li\t.\t%c\t.\tID=%s;NAME=%s\tNOTABS\n", gf22[i].n, gf22[i].t, gf22[i].c[0]+1L,gf22[i].c[1], gf22[i].sd, gf22[i].i, gf22[i].i);
			ki++;
			continue;
		}
		found = 0;
        tint=hashit(gf22[i].i, tsz);
        tsnodd2=stab[tint];
        while( (tsnodd2 != NULL) ) {
#ifdef DBG2
			printf("testing %s vs. %s: %li vs. %li: %li vs. %li\n", gf22[i].i, tsnodd2->blop->qfs, gf22[i].c[0], tsnodd2->blop->q[0], gf22[i].c[1], tsnodd2->blop->q[1]);
#endif
            if(!strcmp(gf22[i].i, tsnodd2->blop->qfs)) {
            // if(!(strcmp(gf22[i].n, tsnodd2->blop->n)) && (tsnodd2->blop->q[0] == gf22[i].c[0])  || && (tsnodd2->blop->q[1] == gf22[i].c[1])) {
				tmpc=strchr(tsnodd2->blop->fs, '_'); // i.e. catch out YXXX_mRNA

				// if((tmpc!=NULL) & ((blop[i].pcti>pcthresh) | (blop[i].eval<evthresh)))
				if((tmpc!=NULL) & (blop[i].pcti>pcthresh))
					printf("%s\tYGD\t%s\t%li\t%li\t.\t%c\t.\tID=%.*s;NAME=%.*s\t%4.4f\n", gf22[i].n, gf22[i].t, gf22[i].c[0]+1L,gf22[i].c[1], gf22[i].sd, (int)(tmpc - tsnodd2->blop->fs), tsnodd2->blop->fs, (int)(tmpc - tsnodd2->blop->fs), tsnodd2->blop->fs, tsnodd2->blop->eval); 
				// else if((tmpc==NULL) & ((blop[i].pcti>pcthresh) | (blop[i].eval<evthresh))) // if not _mRNA, don't know what to do: just it all out.
				else if((tmpc==NULL) & (blop[i].pcti>pcthresh))
					printf("%s\tYGD\t%s\t%li\t%li\t.\t%c\t.\tID=%s;NAME=%s\t%4.4f\n", gf22[i].n, gf22[i].t, gf22[i].c[0]+1L,gf22[i].c[1], gf22[i].sd, tsnodd2->blop->fs, tsnodd2->blop->fs, tsnodd2->blop->eval);
				else {
					printf("%s\tYGD\t%s\t%li\t%li\t.\t%c\t.\tHIT1UNDER%3.1f\n", gf22[i].n, gf22[i].t, gf22[i].c[0]+1L,gf22[i].c[1], gf22[i].sd, pcthresh);
					CONDREALLOC(kk, gbuf2, GBUF, inadarr, int);
					inadarr[kk++]=i;
				}
				numfound2++;
				found = 1;
				break; // gets out of the while, not the for
			}
			tsnodd0=tsnodd2;
            tsnodd2=tsnodd2->n;
		}
		if(!(found)) {
			printf("%s\tYGD\t%s\t%li\t%li\t.\t%c\t.\tNOTFOUND\n", gf22[i].n, gf22[i].t, gf22[i].c[0]+1L,gf22[i].c[1], gf22[i].sd);
			CONDREALLOC(k, gbuf, GBUF, nfiarr, int);
			nfiarr[k++]=i;
		}
	}
	if(k)
		nfiarr=realloc(nfiarr, k*sizeof(int)); // k can be zero
	if(kk)
		inadarr=realloc(inadarr, kk*sizeof(int));

	printf("Of %i, features ignored: %i; features found: %i; found but inadequate: %i; not found at all: %i\n", m7, ki, numfound2, kk, k);

	free(nfiarr);
	free(inadarr);

	return;
}

void mgf222blop(char *fname0, char *fname2, gf22_t *gf22, int m7, blosnod **stab, unsigned tsz, blop_t *blop, int mb) /* match blastoutput to the gf22 ... but only is ABSENTCOLS see alternate function for full gf22 search */
{
	unsigned char found;
	int numfound2=0;
	int i, k=0, kk=0, ki=0;
	char *tmpc=NULL;
    blosnod *tsnodd0=NULL, *tsnodd2=NULL;
    unsigned tint;
	unsigned gbuf=GBUF, gbuf2=GBUF;
	int *nfiarr=malloc(gbuf*sizeof(int)); /* the not found index array */
	int *inadarr=malloc(gbuf2*sizeof(int)); /* inadequate finds: due to pct identity being too low or something */
	float pcthresh=85.0, evthresh=0.00001;
	printf("%s gf22 file with %i records wll not be matched against blast outp file %s with %i records at %2.1f identity:\n", fname0, m7, fname2, mb, pcthresh);
	for(i=0;i<m7;++i) {
		if(gf22[i].fc != ABS) { // in this version of thefunction we're not interested in non-ABS
			printf("%s\tYGD\t%s\t%li\t%li\t.\t%c\t.\tID=%s;NAME=%s\n", gf22[i].n, gf22[i].t, gf22[i].c[0]+1L,gf22[i].c[1], gf22[i].sd, gf22[i].i, gf22[i].i);
			ki++;
			continue;
		}
		found = 0;
        tint=hashit(gf22[i].i, tsz);
        tsnodd2=stab[tint];
        while( (tsnodd2 != NULL) ) {
            if(!strcmp(gf22[i].i, tsnodd2->blop->qfs)) {
            // if(!(strcmp(gf22[i].n, tsnodd2->blop->n)) && (tsnodd2->blop->q[0] == gf22[i].c[0])  || && (tsnodd2->blop->q[1] == gf22[i].c[1])) {
				tmpc=strchr(tsnodd2->blop->fs, '_'); // i.e. catch out YXXX_mRNA

				// if((tmpc!=NULL) & ((blop[i].pcti>pcthresh) | (blop[i].eval<evthresh)))
				if((tmpc!=NULL) & (blop[i].pcti>pcthresh))
					printf("%s\tYGD\t%s\t%li\t%li\t.\t%c\t.\tID=%.*s;NAME=%.*s\n", gf22[i].n, gf22[i].t, gf22[i].c[0]+1L,gf22[i].c[1], gf22[i].sd, (int)(tmpc - tsnodd2->blop->fs), tsnodd2->blop->fs, (int)(tmpc - tsnodd2->blop->fs), tsnodd2->blop->fs);
				// else if((tmpc==NULL) & ((blop[i].pcti>pcthresh) | (blop[i].eval<evthresh))) // if not _mRNA, don't know what to do: just it all out.
				else if((tmpc==NULL) & (blop[i].pcti>pcthresh))
					printf("%s\tYGD\t%s\t%li\t%li\t.\t%c\t.\tID=%s;NAME=%s\n", gf22[i].n, gf22[i].t, gf22[i].c[0]+1L,gf22[i].c[1], gf22[i].sd, tsnodd2->blop->fs, tsnodd2->blop->fs);
				else {
					printf("%s\tYGD\t%s\t%li\t%li\t.\t%c\t.\tHITUND%3.1fID\n", gf22[i].n, gf22[i].t, gf22[i].c[0]+1L,gf22[i].c[1], gf22[i].sd, pcthresh);
					CONDREALLOC(kk, gbuf2, GBUF, inadarr, int);
					inadarr[kk++]=i;
				}
				numfound2++;
				found = 1;
				break; // gets out of the while, not the for
			}
			tsnodd0=tsnodd2;
            tsnodd2=tsnodd2->n;
		}
		if(!(found)) {
			printf("%s\tYGD\t%s\t%li\t%li\t.\t%c\t.\tNOHITATALL\n", gf22[i].n, gf22[i].t, gf22[i].c[0]+1L,gf22[i].c[1], gf22[i].sd);
			CONDREALLOC(k, gbuf, GBUF, nfiarr, int);
			nfiarr[k++]=i;
		}
	}
	if(k)
		nfiarr=realloc(nfiarr, k*sizeof(int)); // k can be zero
	if(kk)
		inadarr=realloc(inadarr, kk*sizeof(int));

	printf("Of %i, features ignored: %i; features found: %i; found but inadequate: %i; not found at all: %i\n", m7, ki, numfound2, kk, k);

	free(nfiarr);
	free(inadarr);

	return;
}

void prtusage()
{
	printf("hblop: hashes a blast output and matches against a GTF annotation file called gf22.\n");
	printf("-b, this is for a blast output file.\n");
	printf("-y: takes a variation of the GFF2 format, a table, column nine gives the Y-name. eg. Output of YG annotator service\n");
	printf("-d: is a flag that merely prints out details of either the -b or -y input file.\n");
	return;
}

int main(int argc, char *argv[])
{
	/* argument accounting */
	if(argc == 1) {
		prtusage();
		exit(EXIT_FAILURE);
	}
	int i, mb, nb, m7, n7;
	opt_t opts={0};
	catchopts(&opts, argc, argv);

	/* declare types to NULL */
	gf22_t *gf22=NULL;
	blop_t *blop=NULL;

	/* declare hashtable parts */
    unsigned htsz, htszblo;
    gf22snod **stab=NULL;
    blosnod **stabblo=NULL;

	if(opts.bstr) {
		blop=processblop(opts.bstr, &mb, &nb);
		htszblo=2*mb/3;
		unsigned char docn =2;
		stabblo=blotochainharr(blop, mb, htszblo, docn);
	}
	if(opts.ystr) { // we're goign to try hashing the ID line of gf22 format */
		gf22=processgf22(opts.ystr, &m7, &n7);
		htsz=2*m7/3; /* our hash table size */
		stab=gf22tochainharr(gf22, m7, htsz); /* we now set up a hash table along side our sequence names from the fasta file */
	}

	/* OK now do stuff, mostly print */
	if((opts.dflg) && (opts.ystr) ) {
		prtgf22(opts.ystr, gf22, m7);
		// prtgf22chaharr(stab, htsz);
	}
	if((opts.dflg) && (opts.bstr) )
		prtblochaharr(stabblo, htszblo);

	/* Now do real stuff */
	if((opts.ystr) && (opts.bstr) ) {
#ifdef DBG
		mgf222blop_d(opts.ystr, opts.bstr, gf22, m7, stabblo, htszblo, blop, mb);
#else
		mgf222blop(opts.ystr, opts.bstr, gf22, m7, stabblo, htszblo, blop, mb);
#endif

	}

final:
	if(opts.ystr) {
		freegf22chainharr(stab, htsz);
		for(i=0;i<m7;++i) {
			free(gf22[i].n);
			free(gf22[i].t);
			free(gf22[i].i);
		}
		free(gf22);
	}
	if(opts.bstr) {
		freeblochainharr(stabblo, htszblo);
		for(i=0;i<mb;++i) {
			free(blop[i].n);
			free(blop[i].tc);
			free(blop[i].fs);
			free(blop[i].qfs);
		}
		free(blop);
	}

	return 0;
}
