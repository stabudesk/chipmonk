/* chipmonk.c: a program to "tackle" BED (genomic features file) files.
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

typedef struct /* onefa */
{
	char *id;
	char *sq;
	unsigned idz, sqz;
} onefa;

typedef struct /* i_s; sequence index and number of symbols */
{
	unsigned idx;
	size_t sylen; /* this is the precise symbol length of the sequence */
	size_t sy[SSZ]; /* used to hold counts of symbols */
	float cgp;
	unsigned ambano[2]; /* number of ambiguous symbols (first), then number of anomalous symbols */
	char *id; // the ID name
	char *sq; // the sequence itself
	unsigned ibf, sbf; // buffers for ID and SQ strings
	unsigned idz, sqz; // actual size  of the ID and SQ strings. Is almost a duplication of sylen, can be removed if memory is a consideration.
} i_s; /* sequence index and number of symbols */

typedef struct /* ia_t integer array type, includes iab the buffer */
{
	int *a;
	unsigned b /* int array buf */, z /* int array size*/;
} ia_t;

typedef struct /* opt_t, a struct for the options */
{
	boole dflg; /* details / information only */
	boole nflg; /* feature names only */
	boole sflg; /* split outout in two files */
	char *istr; /* first bedgraph file, the target of the filtering by the second */
	char *fstr; /* the name of the feature bed file, often converted from GFF3 */
	char *ustr; /* the name of a file with the list of elements to be unified */
	char *pstr; /* depth file name */
	char *qstr; /* 2nddepth file name */
	char *gstr; /* genome file name */
	char *rstr; /* repeatmasker gtf/gff2 file */
	char *ystr; /* the gf22_t type */
	char *astr; /* a fasta file */
} opt_t;

typedef struct /* i4_t */
{
	int sc; /* number of same chromosome (occurences?) */
	float mc; /* min signal value */
	int b1i; /* index of the 1st bgr_t, which satisfies the conditions */
	int lgbi; /* last good bgr_t index */
} i4_t; /* 4 vals of some sort? */

typedef struct /* bgr_t */
{
	char *n;
	size_t nsz; /* size of the name r ID field */
	long c[2]; /* coords: 1) start 2) end */
	float co; /* signal value */
} bgr_t; /* bedgraph row type */

typedef struct /* bgr_t2 */
{
	char *n;
	size_t nsz; /* size of the name r ID field */
	long c[2]; /* coords: 1) start 2) end */
	char *f; /* f for feature .. 4th col */
	size_t fsz; /* size of the feature field*/
} bgr_t2; /* bedgraph row type 2i. column is the feature */

typedef struct /* rmf_t: repeatmasker gff2 file format */
{
	char *n;
	size_t nsz; /* size of the name r ID field */
	long c[2]; /* coords: 1) start 2) cols 4 and 5 */
	char *m; /* the motif string ... 9th column */
	char sd; /* strand + or - */
	size_t msz; /* size of motif string */
} rmf_t;

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
} gf22_t;

typedef struct /* words_t: file with only single words per line */
{
	char *n;
	size_t nsz; /* size of the name r ID field */
} words_t; /* bedgraph row type */

typedef struct /* dpf_t : depth file type ... just chr name, pos and read quant */
{
	char *n;
	size_t nsz; /* size of the name r ID field */
	long p; /* position */
	int d; /* depth reading */
} dpf_t;

typedef struct /* gf_t : genome file type ... just chr name, pos and read quant */
{
	char *n;
	size_t nsz; /* size of the name r ID field */
	long z; /* size of the the chromosome */
} gf_t;

typedef struct /* wseq_t */
{
	size_t *wln;
	size_t wsbuf;
	size_t quan;
	size_t lbuf; /* a buffer for the number of lines */
	size_t numl; /* number of lines, i.e. rows */
	size_t *wpla; /* words per line array: the number of words on each line */
} wseq_t;

struct strchainode /* gf22snod struct */
{
    gf22_t *gf22; /* ptr to a single element, not an array, TODO void* it. */
    struct strchainode *n;
};
typedef struct strchainode gf22snod; /* yes, leave this alone, it's the way a struct can have a ptr ot its own type! */

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

char idlsearch(gf22snod **stab, unsigned tsz, char *line, unsigned lnsz)
{
    char yes=0;
    size_t i;
    gf22snod *tsnod2;

    unsigned tint;
    for(i=0; i<tsz; ++i) {
        tint=hashit(line, tsz);
        if( (stab[tint] == NULL) )
            goto nxt; /* hashtable at that position is empty ... so it's impossible for that string to be there */

        tsnod2=stab[tint];
        while( (tsnod2 != NULL) ) {
            if( (tsnod2->gf22->isz == lnsz) & !(strncmp(tsnod2->gf22->i, line, tsnod2->gf22->isz)) ) {
                yes=1; /* i.e. no=0 */
                goto nxt;
            }
            tsnod2=tsnod2->n;
        }
    }
nxt:    return yes;
}

void prtfa(onefa *fac)
{
	printf(">");
	printf("%s\n", fac->id);
	printf("%s\n", fac->sq);
}

void prtfaf(onefa *fac, FILE *fp)
{
	fprintf(fp, ">");
	fprintf(fp, "%s\n", fac->id);
	fprintf(fp, "%s\n", fac->sq);
}

void prtfa2(onefa *fac)
{
	int i;
	printf("SQZ=%d:", fac->sqz);
	for(i=0;i<3;++i) 
		putchar(fac->sq[i]);
	printf("\n"); 
}

void prtsq(i_s *sqisz, int sz)
{
	int i;
	printf("Number of different sequences=%i\n", sz); 
#ifdef DBG
	for(i=0;i<sz;++i) {
		printf("%s\n", sqisz[i].id);
		printf("%s\n", sqisz[i].sq);
	}
#endif
	return;
}

void prtsqbdg(i_s *sqisz, bgr_t *bgrow, int m, int sz)
{
	int i, j;
	char rangestr[64]={0}; /* generally helpful to say what range is being given */
	for(i=0;i<sz;++i) {
		for(j=0;j<m;++j) 
			if(!strcmp(sqisz[i].id, bgrow[j].n)) {
				sprintf(rangestr, "|range_%li_to_%li", bgrow[j].c[0], bgrow[j].c[1]);
				printf(">%s", sqisz[i].id);
				printf("%s\n", rangestr);
				printf("%.*s\n", (int)(bgrow[j].c[1]-bgrow[j].c[0]), sqisz[i].sq+bgrow[j].c[0]);
				break;
			}
	}
	return;
}

void prti_s(i_s *sqisz, int sz, float *mxcg, float *mncg)
{
	int i;
	char *sqgood;
	*mxcg=.0;
	*mncg=1.;

	size_t tsz;
	for(i=0;i<sz;++i) {
		if(sqisz[i].ambano[1] != 0)
			sqgood="AnoSQ";
		else
			sqgood="SQ";
		tsz = sqisz[i].sy[0] + sqisz[i].sy[1];
		sqisz[i].cgp=(float)sqisz[i].sy[0]/tsz;
		if(sqisz[i].cgp>*mxcg)
			*mxcg=sqisz[i].cgp;
		if(sqisz[i].cgp<*mncg)
			*mncg=sqisz[i].cgp;

		printf("| %s#%i=TOT:%zu CG:%.4f ", sqgood, i, sqisz[i].sylen, sqisz[i].cgp);
	}
	printf("|\n"); 
}

i_s *procfa(char *fname, unsigned *nsq)
{
	FILE *fin;
	char IGLINE, begline;
	size_t lidx, mxsylen, mnsylen;
	unsigned mxamb, mnamb;
	int i, c, sqidx;
	int gbuf;
	i_s *sqisz=NULL;
	int whatint; // a number reflecting the type of symbol read
	unsigned numsq, numano;
	int ididx0=0;

	// OK open the file
	if(!(fin=fopen(fname, "r")) ) { /*should one check the extension of the fasta file ? */
		printf("Error. Cannot open file named \"%s\".\n", fname);
		exit(EXIT_FAILURE);
	}

	IGLINE=0, begline=1;
	lidx=0, mxsylen=0, mnsylen=0XFFFFFFFFFFFFFFFF;
	mxamb=0, mnamb=0xFFFFFFFF;

	sqidx=-1; /* this is slightly dangerous, you need very much to know what you're doing */
	gbuf=GBUF;
	// sqisz=malloc(gbuf*sizeof(i_s));
	sqisz=realloc(sqisz, gbuf*sizeof(i_s));
	for(i=0;i<gbuf;++i) {
		sqisz[i].ibf=GBUF;
		sqisz[i].sbf=GBUF;
		sqisz[i].id=calloc(sqisz[i].ibf, sizeof(char));
		sqisz[i].sq=calloc(sqisz[i].sbf, sizeof(char));
	}
	for(i=gbuf-GBUF;i<gbuf;++i) {
		sqisz[i].ambano[0]=0;
		sqisz[i].ambano[1]=0;
	}
	whatint=0; /* needs explanation */
	ididx0=0;

	while( ( (c = fgetc(fin)) != EOF) ) {
		if(c =='\n') {
			IGLINE=0;
			begline=1;
			lidx++;
		} else if( (begline==1) & (c == '>') ) { /* this condition catches the beginning of a new sequence, and uses it to prepare the nextsequence.*/
			IGLINE =1;
			begline=0; 
			if(sqidx>=0) { /* chancing my arm here ... operating on the past sequence */
				if(sqisz[sqidx].sylen > mxsylen)
					mxsylen = sqisz[sqidx].sylen;
				if(sqisz[sqidx].sylen < mnsylen)
					mnsylen = sqisz[sqidx].sylen;
				if(sqisz[sqidx].ambano[0] > mxamb)
					mxamb = sqisz[sqidx].ambano[0];
				if(sqisz[sqidx].ambano[0] < mnamb)
					mnamb = sqisz[sqidx].ambano[0];

				CONDREALLOC(ididx0, sqisz[sqidx].ibf, GBUF, sqisz[sqidx].id, char);
				sqisz[sqidx].id[ididx0]='\0';
				CONDREALLOC(sqisz[sqidx].sylen, sqisz[sqidx].sbf, GBUF, sqisz[sqidx].sq, char);
				sqisz[sqidx].sq[sqisz[sqidx].sylen]='\0';
				sqisz[sqidx].idz=1+ididx0;
				sqisz[sqidx].sqz=1+sqisz[sqidx].sylen;
			}

			sqidx++;
			if(sqidx==gbuf) {
				gbuf+=GBUF;
				sqisz=realloc(sqisz, gbuf*sizeof(i_s));
				for(i=gbuf-GBUF;i<gbuf;++i) {
					sqisz[i].ibf=GBUF;
					sqisz[i].sbf=GBUF;
					sqisz[i].id=calloc(sqisz[i].ibf, sizeof(char));
					sqisz[i].sq=calloc(sqisz[i].sbf, sizeof(char));
				}
			}
			sqisz[sqidx].idx=sqidx;

			/* resetting stuff */
			sqisz[sqidx].sylen=0;
			ididx0=0;
			for(i=0;i<SSZ;++i)
				sqisz[sqidx].sy[i]=0;
			for(i=0;i<2;++i)
				sqisz[sqidx].ambano[i]=0;
		} else if (IGLINE==1) {
			CONDREALLOC(ididx0, sqisz[sqidx].ibf, GBUF, sqisz[sqidx].id, char);
			sqisz[sqidx].id[ididx0++]=c;
		} else if (IGLINE==0) {
			CONDREALLOC(sqisz[sqidx].sylen, sqisz[sqidx].sbf, GBUF, sqisz[sqidx].sq, char);
			sqisz[sqidx].sq[sqisz[sqidx].sylen]=c;
			sqisz[sqidx].sylen++;
			switch(c) {
				case 'A': case 'a':
					whatint=1; break;
				case 'C': case 'c':
					whatint=2; break;
				case 'G': case 'g':
					whatint=3; break;
				case 'T': case 't':
					whatint=4; break;
				case 'R': case 'r':
					whatint=5; break;
				case 'Y': case 'y':
					whatint=6; break;
				case 'K': case 'k': /* the ketos */
					whatint=7; break;
				case 'M': case 'm': /* the aminoids */
					whatint=8; break;
				case 'S': case 's':
					whatint=9; break;
				case 'W': case 'w':
					whatint=10; break;
				case 'B': case 'b':
					whatint=11; break;
				case 'D': case 'd':
					whatint=12; break;
				case 'H': case 'h':
					whatint=13; break;
				case 'V': case 'v':
					whatint=14; break;
				case 'N': case 'n':
					whatint=15; break;
				case '-':
					whatint=16; break;
				default:
					whatint=17; /* unknown this means your fasta file is naff. */
			}
		}
		if( (whatint == 2) || (whatint == 3) ) {
			sqisz[sqidx].sy[0]++;
			sqisz[sqidx].ambano[1]++;
		} else if (whatint < 5) {
			sqisz[sqidx].sy[1]++;
			sqisz[sqidx].ambano[1]++;
		} else 
			sqisz[sqidx].ambano[0]++;
	}
	fclose(fin);
	/* postprocessing on the final sequence */
	if(sqisz[sqidx].sylen > mxsylen)
		mxsylen = sqisz[sqidx].sylen;
	if(sqisz[sqidx].sylen < mnsylen)
		mnsylen = sqisz[sqidx].sylen;
	if(sqisz[sqidx].ambano[0] > mxamb)
		mxamb = sqisz[sqidx].ambano[0];
	if(sqisz[sqidx].ambano[0] < mnamb)
		mnamb = sqisz[sqidx].ambano[0];

	/* the last sequence */
	CONDREALLOC(ididx0, sqisz[sqidx].ibf, GBUF, sqisz[sqidx].id, char);
	sqisz[sqidx].id[ididx0]='\0';
	CONDREALLOC(sqisz[sqidx].sylen, sqisz[sqidx].sbf, GBUF, sqisz[sqidx].sq, char);
	sqisz[sqidx].sq[sqisz[sqidx].sylen]='\0';
	sqisz[sqidx].idz=1+ididx0;
	sqisz[sqidx].sqz=1+sqisz[sqidx].sylen;

	numsq=sqidx+1, numano=0;
	for(i=0;i<numsq;++i) {
		if(sqisz[i].ambano[1])
			numano++;
	}

	for(i=numsq;i<gbuf;++i) {
		free(sqisz[i].id);
		free(sqisz[i].sq);
	}
	sqisz=realloc(sqisz, numsq*sizeof(i_s));

    *nsq=numsq;
	return sqisz;
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

int catchopts(opt_t *opts, int oargc, char **oargv)
{
	int c;
	opterr = 0;

	while ((c = getopt (oargc, oargv, "dsni:f:u:p:g:r:q:y:a:")) != -1)
		switch (c) {
			case 'd':
				opts->dflg = 1;
				break;
			case 's':
				opts->sflg = 1;
				break;
			case 'n':
				opts->nflg = 1;
				break;
			case 'i':
				opts->istr = optarg;
				break;
			case 'f':
				opts->fstr = optarg;
				break;
			case 'u': /* unify certain bed2 elements into one file */
				opts->ustr = optarg;
				break;
			case 'p': /* depth file */
				opts->pstr = optarg;
				break;
			case 'q': /* 2nd depth file */
				opts->qstr = optarg;
				break;
			case 'g': /* genome file */
				opts->gstr = optarg;
				break;
			case 'r': /* repeatmasker gff2 file */
				opts->rstr = optarg;
				break;
			case 'y': /* based on repeatmasker gff2 file */
				opts->ystr = optarg;
				break;
			case 'a':
				opts->astr = optarg;
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

void free_wseq(wseq_t *wa)
{
	free(wa->wln);
	free(wa->wpla);
	free(wa);
}

words_t *processwordf(char *fname, int *m, int *n)
{
	/* In order to make no assumptions, the file is treated as lines containing the same amount of words each,
	 * except for lines starting with #, which are ignored (i.e. comments). These words are checked to make sure they contain only floating number-type
	 * characters [0123456789+-.] only, one string variable is continually written over and copied into a growing floating point array each time */

	/* declarations */
	FILE *fp=fopen(fname,"r");
	int i;
	size_t couc /*count chars per line */, couw=0 /* count words */, oldcouw = 0;
	int c;
	boole inword=0;
	wseq_t *wa=create_wseq_t(GBUF);
	size_t bwbuf=WBUF;
	char *bufword=calloc(bwbuf, sizeof(char)); /* this is the string we'll keep overwriting. */

	words_t *bedword=malloc(GBUF*sizeof(words_t));

	while( (c=fgetc(fp)) != EOF) { /* grab a char */
		if( (c== '\n') | (c == ' ') | (c == '\t') | (c=='#')) { /* word closing events */
			if( inword==1) { /* first word closing event */
				wa->wln[couw]=couc;
				bufword[couc++]='\0';
				bufword = realloc(bufword, couc*sizeof(char)); /* normalize */
				/* for the struct, we want to know if it's the first word in a line, like so: */
				if(couw==oldcouw) {
					bedword[wa->numl].n=malloc(couc*sizeof(char));
					bedword[wa->numl].nsz=couc;
					strcpy(bedword[wa->numl].n, bufword);
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
					bedword=realloc(bedword, wa->lbuf*sizeof(words_t));
					memset(wa->wpla+(wa->lbuf-WBUF), 0, WBUF*sizeof(size_t));
				}
				wa->wpla[wa->numl] = couw-oldcouw; /* number of words in current line */
				if(couw-oldcouw >4) {
					printf("Error, each row cannot exceed 4 words: revise your input file\n"); 
					/* need to release all memory too */
					free_wseq(wa);
					exit(EXIT_FAILURE);
				}
				oldcouw=couw; /* restart words per line count */
				wa->numl++; /* brand new line coming up */
			}
			inword=0;
		} else if(inword==0) { /* deal with first character of new word, + and - also allowed */
			if(couw == wa->wsbuf-1) {
				wa->wsbuf += GBUF;
				wa->wln=realloc(wa->wln, wa->wsbuf*sizeof(size_t));
				bedword=realloc(bedword, wa->wsbuf*sizeof(words_t));
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
	bedword = realloc(bedword, wa->quan*sizeof(words_t)); /* normalize */
	wa->wpla= realloc(wa->wpla, wa->numl*sizeof(size_t));

	*m= wa->numl;
	int k=wa->wpla[0];
	for(i=1;i<wa->numl;++i)
		if(k != wa->wpla[i])
			printf("Warning: Numcols is not uniform at %i words per line on all lines. This file has one with %zu.\n", k, wa->wpla[i]); 
	*n= k; 
	free_wseq(wa);

	return bedword;
}

bgr_t *processinpf(char *fname, int *m, int *n)
{
	/* In order to make no assumptions, the file is treated as lines containing the same amount of words each,
	 * except for lines starting with #, which are ignored (i.e. comments). These words are checked to make sure they contain only floating number-type
	 * characters [0123456789+-.] only, one string variable is continually written over and copied into a growing floating point array each time */

	/* declarations */
	FILE *fp=fopen(fname,"r");
	int i;
	size_t couc /*count chars per line */, couw=0 /* count words */, oldcouw = 0;
	int c;
	boole inword=0;
	wseq_t *wa=create_wseq_t(GBUF);
	size_t bwbuf=WBUF;
	char *bufword=calloc(bwbuf, sizeof(char)); /* this is the string we'll keep overwriting. */

	bgr_t *bgrow=malloc(GBUF*sizeof(bgr_t));

	while( (c=fgetc(fp)) != EOF) { /* grab a char */
		if( (c== '\n') | (c == ' ') | (c == '\t') | (c=='#')) { /* word closing events */
			if( inword==1) { /* first word closing event */
				wa->wln[couw]=couc;
				bufword[couc++]='\0';
				bufword = realloc(bufword, couc*sizeof(char)); /* normalize */
				/* for the struct, we want to know if it's the first word in a line, like so: */
				if(couw==oldcouw) {
					bgrow[wa->numl].n=malloc(couc*sizeof(char));
					bgrow[wa->numl].nsz=couc;
					strcpy(bgrow[wa->numl].n, bufword);
				} else if((couw-oldcouw)<3) /* it's not the first word, and it's 1st and second col */
					bgrow[wa->numl].c[couw-oldcouw-1]=atol(bufword);
				else if( (couw-oldcouw)==3) { // assume float
					bgrow[wa->numl].co=atof(bufword);
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
					bgrow=realloc(bgrow, wa->lbuf*sizeof(bgr_t));
					memset(wa->wpla+(wa->lbuf-WBUF), 0, WBUF*sizeof(size_t));
				}
				wa->wpla[wa->numl] = couw-oldcouw; /* number of words in current line */
				if(couw-oldcouw >4) {
					printf("Error, each row cannot exceed 4 words: revise your input file\n"); 
					/* need to release all memory too */
					free_wseq(wa);
					exit(EXIT_FAILURE);
				}
				oldcouw=couw; /* restart words per line count */
				wa->numl++; /* brand new line coming up */
			}
			inword=0;
		} else if(inword==0) { /* deal with first character of new word, + and - also allowed */
			if(couw == wa->wsbuf-1) {
				wa->wsbuf += GBUF;
				wa->wln=realloc(wa->wln, wa->wsbuf*sizeof(size_t));
				bgrow=realloc(bgrow, wa->wsbuf*sizeof(bgr_t));
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
	bgrow = realloc(bgrow, wa->quan*sizeof(bgr_t)); /* normalize */
	wa->wpla= realloc(wa->wpla, wa->numl*sizeof(size_t));

	*m= wa->numl;
	int k=wa->wpla[0];
	for(i=1;i<wa->numl;++i)
		if(k != wa->wpla[i])
			printf("Warning: Numcols is not uniform at %i words per line on all lines. This file has one with %zu.\n", k, wa->wpla[i]); 
	*n= k; 
	free_wseq(wa);

	return bgrow;
}

bgr_t2 *processinpf2(char *fname, int *m, int *n) /*fourth column is string, other columns to be ignored */
{
	/* In order to make no assumptions, the file is treated as lines containing the same amount of words each,
	 * except for lines starting with #, which are ignored (i.e. comments). These words are checked to make sure they contain only floating number-type
	 * characters [0123456789+-.] only, one string variable is continually written over and copied into a growing floating point array each time */

	/* declarations */
	FILE *fp=fopen(fname,"r");
	int i;
	size_t couc /*count chars per line */, couw=0 /* count words */, oldcouw = 0;
	int c;
	boole inword=0;
	wseq_t *wa=create_wseq_t(GBUF);
	size_t bwbuf=WBUF;
	char *bufword=calloc(bwbuf, sizeof(char)); /* this is the string we'll keep overwriting. */

	bgr_t2 *bgrow=malloc(GBUF*sizeof(bgr_t2));

	while( (c=fgetc(fp)) != EOF) { /* grab a char */
		if( (c== '\n') | (c == ' ') | (c == '\t') | (c=='#')) { /* word closing events */
			if( inword==1) { /* first word closing event */
				wa->wln[couw]=couc;
				bufword[couc++]='\0';
				bufword = realloc(bufword, couc*sizeof(char)); /* normalize */
				/* for the struct, we want to know if it's the first word in a line, like so: */
				if(couw==oldcouw) {
					bgrow[wa->numl].n=malloc(couc*sizeof(char));
					bgrow[wa->numl].nsz=couc;
					strcpy(bgrow[wa->numl].n, bufword);
				} else if((couw-oldcouw)<3) { /* it's not the first word, and it's 1st and second col */
					bgrow[wa->numl].c[couw-oldcouw-1]=atol(bufword);
				} else if( (couw-oldcouw)==3) { // assume float
					bgrow[wa->numl].f=malloc(couc*sizeof(char));
					bgrow[wa->numl].fsz=couc;
					strcpy(bgrow[wa->numl].f, bufword);
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
					bgrow=realloc(bgrow, wa->lbuf*sizeof(bgr_t2));
					memset(wa->wpla+(wa->lbuf-WBUF), 0, WBUF*sizeof(size_t));
				}
				wa->wpla[wa->numl] = couw-oldcouw; /* number of words in current line */
				oldcouw=couw; /* restart words per line count */
				wa->numl++; /* brand new line coming up */
			}
			inword=0;
		} else if(inword==0) { /* deal with first character of new word, + and - also allowed */
			if(couw == wa->wsbuf-1) {
				wa->wsbuf += GBUF;
				wa->wln=realloc(wa->wln, wa->wsbuf*sizeof(size_t));
				bgrow=realloc(bgrow, wa->wsbuf*sizeof(bgr_t2));
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
	bgrow = realloc(bgrow, wa->quan*sizeof(bgr_t2)); /* normalize */
	wa->wpla= realloc(wa->wpla, wa->numl*sizeof(size_t));

	*m= wa->numl;
	int k=wa->wpla[0];
	for(i=1;i<wa->numl;++i)
		if(k != wa->wpla[i])
			printf("Warning: Numcols is not uniform at %i words per line on all lines. This file has one with %zu.\n", k, wa->wpla[i]); 
	*n= k; 
	free_wseq(wa);

	return bgrow;
}

rmf_t *processrmf(char *fname, int *m, int *n) /*fourth column is string, other columns to be ignored */
{
	/* In order to make no assumptions, the file is treated as lines containing the same amount of words each,
	 * except for lines starting with #, which are ignored (i.e. comments). These words are checked to make sure they contain only floating number-type
	 * characters [0123456789+-.] only, one string variable is continually written over and copied into a growing floating point array each time */

	/* declarations */
	FILE *fp=fopen(fname,"r");
	int i;
	size_t couc /*count chars per line */, couw=0 /* count words */, oldcouw = 0;
	int c;
	boole inword=0;
	wseq_t *wa=create_wseq_t(GBUF);
	size_t bwbuf=WBUF;
	char *bufword=calloc(bwbuf, sizeof(char)); /* this is the string we'll keep overwriting. */

	rmf_t *rmf=malloc(GBUF*sizeof(rmf_t));

	while( (c=fgetc(fp)) != EOF) { /* grab a char */
		if( (c== '\n') | (c == ' ') | (c == '\t') | (c=='#')) { /* word closing events */
			if( inword==1) { /* first word closing event */
				wa->wln[couw]=couc;
				bufword[couc++]='\0';
				bufword = realloc(bufword, couc*sizeof(char)); /* normalize */
				/* for the struct, we want to know if it's the first word in a line, like so: */
				if(couw==oldcouw) {
					rmf[wa->numl].n=malloc(couc*sizeof(char));
					rmf[wa->numl].nsz=couc;
					strcpy(rmf[wa->numl].n, bufword);
				} else if((couw-oldcouw)==3) { /* fourth col */
					rmf[wa->numl].c[couw-oldcouw-3]=atol(bufword)-1L; // change to zero indexing
				} else if((couw-oldcouw)==4) { /* it's not the first word, and it's 1st and second col */
					rmf[wa->numl].c[couw-oldcouw-3]=atol(bufword); // no 0 indexing change required here.
				} else if((couw-oldcouw)==6 )  { /* the strand */
					rmf[wa->numl].sd=bufword[0];
				} else if( (couw-oldcouw)==9) { // the motif string
					rmf[wa->numl].m=malloc(couc*sizeof(char));
					rmf[wa->numl].msz=couc;
					strcpy(rmf[wa->numl].m, bufword);
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
					rmf=realloc(rmf, wa->lbuf*sizeof(rmf_t));
					memset(wa->wpla+(wa->lbuf-WBUF), 0, WBUF*sizeof(size_t));
				}
				wa->wpla[wa->numl] = couw-oldcouw; /* number of words in current line */
				oldcouw=couw; /* restart words per line count */
				wa->numl++; /* brand new line coming up */
			}
			inword=0;
		} else if(inword==0) { /* deal with first character of new word, + and - also allowed */
			if(couw == wa->wsbuf-1) {
				wa->wsbuf += GBUF;
				wa->wln=realloc(wa->wln, wa->wsbuf*sizeof(size_t));
				rmf=realloc(rmf, wa->wsbuf*sizeof(rmf_t));
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
	rmf = realloc(rmf, wa->quan*sizeof(rmf_t)); /* normalize */
	wa->wpla= realloc(wa->wpla, wa->numl*sizeof(size_t));

	*m= wa->numl;
	int k=wa->wpla[0];
	for(i=1;i<wa->numl;++i)
		if(k != wa->wpla[i])
			printf("Warning: Numcols is not uniform at %i words per line on all lines. This file has one with %zu.\n", k, wa->wpla[i]); 
	*n= k; 
	free_wseq(wa);

	return rmf;
}

gf22_t *processgf22(char *fname, int *m, int *n) /* this is dummy nmae .. it's for the special gff format for ;yze annottion */
{
	/* In order to make no assumptions, the file is treated as lines containing the same amount of words each,
	 * except for lines starting with #, which are ignored (i.e. comments). These words are checked to make sure they contain only floating number-type
	 * characters [0123456789+-.] only, one string variable is continually written over and copied into a growing floating point array each time */

	/* declarations */
	FILE *fp=fopen(fname,"r");
	int i;
	size_t couc /*count chars per line */, couw=0 /* count words */, oldcouw = 0;
	int c;
	int idanomals=0; /* this is for testing col 8 ... ID and string often the same, lie to keep it that way too. This will count when not */
	boole inword=0;
	wseq_t *wa=create_wseq_t(GBUF);
	size_t bwbuf=WBUF;
	char *bufword=calloc(bwbuf, sizeof(char)); /* this is the string we'll keep overwriting. */
	char *sm=NULL /*first semicolon marker*/, *fem=NULL /* first equals marker */, *lem=NULL /* last equals marker */;
	int cncols /* canonical numcols ... we will use the first line */;

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
					gf22[wa->numl].isz=12L;
					gf22[wa->numl].i=malloc(gf22[wa->numl].isz*sizeof(char));
					strcpy(gf22[wa->numl].i, "ABSENTIDCOL"); //* absent ID column */
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

dpf_t *processdpf(char *fname, int *m, int *n) /*fourth column is string, other columns to be ignored */
{
	/* In order to make no assumptions, the file is treated as lines containing the same amount of words each,
	 * except for lines starting with #, which are ignored (i.e. comments). These words are checked to make sure they contain only floating number-type
	 * characters [0123456789+-.] only, one string variable is continually written over and copied into a growing floating point array each time */

	/* declarations */
	FILE *fp=fopen(fname,"r");
	int i;
	size_t couc /*count chars per line */, couw=0 /* count words */, oldcouw = 0;
	int c;
	boole inword=0;
	wseq_t *wa=create_wseq_t(GBUF);
	size_t bwbuf=WBUF;
	char *bufword=calloc(bwbuf, sizeof(char)); /* this is the string we'll keep overwriting. */

	dpf_t *dpf=malloc(GBUF*sizeof(dpf_t));

	while( (c=fgetc(fp)) != EOF) { /* grab a char */
		if( (c== '\n') | (c == ' ') | (c == '\t') | (c=='#')) { /* word closing events */
			if( inword==1) { /* first word closing event */
				wa->wln[couw]=couc;
				bufword[couc++]='\0';
				bufword = realloc(bufword, couc*sizeof(char)); /* normalize */
				/* for the struct, we want to know if it's the first word in a line, like so: */
				if(couw==oldcouw) {
					dpf[wa->numl].n=malloc(couc*sizeof(char));
					dpf[wa->numl].nsz=couc;
					strcpy(dpf[wa->numl].n, bufword);
				} else if((couw-oldcouw)==1) /* it's not the first word, and it's 1st and second col */
					dpf[wa->numl].p=atol(bufword)-1L; /* note we convert to o indexed */
				else if((couw-oldcouw)==2) /* it's not the first word, and it's 1st and second col */
					dpf[wa->numl].d=atoi(bufword);
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
					dpf=realloc(dpf, wa->lbuf*sizeof(bgr_t));
					memset(wa->wpla+(wa->lbuf-WBUF), 0, WBUF*sizeof(size_t));
				}
				wa->wpla[wa->numl] = couw-oldcouw; /* number of words in current line */
				oldcouw=couw; /* restart words per line count */
				wa->numl++; /* brand new line coming up */
			}
			inword=0;
		} else if(inword==0) { /* deal with first character of new word, + and - also allowed */
			if(couw == wa->wsbuf-1) {
				wa->wsbuf += GBUF;
				wa->wln=realloc(wa->wln, wa->wsbuf*sizeof(size_t));
				dpf=realloc(dpf, wa->wsbuf*sizeof(dpf_t));
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
	dpf = realloc(dpf, wa->quan*sizeof(dpf_t)); /* normalize */
	wa->wpla= realloc(wa->wpla, wa->numl*sizeof(size_t));

	*m= wa->numl;
	int k=wa->wpla[0];
	for(i=1;i<wa->numl;++i)
		if(k != wa->wpla[i])
			printf("Warning: Numcols is not uniform at %i words per line on all lines. This file has one with %zu.\n", k, wa->wpla[i]); 
	*n= k; 
	free_wseq(wa);

	return dpf;
}

gf_t *processgf(char *fname, int *m, int *n) /* read in a genome file */
{
	/* In order to make no assumptions, the file is treated as lines containing the same amount of words each,
	 * except for lines starting with #, which are ignored (i.e. comments). These words are checked to make sure they contain only floating number-type
	 * characters [0123456789+-.] only, one string variable is continually written over and copied into a growing floating point array each time */

	/* declarations */
	FILE *fp=fopen(fname,"r");
	int i;
	size_t couc /*count chars per line */, couw=0 /* count words */, oldcouw = 0;
	int c;
	boole inword=0;
	wseq_t *wa=create_wseq_t(GBUF);
	size_t bwbuf=WBUF;
	char *bufword=calloc(bwbuf, sizeof(char)); /* this is the string we'll keep overwriting. */

	gf_t *gf=malloc(GBUF*sizeof(gf_t));

	while( (c=fgetc(fp)) != EOF) { /* grab a char */
		if( (c== '\n') | (c == ' ') | (c == '\t') | (c=='#')) { /* word closing events */
			if( inword==1) { /* first word closing event */
				wa->wln[couw]=couc;
				bufword[couc++]='\0';
				bufword = realloc(bufword, couc*sizeof(char)); /* normalize */
				/* for the struct, we want to know if it's the first word in a line, like so: */
				if(couw==oldcouw) {
					gf[wa->numl].n=malloc(couc*sizeof(char));
					gf[wa->numl].nsz=couc;
					strcpy(gf[wa->numl].n, bufword);
				} else if((couw-oldcouw)==1) /* it's not the first word, and it's 1st and second col */
					gf[wa->numl].z=atol(bufword);
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
					gf=realloc(gf, wa->lbuf*sizeof(bgr_t));
					memset(wa->wpla+(wa->lbuf-WBUF), 0, WBUF*sizeof(size_t));
				}
				wa->wpla[wa->numl] = couw-oldcouw; /* number of words in current line */
				oldcouw=couw; /* restart words per line count */
				wa->numl++; /* brand new line coming up */
			}
			inword=0;
		} else if(inword==0) { /* deal with first character of new word, + and - also allowed */
			if(couw == wa->wsbuf-1) {
				wa->wsbuf += GBUF;
				wa->wln=realloc(wa->wln, wa->wsbuf*sizeof(size_t));
				gf=realloc(gf, wa->wsbuf*sizeof(gf_t));
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
	gf = realloc(gf, wa->quan*sizeof(gf_t)); /* normalize */
	wa->wpla= realloc(wa->wpla, wa->numl*sizeof(size_t));

	*m= wa->numl;
	int k=wa->wpla[0];
	for(i=1;i<wa->numl;++i)
		if(k != wa->wpla[i])
			printf("Warning: Numcols is not uniform at %i words per line on all lines. This file has one with %zu.\n", k, wa->wpla[i]); 
	*n= k; 
	free_wseq(wa);

	return gf;
}

void prtbd2ia(bgr_t2 *bed2, int n, ia_t *ia)
{
	int i, j;
	for(i=0;i<ia->z;++i) {
		for(j=0;j<n;++j) {
			if(j==0)
				printf("%s ", bed2[ia->a[i]].n);
			else if(j==3)
				printf("%s ", bed2[ia->a[i]].f);
			else
				printf("%li ", bed2[ia->a[i]].c[j-1]);
			}
			printf("\n"); 
	}
	return;
}

void prtgf22(char *fname, gf22_t *gf22, int m7)
{
	int i;
	for(i=0;i<m7;++i) // note how we cut out the spurious parts of the motif string to leave it pure and raw (slightly weird why two-char deletion is necessary.
		printf("%s\t%li\t%li\t%c\t%s\t%s\n", gf22[i].n, gf22[i].c[0], gf22[i].c[1], gf22[i].sd, gf22[i].t, gf22[i].i);

#ifdef DBG
	printf("You have just seen the %i entries of basic gff2 file called \"%s\".\n", m7, fname); 
#endif
	return;
}

void prtrmf(char *fname, rmf_t *rmf, int m6)
{
	int i;
	for(i=0;i<m6;++i) // note how we cut out the spurious parts of the motif string to leave it pure and raw (slightly weird why two-char deletion is necessary.
		printf("%s\t%li\t%li\t%c\t%.*s\n", rmf[i].n, rmf[i].c[0], rmf[i].c[1], rmf[i].sd, (int)(rmf[i].msz-9), rmf[i].m+7);

	printf("You have just seen the %i entries of repeatmasker gff2 file called \"%s\".\n", m6, fname); 
	return;
}

void bed2in2(char *bed2fn, bgr_t2 *bed2, int m, int n, ia_t *ia) // split into 2 files
{
	int i, j, k=0;
	size_t lfn=strlen(bed2fn);
	char *outfn1=calloc(4+lfn ,sizeof(char));
	char *outfn2=calloc(4+lfn, sizeof(char));
	int rootsz=(int)(strchr(bed2fn, '.')-bed2fn);
	sprintf(outfn1, "%.*s_p1.bed", rootsz, bed2fn);
	sprintf(outfn2, "%.*s_p2.bed", rootsz, bed2fn);
	FILE *of1=fopen(outfn1, "w");
	FILE *of2=fopen(outfn2, "w");
	printf("bgr_t is %i rows by %i columns and is as follows:\n", m, n); 
	for(i=0;i<m;++i) {
		if(i==ia->a[k]){
			for(j=0;j<n;++j) {
				if(j==0)
					fprintf(of2, "%s\t", bed2[i].n);
				else if(j==3)
					fprintf(of2, "%s\n", bed2[i].f);
				else
					fprintf(of2, "%li\t", bed2[i].c[j-1]);
			}
			k++;
		} else {
			for(j=0;j<n;++j) {
				if(j==0)
					fprintf(of1, "%s\t", bed2[i].n);
				else if(j==3)
					fprintf(of1, "%s\n", bed2[i].f);
				else
					fprintf(of1, "%li\t", bed2[i].c[j-1]);
			}
		}
	}
	fclose(of1);
	fclose(of2);
	free(outfn1);
	free(outfn2);
	return;
}

void prtobed(bgr_t *bgrow, int m, int n, float minsig) // print over bed ... a value that is over a certain signal
{
	int i, j;
	printf("bgr_t is %i rows by %i columns and is as follows:\n", m, n); 
	for(i=0;i<m;++i) {
		if(bgrow[i].co >= minsig) {
			for(j=0;j<n;++j) {
				if(j==0)
					printf("%s ", bgrow[i].n);
				else if(j==3)
					printf("%2.6f ", bgrow[i].co);
				else
					printf("%li ", bgrow[i].c[j-1]);
			}
			printf("\n"); 
		}
	}
	return;
}

int *hist_co(bgr_t *bgrow, int m, float mxco, float mnco, int numbuckets)
{
	int i, j;
	float step=(mxco-mnco)/(float)numbuckets;
	float *bucketlimarr=malloc((numbuckets-1)*sizeof(float));
	int *bucketarr=calloc(numbuckets, sizeof(int));
	bucketlimarr[0]=step+mnco;
	for(i=1;i<numbuckets-1;++i) 
		bucketlimarr[i]=bucketlimarr[i-1]+step;

	for(i=0;i<m;++i)
		if(bgrow[i].co>=bucketlimarr[numbuckets-2]) {
			bucketarr[numbuckets-1]++;
			continue;
		} else {
			for(j=0;j<numbuckets-1;++j)
				if(bgrow[i].co < bucketlimarr[j]) {
					bucketarr[j]++;
					break;
				}
		}
	free(bucketlimarr);
	return bucketarr;
}

void prthist(char *histname, int *bucketarr, int numbuckets, int m, float mxco, float mnco)
{
	int i;
	printf("%s value %d-bin hstgrm for: %-24.24s (totels=%04i):\n", histname, numbuckets, histname, m); 
	printf("minval=%4.6f<-", mnco); 
	for(i=0;i<numbuckets;++i) 
		printf("| %i ", bucketarr[i]);
	printf("|->maxval=%4.6f\n", mxco); 
	return;
}

void prtdets(bgr_t *bgrow, int m, int n, char *label)
{
	int i;
	printf("bgr_t is %i rows by %i columns and is as follows:\n", m, n); 
	for(i=0;i<m;++i) {
		if(n==3)
			printf("%s\t%li\t%li\n", bgrow[i].n, bgrow[i].c[0], bgrow[i].c[1]); 
		if(n==4)
			printf("%s\t%li\t%li\t%4.4f\n", bgrow[i].n, bgrow[i].c[0], bgrow[i].c[1], bgrow[i].co); 
	}
	return;
}

void prtdeth(bgr_t *bgrow, int m, int n, char *label) /* Print intensity bedgraph in histogram format */
{
	int i;
	float mxco=.0, mnco=10e20;
	printf("bgr_t is %i rows by %i columns and is as follows:\n", m, n); 
	for(i=0;i<m;++i) {
		if(bgrow[i].co > mxco)
			mxco=bgrow[i].co;
		if(bgrow[i].co < mnco)
			mnco = bgrow[i].co;
	}
	int *hco=hist_co(bgrow, m, mxco, mnco, NUMBUCKETS);
	prthist(label, hco, NUMBUCKETS, m, mxco, mnco);
	free(hco);
	return;
}

void prtdetg(char *fname, gf_t *gf, int m, int n, char *label)
{
	int i;
	printf("%s called \"%s\" is %i rows by %i columns and is as follows:\n", label, fname, m, n); 
	for(i=0;i<m;++i)
		printf("%s\t%li\n", gf[i].n, gf[i].z);

	return;
}

void prtbed2s(bgr_t2 *bed2, int m, int n, words_t *bedword, int m3, int n3, char *label)
{
	/* TODO what you want is a copy of the data structure */
	int i, j, k;
	boole foundifeat;
	printf("Separated feature file %s is %i rows by %i columns and is as follows:\n", label, m, n); 
	for(i=0;i<m;++i) {
		foundifeat=0;
		for(k=0;k<m3;++k) {
			if(!strcmp(bedword[k].n, bed2[i].f) ) {
				foundifeat=1;
				for(j=0;j<n;++j) {
					if(j==0)
						printf("%s ", bed2[i].n);
					else if(j==3)
						printf("%s ", bed2[i].f);
					else
						printf("%li ", bed2[i].c[j-1]);
				}
				printf("\n"); 
			}
			if(foundifeat)
				break;
		}
	}
	return;
}

ia_t *gensplbdx(bgr_t2 *bed2, int m, int n, words_t *bedword, int m3, int n3) /* generate split bed index */
{
	/* TODO what you want is a copy of the data structure:
	 * NOPE! what you want is an array of indices */
	int i, k;
	ia_t *ia=calloc(1, sizeof(ia_t));
	ia->b=GBUF;
	ia->a=calloc(ia->b, sizeof(int));
	boole foundifeat;
	for(i=0;i<m;++i) {
		foundifeat=0;
		for(k=0;k<m3;++k) {
			if(!strcmp(bedword[k].n, bed2[i].f) ) {
				foundifeat=1;
				CONDREALLOC(ia->z, ia->b, GBUF, ia->a, int);
				ia->a[ia->z]=i;
				ia->z++;
			}
			if(foundifeat)
				break;
		}
	}
	ia->a=realloc(ia->a, ia->z*sizeof(int)); /*normalize */
	return ia;
}

void prtbed2fo(char *fname, bgr_t2 *bgrow, int m, int n, char *label) /* print feature beds file features only */
{
	int i;
	printf("%s file called %s is %i rows by %i columns and has following features:\n", label, fname, m, n); 
	printf("You can direct these name into a file and then presient to this program again under the -u option,\n");
	printf("whereupon only those name will be looked at\n");
	for(i=0;i<m;++i)
		printf("%s\n", bgrow[i].f);

	return;
}

void prtmbed(bgr_t **bgra, i4_t *dca, int dcasz, int n) /* the 2D version */
{
	int i, j;
	for(i=0;i<dcasz;++i) {
		for(j=0;j<dca[i].sc;++j) { // we're cycling though all of them, though we're really only interested in the first and last.
			if(j==0) { 
				printf("%s ", bgra[i][j].n);
				printf("%li ", bgra[i][j].c[0]);
			}
			if(j==dca[i].sc-1) { // note this cannot be an else if, because if only one line j = 0 = dca[i]-1.
				printf("%li ", bgra[i][j].c[1]);
				printf("%2.6f ", dca[i].mc);
			}
		}
		printf("\n"); 
	}
	return;
}

void m2beds(bgr_t *bgrow, bgr_t2 *bed2, int m2, int m) /* match up 2 beds */
{
	/* TODO: there could be an issue with intensity l;ines that span the end of one region and the start of another
	 * Need to look into that. this will only introduce a small error though.
	 */
	int i, j;
	int reghits; /* hits for region: number of lines in bed1 which coincide with a region in bed2 */
	int cloci; /* as opposed to hit, catch the number of loci */
	int rangecov=0;
	double assoctval=0;
	int istarthere=0, catchingi=0;
	boole caught;
	for(j=0;j<m2;++j) {
		caught=0;
		reghits=0;
		cloci=0;
		assoctval=0;
		for(i=istarthere;i<m;++i) {
			if( !(strcmp(bgrow[i].n, bed2[j].n)) & (bgrow[i].c[0] >= bed2[j].c[0]) & (bgrow[i].c[1] <= bed2[j].c[1]) ) {
				reghits++;
				rangecov=bgrow[i].c[1] - bgrow[i].c[0]; // range covered by this hit
				cloci+=rangecov;
				assoctval+=rangecov * bgrow[i].co;
				catchingi=i;
				caught=1;
			} else if (caught) { // will catch first untruth after a series of truths.
				caught=2;
				break; // bed1 is ordered so we can forget about trying to match anymore.
			}
		}
		if(caught==2)
			istarthere=catchingi+1;
		printf("Bed2idx %i / name %s / size %li got %i hits from bed1 , being %i loci and total assoc (prob .intensty) val %4.2f\n", j, bed2[j].f, bed2[j].c[1]-bed2[j].c[0], reghits, cloci, assoctval);
		if(istarthere >= m)
			break;
	}
	return;
}

void mgf2bed(char *gfname, char *ffile, gf_t *gf, bgr_t2 *bed2, int m2, int m5) /* match gf to feature bed file */
{
	setlocale(LC_NUMERIC, "");
	int i, j;
	int reghits; /* hits for region: number of lines in bed1 which coincide with a region in bed2 */
	int *acov=calloc(m5, sizeof(int)); /* coverage of this chromosome in the bed file */
	int rangecov=0;
	int istarthere=0, catchingi=0;
	int strmatch;
	boole caught;
	printf("Coverage of \"%s\" (genome size file) by \"%s\" (feature bed file):\n", gfname, ffile); 
	for(j=0;j<m5;++j) {
		caught=0;
		reghits=0;
		for(i=istarthere;i<m2;++i) {
			strmatch=strcmp(gf[j].n, bed2[i].n);
			if( (!strmatch) & (gf[j].z > bed2[i].c[0]) & (gf[j].z >= bed2[i].c[1]) ) {
				reghits++;
				rangecov=bed2[i].c[1] - bed2[i].c[0]; // range covered by this hit
				acov[j] +=rangecov;
				catchingi=i;
				caught=1;
			} else if (caught) { // will catch first untruth after a series of truths.
				caught=2;
				break; // bed1 is ordered so we can forget about trying to match anymore.
			} else if( (!strmatch) & (gf[j].z <= bed2[i].c[0]) & (gf[j].z < bed2[i].c[1]) ) {
				printf("There's a problem with the genome size file ... are you sure it's the right one? Bailing out.\n"); 
				exit(EXIT_FAILURE);
			}
		}
		if(caught==2)
			istarthere=catchingi+1;
		// printf("%s / cov %2.4f got %i hits from bed2\n", gf[j].n, (float)acov[j]/gf[j].z, reghits);
		printf("%s\t%4.2f%%\tof %'li bp\n", gf[j].n, 100.*(float)acov[j]/gf[j].z, gf[j].z);
		if(istarthere >= m2)
			break;
	}
	free(acov);
	return;
}

void mgf2dpf(char *gfname, char *dpffile, gf_t *gf, dpf_t *dpf, int m4, int m5) /* match gf to feature bed file */
{
	setlocale(LC_NUMERIC, "");
	int i, j;
	long *prdp=calloc(m5, sizeof(long)); /* position read depth */
	int istarthere=0, catchingi=0;
	int strmatch;
	boole caught;
	printf("Coverage of \"%s\" (genome size file) by \"%s\" (feature bed file):\n", gfname, dpffile); 
	for(j=0;j<m5;++j) {
		caught=0;
		for(i=istarthere;i<m4;++i) {
			strmatch=strcmp(gf[j].n, dpf[i].n);
			if( (!strmatch) & (gf[j].z > dpf[i].p)) {
				prdp[j] += dpf[i].d;
#ifdef DBG
				printf("prdp: %'li\n", prdp[j]); 
#endif
				catchingi=i;
				caught=1;
			} else if (caught) { // will catch first untruth after a series of truths.
				caught=2;
				break; // bed1 is ordered so we can forget about trying to match anymore.
			} else if( (!strmatch) & (gf[j].z <= dpf[i].p)) {
				printf("There's a problem with the genome size file ... are you sure it's the right one? Bailing out.\n"); 
				exit(EXIT_FAILURE);
			}
		}
		if(caught==2)
			istarthere=catchingi+1;
		printf("%s:\t%'li totreads\t%4.2f avgnumreads\twithin %'li bp\n", gf[j].n, prdp[j], (float)prdp[j]/gf[j].z, gf[j].z);
		if(istarthere >= m4)
			break;
	}
	free(prdp);
	return;
}

void mgf2rmf(char *gfname, char *rmffile, gf_t *gf, rmf_t *rmf, int m6, int m5) /* match gf to feature bed file */
{
	setlocale(LC_NUMERIC, "");
	int i, j;
	int reghits; /* hits for region: number of lines in bed1 which coincide with a region in rmf */
	int *acov=calloc(m5, sizeof(int)); /* coverage of this chromosome in the bed file */
	int rangecov=0;
	int istarthere=0, catchingi=0;
	int strmatch;
	boole caught;
	printf("Coverage of \"%s\" (genome size file) by \"%s\" (feature bed file):\n", gfname, rmffile); 
	for(j=0;j<m5;++j) {
		caught=0;
		reghits=0;
		for(i=istarthere;i<m6;++i) {
			strmatch=strcmp(gf[j].n, rmf[i].n);
			if( (!strmatch) & (gf[j].z > rmf[i].c[0]) & (gf[j].z >= rmf[i].c[1]) ) {
				reghits++;
				rangecov=rmf[i].c[1] - rmf[i].c[0]; // range covered by this hit
				acov[j] +=rangecov;
				catchingi=i;
				caught=1;
			} else if (caught) { // will catch first untruth after a series of truths.
				caught=2;
				break; // ordered so we can forget about trying to match anymore.
			} else if( (!strmatch) & (gf[j].z <= rmf[i].c[0]) & (gf[j].z < rmf[i].c[1]) ) {
				printf("There's a problem with the genome size file ... are you sure it's the right one? Bailing out.\n"); 
				exit(EXIT_FAILURE);
			}
		}
		if(caught==2)
			istarthere=catchingi+1;
		// printf("%s / cov %2.4f got %i hits from rmf\n", gf[j].n, (float)acov[j]/gf[j].z, reghits);
		printf("%s\t%4.2f%%\tof %'li bp\n", gf[j].n, 100.*(float)acov[j]/gf[j].z, gf[j].z);
		if(istarthere >= m6)
			break;
	}
	free(acov);
	return;
}

void md2bedp(dpf_t *dpf, bgr_t2 *bed2, int m2, int m) /* match up a samtools depth file (-d option) and a feature bed file (-f option) and print */
{
	int i, j, min, max;
	int reghits; /* hits for region: number of lines in bed1 which coincide with a region in bed2 */
	int cloci; /* as opposed to hit, catch the number of loci */
	long assoctval=0;
	int istarthere=0, catchingi=0;
	boole caught;
	for(j=0;j<m2;++j) {
		caught=0;
		reghits=0;
		cloci=0;
		assoctval=0;
		min=9999999;
		max=0;
		for(i=istarthere;i<m;++i) {
			if( !(strcmp(dpf[i].n, bed2[j].n)) & (dpf[i].p >= bed2[j].c[0]) & (dpf[i].p < bed2[j].c[1]) ) {
				reghits++;
				cloci++;
				assoctval+= dpf[i].d;
				catchingi=i;
				if(dpf[i].d<min)
					min=dpf[i].d;
				if(dpf[i].d>max)
					max=dpf[i].d;
				caught=1;
			} else if (caught) { // will catch first untruth after a series of truths.
				caught=2;
				break; // bed1 is ordered so we can forget about trying to match anymore.
			}
		}
		if(caught==2)
			istarthere=catchingi+1;
		// printf("Bed2idx %i / name %s / size %li got %i hits from dpf , being %i loci and accumulated depth val of %lu\n", j, bed2[j].f, bed2[j].c[1]-bed2[j].c[0], reghits, cloci, assoctval);
		printf("%s\t%li\t%li\t%s\t%i\t%i\t%li\t%4.4f\n", bed2[j].n, bed2[j].c[0], bed2[j].c[1], bed2[j].f, min, max, assoctval, (float)assoctval/cloci);

		if(istarthere >= m)
			break;
	}
	return;
}

i4_t *difca(bgr_t *bgrow, int m, int *dcasz, float minsig) /* An temmpt to merge bgraph quickly, no hope */
{
	int i, goodi=0 /* the last i at which minsig was satisfied */;
	boole seenminsig=0;
	/* how many different chromosomes are there? the dc (different chromsosome array */
	int dcbf=GBUF, dci=0;
	i4_t *dca=calloc(dcbf, sizeof(i4_t));
	char *tstr=NULL;
	/* find first bgrow element which is over the minimum coverage */
	for(i=0;i<m;++i)
		if(bgrow[i].co >= minsig) {
			tstr=malloc(bgrow[i].nsz*sizeof(char)); /* tmp string */
			strcpy(tstr, bgrow[i].n);
			dca[dci].sc++;
			dca[dci].mc=bgrow[i].co;
			dca[dci].b1i=i;
			dca[dci].lgbi=i;
			seenminsig=1;
			goodi=i;
			break;
		}
	if(!seenminsig) {
		printf("Error. No bedgraph element was able to satisfy the minimum signal value that was specified: abandoning ship.\n");
		exit(EXIT_FAILURE);
	}

	for(i=goodi+1;i<m;++i) {
		/* the same now means same name and contiguous */
		if( (!strcmp(tstr, bgrow[i].n)) & (bgrow[i].c[0] == bgrow[dca[dci].lgbi].c[1]) & (bgrow[i].co >= minsig) ) {
			dca[dci].sc++;
			dca[dci].lgbi=i;
			if(bgrow[i].co<dca[dci].mc)
				dca[dci].mc=bgrow[i].co;
		} else if (bgrow[i].co >= minsig) {
			CONDREALLOC(dci, dcbf, GBUF, dca, i4_t);
			dci++;
			dca[dci].sc++;
			dca[dci].mc=bgrow[i].co;
			dca[dci].b1i=i;
			dca[dci].lgbi=i;
			/* new string could be differnt length*/
			tstr=realloc(tstr, bgrow[i].nsz*sizeof(char)); /* tmp string */
			strcpy(tstr, bgrow[i].n);
		}
	}
	dca=realloc(dca, (dci+1)*sizeof(i4_t));
#ifdef DBG
	printf("Num of different chromcontigs=%i. How many of each? Let's see:\n", dci+1); 
	printf("dcbf=%i. 4-tupe is sc/mc/b1i/lgbi\n", dcbf); 
	for(i=0;i<=dci;++i) 
		printf("%i/%4.2f/%i/%i ",dca[i].sc, dca[i].mc, dca[i].b1i, dca[i].lgbi); 
	printf("\n"); 
#endif
	*dcasz=dci+1;
	free(tstr);
	return dca;
}

void prtusage()
{
	printf("chipmonk: this takes a bedgraph file, specified by -i, probably the bedgraph from a MACS2 intensity signal,\n");
	printf("and another bedgraph file, specified by -f, and merges the first into lines defined by the second.\n");
	printf("Before filtering however, please run with the -d (details) option. This will showi a rough spread of the values,\n");
	printf("so you can run a second time choosing filtering value (-f) more easily.\n");
	return;
}

int main(int argc, char *argv[])
{
	/* argument accounting */
	if(argc == 1) {
		prtusage();
		exit(EXIT_FAILURE);
	}
	int i, m, n /*rows,cols for bgr_t*/, m2, n2, m3, n3, m4, n4, m4b, n4b, m5, n5, m6, n6, m7, n7 /*gf22 dims */;
    unsigned numsq;
	opt_t opts={0};
	catchopts(&opts, argc, argv);

	/* Read in files according to what's defined in options */
	bgr_t *bgrow=NULL; /* usually macs signal */
	bgr_t2 *bed2=NULL; /* usually bed file from gff */
	words_t *bedword=NULL; /* usually feature names of interest */
	dpf_t *dpf=NULL; /* usually feature names of interest */
	dpf_t *dpf2=NULL; /* usually feature names of interest */
	gf_t *gf=NULL; /* usually genome size file */
	rmf_t *rmf=NULL; /* usually genome size file */
	gf22_t *gf22=NULL;
	i_s *sqisz=NULL;

	/* for the ystr, gf22 hash handling */
    unsigned htsz;
    gf22snod **stab=NULL;

	if(opts.istr)
		bgrow=processinpf(opts.istr, &m, &n);
	if(opts.fstr)
		bed2=processinpf2(opts.fstr, &m2, &n2);
	if(opts.ustr)
		bedword=processwordf(opts.ustr, &m3, &n3);
	if(opts.pstr)
		dpf=processdpf(opts.pstr, &m4, &n4);
	if(opts.qstr)
		dpf2=processdpf(opts.pstr, &m4b, &n4b);
	if(opts.gstr)
		gf=processgf(opts.gstr, &m5, &n5);
	if(opts.rstr)
		rmf=processrmf(opts.rstr, &m6, &n6);
	if(opts.ystr) { // we're goign to try hashing the ID line of gf22 format */
		gf22=processgf22(opts.ystr, &m7, &n7);
		htsz=2*m7/3; /* our hash table size */
		stab=gf22tochainharr(gf22, m7, htsz); /* we now set up a hash table along side our sequence names from the fasta file */
	}
	if(opts.astr)
		sqisz=procfa(opts.astr, &numsq);

	/* conditional execution of certain functions depending on the options */
	if((opts.dflg) && (opts.istr)) {
		prtdets(bgrow, m, n, "Target bedgraph (1st) file");
		goto final;
	}
	if((opts.dflg) && (opts.gstr)) {
		prtdetg(opts.gstr, gf, m5, n5, "Size file");
		goto final;
	}
	if((opts.nflg) && (opts.fstr)) {
		prtbed2fo(opts.fstr, bed2, m2, n2, "Feature (bed2)");
		goto final;
	}
	// prtbed2(bed2, m2, MXCOL2VIEW);
	if((opts.istr) && (opts.fstr))
		m2beds(bgrow, bed2, m2, m);

	if((opts.ustr) && (opts.fstr) && (!opts.sflg)) {
		printf("bedwords:\n"); 
		for(i=0;i<m3;++i)
			printf("%s\n", bedword[i].n);
	}

	if((opts.pstr) && (opts.fstr) )
		md2bedp(dpf, bed2, m2, m4);

	if((opts.dflg) && (opts.rstr) )
		prtrmf(opts.rstr, rmf, m6);

	if((opts.dflg) && (opts.ystr) ) {
		prtgf22(opts.ystr, gf22, m7);
		prtgf22chaharr(stab, htsz);
	}

	if((opts.dflg) && (opts.astr) ) 
		prtsq(sqisz, numsq);

	/* bedgraph and fasta file */
	if((opts.istr) && (opts.astr) )
		prtsqbdg(sqisz, bgrow, m, numsq);

	if((opts.gstr) && (opts.rstr) )
		mgf2rmf(opts.gstr, opts.rstr, gf, rmf, m6, m5);

	if((opts.gstr) && (opts.pstr) )
		mgf2dpf(opts.gstr, opts.pstr, gf, dpf, m4, m5);

	if((opts.gstr) && (opts.fstr) )
		mgf2bed(opts.gstr, opts.fstr, gf, bed2, m2, m5);
	// if((opts.ustr) && (opts.fstr) && opts.sflg)
	// 	prtbed2s(bed2, m2, MXCOL2VIEW, bedword, m3, n3, "bed2 features that are in interesting-feature-file");

	ia_t *ia=NULL;
	if((opts.ustr) && (opts.fstr) && opts.sflg) {
		ia=gensplbdx(bed2, m2, n2, bedword, m3, n3);
		bed2in2(opts.fstr, bed2, m2, n2, ia);
	}

final:
	if(opts.pstr) {
		for(i=0;i<m4;++i)
			free(dpf[i].n);
		free(dpf);
	}
	if(opts.qstr) {
		for(i=0;i<m4b;++i)
			free(dpf2[i].n);
		free(dpf2);
	}
	if(opts.istr) {
		for(i=0;i<m;++i)
			free(bgrow[i].n);
		free(bgrow);
	}
	if(opts.fstr) {
		for(i=0;i<m2;++i) {
			free(bed2[i].n);
			free(bed2[i].f);
		}
		free(bed2);
	}
	if(opts.rstr) {
		for(i=0;i<m6;++i) {
			free(rmf[i].n);
			free(rmf[i].m);
		}
		free(rmf);
	}
	if(opts.gstr) {
		for(i=0;i<m5;++i)
			free(gf[i].n);
		free(gf);
	}
	if(opts.ustr) {
		for(i=0;i<m3;++i)
			free(bedword[i].n);
		free(bedword);
	}
	if(opts.ystr) {
		for(i=0;i<m7;++i) {
			free(gf22[i].n);
			free(gf22[i].t);
			free(gf22[i].i);
		}
		free(gf22);
	}
	if(opts.astr) {
		freegf22chainharr(stab, htsz);
		for(i=0;i<numsq;++i) {
			free(sqisz[i].id);
			free(sqisz[i].sq);
		}
		free(sqisz);
	}

	return 0;
}
