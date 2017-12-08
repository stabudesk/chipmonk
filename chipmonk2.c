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

typedef enum { /* type category: CDS, gene or whitever */
	GNE=1, /* gene */
	CDS=2, /* cds */
	MRN=3, /* cds */
	ARS=4, /* cds */
	TEL=5, /* cds */
	INR=6, /* cds */
	NCE=7, /* uncat */
	TRG=8, /* uncat */
	NRG=9, /* uncat */
	RRG=10, /* uncat */
	ORR=11, /* uncat */
	TLR=12,
	XEL=13,
	XEC=14,
	CE1=15,
	CE2=16,
	CE3=17,
	ACS=18,
	SRG=19,
	P1T=20,
	LTR=21,
	TEG=22,
	LT0=23,
	YPE=24,
	RGN=25,
	MAS=26,
	PSU=27,
	UI5=28,
	CE0=29,
	ETS=30,
	ITS=31,
	NTR=32,
	SNR=33,
	BRF=34,
	TR0=35,
	SMT=36,
	WRG=37,
	XRG=38,
	YRG=39,
	Z1R=40,
	Z2R=41,
	MTR=42,
	IER=43,
	UCT=44 /* uncat */
} tcat;
#define TCQUAN 44
char *tcnames[TCQUAN]={"gene", "CDS", "mRNA", "ARS", "telomere", "intron", "noncoding_exon", "tRNA_gene", "ncRNA_gene", "rRNA_gene", "origin_of_replication", "telomeric_repeat", "X_element", "X_element_combinatorial_repeat", "centromere_DNA_Element_I", "centromere_DNA_Element_II", "centromere_DNA_Element_III", "ARS_consensus_sequence", "snoRNA_gene", "plus_1_translational_frameshift", "LTR_retrotransposon", "transposable_element_gene", "long_terminal_repeat",  "Y_prime_element", "region", "matrix_attachment_site", "pseudogene", "five_prime_UTR_intron", "centromere", "external_transcribed_spacer_region", "internal_transcribed_spacer_region", "non_transcribed_region", "snRNA_gene", "blocked_reading_frame", "telomerase_RNA_gene", "silent_mating_type_cassette_array", "W_region", "X_region", "Y_region", "Z1_region", "Z2_region", "mating_type_region", "intein_encoding_region", "uncat"};
// const char *importantf="gene";

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
	char *lstr; /* YG list with IDs */
	char *pstr; /* depth file name */
	char *qstr; /* 2nddepth file name */
	char *gstr; /* genome file name */
	char *rstr; /* repeatmasker gtf/gff2 file */
	char *ystr; /* the gf22_t type */
	char *hstr; /* the gf23_t type */
	char *astr; /* a fasta file */
	char *bstr; /* blast output format 7 */
	char *xstr; /* in gfmatchup this was zstr ... it a special add on to gf22 ... it's the fcat enum so you can select types inside gf22 */
	char *zstr; /* the tcname of interest */
} opt_t;

typedef struct /* i4_t */
{
	int sc; /* number of same chromosome (occurences?) */
	float mc; /* min signal value */
	int b1i; /* index of the 1st bgr_t, which satisfies the conditions */
	int lgbi; /* last good bgr_t index */
} i4_t; /* 4 vals of some sort? */

typedef struct /* bgr_t: bedgraph file chrom|start|exend|float */
{
	char *n;
	size_t nsz; /* size of the name r ID field */
	long c[2]; /* coords: 1) start 2) end */
	float co; /* signal value */
} bgr_t; /* bedgraph row type */

typedef struct /* ygl_t */
{
	char *s /* S-id */, *y /* Y-id */, *g /* Gene-id, if any */, *d /* Desc, if any */;
	size_t ssz, ysz, gsz, dsz;
} ygl_t; /* Listing from Yeast Genome */

typedef struct /* bgr_t2 */
{
	char *n;
	size_t nsz; /* size of the name r ID field */
	long c[2]; /* coords: 1) start 2) end */
	char *f, *s, *t; /* f for feature .. 4th col */
	size_t fsz, ssz, tsz; /* size of the feature field*/
	char sd; /* strand + or - */
	tcat tc;
} bgr_t2; /* bedgraph row type 2i. column is the feature */

typedef struct /* blop_t: the b option, bast output format 7 */
{
	char *n;
	size_t nsz; /* size of the name r ID field */
	char *tc;
	size_t tcsz; /* target string ..yep */
	char *qfs; /* query Feature substring */
	char *fs; /* target feature substring */
	size_t fssz; /* size of fs */
	size_t qfssz; /* size of fs */
	int al; /* alignment length */
	float eval;
	float pcti;
	long c[2]; /* coords on the target : 1) start 2) cols 8 and 9 */
	unsigned qfnoc; /* query feature name number of occurrences */
} blop_t;

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
	fcat fc;
	unsigned char blofound; // relates to specific -b -y operation ... did this get a hit? if so skip othe hits; zero if not.
} gf22_t;

typedef struct /* gf23_t: the h option, based on rmf_t */
{
	char *n;
	size_t nsz; /* size of the name r ID field */
	long c[2]; /* coords: 1) start 2) cols 4 and 5 */
	char *t; /* the typestring ... 3rd column */
	char *i; /* the ID string ... last column */
	char sd; /* strand + or - */
	size_t isz; /* size of iD string */
	size_t tsz; /* size of type stringstring */
} gf23_t;

typedef struct /* words_t: file with only single words per line */
{
	char *n;
	size_t nsz; /* size of the name r ID field */
} words_t;

typedef struct /* dpf_t : depth file type ... just chr name, pos and read quant: so, single base-locus positions, yes? */
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

struct strchainode2 /* bt2snod struct ..e for the -f option */
{
    bgr_t2 *bed2; /* ptr to a single element, not an array, TODO void* it .. oh yeah sure. */
    struct strchainode2 *n;
};
typedef struct strchainode2 bt2snod; /* yes, leave this alone, it's the way a struct can have a ptr ot its own type! */

struct strchainode3 /* gf23snod struct */
{
    gf23_t *gf23; /* ptr to a single element, not an array, TODO void* it. */
    struct strchainode3 *n;
};
typedef struct strchainode3 gf23snod; /* yes, leave this alone, it's the way a struct can have a ptr ot its own type! */

struct strchainodeyg /* ygsnod struct */
{
    ygl_t *yglst; /* ptr to a single element, not an array, TODO void* it. */
    struct strchainodeyg *n;
};
typedef struct strchainodeyg ygsnod; /* yes, leave this alone, it's the way a struct can have a ptr ot its own type! */

struct strchainodeblo /* blosnod struct */
{
    blop_t *blop; /* ptr to a single element, not an array, TODO void* it. */
    struct strchainodeblo *n;
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
			case 's':
				opts->sflg = 1;
				break;
			case 'n':
				opts->nflg = 1;
				break;
			case 'l':
				opts->lstr = optarg;
				break;
			case 'b':
				opts->bstr = optarg;
				break;
			case 'x':
				opts->xstr = optarg;
				break;
			case 'z':
				opts->zstr = optarg;
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
			case 'h': /* w303sgd: based on repeatmasker gff2 file */
				opts->hstr = optarg;
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

bt2snod **bt2tochainharr(bgr_t2 *bed2, unsigned numsq, unsigned tsz)
{
    unsigned i;

    bt2snod **stab=malloc(tsz*sizeof(bt2snod *));
    for(i=0;i<tsz;++i) 
        stab[i]=NULL; /* _is_ a valid ptr, but it's unallocated. Initialization is possible though. */
    bt2snod *tsnod0, *tsnod2;

    unsigned tint;
    for(i=0; i<numsq; ++i) {
        tint=hashit(bed2[i].f, tsz);
        if( (stab[tint] == NULL) ) {
            stab[tint]=malloc(sizeof(bt2snod));
            stab[tint]->bed2=bed2+i;
            stab[tint]->n=NULL;
            continue;
        }
        tsnod2=stab[tint];
        while( (tsnod2 != NULL) ){
            if(!strcmp(tsnod2->bed2->f, bed2[i].f)) {
                goto nxt;
            }
            tsnod0=tsnod2;
            tsnod2=tsnod2->n;
        }
        tsnod0->n=malloc(sizeof(bt2snod));
        tsnod0->n->bed2=bed2+i;
        tsnod0->n->n=NULL;
nxt:        continue;
    }
    return stab;
}

gf23snod **gf23tochainharr(gf23_t *gf23, unsigned numsq, unsigned tsz)
{
    unsigned i;

    gf23snod **stab=malloc(tsz*sizeof(gf23snod *));
    for(i=0;i<tsz;++i) 
        stab[i]=NULL; /* _is_ a valid ptr, but it's unallocated. Initialization is possible though. */
    gf23snod *tsnod0, *tsnod2;

    unsigned tint;
    for(i=0; i<numsq; ++i) {
        tint=hashit(gf23[i].i, tsz);
        if( (stab[tint] == NULL) ) {
            stab[tint]=malloc(sizeof(gf23snod));
            stab[tint]->gf23=gf23+i;
            stab[tint]->n=NULL;
            continue;
        }
        tsnod2=stab[tint];
        while( (tsnod2 != NULL) ){
            if(!strcmp(tsnod2->gf23->i, gf23[i].i)) {
                goto nxt;
            }
            tsnod0=tsnod2;
            tsnod2=tsnod2->n;
        }
        tsnod0->n=malloc(sizeof(gf23snod));
        tsnod0->n->gf23=gf23+i;
        tsnod0->n->n=NULL;
nxt:        continue;
    }
    return stab;
}

ygsnod **ygtochainharr(ygl_t *yglst, unsigned numsq, unsigned tsz)
{
    unsigned i;

    ygsnod **stab=malloc(tsz*sizeof(ygsnod *));
    for(i=0;i<tsz;++i) 
        stab[i]=NULL; /* _is_ a valid ptr, but it's unallocated. Initialization is possible though. */
    ygsnod *tsnod0, *tsnod2;

    unsigned tint;
    for(i=0; i<numsq; ++i) {
        tint=hashit(yglst[i].y, tsz);
        if( (stab[tint] == NULL) ) {
            stab[tint]=malloc(sizeof(ygsnod));
            stab[tint]->yglst=yglst+i;
            stab[tint]->n=NULL;
            continue;
        }
        tsnod2=stab[tint];
        while( (tsnod2 != NULL) ){
            if(!strcmp(tsnod2->yglst->y, yglst[i].y)) {
                goto nxt;
            }
            tsnod0=tsnod2;
            tsnod2=tsnod2->n;
        }
        tsnod0->n=malloc(sizeof(ygsnod));
        tsnod0->n->yglst=yglst+i;
        tsnod0->n->n=NULL;
nxt:        continue;
    }
    return stab;
}

blosnod **blotochainharr(blop_t *blop, unsigned numsq, unsigned tsz)
{
    unsigned i;

    blosnod **stab=malloc(tsz*sizeof(blosnod *));
    for(i=0;i<tsz;++i) 
        stab[i]=NULL; /* _is_ a valid ptr, but it's unallocated. Initialization is possible though. */
    blosnod *tsnod0, *tsnod2;

    unsigned tint;
    for(i=0; i<numsq; ++i) {
        tint=hashit(blop[i].qfs, tsz);
        if( (stab[tint] == NULL) ) {
            stab[tint]=malloc(sizeof(blosnod));
            stab[tint]->blop=blop+i;
            stab[tint]->n=NULL;
            continue;
        }
        tsnod2=stab[tint];
        while( (tsnod2 != NULL) ){
            if(!strcmp(tsnod2->blop->qfs, blop[i].qfs)) {
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

void prtbt2chaharr(bt2snod **stab, unsigned tsz)
{
    unsigned i;
    bt2snod *tsnod2;
    for(i=0;i<tsz;++i) {
        printf("Tablepos %i: ", i); 
        tsnod2=stab[i];
        while(tsnod2) {
            printf("\"%s\" ", tsnod2->bed2->f); 
            tsnod2=tsnod2->n;
        }
        printf("\n"); 
    }
    return;
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

void prtgf23chaharr(gf23snod **stab, unsigned tsz)
{
    unsigned i;
    gf23snod *tsnod2;
    for(i=0;i<tsz;++i) {
        printf("Tablepos %i: ", i); 
        tsnod2=stab[i];
        while(tsnod2) {
            printf("\"%s\" ", tsnod2->gf23->i); 
            tsnod2=tsnod2->n;
        }
        printf("\n"); 
    }
    return;
}

void prtygchaharr(ygsnod **stab, unsigned tsz)
{
    unsigned i;
    ygsnod *tsnod2;
    for(i=0;i<tsz;++i) {
        printf("Tablepos %i: ", i); 
        tsnod2=stab[i];
        while(tsnod2) {
            printf("\"%s\" ", tsnod2->yglst->y); 
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

void freegf23chainharr(gf23snod **stab, size_t tsz)
{
    int i;
    gf23snod *tsnod0, *tsnod2;
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

void freebt2chainharr(bt2snod **stab, size_t tsz)
{
    int i;
    bt2snod *tsnod0, *tsnod2;
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

void freeygchainharr(ygsnod **stab, size_t tsz)
{
    int i;
    ygsnod *tsnod0, *tsnod2;
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

char idlsearch3(gf23snod **stab, unsigned tsz, char *line, unsigned lnsz)
{
    char yes=0;
    size_t i;
    gf23snod *tsnod2;

    unsigned tint;
    for(i=0; i<tsz; ++i) {
        tint=hashit(line, tsz);
        if( (stab[tint] == NULL) )
            goto nxt; /* hashtable at that position is empty ... so it's impossible for that string to be there */

        tsnod2=stab[tint];
        while( (tsnod2 != NULL) ) {
            if( (tsnod2->gf23->isz == lnsz) & !(strncmp(tsnod2->gf23->i, line, tsnod2->gf23->isz)) ) {
                yes=1; /* i.e. no=0 */
                goto nxt;
            }
            tsnod2=tsnod2->n;
        }
    }
nxt:    return yes;
}

void prtblop(char *fname, blop_t *blop, int mb)
{
	int i;
	for(i=0;i<mb;++i) // note how we cut out the spurious parts of the motif string to leave it pure and raw (slightly weird why two-char deletion is necessary.
		printf("%i\t%s\t%i\t%s\tqfs:%s\ttfs:%s\t%li\t%li\t%4.2f\t%4.2f\n", blop[i].qfnoc, blop[i].n, blop[i].al, blop[i].tc, blop[i].qfs, blop[i].fs, blop[i].c[0], blop[i].c[1], blop[i].pcti, blop[i].eval);

	printf("You have just seen the %i entries of blast output file called \"%s\".\n", mb, fname); 
	return;
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
	printf("Number of different sequences=%i\n", sz); 
#ifdef DBG
	int i;
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
	char rangestr[128]={0}; /* generally helpful to say what range is being given */
	for(j=0;j<m;++j) {
		for(i=0;i<sz;++i) {
			if(!strcmp(sqisz[i].id, bgrow[j].n)) {
				sprintf(rangestr, "|range_%li_to_%li", bgrow[j].c[0], bgrow[j].c[1]);
				printf(">%s", sqisz[i].id);
				printf("%s\n", rangestr);
				printf("%.*s\n", (int)(bgrow[j].c[1]-bgrow[j].c[0]), sqisz[i].sq+bgrow[j].c[0]);
				break;
			}
		}
	}
	return;
}

void prtsqbdgf(i_s *sqisz, bgr_t2 *bgrow, int m, int sz)
{
	int i, j;
	char rangestr[128]={0}; /* generally helpful to say what range is being given */
	for(j=0;j<m;++j) {
		for(i=0;i<sz;++i) {
			if(!strcmp(sqisz[i].id, bgrow[j].n)) {
				sprintf(rangestr, "|range_%li_to_%li|%s", bgrow[j].c[0], bgrow[j].c[1], bgrow[j].f);
				printf(">%s", sqisz[i].id);
				printf("%s\n", rangestr);
				printf("%.*s\n", (int)(bgrow[j].c[1]-bgrow[j].c[0]), sqisz[i].sq+bgrow[j].c[0]);
				break;
			}
		}
	}
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

void prtmbed2wof(i_s *sqisz, int sz, bgr_t2 *bgrow, int m2, bt2snod **stab2, unsigned tsz, words_t *bedword, int m3) /* search beds2 feature names from file of words */
{
	unsigned char found;
	int numfound2=0;
	int i, j, k=0;
	long kk;
	char rangestr[128]={0}; /* generally helpful to say what range is being given */
    bt2snod *tsnodd0=NULL, *tsnodd2=NULL;
    unsigned tint;
	unsigned gbuf=GBUF;
	int *nfiarr=malloc(gbuf*sizeof(int)); /* the not found index array */
	for(i=0;i<m3;++i) {
		found = 0;
        tint=hashit(bedword[i].n, tsz);
        tsnodd2=stab2[tint];
        while( (tsnodd2 != NULL) ) {
            if(!strcmp(tsnodd2->bed2->f, bedword[i].n)) {
				for(j=0;j<sz;++j) { /* have to go through each chromosome */
					if(!strcmp(sqisz[j].id, tsnodd2->bed2->n)) {
						sprintf(rangestr, "|range_%li_to_%li|%s", tsnodd2->bed2->c[0], tsnodd2->bed2->c[1], tsnodd2->bed2->f);
						printf(">%s", sqisz[j].id);
						printf("%s\n", rangestr);
						if(tsnodd2->bed2->sd =='+')
							printf("%.*s\n", (int)(tsnodd2->bed2->c[1] - tsnodd2->bed2->c[0]), sqisz[j].sq + tsnodd2->bed2->c[0]);
						else if(tsnodd2->bed2->sd =='-') {
							kk=tsnodd2->bed2->c[1]-1;
							while(kk>=tsnodd2->bed2->c[0]) { switch (sqisz[j].sq[kk--]) { case 'A': putchar('T'); break; case 'C': putchar('G'); break; case 'G': putchar('C'); break; case 'T': putchar('A'); break; default: break; } }
							putchar('\n');
						}
						break;
					}
				}
				numfound2++;
				found = 1;
				break;
			}
			tsnodd0=tsnodd2;
            tsnodd2=tsnodd2->n;
		}
		if(!found) { // we're counting non-found items.
			CONDREALLOC(k, gbuf, GBUF, nfiarr, int);
			nfiarr[k++]=i;
		}
	}
	if(k) {
		nfiarr=realloc(nfiarr, k*sizeof(int)); // k can be zero
		fprintf(stderr, "%i names found. %i names not found.\n", numfound2, k);
	} else
		fprintf(stderr, "All supplied %i names were found.\n", numfound2);

	free(nfiarr);

	return;
}

void prtsqbdgf2(i_s *sqisz, bgr_t2 *bgrow, int m, int sz, tcat ztcat)
{
	int i, j;
	char rangestr[128]={0}; /* generally helpful to say what range is being given */
	for(j=0;j<m;++j) {
		for(i=0;i<sz;++i) {
			if(bgrow[j].tc != ztcat) // only this feature type will figure
				continue;
			if(!strcmp(sqisz[i].id, bgrow[j].n)) {
				sprintf(rangestr, "|range_%li_to_%li|%s", bgrow[j].c[0], bgrow[j].c[1], bgrow[j].f);
				printf(">%s", sqisz[i].id);
				printf("%s\n", rangestr);
				printf("%.*s\n", (int)(bgrow[j].c[1]-bgrow[j].c[0]), sqisz[i].sq+bgrow[j].c[0]);
				break;
			}
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

ygl_t *processygl(char *fname, int *m, int *n)
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

	ygl_t *yglst=malloc(GBUF*sizeof(ygl_t));

	while( (c=fgetc(fp)) != EOF) { /* grab a char */
		if( (c== '\n') | (c == '\t') | (c=='#')) { /* word closing events */
			if( inword==1) { /* first word closing event */
				wa->wln[couw]=couc;
				bufword[couc++]='\0';
				bufword = realloc(bufword, couc*sizeof(char)); /* normalize */
				/* for the struct, we want to know if it's the first word in a line, like so: */
				if(couw==oldcouw) {
					yglst[wa->numl].s=malloc(couc*sizeof(char));
					yglst[wa->numl].ssz=couc;
					strcpy(yglst[wa->numl].s, bufword);
				} else if((couw-oldcouw)== 1) {
					yglst[wa->numl].y=malloc(couc*sizeof(char));
					yglst[wa->numl].ysz=couc;
					strcpy(yglst[wa->numl].y, bufword);
				} else if((couw-oldcouw)== 3) {
					yglst[wa->numl].g=malloc(couc*sizeof(char));
					yglst[wa->numl].gsz=couc;
					strcpy(yglst[wa->numl].g, bufword);
				} else if((couw-oldcouw)== 4) {
					yglst[wa->numl].d=malloc(couc*sizeof(char));
					yglst[wa->numl].dsz=couc;
					strcpy(yglst[wa->numl].d, bufword);
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
					yglst=realloc(yglst, wa->lbuf*sizeof(ygl_t));
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
				yglst=realloc(yglst, wa->wsbuf*sizeof(ygl_t));
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
	yglst = realloc(yglst, wa->quan*sizeof(ygl_t)); /* normalize */
	wa->wpla= realloc(wa->wpla, wa->numl*sizeof(size_t));

	*m= wa->numl;
	int k=wa->wpla[0];
	for(i=1;i<wa->numl;++i)
		if(k != wa->wpla[i])
			printf("Warning: Numcols is not uniform at %i words per line on all lines. This file has one with %zu.\n", k, wa->wpla[i]); 
	*n= k; 
	free_wseq(wa);

	return yglst;
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
	boole inword=0, foundtype;
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
				} else if( (couw-oldcouw)==3) {
					bgrow[wa->numl].f=malloc(couc*sizeof(char));
					bgrow[wa->numl].fsz=couc;
					strcpy(bgrow[wa->numl].f, bufword);
				} else if((couw-oldcouw)==5) { /* it's not the first word, and it's 1st and second col */
					bgrow[wa->numl].sd=bufword[0];
				} else if( (couw-oldcouw)==6) {
					bgrow[wa->numl].s=malloc(couc*sizeof(char));
					bgrow[wa->numl].ssz=couc;
					strcpy(bgrow[wa->numl].s, bufword);
				} else if( (couw-oldcouw)==7) {
					bgrow[wa->numl].t=malloc(couc*sizeof(char));
					bgrow[wa->numl].tsz=couc;
					strcpy(bgrow[wa->numl].t, bufword);
					foundtype=0;
					for(i=0;i<TCQUAN-1;++i) 
						if(!strcmp(bufword, tcnames[i])) {
							bgrow[wa->numl].tc=i+1;
							foundtype=1;
							break;
						}
					if(!foundtype)
						bgrow[wa->numl].tc=UCT;
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
	char *tmpc=NULL, *tmpc2=NULL;
	char *oldqfs=calloc(16, sizeof(char));
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
					tmpc2=strchr(bufword, 'r'); // hacky ... avoiding strtok by grabbing r after 2nd |
					if( (tmpc == NULL) | (tmpc2==NULL) | (tmpc2 <= tmpc))
						printf("query name did not have a : or r or r appeared before : ... you'll get a segfault.\n");
					blop[wa->numl].n=malloc(couc*sizeof(char));
					blop[wa->numl].nsz=couc;
					strcpy(blop[wa->numl].n, bufword);
					qfpresz=(size_t)(tmpc2-tmpc)-1UL;
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

	if(idanomals)
		printf("Warning, idanomals was not zero, it was %i, so some IDs and are not exactly the same as NAMEs.\n", idanomals); 

	return blop;
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

gf23_t *processgf23(char *fname, int *m, int *n) /* this is dummy nmae .. it's for the special gff format for ;yze annottion */
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
	char *sm=NULL /*first semicolon marker*/, *fem=NULL /* first equals marker */;
	int cncols /* canonical numcols ... we will use the first line */;

	gf23_t *gf23=malloc(GBUF*sizeof(gf23_t));

	while( (c=fgetc(fp)) != EOF) { /* grab a char */
		if( (c== '\n') | (c == ' ') | (c == '\t') | (c=='#')) { /* word closing events */
			if( inword==1) { /* first word closing event */
				wa->wln[couw]=couc;
				bufword[couc++]='\0';
				bufword = realloc(bufword, couc*sizeof(char)); /* normalize */
				/* for the struct, we want to know if it's the first word in a line, like so: */
				if(couw==oldcouw) {
					gf23[wa->numl].n=malloc(couc*sizeof(char));
					gf23[wa->numl].nsz=couc;
					strcpy(gf23[wa->numl].n, bufword);
				} else if((couw-oldcouw)==3) { /* fourth col */
					gf23[wa->numl].c[couw-oldcouw-3]=atol(bufword)-1L; // change to zero indexing
				} else if((couw-oldcouw)==4) { /* it's not the first word, and it's 1st and second col */
					gf23[wa->numl].c[couw-oldcouw-3]=atol(bufword); // no 0 indexing change required here.
				} else if((couw-oldcouw)==6 )  { /* the strand: simply grab first character, it will be + or - */
					gf23[wa->numl].sd=bufword[0];
				} else if( (couw-oldcouw)==2) { // the type string
					gf23[wa->numl].t=malloc(couc*sizeof(char));
					gf23[wa->numl].tsz=couc;
					strcpy(gf23[wa->numl].t, bufword);
				} else if( (couw-oldcouw)==9) { // the ID string at the end, this time at col 10
					fem=NULL;
					sm=NULL;
					fem=strchr(bufword, '=');
					sm=strchr(bufword, ';');
					if(fem) {
						if(sm) {
							gf23[wa->numl].isz=(size_t)(sm-fem);
							gf23[wa->numl].i=malloc(gf23[wa->numl].isz*sizeof(char));
							memcpy(gf23[wa->numl].i, fem+1, (gf23[wa->numl].isz-1)*sizeof(char)); // strncpy writes an extra bit for \0
							gf23[wa->numl].i[gf23[wa->numl].isz-1]='\0'; // null terminate
						} else {
							gf23[wa->numl].isz=1+strlen(fem+1);
							gf23[wa->numl].i=malloc(gf23[wa->numl].isz*sizeof(char));
							memcpy(gf23[wa->numl].i, fem+1, (gf23[wa->numl].isz-1)*sizeof(char)); // strncpy writes an extra bit for \0
							gf23[wa->numl].i[gf23[wa->numl].isz-1]='\0'; // null terminate
						}
					} else {
						gf23[wa->numl].i=malloc(couc*sizeof(char));
						gf23[wa->numl].isz=couc;
						strcpy(gf23[wa->numl].i, bufword);
					}
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
					gf23=realloc(gf23, wa->lbuf*sizeof(gf23_t));
					memset(wa->wpla+(wa->lbuf-WBUF), 0, WBUF*sizeof(size_t));
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
				gf23=realloc(gf23, wa->wsbuf*sizeof(gf23_t));
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
	gf23 = realloc(gf23, wa->quan*sizeof(gf23_t)); /* normalize */
	wa->wpla= realloc(wa->wpla, wa->numl*sizeof(size_t));

	*m= wa->numl;
	*n= cncols; 
	free_wseq(wa);

	return gf23;
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
		printf("%s\tYGD\t%s\t%li\t%li\t.\t%c\t.\tID=%s;NAME=%s\t%i\n", gf22[i].n, gf22[i].t, gf22[i].c[0], gf22[i].c[1], gf22[i].sd, gf22[i].i, gf22[i].i, gf22[i].fc);

#ifdef DBG
	printf("You have just seen the %i entries of basic gff2 file called \"%s\".\n", m7, fname); 
#endif
	return;
}

void prtgf23(char *fname, gf23_t *gf23, int m)
{
	int i;
	for(i=0;i<m;++i) // note how we cut out the spurious parts of the motif string to leave it pure and raw (slightly weird why two-char deletion is necessary.
		printf("%s\t%li\t%li\t%c\t%s\t%s\n", gf23[i].n, gf23[i].c[0], gf23[i].c[1], gf23[i].sd, gf23[i].t, gf23[i].i);

#ifdef DBG
	printf("You have just seen the %i entries of basic gff2 file called \"%s\".\n", m, fname); 
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

void prtygl(char *fname, ygl_t *yglst, int myg)
{
	int i;
	for(i=0;i<myg;++i) 
		printf("%s\t%s\t%s\t%s\n", yglst[i].s, yglst[i].y, yglst[i].g, yglst[i].d);

	printf("You have just seen the %i entries of yeast genome lst file called \"%s\".\n", myg, fname); 
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

void prtbed2fo_(char *fname, bgr_t2 *bgrow, int m, int n, char *label) /* will print all -f feature bed file */
{
	int i;
	printf("%s file called %s is %i rows by %i columns and has following features:\n", label, fname, m, n); 
	printf("Type quantities:\n"); 
	tcat tcarr[TCQUAN]={0};
	for(i=0;i<m;++i) {
		tcarr[bgrow[i].tc-1]++;
		if(bgrow[i].tc == UCT) 
			printf("%i|%s ", i, bgrow[i].t);
		printf("\n"); 
	}
	// for(i=0;i<TCQUAN;++i)
	// 	printf((i==TCQUAN-1)?"%s\n":"%s\t", tcnames[i]); 
	// for(i=0;i<TCQUAN;++i)
	// 	printf((i==TCQUAN-1)?"%i\n":"%i\t", tcarr[i]); 
	for(i=0;i<TCQUAN;++i)
		printf("%s:%i\n", tcnames[i], tcarr[i]); 

	return;
}

void prtbed2fo(char *fname, bgr_t2 *bgrow, int m, int n, char *label) /* will print all -f feature bed file */
{
	int i;
	printf("%s file called %s is %i rows by %i columns and has following features:\n", label, fname, m, n); 
	printf("You can direct these names into a file and then present to this program again under the -u option,\n");
	printf("whereupon only those names will be looked at\n");
	for(i=0;i<m;++i) {
		printf("%s\t%c\t%s\t%s\n", bgrow[i].f, bgrow[i].sd, bgrow[i].s, bgrow[i].t);
	}

	return;
}

ia_t *mablopbed2(blop_t *blop, bgr_t2 *bed2, int mb, int m2) /* match up blop and bed2: blop is usually un ordered! */
{
	int i, j, k;
	int reghits; /* hits for region: number of lines in bed1 which coincide with a region in bed2 */
	int rbf /*reghit buf */;
	int *ra=NULL; /* array for per-outer-loop rhits */
	int cloci; /* as opposed to hit, catch the number of loci */
	float cpct; /* cutoff pct */
	int rangecov=0;
	long rbeg, rend; /* real start, real end */
	// int istarthere=0, catchingi=0;
	boole startcaught, endcaught; // final two imply other end is not caught
	ia_t *iabed2=malloc(sizeof(ia_t));
	iabed2->b=GBUF;
	iabed2->z=0;
	iabed2->a=calloc(iabed2->b, sizeof(int));

	/* outloop governed by bed2 ... it is more likely to have whole repat sections inside it */
	for(j=0;j<m2;++j) {
		startcaught=0;
		endcaught=0;
		reghits=0;
		rbf=GBUF;
		ra=calloc(rbf, sizeof(int));
		cloci=0;
		for(i=0;i<mb;++i) {
			if( !strcmp(blop[i].tc, bed2[j].n) ) {
				if((blop[i].c[0] >= bed2[j].c[0]) & (blop[i].c[1] <= bed2[j].c[1]) ) {
					CONDREALLOC(reghits, rbf, GBUF, ra, int);
					ra[reghits]=i;
					reghits++;
					CONDREALLOC(iabed2->z, iabed2->b, GBUF, iabed2->a, int);
					iabed2->a[iabed2->z]=j;
					iabed2->z++;
					rend=blop[i].c[1];
				   	rbeg=blop[i].c[0]; // range covered by this hit
					rangecov=rend-rbeg; // range covered by this hit
#ifdef DBG
					printf("r:%li-%li\n", rend, rbeg);
#endif
					cloci+=rangecov;
					startcaught=1;
					endcaught=1;
				} else if((blop[i].c[0] >= bed2[j].c[0]) & (blop[i].c[0] < bed2[j].c[1]) & (blop[i].c[1] > bed2[j].c[1]) ) {
					CONDREALLOC(reghits, rbf, GBUF, ra, int);
					ra[reghits]=i;
					reghits++;
					rend=bed2[j].c[1];
					rbeg=blop[i].c[0]; // range covered by this hit
					rangecov=rend-rbeg; // range covered by this hit
#ifdef DBG
					printf("r:%li-%li\n", rend, rbeg);
#endif
					cloci+=rangecov;
					startcaught=1;
				} else if((blop[i].c[0] < bed2[j].c[0]) & (blop[i].c[1] >= bed2[j].c[0]) & (blop[i].c[1] <= bed2[j].c[1]) ) {
					CONDREALLOC(reghits, rbf, GBUF, ra, int);
					ra[reghits]=i;
					reghits++;
					rend=blop[i].c[1];
					rbeg=bed2[j].c[0];
					rangecov=rend-rbeg; // range covered by this hit
#ifdef DBG
					printf("r:%li-%li\n", rend, rbeg);
#endif
					cloci+=rangecov;
					endcaught=1;
				}
			// } else if((startcaught) | (endcaught)) { // because it's an else to the above if, will catch first untruth after a series of truths.
			// 	startcaught=2;
			// 	endcaught=2;
			// 	break; // bed1 is ordered so we can forget about trying to match anymore.
			}
		}
		// if((startcaught==2) & (endcaught==2))
		// 	istarthere=catchingi+1;
		if(startcaught & endcaught) { // only print bed2 lines if they whollycaught anything.
#ifdef DBG
			printf("bed2idx %i in chr:%s w/ feature %s / xtnt %li got %i hits from blop being %i loci\n", j, bed2[j].n, bed2[j].f, bed2[j].c[1]-bed2[j].c[0], reghits, cloci);
#else
			cpct=(100.*cloci)/(bed2[j].c[1]-bed2[j].c[0]);
			if(cpct>CUTOFFPCT) {
				printf("%s\t%li\t%li\t%s\t%i\t%4.4f%%\t", bed2[j].n, bed2[j].c[0], bed2[j].c[1], bed2[j].f, cloci, cpct);
				for(k=0;k<reghits;++k) 
					printf("%s|", blop[k].n);
				printf("\n"); 
			}
#endif
		}
		free(ra);
		// if(istarthere >= mb)
		// 	break;
	}
	printf("Postheader:\n"); 
	printf("%s\t%s\t%s\t%s\t%s\n", "GF22NAM", "BEGGF22", "ENDGFF2", "FEATIDNAM", "PCTCOVBYRMF");
	return iabed2;
}

void prtsqbdggf22(i_s *sqisz, gf22_t *gf22, int m, int sz, fcat fc)
{
	int i, j;
	long kk;
	char rangestr[128]={0}; /* generally helpful to say what range is being given */
	if(fc!=ALL) {
		for(j=0;j<m;++j) 
			for(i=0;i<sz;++i) {
				if((gf22[j].fc == fc) && !(strcmp(sqisz[i].id, gf22[j].n))) {
					sprintf(rangestr, "|ID:%s|range_%li:%li|%c", gf22[j].i, gf22[j].c[0], gf22[j].c[1], gf22[j].sd);
					printf(">%s", sqisz[i].id);
					printf("%s\n", rangestr);
					if(gf22[j].sd =='+')
						printf("%.*s\n", (int)(gf22[j].c[1]-gf22[j].c[0]), sqisz[i].sq+gf22[j].c[0]);
					else if(gf22[j].sd =='-') {
						kk=gf22[j].c[1]-1;
						while(kk>=gf22[j].c[0]) { switch (sqisz[i].sq[kk--]) { case 'A': putchar('T'); break; case 'C': putchar('G'); break; case 'G': putchar('C'); break; case 'T': putchar('A'); break; default: break; } }
						putchar('\n');
					}
					break;
				}
		}
	} else { // will print out all annotated sequences
		for(j=0;j<m;++j) 
			for(i=0;i<sz;++i) {
				if(!strcmp(sqisz[i].id, gf22[j].n)) {
					if(gf22[j].sd =='+') {
						sprintf(rangestr, "|ID:%s|range_%li:%li|%c", gf22[j].i, gf22[j].c[0], gf22[j].c[1], gf22[j].sd);
						printf(">%s", sqisz[i].id);
						printf("%s\n", rangestr);
						printf("%.*s\n", (int)(gf22[j].c[1]-gf22[j].c[0]), sqisz[i].sq+gf22[j].c[0]);
					} else if(gf22[j].sd =='-') {
						sprintf(rangestr, "|ID:%s|range_%li:%li|%c", gf22[j].i, gf22[j].c[0], gf22[j].c[1], gf22[j].sd);
						printf(">%s", sqisz[i].id);
						printf("%s\n", rangestr);
						kk=gf22[j].c[1]-1;
						while(kk>=gf22[j].c[0]) { switch (sqisz[i].sq[kk--]) { case 'A': putchar('T'); break; case 'C': putchar('G'); break; case 'G': putchar('C'); break; case 'T': putchar('A'); break; default: break; } }
						putchar('\n');
					} else if(gf22[j].sd =='.') {
						sprintf(rangestr, "|ID:%s|range_%li:%li|.+", gf22[j].i, gf22[j].c[0], gf22[j].c[1]);
						printf(">%s", sqisz[i].id);
						printf("%s\n", rangestr);
						printf("%.*s\n", (int)(gf22[j].c[1]-gf22[j].c[0]), sqisz[i].sq+gf22[j].c[0]);
						sprintf(rangestr, "|ID:%s|range_%li:%li|.-", gf22[j].i, gf22[j].c[0], gf22[j].c[1]);
						printf(">%s", sqisz[i].id);
						printf("%s\n", rangestr);
						kk=gf22[j].c[1]-1;
						while(kk>=gf22[j].c[0]) { switch (sqisz[i].sq[kk--]) { case 'A': putchar('T'); break; case 'C': putchar('G'); break; case 'G': putchar('C'); break; case 'T': putchar('A'); break; default: break; } }
						putchar('\n');
					}
					break;
				}
		}
	}
	return;
}

void mbed2yg(char *fname0, char *fname2, bgr_t2 *bgrow, int m2, ygl_t *yglst, int myg, ygsnod **stab, unsigned tsz)
{
	unsigned char found;
	int i, k=0;
    ygsnod *tsnod0=NULL, *tsnod2=NULL;
    unsigned tint;
	unsigned gbuf=GBUF;
	int *nfiarr=malloc(gbuf*sizeof(int)); /* the not found index array */
	printf("%s bed file with %i records matched again YG list file %s with %i records:\n", fname0, m2, fname2, myg); 
	for(i=0;i<m2;++i) {
		found = 0;
        tint=hashit(bgrow[i].f, tsz);
        tsnod2=stab[tint];
        while( (tsnod2 != NULL) ){
            if(!strcmp(tsnod2->yglst->y, bgrow[i].f)) {
				printf("%s\t%s\t%s matches %s\t%s\t%s\n", bgrow[i].f, bgrow[i].s, bgrow[i].t, tsnod2->yglst->s, tsnod2->yglst->g, tsnod2->yglst->d);
				found = 1;
				break;
			}
            tsnod0=tsnod2;
            tsnod2=tsnod2->n;
		}
		if(!found) {
			CONDREALLOC(k, gbuf, GBUF, nfiarr, int);
			nfiarr[k++]=i;
		}
	}
	nfiarr=realloc(nfiarr, k*sizeof(int));
	printf("features from %s which did not match:\n", fname0);
	for(i=0;i<k;++i) 
		printf("%i|%s ", nfiarr[i], bgrow[nfiarr[i]].f); 
	printf("\n"); 
	free(nfiarr);

	return;
}

void mbed2ya(char *fname0, char *fname2, bgr_t2 *bgrow, int m2, bt2snod **stab2, unsigned htsz2, gf22_t *gf22, int m7, gf22snod **stab, unsigned tsz)
{
	unsigned char found;
	int numfound=0i, numfound2=0;
	int i, k=0;
    gf22snod *tsnod0=NULL, *tsnod2=NULL;
    bt2snod *tsnodd0=NULL, *tsnodd2=NULL;
    unsigned tint;
	unsigned gbuf=GBUF;
	int *nfiarr=malloc(gbuf*sizeof(int)); /* the not found index array */
	printf("%s bed file with %i records matched against YA file %s with %i records:\n", fname0, m2, fname2, m7); 
	for(i=0;i<m2;++i) {
		found = 0;
        tint=hashit(bgrow[i].f, tsz);
        tsnod2=stab[tint];
        while( (tsnod2 != NULL) ){
            if(!strcmp(tsnod2->gf22->i, bgrow[i].f)) {
				printf("%s\t%s\t%s matches.\n", bgrow[i].f, bgrow[i].s, bgrow[i].t);
				numfound++;
				found = 1;
				break;
			}
            tsnod0=tsnod2;
            tsnod2=tsnod2->n;
		}
		if(!found) {
			CONDREALLOC(k, gbuf, GBUF, nfiarr, int);
			nfiarr[k++]=i;
		}
	}
	nfiarr=realloc(nfiarr, k*sizeof(int));
	printf("%i features found, though %i features in \"%s\" did not match \"%s\".\n", numfound, k, fname0, fname2);
#if DBG
	printf("features from %s which did not match:\n", fname0);
	for(i=0;i<k;++i) 
		printf("%i|%s ", nfiarr[i], bgrow[nfiarr[i]].f); 
	printf("\n"); 
#endif
	free(nfiarr);
	gbuf=GBUF;
	k=0;
	nfiarr=malloc(gbuf*sizeof(int));
	for(i=0;i<m7;++i) {
		found = 0;
        tint=hashit(gf22[i].i, htsz2);
        tsnodd2=stab2[tint];
        while( (tsnodd2 != NULL) ) {
			printf("%s vs. %s\n", tsnodd2->bed2->f, gf22[i].i);
            if(!strcmp(tsnodd2->bed2->f, gf22[i].i)) {
				printf("%s in \"%s\" matches at %s, %li to %li.\n", gf22[i].i, fname2, tsnodd2->bed2->n, tsnodd2->bed2->c[0], tsnodd2->bed2->c[1]);
				numfound2++;
				found = 1;
				break;
			}
            tsnodd0=tsnodd2;
            tsnodd2=tsnodd2->n;
		}
		if(!found) {
			CONDREALLOC(k, gbuf, GBUF, nfiarr, int);
			nfiarr[k++]=i;
		}
	}
	nfiarr=realloc(nfiarr, k*sizeof(int));
	printf("%i features found, though %i features in \"%s\" did not match \"%s\".\n", numfound2, k, fname2, fname0);
	free(nfiarr);

	return;
}

void mbed2ya1(char *fname0, char *fname2, bgr_t2 *bgrow, int m2, bt2snod **stab2, unsigned htsz2, gf22_t *gf22, int m7, gf22snod **stab, unsigned tsz, tcat ztcat) /* S288 ref is very big .. we only want genes */
{
	/* this version of mbed2ya only check for 1 feature type, usually gene. */
	unsigned char found;
	int numfound=0, numfound2=0;
	int i, k=0;
    gf22snod *tsnod0=NULL, *tsnod2=NULL;
    bt2snod *tsnodd0=NULL, *tsnodd2=NULL;
    unsigned tint;
	unsigned gbuf=GBUF;
	int *nfiarr=malloc(gbuf*sizeof(int)); /* the not found index array */
	printf("%s bed file with %i records matched against YA file %s with %i records:\n", fname0, m2, fname2, m7); 
	for(i=0;i<m2;++i) {
		if(bgrow[i].tc != ztcat)
			continue;
		found = 0;
        tint=hashit(bgrow[i].f, tsz);
        tsnod2=stab[tint];
        while( (tsnod2 != NULL) ){
            if(!strcmp(tsnod2->gf22->i, bgrow[i].f)) {
#ifdef DBG
				printf("%s\t%s\t%s matches.\n", bgrow[i].f, bgrow[i].s, bgrow[i].t);
#endif
				numfound++;
				found = 1;
				break;
			}
            tsnod0=tsnod2;
            tsnod2=tsnod2->n;
		}
		if(!found) {
			CONDREALLOC(k, gbuf, GBUF, nfiarr, int);
			nfiarr[k++]=i;
		}
	}
	nfiarr=realloc(nfiarr, k*sizeof(int));
	printf("%i features found, though %i features in \"%s\" did not match \"%s\".\n", numfound, k, fname0, fname2);
#if DBG
	printf("features from %s which did not match:\n", fname0);
	for(i=0;i<k;++i) 
		printf("%i|%s ", nfiarr[i], bgrow[nfiarr[i]].f); 
	printf("\n"); 
#endif
	free(nfiarr);
	gbuf=GBUF;
	k=0;
	nfiarr=malloc(gbuf*sizeof(int));
	for(i=0;i<m7;++i) {
		found = 0;
        tint=hashit(gf22[i].i, htsz2); // be careful, had copy-pasted and forgot about changing tsz to htsz2
        tsnodd2=stab2[tint];
        while( (tsnodd2 != NULL) ) {
			if(tsnodd2->bed2->tc != ztcat)
				goto keepmoving; // not much point really, now we're on the row of the hashtable
#ifdef DBG2
			printf("%s vs. %s\n", tsnodd2->bed2->f, gf22[i].i);
#endif
            if(!strcmp(tsnodd2->bed2->f, gf22[i].i)) {
#ifdef DBG
				printf("%s in \"%s\" matches at %s, %li to %li.\n", gf22[i].i, fname2, tsnodd2->bed2->n, tsnodd2->bed2->c[0], tsnodd2->bed2->c[1]);
#endif
				numfound2++;
				found = 1;
				break;
			}
keepmoving: tsnodd0=tsnodd2;
            tsnodd2=tsnodd2->n;
		}
		if(!found) {
			CONDREALLOC(k, gbuf, GBUF, nfiarr, int);
			nfiarr[k++]=i;
		}
	}
	nfiarr=realloc(nfiarr, k*sizeof(int));
	printf("%i features found, though %i features in \"%s\" did not match \"%s\".\n", numfound2, k, fname2, fname0);
	free(nfiarr);

	return;
}

void mbed2wof(char *fname0, char *fname2, bgr_t2 *bgrow, int m2, bt2snod **stab2, unsigned tsz, words_t *bedword, int m3) /* search beds2 feature names from file of words */
{
	unsigned char found;
	int numfound2=0;
	int i, k=0;
    bt2snod *tsnodd0=NULL, *tsnodd2=NULL;
    unsigned tint;
	unsigned gbuf=GBUF;
	int *nfiarr=malloc(gbuf*sizeof(int)); /* the not found index array */
	printf("%s bed file with %i records matched against words file file %s with %i records:\n", fname0, m2, fname2, m3); 
	for(i=0;i<m3;++i) {
		found = 0;
        tint=hashit(bedword[i].n, tsz);
        tsnodd2=stab2[tint];
        while( (tsnodd2 != NULL) ) {
            if(!strcmp(tsnodd2->bed2->f, bedword[i].n)) {
				printf("%s in \"%s\" matches at %s, %li to %li.\n", bedword[i].n, fname2, tsnodd2->bed2->n, tsnodd2->bed2->c[0], tsnodd2->bed2->c[1]);
				numfound2++;
				found = 1;
				break;
			}
keepmoving: tsnodd0=tsnodd2;
            tsnodd2=tsnodd2->n;
		}
		if(!found) { // we're counting non-found items.
			CONDREALLOC(k, gbuf, GBUF, nfiarr, int);
			nfiarr[k++]=i;
		}
	}
	if(k) {
		nfiarr=realloc(nfiarr, k*sizeof(int)); // k can be zero
		printf("%i features found, though %i features in \"%s\" did not match \"%s\".\n", numfound2, k, fname2, fname0);
	} else
		printf("%i features found, which were all in \"%s\"\n", numfound2, fname2);

	free(nfiarr);

	return;
}

void mblop2gf22(char *fname0, char *fname2, gf22_t *gf22, int m7, gf22snod **stab, unsigned tsz, blop_t *blop, int mb) /* match blastoutput to the gf22 ... but only is ABSENTCOLS see alternate function for full gf22 search ... messed up!*/
{
	unsigned char found;
	int numfound2=0;
	int i, k=0;
	char *tmpc=NULL;
    gf22snod *tsnodd0=NULL, *tsnodd2=NULL;
    unsigned tint;
	unsigned gbuf=GBUF;
	int *nfiarr=malloc(gbuf*sizeof(int)); /* the not found index array */
	float thresh=0.95;
	printf("%s gf22 file with %i records matched against blast outp file %s with %i records at %2.1f identity:\n", fname0, m7, fname2, mb, thresh);
	for(i=0;i<mb;++i) {
		found = 0;
        tint=hashit(blop[i].qfs, tsz);
        tsnodd2=stab[tint];
        while( (tsnodd2 != NULL) ) {
			/* note how purposefully we match on the query feature, but are really interested in the target feature! */
            if(!strcmp(tsnodd2->gf22->i, blop[i].qfs)) {
				if(tsnodd2->gf22->fc != ABS) { //found, but we're not interested in this type
					if(tsnodd2->gf22->blofound==1)
						tsnodd2->gf22->blofound=2; // second time found
					else if(tsnodd2->gf22->blofound==0)
						tsnodd2->gf22->blofound=1;
					break;
				}
				if(tsnodd2->gf22->blofound)
					break; //this gff id already found .. i.e. top blast result
				tsnodd2->gf22->blofound=1;
				tmpc=strchr(blop[i].fs, '_');
				//orig:
				// printf("%s in \"%s\" matches queryname at %s, %li to %li with pcti %4.2f and e-val %4.6f\n", blop[i].qfs, fname2, tsnodd2->gf22->n, tsnodd2->gf22->c[0], tsnodd2->gf22->c[1], blop[i].pcti, blop[i].eval);
				// plain, no pcti nor evals
				if((tmpc!=NULL) & (blop[i].pcti>95.0))// if there's a _ it's an _mRNA
					printf("%s\tYGD\t%s\t%li\t%li\t.\t%c\t.\tID=%.*s;NAME=%.*s\n", tsnodd2->gf22->n, tsnodd2->gf22->t, tsnodd2->gf22->c[0]+1L,tsnodd2->gf22->c[1], tsnodd2->gf22->sd, (int)(tmpc-blop[i].fs), blop[i].fs, (int)(tmpc-blop[i].fs), blop[i].fs);
				else if((tmpc==NULL) & (blop[i].pcti>95.0)) // if not _mRNA, don't know what to do: just it all out.
					printf("%s\tYGD\t%s\t%li\t%li\t.\t%c\t.\tID=%s;NAME=%s\n", tsnodd2->gf22->n, tsnodd2->gf22->t, tsnodd2->gf22->c[0]+1L,tsnodd2->gf22->c[1], tsnodd2->gf22->sd, blop[i].fs, blop[i].fs);
				numfound2++;
				found = 1;
				break;
			}
			tsnodd0=tsnodd2;
            tsnodd2=tsnodd2->n;
		}
		if(!(found) & (tsnodd2->gf22->fc==ABS)) {
			printf("%s\tYGD\t%s\t%li\t%li\t.\t%c\t.\tNOHIT>%2.1fID\n", tsnodd2->gf22->n, tsnodd2->gf22->t, tsnodd2->gf22->c[0]+1L,tsnodd2->gf22->c[1], tsnodd2->gf22->sd, thresh);
			CONDREALLOC(k, gbuf, GBUF, nfiarr, int);
			nfiarr[k++]=i;
		} else if((tsnodd2->gf22->blofound==1) & (tsnodd2->gf22->fc!=ABS)) { // print out normal gf22
			printf("%s\tYGD\t%s\t%li\t%li\t.\t%c\t.\tID=%s;NAME=%s\n", tsnodd2->gf22->n, tsnodd2->gf22->t, tsnodd2->gf22->c[0]+1L,tsnodd2->gf22->c[1], tsnodd2->gf22->sd, tsnodd2->gf22->i, tsnodd2->gf22->i);
		}
	}
	if(k) {
		nfiarr=realloc(nfiarr, k*sizeof(int)); // k can be zero
		// printf("%i features found, though %i features in \"%s\" did not match \"%s\".\n", numfound2, k, fname2, fname0);
	} else
		// printf("%i features found, which were all in \"%s\"\n", numfound2, fname2);

	free(nfiarr);

	return;
}

void mgf222blop(char *fname0, char *fname2, gf22_t *gf22, int m7, blosnod **stab, unsigned tsz, blop_t *blop, int mb) /* match blastoutput to the gf22 ... but only is ABSENTCOLS see alternate function for full gf22 search */
{
	unsigned char found;
	int numfound2=0;
	int i, k=0;
	char *tmpc=NULL;
    blosnod *tsnodd0=NULL, *tsnodd2=NULL;
    unsigned tint;
	unsigned gbuf=GBUF;
	int *nfiarr=malloc(gbuf*sizeof(int)); /* the not found index array */
	float thresh=0.8;
	printf("%s gf22 file with %i records matched against blast outp file %s with %i records at %2.1f identity:\n", fname0, m7, fname2, mb, thresh);
	for(i=0;i<m7;++i) {
		if(gf22[i].fc != ABS) { // in this version of thefunction we're not interested in non-ABS
			printf("%s\tYGD\t%s\t%li\t%li\t.\t%c\t.\tID=%s;NAME=%s\n", gf22[i].n, gf22[i].t, gf22[i].c[0]+1L,gf22[i].c[1], gf22[i].sd, gf22[i].i, gf22[i].i);
			continue;
		}
		found = 0;
        tint=hashit(gf22[i].i, tsz);
        tsnodd2=stab[tint];
        while( (tsnodd2 != NULL) ) {
            if(!strcmp(tsnodd2->blop->qfs, gf22[i].i)) {
				tmpc=strchr(tsnodd2->blop->fs, '_');
				// plain, no pcti nor evals
				if((tmpc!=NULL) & (blop[i].pcti>95.0))// if there's a _ it's an _mRNA
					printf("%s\tYGD\t%s\t%li\t%li\t.\t%c\t.\tID=%.*s;NAME=%.*s\n", gf22[i].n, gf22[i].t, gf22[i].c[0]+1L,gf22[i].c[1], gf22[i].sd, (int)(tmpc-tsnodd2->blop->fs), tsnodd2->blop->fs, (int)(tmpc-tsnodd2->blop->fs), tsnodd2->blop->fs);
				else if((tmpc==NULL) & (blop[i].pcti>95.0)) // if not _mRNA, don't know what to do: just it all out.
					printf("%s\tYGD\t%s\t%li\t%li\t.\t%c\t.\tID=%s;NAME=%s\n", gf22[i].n, gf22[i].t, gf22[i].c[0]+1L,gf22[i].c[1], gf22[i].sd, tsnodd2->blop->fs, tsnodd2->blop->fs);
				else
					printf("%s\tYGD\t%s\t%li\t%li\t.\t%c\t.\tHIT1UNDER%3.1f\n", gf22[i].n, gf22[i].t, gf22[i].c[0]+1L,gf22[i].c[1], gf22[i].sd, thresh);
				numfound2++;
				found = 1;
				break;
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
	if(k) {
		nfiarr=realloc(nfiarr, k*sizeof(int)); // k can be zero
		// printf("%i features found, though %i features in \"%s\" did not match \"%s\".\n", numfound2, k, fname2, fname0);
	} else
		// printf("%i features found, which were all in \"%s\"\n", numfound2, fname2);

	free(nfiarr);

	return;
}

void prtbed2fo2(char *fname, bgr_t2 *bgrow, int m, int n, char *label) /* all -f bed file features except dots and _mRNA ending ones */
{
	int i;
	char *usc;
	printf("%s file called %s is %i rows by %i columns and has following features:\n", label, fname, m, n); 
	printf("You can direct these name into a file and then present to this program again under the -u option,\n");
	printf("whereupon only those names will be looked at\n");
	printf("Note: only feature names without dot and underscore will be printed\n");
	for(i=0;i<m;++i) {
		usc=NULL;
		if(!strcmp(bgrow[i].f, "."))
			continue;
		else {
			usc=strchr(bgrow[i].f, '_');
			if(usc)
				continue;
			printf("%s\n", bgrow[i].f);
		}
	}

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
	printf("chipmonk: takes a variety of genomic features ifles (principally yeast, but can be modified) with following options:\n");
	printf("-i, this is for a bedgraph type file.\n");
	printf("-f, this is for a feature file in bed format, probably converted from gff.\n");
	printf("-d, this is a flag, it stands of details, to be used with a single filename specification option, merely allows you see details of the file.\n");
	printf("-a: takes a fasta file.\n");
	printf("-h: takes a GFF3 format, but rejects . and mRNA and grabs ID on 10th not 9th column.\n");
	printf("-y: takes a variation of the GFF2 format, a table, column nine gives the Y-name. eg. Output of YG annotator service\n");
	printf("-l: takes a Yeast Genome listing file ... obtained from their website.\n");
	printf("-n: the names flag ... to only print out feature names.\n");
	return;
}

int main(int argc, char *argv[])
{
	/* argument accounting */
	if(argc == 1) {
		prtusage();
		exit(EXIT_FAILURE);
	}
	int i, m, n /*rows,cols for bgr_t*/, m2, n2 /* for the bed2 */, m3, n3 /*for the words_t file */, m4, n4, m4b, n4b, m5, n5, m6, n6, m7, n7 /*gf22 dims */, m8, n8 /* gf23 dims */, myg, nyg, mb, nb /* for the blop */;
    unsigned numsq; /* number of sequences in the -a (fasta) option */
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
	gf23_t *gf23=NULL;
	ygl_t *yglst=NULL;
	i_s *sqisz=NULL;

	/* for the ystr, gf22 hash handling */
    unsigned htsz, htsz2, htsz3, htszyg, htszblo;
    gf22snod **stab=NULL;
    bt2snod **stab2=NULL;
    gf23snod **stab3=NULL;
    ygsnod **stabyg=NULL;
    blosnod **stabblo=NULL;
	blop_t *blop=NULL;
	tcat ztcat;
	boole zfound=0;

	if(opts.lstr) {
		yglst=processygl(opts.lstr, &myg, &nyg);
		htszyg=2*myg/3;
		stabyg=ygtochainharr(yglst, myg, htszyg);
	}
	if(opts.bstr) {
		blop=processblop(opts.bstr, &mb, &nb);
		htszblo=2*mb/3;
		stabblo=blotochainharr(blop, mb, htszblo);
	}
	if(opts.istr)
		bgrow=processinpf(opts.istr, &m, &n);
	if(opts.fstr) {
		bed2=processinpf2(opts.fstr, &m2, &n2);
		htsz2=2*m2/3; /* our hash table size */
		stab2=bt2tochainharr(bed2, m2, htsz2);
	}
	if(opts.zstr){
		for(i=0;i<TCQUAN-1;++i) 
			if(!strcmp(tcnames[i], opts.zstr)) {
				zfound=1;
				ztcat=i+1;
				break;
			}
		if(!zfound) {
			printf("Error. That string you entered with -z is not valid, check the tcnames global\n"); 
			goto final;
		}
	}
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
	if(opts.hstr) { // hacky form of gff3 reading ... */
		gf23=processgf23(opts.hstr, &m8, &n8);
		htsz3=2*m8/3; /* our hash table size */
		stab3=gf23tochainharr(gf23, m8, htsz3); /* we now set up a hash table along side our sequence names from the fasta file */
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

	/* combo: df!sa */
	if((opts.dflg) && (opts.fstr) && !opts.sflg && !opts.astr && !opts.nflg && !opts.istr) {
		prtbed2fo(opts.fstr, bed2, m2, n2, "Feature (bed2)");
		// prtbt2chaharr(stab2, htsz2);
		goto final;
	}

	/* combo: nf */
	if((opts.nflg) && (opts.fstr)) {
		// prtbed2fo(opts.fstr, bed2, m2, n2, "Feature (bed2)");
		prtbed2fo2(opts.fstr, bed2, m2, n2, "Feature (bed2)"); // alterantive skipping . and _
		goto final;
	}
	// prtbed2(bed2, m2, MXCOL2VIEW);
	/* combo: if */
	if((opts.istr) && (opts.fstr))
		m2beds(bgrow, bed2, m2, m);

	if((opts.ustr) && !(opts.fstr) && (!opts.sflg) && opts.dflg) {
		printf("bedwords: "); 
		for(i=0;i<m3;++i)
			printf("%s ", bedword[i].n);
		printf("\n"); 
	}

	/* -d -u and -f, bit stealth here, looking for words in -f's features */
	/* use as ./chipmonk -d -f sacchmys_annotsnochrline.bed -u gna.txt, where your gna.txt has Y-names for example */
	if((opts.ustr) && (opts.fstr) && opts.dflg && !(opts.astr))
		mbed2wof(opts.fstr, opts.ustr, bed2, m2, stab2, htsz2, bedword, m3);

	if(opts.ustr && opts.fstr && !(opts.dflg) && opts.astr) {
		prtmbed2wof(sqisz, numsq, bed2, m2, stab2, htsz2, bedword, m3);
	}

	if((opts.pstr) && (opts.fstr) )
		md2bedp(dpf, bed2, m2, m4);

	if((opts.dflg) && (opts.lstr) ) {
		// ./chipmonk -d -l yg0.lst
		prtygl(opts.lstr, yglst, myg);
		prtygchaharr(stabyg, htszyg);
	}

	if((opts.fstr) && (opts.lstr) )
		// ./chipmonk -f sacchmys_annotsnochrline.bed -l ygallg.tsv
		// with his yo can match up the S228 annotation witht he YG list
		mbed2yg(opts.fstr, opts.lstr, bed2, m2, yglst, myg, stabyg, htszyg);

	if((opts.fstr) && (opts.ystr) && opts.zstr)
		mbed2ya1(opts.fstr, opts.ystr, bed2, m2, stab2, htsz2, gf22, m7, stab, htsz, ztcat);

	if((opts.dflg) && (opts.rstr) )
		prtrmf(opts.rstr, rmf, m6);

	/* ystr handling */
	if((opts.dflg) && (opts.ystr) ) {
		prtgf22(opts.ystr, gf22, m7);
		// prtgf22chaharr(stab, htsz);
	}
	fcat fc;
	if((opts.astr) && (opts.ystr) && (opts.xstr)) {
		// prtgf22(opts.ystr, gf22, m7);
		fc=getfc(opts.xstr);
		prtsqbdggf22(sqisz, gf22, m7, numsq, fc);
		// prtgf22chaharr(stab, htsz);
	} else if((opts.astr) && (opts.ystr) && !(opts.xstr)) {
		printf("I'm sorry but you must specify a gf22 file type via the -x option: ... try ALL\n"); 
		goto final;
	}

	if((opts.dflg) && (opts.hstr) ) {
		prtgf23(opts.hstr, gf23, m8);
		prtgf23chaharr(stab3, htsz3);
	}

	if((opts.dflg) && (opts.astr) ) 
		prtsq(sqisz, numsq);

	/* bedgraph and fasta file */
	if((opts.istr) && (opts.astr) )
		prtsqbdg(sqisz, bgrow, m, numsq);

	/* bed2 feature and fasta file */
	if((opts.fstr) && (opts.astr) && !opts.zstr && !opts.ustr)
		prtsqbdgf(sqisz, bed2, m2, numsq);

	if((opts.fstr) && (opts.astr) && opts.zstr)
		prtsqbdgf2(sqisz, bed2, m2, numsq, ztcat);

	if((opts.gstr) && (opts.rstr) )
		mgf2rmf(opts.gstr, opts.rstr, gf, rmf, m6, m5);

	if((opts.gstr) && (opts.pstr) )
		mgf2dpf(opts.gstr, opts.pstr, gf, dpf, m4, m5);

	if((opts.gstr) && (opts.fstr) )
		// i.e. ./chipmonk_d -f sacchmys_annotsnochrline.bed -g S288_maniid.sizes
		mgf2bed(opts.gstr, opts.fstr, gf, bed2, m2, m5);

	// if((opts.ustr) && (opts.fstr) && opts.sflg)
	// 	prtbed2s(bed2, m2, MXCOL2VIEW, bedword, m3, n3, "bed2 features that are in interesting-feature-file");

	ia_t *ia=NULL;
	if((opts.ustr) && (opts.fstr) && opts.sflg) {
		ia=gensplbdx(bed2, m2, n2, bedword, m3, n3);
		bed2in2(opts.fstr, bed2, m2, n2, ia);
	}
	
	/** bstr, dealing with **/
	if((opts.dflg) && (opts.bstr) )
		prtblochaharr(stabblo, htszblo);
		// prtblop(opts.bstr, blop, mb);
	/*routine to discover where blastoutput and annotation file match */
	ia_t *iabed2=NULL;
	if((opts.bstr) && (opts.fstr)) {
		iabed2=mablopbed2(blop, bed2, mb, m2);
		printf("%i indices where %s matched by %s:\n", iabed2->z, opts.ystr, opts.rstr); 
		for(i=0;i<iabed2->z;++i) 
			printf("%s ", bed2[iabed2->a[i]].f); 
		printf("\n"); 
		free(iabed2->a);
		free(iabed2);
	}
	if((opts.ystr) && (opts.bstr) )
		mgf222blop(opts.ystr, opts.bstr, gf22, m7, stabblo, htszblo, blop, mb);

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
		freebt2chainharr(stab2, htsz2);
		for(i=0;i<m2;++i) {
			free(bed2[i].n);
			free(bed2[i].f);
			free(bed2[i].s);
			free(bed2[i].t);
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
	if(opts.lstr) {
		freeygchainharr(stabyg, htszyg);
		for(i=0;i<myg;++i) {
			free(yglst[i].s);
			free(yglst[i].y);
			free(yglst[i].g);
			free(yglst[i].d);
		}
		free(yglst);
	}
	if(opts.ystr) {
		freegf22chainharr(stab, htsz);
		for(i=0;i<m7;++i) {
			free(gf22[i].n);
			free(gf22[i].t);
			free(gf22[i].i);
		}
		free(gf22);
	}
	if(opts.hstr) {
		freegf23chainharr(stab3, htsz3);
		for(i=0;i<m8;++i) {
			free(gf23[i].n);
			free(gf23[i].t);
			free(gf23[i].i);
		}
		free(gf23);
	}
	if(opts.astr) {
		for(i=0;i<numsq;++i) {
			free(sqisz[i].id);
			free(sqisz[i].sq);
		}
		free(sqisz);
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
