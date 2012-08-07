#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include "types.h"
#include "debugdef.h"
#include "errordef.h"
#include "protodef.h"
#include "spacedef.h"
#include "streedef.h"
#include "maxmatdef.h"




// 4 using pthread parallelization
void *procmaxmatches_pt(void *num);
void printbaseprofile();
void printcontig();
void printcontig2();
int searchrefloc(int, int*);




typedef struct {

  Suffixtree stree;            // the suffix tree of the subject-sequence

  Multiseq *subjectmultiseq,   // reference to multiseq of subject
           querymultiseq;      // the Multiseq record of the queries

  Uint minmatchlength,         // minimum length of a match
       maxdesclength,          // maximum length of a description
       currentquerylen,        // length of the current query sequence
       maxmismatch,            // # of maximum mismatches
       minsimpct,              // minimum % similarity
       strand;                 // current strand

  Uint numqseqs;

  Uchar **queryseqs;
  Uchar **queryids;
  Ushort *querylens;
  
  // store temporary top hits for each read
  // 0-10, indexed by number of mismatches
  // char * represents the output string
  char * tophits[11];

  // each mapping between read and ref may have multiple seeds, separated by mutations
  // we only need one seed to anchor the read to a particular ref loc
  int hitsloc[10000];

  // start and end loc of ref seqs in the whole refseq string.
  // used to convert locations
  int refloc[MAXREFNUM*2];

  // number of ref seqs
  int reflocnum;

  // output file for mapping info
  FILE *outmap;


} Matchprocessinfo;

