#ifndef MAXMATDEF_H
#define MAXMATDEF_H

#include <limits.h>
#include "multidef.h"


// Definition of some constants and types for computing maximal matches


#define SEPARATOR                UCHAR_MAX    // separator symbol in multiple seq
#define MMREPLACEMENTCHARSUBJECT        'x'   // replace illegal chars in reference seq
#define MMREPLACEMENTCHARQUERY          'n'   // replace illegal chars in query reads


#define MAXNUMOFQUERYFILES  32
#define MAXREFNUM           1000000    // maximum number of reference seqs


// scores for Smith-Waterman alignment
#define GAPOPEN          -5
#define GAPEXTEND        -2
#define BASEMATCH         1
#define BASEMISMATCH     -2


// parameters from command line arguments
typedef struct {


  Uint minmatchlength,          // seed length (exact match)
       numofqueryfiles,         // number of query files
       nummaxmis,               // max # of mismatches
       minpct,                  // min % similarity
       nump;                    // number of processors
  
  char program[PATH_MAX+1],     // the path of the program
       subjectfile[PATH_MAX+1], // filename of reference sequence
       queryfilelist[MAXNUMOFQUERYFILES][PATH_MAX+1],
                                // filenames of query sequences
       outprefix[200];          // output prefix
  
} MMcallinfo;


#endif

