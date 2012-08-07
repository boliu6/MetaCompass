/*
  Part of the code are borrow from Stefan Kurtz

  This file contains functions to appropriately call the function findmumcandidates()
  and findmaxmatches(), and to process their result according to the options given by the user.
*/

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
#include "procmaxmat.h"




// from findmaxmat.c
Sint findmaxmatches(Suffixtree *stree,
                    Uint minmatchlength,
                    void *processinfo,
                    Uchar *query,
                    Uint querylen,
                    Uint seqnum);




// from maxmatinp.c
Sint scanmultiplefastafile(Multiseq *multiseq,
			   char *filename,
			   Uchar replacewildcardchar,
			   Uchar *input,
			   Uint inputlen);




// complementary dna base
#define ASSIGNMAXMATCOMPLEMENT(VL,VR)			\
  switch (VR)						\
    {							\
    case 'a':						\
      VL = (Uchar) 't';					\
      break;						\
    case 'c':						\
      VL = (Uchar) 'g';					\
      break;						\
    case 'g':						\
      VL = (Uchar) 'c';					\
      break;						\
    case 't':						\
      VL = (Uchar) 'a';					\
      break;						\
    default:                        /* anything */	\
      VL = (Uchar) 'n';					\
      break;						\
    }							\



// reverse complementary
static void revcom(Uchar *seq, Uint seqlen) {
  Uchar *front, *back, tmp = 0;
  for (front = seq, back = seq + seqlen - 1; front <= back; front++, back--) {
    
    ASSIGNMAXMATCOMPLEMENT (tmp, *front);
    ASSIGNMAXMATCOMPLEMENT (*front, *back);
    *back = tmp;
  }
}



// max length of fasta title
static Sint getmaxdesclen(Multiseq *multiseq)
{
  Uint desclen, maxdesclen, seqnum;
  if(multiseq->numofsequences == 0) {
      ERROR0("multiple sequence contains 0 sequences");
      return -1;
  }
  maxdesclen = DESCRIPTIONLENGTH(multiseq,0);
  for(seqnum = UintConst(1); seqnum < multiseq->numofsequences; seqnum++) {
    desclen = DESCRIPTIONLENGTH(multiseq,seqnum);
    if(desclen > maxdesclen) {
	maxdesclen = desclen;
    }
  }
  return (Sint) maxdesclen;
}




// these global variables are initialized in the following function,
// and then used by function procmaxmatches_pt() for multiple threads;
int nump = 1;
MMcallinfo *mmcallinfo;
Multiseq *subjectmultiseq;
Matchprocessinfo matchprocessinfotmp;
Suffixtree streetmp;



// 1. Constructs the suffix tree,
// 2. Initializes the Matchprocessinfo struct appropriately,
// 3. over all sequences in querymultiseq.
// 4. Finally, the space allocated for the suffix tree is freed.

Sint procmaxmatches(MMcallinfo *mmcallinfoloc,Multiseq *subjectmultiseqloc) {


  // the following 3 lines assign info. to global variables
  
  mmcallinfo = mmcallinfoloc;
  subjectmultiseq = subjectmultiseqloc;
  nump = mmcallinfoloc -> nump;
  /**********************************************************************************/
  

  

  // Build suffix tree, store it in streetmp
  fprintf(stderr,"# construct suffix tree ... ");
  time_t seconds = time(NULL);

  if (constructstree(&streetmp,
		     subjectmultiseq->sequence,
		     subjectmultiseq->totallength) != 0) {
    return -1;
  }

  fprintf(stderr,"finished. (%d seconds)\n", time(NULL) - seconds);
  /**********************************************************************************/



  
  // allocate memory to store reference sequences
  fprintf(stderr,"# read query sequence file ... ");
  seconds = time(NULL);
  Uint reffilelen;
  Uchar *reffilecontent = CREATEMEMORYMAP (mmcallinfo->subjectfile, True, &reffilelen);
  if (reffilecontent == NULL || reffilelen == 0) {
    ERROR2("cannot open file \"%s\" or file \"%s\" is empty",mmcallinfo->subjectfile,
                                                             mmcallinfo->subjectfile);
    return -1;
  }
  /**********************************************************************************/



  
  // store the start and end loc of each ref seq in the concatenated string,
  // which is used to convert mapping locations
  
  matchprocessinfotmp.reflocnum = 0;                   // # of ref seqs
  unsigned sloc = 0, eloc = 0, tag = 0, j = 0, i = 0;

  for (i = 0; i < MAXREFNUM*2; i++)
    matchprocessinfotmp.refloc[i] = 0;                 // initialize the array
  
  for (i = 0; i < reffilelen; i++) {                   // read ref seq file char by char

    // start of a new fasta, store loc info. of previous record
    if ('>' == *(reffilecontent+i)) {           
      if (tag == 1) {
	matchprocessinfotmp.refloc[j++] = sloc;
	matchprocessinfotmp.refloc[j++] = eloc;
	matchprocessinfotmp.reflocnum++;
	eloc++;
	tag = 0;                                       // title line
      }
    }
    else if (*(reffilecontent+i) == '\n') {
      if (tag == 0) {                                  // end of title line
	sloc = eloc;  
	tag = 1;                                       // will be in seq lines from now on
      }
    }
    else if (tag == 0) {}                              // inside of title line
    else if (tag == 1) {                               // sequence line
      eloc++;
    }
  }
  matchprocessinfotmp.reflocnum++;                     // store information of last record
  matchprocessinfotmp.refloc[j++] = sloc;
  matchprocessinfotmp.refloc[j++] = eloc;
  DELETEMEMORYMAP(reffilecontent);
  /***********************************************************************************************/

  


  // load in query reads

  Uint filenum = 0;

  // compute number of sequences
  
  Uchar eachline[10000];
  Uint numseqs = 0;
  FILE *fp = fopen(mmcallinfo->queryfilelist[filenum], "r");
  while (fgets(eachline, 10000, fp)) {
    if (eachline[0] == '>') {
      ++numseqs;
    }
  }

  matchprocessinfotmp.numqseqs = numseqs;
  matchprocessinfotmp.queryseqs = malloc(sizeof(Uchar *) * numseqs);
  matchprocessinfotmp.queryids = malloc(sizeof(Uchar *) * numseqs);
  matchprocessinfotmp.querylens = malloc(sizeof(Ushort) * numseqs);

  rewind(fp);
  Uint seqnum = 0;
  Uchar id[1000];
  Uchar dnaseq[10000];
  Uchar first = 1;
  Uchar seqid[1000];
  while (fgets(eachline, 10000, fp)) {

    if (eachline[0] == '>') {


      if (first) {
	first = 0;
      }
      else {
	strcpy(matchprocessinfotmp.queryids[seqnum], seqid);
	matchprocessinfotmp.queryseqs[seqnum] = malloc(sizeof(Uchar) * strlen(dnaseq) + 1);
	strcpy(matchprocessinfotmp.queryseqs[seqnum], dnaseq);
	matchprocessinfotmp.querylens[seqnum] = strlen(dnaseq);
	//printf("%d\t%d\t%s\n", seqnum, strlen(dnaseq), dnaseq);
	seqnum++;
      }
      
      eachline[strlen(eachline)-1] = '\0';
      dnaseq[0] = '\0';

      matchprocessinfotmp.queryids[seqnum] = malloc(sizeof(Uchar) * strlen(eachline));
      Ushort chari = 0;
      for (chari = 1; chari < strlen(eachline); ++chari) {
	if (eachline[chari] == ' ') {
	  eachline[chari] = '\0';
	  break;
	}

      }
      //strcpy(matchprocessinfotmp.queryids[seqnum], eachline+1);
      strcpy(seqid, eachline+1);

    }
    
    else if (eachline[0] == '\n' || eachline[0] == ' ')
      continue;
    
    else {

      Ushort chari = 0;
      for (chari = 0; chari < strlen(eachline) - 1; ++chari) {
	eachline[chari] = tolower(eachline[chari]);
      }
      eachline[chari] = '\0';
      strcat(dnaseq, eachline);

      /*
      eachline[strlen(eachline)-1] = '\0';
      matchprocessinfotmp.querylens[seqnum] = strlen(eachline);
      //printf("%d\n", seqnum);
      matchprocessinfotmp.queryseqs[seqnum] = malloc(sizeof(Uchar) * strlen(eachline) + 1);

      Ushort chari = 0;
      for (chari = 0; chari < strlen(eachline); ++chari) {
	eachline[chari] = tolower(eachline[chari]);
      }

      strcpy(matchprocessinfotmp.queryseqs[seqnum], eachline);
      */

      
      
      //++seqnum;
    }
  }

  strcpy(matchprocessinfotmp.queryids[seqnum], seqid);
  matchprocessinfotmp.queryseqs[seqnum] = malloc(sizeof(Uchar) * strlen(dnaseq) + 1);
  strcpy(matchprocessinfotmp.queryseqs[seqnum], dnaseq);
  matchprocessinfotmp.querylens[seqnum] = strlen(dnaseq);
  //printf("%d\t%d\t%s\n", seqnum, strlen(dnaseq), dnaseq);  
  
  fclose(fp);
  
  fprintf(stderr,"finished. (%d seconds)\n", time(NULL) - seconds);
  //DELETEMEMORYMAP(filecontent);
  /***********************************************************************************************/


  

  // output the concatenated reference sequences
  /*
  char outrefseqfile[200];
  strcpy(outrefseqfile, mmcallinfo->outprefix);
  FILE *outrefseq = fopen(strcat(outrefseqfile, ".ref"), "w");
  fprintf(outrefseq, "%s\n", subjectmultiseq->sequence);
  fclose(outrefseq);
  */
  /***********************************************************************************************/

  


  // create file for storing mapping information

  char outmapfile[200];
  strcpy(outmapfile, mmcallinfo->outprefix);
  matchprocessinfotmp.outmap = fopen(strcat(outmapfile, ".map"), "w");
  if (matchprocessinfotmp.outmap == NULL) {
    fprintf(stderr,"Cannot open output file %s.\n", mmcallinfo->outprefix);
    return 0;
  }
  /***********************************************************************************************/
  
  


  // everything is ready, now create threads and do the mapping

  fprintf(stderr,"# map query sequences to references ... ");
  seconds = time(NULL);
  const int NUMP = mmcallinfoloc -> nump; // # processors
  pthread_t threads[NUMP];                // create a set of threads
  int thread_args[NUMP];                  // argument for each thread
  int rc;
  int k = 0;

  for (k=0; k<NUMP; ++k) {
    thread_args[k] = k;                   // argument is thread number
    rc = pthread_create(&threads[k], NULL, procmaxmatches_pt, (void*) &thread_args[k]);
  }

  for (k=0; k<NUMP; ++k) {
    rc = pthread_join(threads[k], NULL);
    assert(0 == rc);
  }
  fprintf(stderr,"finished. (%d seconds)\n", time(NULL) - seconds);
  /***********************************************************************************************/


  

  fclose(matchprocessinfotmp.outmap);
  freemultiseq(&matchprocessinfotmp.querymultiseq);
  freestree(&streetmp);

  return 0;
}


static int iindexall = 0;
static pthread_mutex_t plock;

void *procmaxmatches_pt(void *num) {

  const int pnum = *((int *) num); // thread number
  Sint retcode;                    // code returned from calling functions

  


  // initialize match process info, every process gets its own copy
  /******************************************************************************/
  Matchprocessinfo matchprocessinfo;
  matchprocessinfo.stree = streetmp;
  matchprocessinfo.subjectmultiseq = subjectmultiseq;
  matchprocessinfo.minmatchlength = mmcallinfo->minmatchlength;
  matchprocessinfo.maxmismatch = mmcallinfo->nummaxmis;
  matchprocessinfo.minsimpct = mmcallinfo->minpct;
  matchprocessinfo.querymultiseq = matchprocessinfotmp.querymultiseq;
  matchprocessinfo.outmap = matchprocessinfotmp.outmap;
  matchprocessinfo.reflocnum = matchprocessinfotmp.reflocnum;

  matchprocessinfo.queryseqs = matchprocessinfotmp.queryseqs;
  matchprocessinfo.queryids = matchprocessinfotmp.queryids;
  matchprocessinfo.querylens = matchprocessinfotmp.querylens;
  matchprocessinfo.numqseqs = matchprocessinfotmp.numqseqs;


  
  // avoid redudant seeds
  Uint i = 0;
  for (i = 0; i < 10000; i++)
    matchprocessinfo.hitsloc[i] = -1;
  
  // max length fasta title, including ">"
  retcode = getmaxdesclen(subjectmultiseq);
  if(retcode < 0)
    exit(-1);
  matchprocessinfo.maxdesclength = (Uint) retcode;

  // ref loc info
  for (i = 0; i < MAXREFNUM*2; ++i)
    matchprocessinfo.refloc[i] = matchprocessinfotmp.refloc[i];
  
  // stores information of top hits
  for (i = 0; i < 11; ++i)
    matchprocessinfo.tophits[i] = (char*) malloc(sizeof(char)*1000000);
  /******************************************************************************/




  // now do the mapping
  /******************************************************************************/
  Uchar *seq, *start, *end;
  Multiseq *multiseq = &matchprocessinfo.querymultiseq;
  seq = multiseq->sequence;

  Uint numseq = matchprocessinfo.numqseqs;
  Uint iindex = 0;

  //printf("%s\t%d\n", matchprocessinfo.queryseqs[0], matchprocessinfo.querylens[0]);
  

  while(1) {
      
    // iindexall keeps track of # sequences that have been processed by all processors
    // iindex is the # of sequence currently under processing in this thread
    pthread_mutex_lock(&plock);
    iindex = iindexall++;
    pthread_mutex_unlock(&plock);
    
    // all sequences have been analyzed
    if (iindex >= numseq) 
      break;

    matchprocessinfo.strand = 1;
    findmaxmatches(&matchprocessinfo.stree,
		   matchprocessinfo.minmatchlength,
		   &matchprocessinfo,
		   matchprocessinfo.queryseqs[iindex],
		   matchprocessinfo.querylens[iindex],
		   iindex);
    
    matchprocessinfo.strand = -1;
    revcom(matchprocessinfo.queryseqs[iindex], matchprocessinfo.querylens[iindex]);
    findmaxmatches(&matchprocessinfo.stree,
		   matchprocessinfo.minmatchlength,
		   &matchprocessinfo,
		   matchprocessinfo.queryseqs[iindex],
		   matchprocessinfo.querylens[iindex],
		   iindex);
    
  }
  /******************************************************************************/


  // release memory
  for (i = 0; i < 11; i++)
    free(matchprocessinfo.tophits[i]);
  
}

