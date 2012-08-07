#include <stdio.h>
#include <ctype.h>
#include <string.h>

#include "streedef.h"
#include "spacedef.h"

#include "maxmatdef.h"
// SEPARATOR

#include "procmaxmat.h"

#include "protodef.h"


/*
  This file contains functions to compute maximal matches of some
  minimum length between the subject-sequence and the query-sequence.
  This is done by a traversal of the suffix tree of the subject-sequence.

  For each suffix, say (s), of the query, the location ploc
  of some prefix (p) of (s) is determined. 
  Let pmax be the longest prefix
  of (s) that occurs as a substring of the subject-sequence.
  (p) is determined as follows.

  1, If the length of pmax is <= minmatchlength, 
  then p=pmax.

  2, If the length of pmax is > minmatchlength, then
  (p) is the prefix of pmax of length minmatchlength.

  3, Given ploc, the location maxloc of pmax is 
  determined, thereby keeping track of the branching nodes 
  visited during this matching step. Finally, the suffix tree below location
  ploc is traversed in a depth first strategy.
  This delivers the set of suffixes representing
  a match of length at least minmatchlength against a prefix of (s).

  4, Using a stack to keep track of the length of the longest common prefix 
  of each encountered suffix and (s).

  5, Finally, it is checked whether the
  match between the match is maximal by comparing the characters
  to the left and to the right of the two instances of the match.
*/


/*
  For the depth first traversal we need a stack containing elements
  of the following type. Each stack element stores information for
  a branching node of the suffix tree. The top of the stack corresponds
  to the branching node, say (b), currently visited.
  querycommondepth is the length of the longest common prefix
  of (s) and all suffixes represented by leaves in the subtree below 
  (b). onmaxpath is true iff the sequence 
  corresponding to (b) is a prefix of the query.
*/




typedef struct
{
  Uint querycommondepth;
  BOOL onmaxpath;
} Nodeinfo;




// The stack is represented by a dynamic array of elements of type Nodeinfo.
DECLAREARRAYSTRUCT(Nodeinfo);




typedef struct {
  
  Suffixtree *stree;              // reference to suffix tree of subject-seq

  ArrayNodeinfo commondepthstack; // stack to store depth values

  ArrayPathinfo matchpath;        // path of br. nodes from ploc to maxloc

  Location maxloc;                // location of \texttt{pmax}

  Uchar *query,                   // the query string
        *querysuffix;             // current suffix of query

  Uint querylen,                  // length of the current query
       queryseqnum,               // number of query sequence
       minmatchlength,            // min length of a match to be reported
       depthofpreviousmaxloc;     // the depth of the previous maxloc

  Matchprocessinfo *processinfo;  // first arg. when calling previous function
  
} Maxmatchinfo;




/*
  The following function is applied to each leaf visited during
  the depth first traversal. It checks if corresponding match is
  left maximal. In this case, the length of the longest common
  prefix of this suffix and s is computed.
*/
static pthread_mutex_t plockout;
static Sint processleaf(Uint seedlocinr,/*@unused@*/ Bref lcpnode,void *info) {
  

  Maxmatchinfo *maxmatchinfo = (Maxmatchinfo *) info;



  // test if this is a MUM, and can be used as a seed

  if(seedlocinr == 0 ||
     maxmatchinfo->query == maxmatchinfo->querysuffix ||
     maxmatchinfo->stree->text[seedlocinr - 1] != 
     *(maxmatchinfo->querysuffix-1)) {}
  else {
    return 0;
  }
  /*******************************************************************************/
  


  
  // compute length of seed

  Uint seedlen;
  if(maxmatchinfo->commondepthstack.nextfreeNodeinfo == 0)
    seedlen = maxmatchinfo->maxloc.locstring.length;
  else {
    Nodeinfo *father = maxmatchinfo->commondepthstack.spaceNodeinfo + 
      maxmatchinfo->commondepthstack.nextfreeNodeinfo-1;
    if(father->onmaxpath &&
       LEAFADDR2NUM(maxmatchinfo->stree,
		    maxmatchinfo->maxloc.nextnode.address) == seedlocinr)
      seedlen = maxmatchinfo->maxloc.locstring.length;
    else
      seedlen = father->querycommondepth;
  }
  /*******************************************************************************/


  
  
  Matchprocessinfo *processinfo = maxmatchinfo->processinfo;
  Uint seqnum = maxmatchinfo->queryseqnum;
  Ushort readlen = processinfo->querylens[seqnum];

  
  Uint seedlocinq = (Uint) (maxmatchinfo->querysuffix - maxmatchinfo->query);   // start of seed in query read
  int readstartinr = seedlocinr - seedlocinq;                                   // start of read in ref string

  // readstartinr will be changed because of insertions and delections
  int readstartinrold = readstartinr;                        
  int readendinr = readstartinr + readlen - 1;                                // end of read in ref

  
  // check if end of read is out of the reference range
  if (readstartinr < 3 || readendinr >= processinfo->subjectmultiseq->totallength) { return 0;}

  int i = 0, j = 0;
  Uint maxmismatch = processinfo->maxmismatch;                 // cutoff for # mismatches
  Uint minsimpct = processinfo->minsimpct;                     // cutoff for % similarity
  
  


  // one mapping between one read and the reference may have more than one seed, and we only need one
  // if there are insertions, they do not match up then we'll check later after dynamic programming
  for (i = 0; i < 10000; i++) {

    //  all hitsloc has been scanned, this seed is not redundant
    if (processinfo->hitsloc[i] == -1) {
      processinfo->hitsloc[i] = readstartinr;         // read start loc in ref string
      break;
    }
    else if (processinfo->hitsloc[i] == readstartinr) // this seed is redundant
      return 0;
    
  }
  if (i == 10000) { return 0;}
  /*******************************************************************************/
  

  //printf("%d\t%d\n", seedlocinq, seedlocinr);

  // align surrounding regions of seed between query read and ref
  Uint mism = 0;   // number of mismatches
  Uint misloc[10]; // locs of mismatches
  char snps[20];  // mismatches, e.g., a->t


  
  /***********************************************************************************/
  /* align right region with dynamic programming */
  /***********************************************************************************/
  char *query = processinfo->queryseqs[seqnum] + seedlocinq + seedlen; // start of right region of query
  char *ref = processinfo->subjectmultiseq->sequence + seedlocinr + seedlen;                  // start of right region of ref
  Ushort qlen = readlen - seedlocinq - seedlen;                // length of right region of query
  Ushort rlen = qlen + 3;                                      // allow max 3 insertions, so ref is 3bp longer
  //fprintf(stderr, "%d\n", seqnum);

  // store dynamic programming
  if (rlen >= 200 || qlen >= 200) {
    return 0;
  }
  int dyntabright[rlen+1][qlen+1];                // alignment scores
  unsigned short dynstateright[rlen+1][qlen+1];   // paths: 1 left, 3 above, 0 diagnal match, 2 diagnal mismatch


  // allow max 3 gaps
  // if 3 gaps are in ref, then ref is 6bp longer than query

  //printf("initialize dp table\n");
  /***********************************************************************************/
  /* initialize dp table */
  /***********************************************************************************/
  dyntabright[0][0] = 0;
  for (j = 1; j <= qlen; ++j) {
    // row 0, columns 1, 2, 3
    if (j <= 3) {
      //printf("%d\t%d\n", seqnum, qlen);
      //fprintf(stderr, "%d\t%d\t%d\n", seqnum, j, qlen);
      dyntabright[0][j] = GAPOPEN + (j-1)*GAPEXTEND;
      dynstateright[0][j] = 1;
    }
    // 3bp off diagnal
    else {
      dyntabright[j-4][j] = -1000;
      dynstateright[j-4][j] = 2;
    }
  }
  
  for (i = 1; i <= rlen; ++i) {
    // column 0, row 1, 2, 3
    if (i <= 3) {
      dyntabright[i][0] = GAPOPEN + (i-1)*GAPEXTEND;
      dynstateright[i][0] = 3;
    }
    // 3bp off diagnal
    else {
      dyntabright[i][i-4] = -1000;
      dynstateright[i][i-4] = 2;
    }
    //printf("%d\t%d\n", i, i-4);
    //printf("qlen rlen%d\t%d\n", qlen, rlen);
  }
  
  //printf("fill dp table\n");
  //printf("%d\t%d\n", qlen, rlen);
  /*******************************************************************************/
  /* fill dp table, only 3bp off diagnal */
  /*******************************************************************************/
  for (j = 1; j <= qlen; ++j) {  // for each column

    /*
    if ((char) ref[j] == (char) SEPARATOR) { // crossed the boundary of ref seq
      //return 0;
    }
    */

    int i = j-3 > 1 ? j-3 : 1;
    for (; i <= j+3 && i <= rlen; ++i) {
      
      // 3 possible scores
      int s1 = (dynstateright[i][j-1] != 1 ? GAPOPEN : GAPEXTEND) + dyntabright[i][j-1];    // from left
      int s2 = (ref[i-1] == query[j-1] ? BASEMATCH : BASEMISMATCH) + dyntabright[i-1][j-1]; // from diagnal
      int s3 = (dynstateright[i-1][j] != 3 ? GAPOPEN : GAPEXTEND) + dyntabright[i-1][j];    // from above
      
      // find the best
      int bests = 0;
      if (s2 > s1) {
	if (s2 > s3) {  // diagnal is the best
	  bests = s2;
	  dynstateright[i][j] = ref[i-1] == query[j-1] ? 0 : 2;
	}
	else {
	  bests = s3;    // above is the best
	  dynstateright[i][j] = 3;
	}
      }
      else {
	if (s1 >= s3) {  // left is the best
	  bests = s1;
	  dynstateright[i][j] = 1;
	}
	else {           // above is the best
	  bests = s3;
	  dynstateright[i][j] = 3;
	}
      }
      dyntabright[i][j] = bests;
      //printf("%c\t%c\n", ref[i-1], query[j-1]);
      //printf("%d\t%d\t%d\n", i, j, bests);
      //printf("%d\t%d\t%d\n\n", s1, s2, s3);
    }
  }
  /*******************************************************************************/

  

  //printf("find best score\n");
  //printf("%d\t%d\n", qlen, rlen);
  /*******************************************************************************/
  /* find the best alignment score, and loc in the last column */
  /*******************************************************************************/
  // locus of best score
  i = rlen;
  j = qlen;
  int alnscore = -1000;
  int k = j-3 > 0 ? j-3: 0;
  for (; k <= rlen; ++k) {
    if (dyntabright[k][qlen] >= alnscore) {
      alnscore = dyntabright[k][qlen];
      i = k;
    }
  }
  /*******************************************************************************/

  

  //printf("trace back\n");
  /*******************************************************************************/
  /* trace back to build alignment */
  /*******************************************************************************/
  //printf("%d\t%d\n", i, j);
  while (i > 0 || j > 0) {
    //printf("%d\t%d\n", j, seedlen);
    // from diagnal, match
    if (dynstateright[i][j] == 0) {
      --i;
      --j;
    }

    // from diagnal, mismatch
    else if (dynstateright[i][j] == 2) {
      misloc[mism] = j - 1 + seedlocinq + seedlen;
      snps[mism*2] = query[j-1];
      snps[mism*2+1] = ref[i-1];
      mism++;
      --i;
      --j;
    }

    // from left, insertion in query, gap in ref
    else if (dynstateright[i][j] == 1) {
      misloc[mism] = j - 1 + seedlocinq + seedlen;
      snps[mism*2] = query[j-1];
      snps[mism*2+1] = '-';
      mism++;
      --j;
      --readendinr;
    }

    // from above
    else if (dynstateright[i][j] == 3) {
      misloc[mism] = j - 1 + seedlocinq + seedlen;
      snps[mism*2] = '-';
      snps[mism*2+1] = ref[i-1];
      mism++;
      --i;
      ++readendinr;
    }

    else {
      fprintf(stderr, "something wrong right %d\t%d\t%d\n", i, j, readstartinr);
      exit(1);
    }

    // test if within criteria
    if (mism > maxmismatch || (mism*1.0 / readlen) * 100.0 > 100.0 - minsimpct)
      return 0;

  }
  // end of align right region




  /* processing left region */
  /***********************************************************************************/
  query = processinfo->queryseqs[seqnum];
  ref = processinfo->subjectmultiseq->sequence + seedlocinr - seedlocinq;
  qlen = seedlocinq;
  rlen = qlen;
  int offset = 3;
  if (readstartinr >= 3) {
    offset = 6;
    rlen += 3;
    ref -= 3;
    // if crossed boundary, do not include additional bases
    /*
    if ((char) *(ref-1) != (char) SEPARATOR &&
	(char) *(ref-2) != (char) SEPARATOR &&
	(char) *(ref-3) != (char) SEPARATOR) {

      rlen += 3;
      ref -= 3;
      
    }
    */

  }

  if (rlen >= 200 || qlen >= 200) {
    return 0;
  }
  int dyntableft[rlen+1][qlen+1];
  unsigned short dynstateleft[rlen+1][qlen+1]; // 1 left, 3 above, 0 diagnal match, 2 diagnal mismatch
  
  
  /* initialize dp table, left */
  /***********************************************************************************/
  dyntableft[0][0] = 0;
  for (j = 1; j <= qlen; ++j) {

    dyntableft[j-1][j] = -1000;
    dynstateleft[j-1][j] = 2;
    
  }
  
  for (i = 1; i <= rlen; ++i) {

    // column 0, row 1, 2, 3, 4, 5, 6
    if (i <= offset) {
      dyntableft[i][0] = 0;
      dynstateleft[i][0] = 3;
    }
    
    // 3bp off diagnal
    else {
      dyntableft[i][i-7] = -1000;
      dynstateleft[i][i-7] = 2;
    }
  }
  

  
  /*******************************************************************************/
  /* fill dp table, only 3bp off diagnal, left */
  /*******************************************************************************/
  for (j = 1; j <= qlen; ++j) {  // for each column

    int i = j;
    for (; i <= j+6 && i <= rlen; ++i) {

      // 3 possible scores
      int s1 = (dynstateleft[i][j-1] != 1 ? GAPOPEN : GAPEXTEND) + dyntableft[i][j-1];     // from left
      int s2 = (ref[i-1] == query[j-1] ? BASEMATCH : BASEMISMATCH) + dyntableft[i-1][j-1]; // from diagnal
      int s3 = (dynstateleft[i-1][j] != 3 ? GAPOPEN : GAPEXTEND) + dyntableft[i-1][j];     // from above

      // find the best
      int bests = 0;
      if (s2 > s1) {
	if (s2 > s3) {  // diagnal is the best
	  bests = s2;
	  dynstateleft[i][j] = ref[i-1] == query[j-1] ? 0 : 2;
	}
	else {
	  bests = s3;    // above is the best
	  dynstateleft[i][j] = 3;
	}
      }
      else {
	if (s1 >= s3) {  // left is the best
	  bests = s1;
	  dynstateleft[i][j] = 1;
	}
	else {           // above is the best
	  bests = s3;
	  dynstateleft[i][j] = 3;
	}
      }
      dyntableft[i][j] = bests;
      //printf("%c\t%c\n", ref[i-1], query[j-1]);
      //printf("%d\t%d\t%d\n", i, j, bests);
      //printf("%d\t%d\t%d\n\n", s1, s2, s3);
    }
  }
  /*******************************************************************************/


  /*
  for (i = 0; i <= rlen; ++i) {
    for (j = 0; j <= qlen; ++j) {
      printf("%d\t", dyntableft[i][j]);
    }
    printf("\n");
  }

  for (i = 0; i <= rlen; ++i) {
    for (j = 0; j <= qlen; ++j) {
      printf("%d\t", dynstateleft[i][j]);
    }
    printf("\n");
  }
  */

  
  /* find the best alignment score, and loc in the last column, left */
  /*******************************************************************************/
  // locus of best score
  i = rlen;
  j = qlen;
  /*******************************************************************************/
  //printf("best %d\t%d\t%d\n\n", i, j, 1);
  

  /*******************************************************************************/
  /* trace back to build alignment, left */
  /*******************************************************************************/
  while (j != 0) {

    //if ((char) ref[i] == (char) SEPARATOR) { // crossed the boundary of ref seq
      //return 0;
    //}

    //printf("%d\t%d\t%d\t%d\n\n", i, j, dynstateleft[i][j], dyntableft[i][j]);
    // from diagnal, match
    if (dynstateleft[i][j] == 0) {
      --i;
      --j;
    }

    // from diagnal, mismatch
    else if (dynstateleft[i][j] == 2) {
      misloc[mism] = j - 1;
      snps[mism*2] = query[j-1];
      snps[mism*2+1] = ref[i-1];
      mism++;
      --i;
      --j;
    }

    // from left, insertion in query, gap in ref
    else if (dynstateleft[i][j] == 1) {
      ++readstartinr;
      misloc[mism] = j - 1;
      snps[mism*2] = query[j-1];
      snps[mism*2+1] = '-';
      mism++;
      --j;
    }

    // from right
    else if (dynstateleft[i][j] == 3) {
      --readstartinr;
      misloc[mism] = j - 1;
      snps[mism*2] = '-';
      snps[mism*2+1] = ref[i-1];
      mism++;
      --i;
    }

    else {
      fprintf(stderr, "something wrong right %d\t%d\t%d\n", i, j, readstartinr);
      exit(1);
    }

    // test if within criteria
    if (mism > maxmismatch || (mism*1.0 / readlen) * 100.0 > 100.0 - minsimpct)
      return 0;
  }
  // end of align right region
  // end of aligning left region

  

  
  /*
    after dynamic programming readstartinr may have changed
  */
  if (readstartinr != readstartinrold) {
    for (i = 0; i < 10000; i++) {
      if (processinfo->hitsloc[i] == -1) {          // this seed is not redundant
	processinfo->hitsloc[i] = readstartinr;         // read start loc in ref string
	break;
      }
      else if (processinfo->hitsloc[i] == readstartinr) // this seed is redundant
	return 0;
    }
  }
  /*******************************************************************************/

  

  
  /*
    if mappings with fewer mismatches have been observed, ignore this one
  */
  for (i = 0; i < mism; i++) {
    if (*processinfo->tophits[i] != '\0') {
      return 0;
    }
  }
  /*******************************************************************************/

  


  /*
    array to store output strings for each mismatches
  */
  char *outtmp = processinfo->tophits[mism];
  int refseqnum = 0;
  int *locptr = processinfo->refloc;
  int reflocnum = processinfo->reflocnum;
  int tmps = 0;


  // binary search
  int locs = readstartinr;
  int loce = readendinr;
  refseqnum = -2;
  int ranges = 0;
  int rangee = (reflocnum - 1)*2;

  while (refseqnum == -2) {

    if (rangee == ranges) {                              // only one element left
	  
      int refs = locptr[ranges];
      int refe = locptr[ranges+1];
      if (locs >= refs && locs <= refe) {           // head is within
	if (loce <= refe - 1) {                         // tail is also within
	  tmps = refs;
	  refseqnum = ranges/2;
	}
	else {                                           // tail is NOT within
	  refseqnum = -1;
	}
      }
      else {
	refseqnum = -1;
      }
	  
    }

    else {
      int mid = (rangee + ranges) / 2;
      if (mid % 2 == 1) {                                 // odd number of elements
	mid--;
      }

      int refs = locptr[mid];
      int refe = locptr[mid+1];
      if (locs >= refs && locs <= refe) {           // head is within
	if (loce <= refe - 1) {                         // tail is also within
	  tmps = refs;
	  refseqnum = mid/2;
	}
	else {                                           // tail is NOT within
	  refseqnum = -1;
	}
      }
      else if (locs < locptr[mid]) {
	rangee = mid - 2;
      }
      else {
	ranges = mid + 2;
      }
    }
	
  }
  if (refseqnum == -1)
    return 0;
  /***************************************************************************************/




  // print query read ids
  sprintf(outtmp + strlen(outtmp), "%s\t", processinfo->queryids[seqnum]);

  // print reference read ids
  Uint descs = processinfo->subjectmultiseq->startdesc[refseqnum];
  Uint desce = processinfo->subjectmultiseq->startdesc[refseqnum+1];
  Uchar *descref = processinfo->subjectmultiseq->descspace.spaceUchar + descs;
  for(i=0; i< desce - descs + 1; i++) {
    if(isspace((int) descref[i]))
      break;
    sprintf(outtmp + strlen(outtmp), "%c", descref[i]);
  }
  strcat(outtmp, "\t");
      

  // strand
  if (processinfo->strand == 1)
    strcat(outtmp, "+\t");
  else 
    strcat(outtmp, "-\t");


  // print out additional information
  sprintf(outtmp + strlen(outtmp), "%d\t%.2f\t%d\t%d\t", readstartinr - tmps +1, 100.0 - mism*100.0/readlen*1.0, mism, readlen);
  for (i = mism - 1; i >= 0; i--) {
            
    switch (snps[i*2]) {
    case 'a': break;
    case 't': break;
    case 'g': break;
    case 'c': break;
    case '-': break;
    case (char) SEPARATOR:
      snps[i*2] = '-';
      break;
    default:
      snps[i*2] = 'n';
      break;
    }

    switch (snps[i*2+1]) {
    case 'a': break;
    case 't': break;
    case 'g': break;
    case 'c': break;
    case '-': break;
    case (char) SEPARATOR:
      snps[i*2+1] = '-';
      break;
    default:
      snps[i*2+1] = 'n';
      break;
    }

    sprintf(outtmp + strlen(outtmp), "%d%c%c;", misloc[i]+1, snps[i*2], snps[i*2+1]);
  }
  sprintf(outtmp + strlen(outtmp), "\n");
  return 0;
}




/*
  The following functions inherits information from the maximal matchpath
  to the stack used in the depth first traversal.
*/
static void inheritfrompath(ArrayPathinfo *matchpath,Location *maxloc,
                            Nodeinfo *stacktop,Bref nodeptr,
                            Uint accessindex,
                            Uint inheritdepth)
{
  if(accessindex > matchpath->nextfreePathinfo)
  {
    stacktop->onmaxpath = False;
    stacktop->querycommondepth = inheritdepth;
  } else
  {
    if(accessindex == matchpath->nextfreePathinfo)
    {
      if(maxloc->remain == 0)
      {
        if(maxloc->nextnode.address == nodeptr)
        {
          stacktop->onmaxpath = True;
          stacktop->querycommondepth = maxloc->locstring.length;
        } else
        {
          stacktop->onmaxpath = False;
          stacktop->querycommondepth = inheritdepth;
        }
      } else
      {
        stacktop->onmaxpath = False;
        if(maxloc->nextnode.address == nodeptr)
        {
          stacktop->querycommondepth = maxloc->locstring.length;
        } else
        {
          stacktop->querycommondepth = inheritdepth;
        }
      }
    } else
    {
      if(matchpath->spacePathinfo[accessindex].ref == nodeptr)
      {
        stacktop->onmaxpath = True;
        stacktop->querycommondepth 
          = matchpath->spacePathinfo[accessindex].depth;
      } else
      {
        stacktop->onmaxpath = False;
        stacktop->querycommondepth = inheritdepth;
      }
    }
  }
}

/*
  The following function is called whenever during a depth first traversal
  of a subtree of the suffix tree, each time
  a branching node is visited for the first time.
  The arguments are the branching node \texttt{nodeptr} and 
  a pointer \texttt{info} to some information passed to the function 
  \texttt{depthfirststree}. If the \texttt{commondepthstack}
  is empty or the father of the current node is on the maximal path,
  then the commondepthstack inherits information from the appropriate
  value of the maximal match path. Otherwise, the information about
  the maximal matching prefix length is propagated down the stack.
  The function always return \texttt{True} and thus the depth first
  traversal continues.
*/
static BOOL processbranch1(Bref nodeptr,void *info)
{
  Maxmatchinfo *maxmatchinfo = (Maxmatchinfo *) info;
  Nodeinfo *stacktop, 
           *father;
  GETNEXTFREEINARRAY(stacktop,&maxmatchinfo->commondepthstack,Nodeinfo,32);
  if(stacktop == maxmatchinfo->commondepthstack.spaceNodeinfo)
  {
    inheritfrompath(&maxmatchinfo->matchpath,
                    &maxmatchinfo->maxloc,
                    stacktop,
                    nodeptr,
                    0,
                    maxmatchinfo->minmatchlength);
  } else
  {
    father = stacktop-1;
    if(father->onmaxpath)
    {
      inheritfrompath(&maxmatchinfo->matchpath,
                      &maxmatchinfo->maxloc,
                      stacktop,
                      nodeptr,
                      maxmatchinfo->commondepthstack.nextfreeNodeinfo-1,
                      father->querycommondepth);
    } else
    {
      stacktop->onmaxpath = False;
      stacktop->querycommondepth = father->querycommondepth;
    }
  }
  return True;
}

/*
  The following function is called whenever a branching node 
  \texttt{nodeptr}
  is visited for the second time (i.e.\ the entire subtree below 
  \texttt{nodeptr} has been processed).
  \texttt{info} is a pointer to some information passed to the function
  \texttt{depthfirststree}.
*/

static Sint processbranch2(/*@unused@*/ Bref nodeptr,void *info)
{
  Maxmatchinfo *maxmatchinfo = (Maxmatchinfo *) info;

  maxmatchinfo->commondepthstack.nextfreeNodeinfo--;
  return 0;
}


/*
  This function computes the maximal matches below location ploc.
  All global information is passed via the Maxmatchinfo struct.
  At first the number rescanprefixlength is determined. This is the length of the 
  the current querysuffix that definitely match some suffix 
  of the subject sequence. Then the suffix tree is scanned starting at
  ploc to find maxloc. During this matching phase,
  all branching nodes visited are stored in matchpath.
  If ploc is a leaf location, the the corresponding leaf
  is directly processed by the function processleaf.
  Otherwise, the nextnode component of ploc
  is pushed on a stack and a depth first traveral of the 
  the suffix tree below node nextnode is performed.
*/
static Sint enumeratemaxmatches (Maxmatchinfo *maxmatchinfo,
                                 Location *ploc) {

  //printf("hi\n");
  Uint rescanprefixlength;
  maxmatchinfo->matchpath.nextfreePathinfo = 0;


  // has been scanned by previous k-mer seed before
  // this seed will be scanned by the length - 1
  if(maxmatchinfo->depthofpreviousmaxloc > UintConst(1))
    rescanprefixlength = maxmatchinfo->depthofpreviousmaxloc - 1;
  else
    rescanprefixlength = 0;

  
  // finds the max match
  (void) findprefixpathstree (maxmatchinfo->stree, 
                              &maxmatchinfo->matchpath,
                              &maxmatchinfo->maxloc, 
                              ploc,
                              maxmatchinfo->querysuffix
                                  + maxmatchinfo->minmatchlength,
                              maxmatchinfo->query 
			          + maxmatchinfo->querylen - 1,
                              rescanprefixlength);
  maxmatchinfo->depthofpreviousmaxloc = maxmatchinfo->maxloc.locstring.length;
  maxmatchinfo->commondepthstack.nextfreeNodeinfo = 0;
  
  if(ploc->nextnode.toleaf) {
    if(processleaf(LEAFADDR2NUM(maxmatchinfo->stree,ploc->nextnode.address),
                   NULL,(void *) maxmatchinfo) != 0) {
      //printf("hello\n");
      return -1;
    }
  }
  else {
    (void) processbranch1(ploc->nextnode.address,(void *) maxmatchinfo);
    if(depthfirststree(maxmatchinfo->stree,
		       &ploc->nextnode,
		       processleaf,
		       processbranch1,
		       processbranch2,
		       NULL,
		       NULL,
		       (void *) maxmatchinfo) != 0)
      return -2;
  }
  return 0;
}


/*
  This function finds all maximal matches between the 
  subject sequence and the query sequence of length at least
  minmatchlength. Initially, the function appropriately intializes the
  maxmatchinfo struct. It then scans query to find
  ploc for the longest suffix of query.
  The depth of ploc is stored in depthofpreviousmaxloc.
  In the for-loop each instance of ploc is determined
  and processed further by enumeratemaxmatches whenever its 
  depth is longer than the minimum match length.
*/
Sint findmaxmatches(Suffixtree *stree,             // suffix tree of reference 
                    Uint minmatchlength,           // minimun match length
                    Matchprocessinfo *processinfo, // how to perform the match
                    Uchar *query,                  // start loc of this read in the whole query seq string
                    Uint querylen,
                    Uint queryseqnum) {

  
  Uchar *querysubstringend;  // query + minmatchlength - 1
  Location ploc;
  Maxmatchinfo maxmatchinfo;


  //printf("%s\t%d\n", query, querylen);
  

  if(querylen < minmatchlength)
    return 0;
  

  
  // initialize maxmatchinfo
  maxmatchinfo.stree = stree;
  INITARRAY(&maxmatchinfo.commondepthstack,Nodeinfo);
  INITARRAY(&maxmatchinfo.matchpath,Pathinfo);
  maxmatchinfo.querysuffix = maxmatchinfo.query = query;
  maxmatchinfo.querylen = querylen;
  maxmatchinfo.minmatchlength = minmatchlength;
  maxmatchinfo.queryseqnum = queryseqnum;
  maxmatchinfo.processinfo = processinfo;
  querysubstringend = query + minmatchlength - 1;
  


  
  // hitsloc: start locs of reads on ref
  // used to avoid aligning redundant seeds
  int i = 0;
  for (i = 0; i < 10000; i++) {
    if (maxmatchinfo.processinfo->hitsloc[i] != -1)
      maxmatchinfo.processinfo->hitsloc[i] = -1;
    else
      break;
  }



  
  // initialize the mapping information for a read
  if (maxmatchinfo.processinfo->strand == 1) {
    int i = 0;
    for (i = 0; i < 11; i++) {
      *maxmatchinfo.processinfo->tophits[i] = '\0';
    }
  }

  

  // finds longest suffix, store info. in ploc
  // see if the first k mers in query match ref
  // if yes then, ploc.locstring.length is 18
  // defined in streesrc/streefiledoc.c
  (void) scanprefixfromnodestree(stree, &ploc, ROOT(stree), 
                                  query, querysubstringend,0);
  maxmatchinfo.depthofpreviousmaxloc = ploc.locstring.length;

  for (;
       querysubstringend < query + querylen - 1; 
       maxmatchinfo.querysuffix++, querysubstringend++) {


    
    // if previous 18mer match, then enumerate
    //printf("%d\t%d\n", ploc.locstring.length, querylen);
    if(ploc.locstring.length >= minmatchlength &&
       enumeratemaxmatches(&maxmatchinfo,&ploc) != 0) {

      //printf("1\n");
      return -1;
      
    }
    
    
    // ploc.locstring.length == 0
    // previous 18mer does not match at all
    // then search from root
    if (ROOTLOCATION (&ploc))
      (void) scanprefixfromnodestree (stree, &ploc, ROOT (stree), 
				      maxmatchinfo.querysuffix+1, 
				      querysubstringend+1,0);
    // previous 18mer partial match
    else {
      linklocstree (stree, &ploc, &ploc);
      // defined in streesrc/streefiledoc.c
      (void) scanprefixstree (stree, &ploc, &ploc,
			      maxmatchinfo.querysuffix
			      + ploc.locstring.length+1,
			      querysubstringend+1,0);
    }
  }
  
  
  // why this?
  while (!ROOTLOCATION (&ploc) && ploc.locstring.length >= minmatchlength) {

    //printf("hellono\n");
    if(enumeratemaxmatches (&maxmatchinfo,&ploc) != 0) {
      return -2;
    }
    linklocstree (stree, &ploc, &ploc);
    maxmatchinfo.querysuffix++;
  }
  

  // search strand 1 first and -1 second
  // so when we get to -1, we print out all the info.
  if (maxmatchinfo.processinfo->strand == -1) {
    int i = 0;
    for (i = 0; i < 11; i++) {

      // print out only the top best hit
      if (*maxmatchinfo.processinfo->tophits[i] != '\0') {

	// print out mapping information
	pthread_mutex_lock(&plockout);
	fprintf(maxmatchinfo.processinfo->outmap, "%s", maxmatchinfo.processinfo->tophits[i]);
	pthread_mutex_unlock(&plockout);
	break; // print out only the top best hit
	
      }
      
    }
    
  }


  FREEARRAY(&maxmatchinfo.commondepthstack,Nodeinfo);
  FREEARRAY(&maxmatchinfo.matchpath,Pathinfo);

  //printf("%s\t%d\n", query, querylen);
  return 0;
}
