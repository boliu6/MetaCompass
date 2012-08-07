#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "spacedef.h"
#include "protodef.h"
#include "debugdef.h"
#include "maxmatdef.h"
#include "fasta.h"


/*
  For each sequence in the Multiseq struct, the dynamic
  array startdesc stores the positions in the dynamic array
  descspace where the sequence description starts. The 
  following macro checks if enough memory has been allocated
  for startdesc. If not, then this is done by incrementing
  the size of the array by 128 entries. Finally, the appropriate
  entry in startdesc is assigned the correct value.
*/
#define STORESTARTDESC\
        if(multiseq->numofsequences >= allocatedstartdesc)\
        {\
          allocatedstartdesc += 128;\
          multiseq->startdesc\
            = ALLOCSPACE(multiseq->startdesc,unsigned long,allocatedstartdesc);\
        }\
        multiseq->startdesc[multiseq->numofsequences]\
          = multiseq->descspace.nextfreeUchar


/*
  The following function scans a string containing the content of
  a multiple fasta formatted file. The parameter are as follows:
  1, multiseq is the Multiseq-record to store the scanned information.
  2, filename is the information from which the file contents was read.
  3, replacewildcardchar is the character used to
     replace a wildcard (then it should be different from the
     characters occuring in DNA sequences) or 0 if wildcards are not replaced. 
  4, input points to the inputstring to be scanned,
  5, inputlen is the length of the input.
  
  Each sequence description begins with the symbol >.
  If it does, then this symbol is skipped. The rest of the 
  line up to the first white space character is stored in descspace.
  Otherwise, the rest of the line is discarded. 
  The remaining lines (until the next symbol > or the end of the input string)
  are scanned for alphanumeric characters which make up the sequence.
  White spaces are ignored. Upper case characters are transformed to lower case.
  The input string must contain at least one sequence.
  If an error, an negative error code is returned.
  If success, the return code is 0.

  input points to the memory address that stores the raw fasta file data.
  This memory space is then used to store the concatenated reference sequences.
*/
signed long scanmultiplefastafile (Multiseq *multiseq,
				   char *filename,
				   Uchar replacewildcardchar,
				   Uchar *input,
				   unsigned long inputlen) {

  
  Uchar *inputptr,              // read each character of input
        *newptr,                // points to the transformed
        tmpchar;                // temporary character
  Uint allocatedstartdesc = 0;  // # characters allocated for startdesc
  BOOL indesc = False,          // inside description part of sequence
       copydescription = False; // currently copying the description 


  initmultiseq (multiseq);
  multiseq->originalsequence = NULL;

  newptr = multiseq->sequence = input;

  // read each character
  for (inputptr = input; inputptr < input + inputlen; inputptr++) {

    // status: insisde of fasta title line
    if (indesc) {

      // status: copying the fasta ID
      if(copydescription) {

	// status: current character is space
        if(isspace((int) *inputptr)) {

	  // stop copying fasta ID
	  copydescription = False;
          STOREINARRAY (&multiseq->descspace, 
                        Uchar,
                        4096,
                        (Uchar) '\n');
        }

	// status: current character is NOT space
	// then copy this character into descspace
	else {
          STOREINARRAY (&multiseq->descspace, 
                        Uchar,
                        4096,
                        *inputptr);
        }
      }
      
      // status: end of fasta title line
      if (*inputptr == '\n')
        indesc = False;
    }

    // status: outside of fasta title line
    else {

      // start of fasta title
      if ('>' == *inputptr) {
	
	// store the start index of each stored fasta ID
        STORESTARTDESC;
        if (multiseq->numofsequences > 0) {
	  
          STOREINARRAY (&multiseq->markpos, 
                        Uint,
                        128,
                        (unsigned long) (newptr - multiseq->sequence));

          *newptr++ = SEPARATOR; // separator is UCHAR_MAX
        }
        multiseq->numofsequences++;
        indesc = True;
        copydescription = True;
      }

      // status: fasta sequences
      else {

	// current character
        tmpchar = *inputptr;

	if (!isspace((int) tmpchar)) { 	// test if it's space
	  
          tmpchar = (Uchar) tolower((int) tmpchar); // convert to lower case

	  // replace unknown characters or not
          if (replacewildcardchar != 0) { 
            switch (tmpchar) {
              case 'a':
              case 'c':
              case 'g':
              case 't':
                break;
              default:
                tmpchar = replacewildcardchar; // set other character to 'x'
            }
          }
          *newptr++ = tmpchar; // store this character
        }
      }
    }
  } // end of reading file
  
  STORESTARTDESC;
  if (multiseq->numofsequences == 0)
  {
    ERROR0 ("no sequences in multiple fasta file");
    return -2;
  }
  multiseq->totallength = (unsigned long) (newptr - multiseq->sequence);
  if(multiseq->totallength == 0)
  {
    ERROR0 ("empty sequence in multiple fasta file");
    return -3;
  }

  *newptr = '\0'; // set last character to nothing
  
  return 0;
}


/*
  Load in the subject sequence file.
  Deliver the parsed multiple sequences in the corresponding
  Multiseq object. The files are read via memory mapping.
  If an error occurs, then the function delivers a negative error code.
  Otherwise the error code is 0.
*/
signed long loadinrefseq (Multiseq *subjectmultiseq,
			    char *subjectfile)
{
  unsigned long filelen;
  Uchar *filecontent;

  filecontent = CREATEMEMORYMAP (subjectfile, True, &filelen);
  if (filecontent == NULL || filelen == 0)
  {
    ERROR2("cannot open file \"%s\" or file \"%s\" is empty", subjectfile, subjectfile);
    return -1;
  }

  if (scanmultiplefastafile (subjectmultiseq, subjectfile,
                             MMREPLACEMENTCHARSUBJECT,     // replace other characters with 'x'
                             filecontent, 
                             filelen) != 0)
  {
    return -2;
  }
  return 0;
}
