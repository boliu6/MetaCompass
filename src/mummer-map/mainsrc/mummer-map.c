#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "types.h"
#include "debugdef.h"
#include "errordef.h"
#include "protodef.h"
#include "spacedef.h"
#include "megabytes.h"
#include "maxmatdef.h"
#include "fasta.h"


/*
  imported from maxmatopt.c
*/
signed long parsemaxmatoptions(MMcallinfo *maxmatcallinfo,
                         Argctype argc,
                         char **argv);


int main(int argc, char *argv[])
{
  signed long retcode;
  MMcallinfo mmcallinfo;
  Multiseq subjectmultiseq;

  //initclock();
  retcode = parsemaxmatoptions (&mmcallinfo, argc, argv);
  if (retcode < 0)
  {
    STANDARDMESSAGE;  // return error code and show message
  }
  if (retcode == 1)   // program was called with option -help
  {
    checkspaceleak ();
    mmcheckspaceleak ();
    return EXIT_SUCCESS;
  }


  // keep track of the running time of different parts of the code
  time_t seconds = time(NULL); 

  
  /*
    Read reference sequences file,
    concatenate them into one single reference string
  */
  fprintf(stderr, "# read reference sequences file ... ");
  if (loadinrefseq(&subjectmultiseq,
		   &mmcallinfo.subjectfile[0]) != 0)
  {
    fprintf(stderr,"%s: %s\n",argv[0],messagespace());
    return EXIT_FAILURE;
  }
  fprintf(stderr, "finished. (%d seconds)\n", time(NULL) - seconds);
  //printf("%s\n", subjectmultiseq.sequence);
  /**************************************************/


  /*
   this function does everything
  */
  if(procmaxmatches(&mmcallinfo,&subjectmultiseq) != 0) {
    fprintf(stderr,"%s: %s\n",argv[0],messagespace());
    return EXIT_FAILURE;
  }
  /**************************************************/


  /*
   delete momery allocations
  */
  freemultiseq(&subjectmultiseq);
  checkspaceleak();
  mmcheckspaceleak();
  /**************************************************/
  
  fprintf(stderr, "# total running time: %d seconds.\n", time(NULL) - seconds);
  return EXIT_SUCCESS;
}
