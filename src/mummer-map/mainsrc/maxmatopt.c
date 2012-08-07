#include <stdio.h>
#include <stdlib.h>

#include <string.h>
// strcpy function

#include "types.h"
#include "optdesc.h"
#include "debugdef.h"
#include "errordef.h"
#include "protodef.h"
#include "maxmatdef.h"


// read command line options and store them in mmcallinfo
// The default value for the minimal unique match length.
#define DEFAULTMINUNIQUEMATCHLEN 18




typedef enum {
  OPTLEASTLENGTH,
  OPTMAXMIS,
  OPTMINPCT,
  OPTPREFIX,
  OPTNUMP,
  OPTH,
  OPTHELP,
  NUMOFOPTIONS
} Optionnumber;





// help text for the option -l
static void makeleastlengthtext(char *spacefortext) {
  
  sprintf(spacefortext,
	  "set the minimum length of seeds (default: %lubp)",
	  (Showuint) DEFAULTMINUNIQUEMATCHLEN);
}




// display help information
static void showusage(char *program,
		      OptionDescription *options,
                      Uint numofoptions) {

  fprintf(stderr, "\nUsage: mummer-map [options] <reference-file> <query-file>\n\n"
	  "Map short reads in <query-file> to reference sequences in <reference-file>.\n"
	  "Exact seeds (maximal matches) are found by MUMmer,and are extended.\n"
	  "\n");
  fprintf(stderr, "Options:\n");
  showoptions(stderr, program, options, numofoptions);
  fprintf(stderr, "\n");
  fprintf(stderr, "Version: 0.8.6\n");
  fprintf(stderr, "Contact: Bo Liu <boliu@umiacs.umd.edu>\n\n");
  
}




// parse command line options. If everything is okay, 0 is returned and the 
// mmcallinfo is correctly initialized. Otherwise, a negative value is returned.
Sint parsemaxmatoptions(MMcallinfo *mmcallinfo,
			Argctype argc,
			char **argv) {
  

  OptionDescription options[NUMOFOPTIONS];   // store the options
  Sint optval;                               // neg. return val. if error, otherwise option number
  Uint argnum;                               // pointer to argv
  Scaninteger readint;                       // temporary integer to read value from string
  char readchar[200];
  char leastlengthtext[128+1];


  initoptions(&options[0],(Uint) NUMOFOPTIONS); // initialization


  makeleastlengthtext(&leastlengthtext[0]);
  ADDOPTION(OPTLEASTLENGTH,"-l",&leastlengthtext[0]);
  ADDOPTION(OPTNUMP,"-p",
	    "number of processors (default: 1)");
  ADDOPTION(OPTMAXMIS,"-m",
	    "maximum number of mismatches including gaps (default: 10; max: 10)");
  ADDOPTION(OPTMINPCT,"-s",
	    "minimum % of simialrity (default: 90)");
  ADDOPTION(OPTPREFIX,"-o",
	    "output prefix (default: output)");
  ADDOPTION(OPTH,"-h",
	    "show possible options");
  ADDOPTION(OPTHELP,"-help",
            "show possible options");

  
  mmcallinfo->nump = 1;
  mmcallinfo->nummaxmis = 10;
  mmcallinfo->minpct = 90;
  strcpy(mmcallinfo->outprefix, "output");
  mmcallinfo->minmatchlength = (Uint) DEFAULTMINUNIQUEMATCHLEN;


  if(argc == 1) {
    showusage(argv[0],&options[0],(Uint) NUMOFOPTIONS);
    return 1;
  }
  

  for(argnum = UintConst(1);
      argnum < (Uint) argc && argv[argnum][0] == '-'; 
      argnum++) {

    optval = procoption(options,(Uint) NUMOFOPTIONS,argv[argnum]);

    if(optval < 0)
      return -1;
    
    switch(optval) {

    case OPTLEASTLENGTH:  // seed length
      argnum++;
      if(argnum > (Uint) (argc-2)) {
	ERROR1("missing argument for option %s",
	       options[OPTLEASTLENGTH].optname);
	return -2;
      }
      if(sscanf(argv[argnum],"%ld",&readint) != 1 || readint <= 0) {
	ERROR2("argument %s for option %s is not a positive integer",
	       argv[argnum],options[OPTLEASTLENGTH].optname);
	return -3;
      }
      mmcallinfo->minmatchlength = (Uint) readint;
      break;
      

    case OPTNUMP:  // # of processors
      argnum++;
      if(argnum > (Uint) (argc-2)) {
	ERROR1("missing argument for option %s",
	       options[OPTNUMP].optname);
	return -2;
      }
      if(sscanf(argv[argnum],"%ld",&readint) != 1 || readint <= 0) {
	ERROR2("argument %s for option %s is not a positive integer",
	       argv[argnum],options[OPTNUMP].optname);
	return -3;
      }
      mmcallinfo->nump = (Uint) readint;
      break;
      

    case OPTMAXMIS:  // # of maximum mismatches
      argnum++;
      if(argnum > (Uint) (argc-2)) {
	ERROR1("missing argument for option %s",
	       options[OPTMAXMIS].optname);
	return -2;
      }
      if(sscanf(argv[argnum],"%ld",&readint) != 1 || readint <= 0) {
	ERROR2("argument %s for option %s is not a positive integer",
	       argv[argnum],options[OPTMAXMIS].optname);
	return -3;
      }
      mmcallinfo->nummaxmis = (Uint) readint;
      if (mmcallinfo->nummaxmis > 10) return -3;
      break;
      

    case OPTMINPCT:  // # of maximum mismatches
      argnum++;
      if(argnum > (Uint) (argc-2)) {
	ERROR1("missing argument for option %s",
	       options[OPTMINPCT].optname);
	return -2;
      }
      if(sscanf(argv[argnum],"%ld",&readint) != 1 || readint <= 0 || readint > 100) {
	ERROR2("argument %s for option %s is not between 0-100",
	       argv[argnum],options[OPTMINPCT].optname);
	return -3;
      }
      mmcallinfo->minpct = (Uint) readint;
      break;
      

    case OPTPREFIX:  // output prefix
      argnum++;
      if(argnum > (Uint) (argc-2)) {
	ERROR1("missing argument for option %s",
	       options[OPTPREFIX].optname);
	return -2;
      }
      if(sscanf(argv[argnum],"%s",&readchar) != 1 || readchar == '\0') {
	ERROR2("argument %s for option %s is not valid",
	       argv[argnum],options[OPTPREFIX].optname);
	return -3;
      }
      strcpy(mmcallinfo->outprefix, readchar);
      break;


    case OPTH:
    case OPTHELP:
      showusage(argv[0],&options[0],(Uint) NUMOFOPTIONS);
      return 1;
    }
  }

  
  if(argnum > (Uint) (argc-2)) {
    ERROR0("missing file arguments");
    return -4;
  }
  if(safestringcopy(&mmcallinfo->program[0],argv[0],PATH_MAX) != 0) {
    return -5;
  }
  if(safestringcopy(&mmcallinfo->subjectfile[0],argv[argnum],PATH_MAX) != 0) {
    return -6;
  }

  for(argnum++, mmcallinfo->numofqueryfiles = 0; 
      argnum < (Uint) argc;
      mmcallinfo->numofqueryfiles++, argnum++) {
    if(mmcallinfo->numofqueryfiles >= (Uint) MAXNUMOFQUERYFILES) {
      ERROR1("too many query files, maximal number is %lu",
	     (Showuint) MAXNUMOFQUERYFILES);
      return -7;
    }
    if(safestringcopy(&mmcallinfo->queryfilelist
		      [mmcallinfo->numofqueryfiles][0],
		      argv[argnum],PATH_MAX) != 0) {
      return -8;
    }
  }

  return 0;
}
