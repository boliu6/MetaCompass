#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "spacedef.h"
#include "protodef.h"
#include "debugdef.h"
#include "maxmatdef.h"



signed long scanmultiplefastafile(Multiseq *multiseq,
				  char *filename,
				  Uchar replacewildcardchar,
				  Uchar *input,
				  unsigned long inputlen);


signed long loadinrefseq(Multiseq *subjectmultiseq,
			 char *subjectfile);
