#ifndef ARRAYDEF_H
#define ARRAYDEF_H


#include "types.h"
#include "errordef.h"
#include "spacedef.h"



/*
  This file defines macros to conveniently declare and 
  manipulate dynamic arrays whose size grow on demand. Each dynamic 
  array over some type T is implemented by a structure consisting of three components:
  1, space##T is a pointer to the space block of type T allocated for the array.
  2, allocated##T is an Uint storing the number of entries.
  3, nextfree##T holds the smallest index of the array where no value is stored.
  Here ## is the concatenation operator of the C-preprocessor.
  The following macro expands to a corresponding type definition over 
  some given TYPE.
*/
#define DECLAREARRAYSTRUCT(TYPE)\
        typedef struct\
        {\
          TYPE *space##TYPE;\
          Uint allocated##TYPE, nextfree##TYPE;\
        } Array##TYPE


// initializes an empty dynamic array.
#define INITARRAY(A,TYPE)\
        (A)->space##TYPE = NULL;\
        (A)->allocated##TYPE = (A)->nextfree##TYPE = 0


/*
  check if the integer nextfree##T points to an index for which
  the space is not allocated yet. If this is the case,
  the number of cells allocated is incremented by L.
*/
#define CHECKARRAYSPACE(A,TYPE,L)\
        if((A)->nextfree##TYPE >= (A)->allocated##TYPE)\
        {\
          (A)->allocated##TYPE += L;\
          (A)->space##TYPE\
             = (TYPE *) allocandusespaceviaptr(__FILE__,(Uint) __LINE__,\
                                               (A)->space##TYPE,\
                                               (Uint) sizeof(TYPE),\
                                               (A)->allocated##TYPE);\
        }\
        NOTSUPPOSEDTOBENULL((A)->space##TYPE)


/*
  check the space and deliver a pointer P 
  to the next free element in the array.
*/
#define GETNEXTFREEINARRAY(P,A,TYPE,L)\
        CHECKARRAYSPACE(A,TYPE,L);\
        P = (A)->space##TYPE + (A)->nextfree##TYPE++;


/*
  check the space and store V in the 
  nextfree component of the array. nextfree is incremented.
*/
#define STOREINARRAY(A,TYPE,L,VAL)\
        CHECKARRAYSPACE(A,TYPE,L);\
        (A)->space##TYPE[(A)->nextfree##TYPE++] = VAL


// free the space for an array if it is not NULL.
#define FREEARRAY(A,TYPE)\
        if((A)->space##TYPE != NULL)\
        {\
          FREESPACE((A)->space##TYPE);\
        }


/* 
  Some declarations for the most common array types.
*/
DECLAREARRAYSTRUCT(Uchar);
DECLAREARRAYSTRUCT(Ushort);
DECLAREARRAYSTRUCT(char);
DECLAREARRAYSTRUCT(Uint);
DECLAREARRAYSTRUCT(Sint);
DECLAREARRAYSTRUCT(PairUint);
DECLAREARRAYSTRUCT(ThreeUint);


typedef ArrayUint  ArrayPosition;
typedef ArrayUchar ArrayCharacters;


/*
  The following array type has some extra components. However, it can be 
  manipulated by the macros above since the record-components
  spaceStrings, nextfreeStrings, and 
  allocatedStrings is declared appropriately.
*/
typedef struct
{
  Stringtype *spaceStrings;
  Uchar *stringbuffer;
  Uint stringbufferlength, nextfreeStrings, allocatedStrings;
} ArrayStrings;   


#endif

