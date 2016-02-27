/*
 This file is a part of MuGI software distributed under GNU GPL 2 licence.
 The homepage of the MuGI project is http://sun.aei.polsl.pl/mugi/
 
 Authors: Agnieszka Danek, Sebastian Deorowicz
 
 Version: 1.0
 Date   : 2015-Mar
 */

#ifndef INDEX_v2_common_h
#define INDEX_v2_common_h

#define LUT_ENCODED_CHARS 12

#define LUT_MAX  ((1U<<(LUT_ENCODED_CHARS*2)))//16777216
//one for end_of_array guard


#define VT_LUT_STEP 1000






#include <errno.h>
#include <cstring>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include "coding.h"

#define BIT_VECTORS_INPUT

#define SNP 0
#define INS 1
#define DEL 2
#define SV  3
#define END  0xFF

#define NO_BITS_UINT 32U // (uint32_t) 8*sizeof(uint32_t)
#define NO_BITS_UINT_MINUS_ONE 31U

typedef unsigned char uchar;

#define MAX_READ_SIZE 1000


struct variant
{
    uint32_t prePos;  //position 'one before variant is found'
    uint32_t alt;  //alt sign for SNP or delLen for DEL or idx in ins_sv_Arr for INS/SV
};

struct singleDelFromRight
{
    uint32_t nextRefPos;
    uint32_t vtNo;
};

extern unsigned char* reference;
extern uint32_t refSize;

using namespace std;

extern uint32_t parts, blocking, cnt_uniq_parts, cnt_uniq_var, cnt_var;
extern uint32_t *all_var;
extern uint32_t *uniq_var;
extern uint64_t *uniq_parts;



extern unsigned char *s_ch;

extern uint32_t kl; //k-mer length

extern uint32_t maxDel, maxSVdel;

extern uint32_t * listSV, noSV;

extern uint32_t lastKmerNoVar;
extern uint32_t noVar, filtered;
extern int32_t inBitV;
extern char  **inputBuffer;
//extern char *inputBuffer_t;
extern char **indNames;
//uint32_t m = 0;

extern uint32_t maxMeaningfulVTinDelReg ;

extern uchar * vtType;


extern unsigned char * ins_sv_arr;




extern singleDelFromRight * delFromRight;


extern variant * vt;
extern uint32_t * prevSV, *prevDeletion;


extern uint64_t **vtVector;



extern uint64_t *resV;
extern uint32_t * vtList;
extern uint32_t max_vtToCheck, vtListSize;
extern uint32_t * vtToCheck;
extern uint64_t resSize;



#endif
