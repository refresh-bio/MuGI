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
//#include "fasta.h"
#include "defs.h"
#include "asmlib.h"
#include "coding.h"


#include <errno.h>
#include <cstring>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>


#define SNP 0
#define INS 1
#define DEL 2
#define SV  3
#define END  0xFF

#define NO_BITS_UINT 32U // (uint32_t) 8*sizeof(uint32_t)
#define NO_BITS_UINT_MINUS_ONE 31U

#define MAX_READ_SIZE 1000
using namespace std;

extern uint32_t parts, blocking, cnt_uniq_parts, cnt_uniq_var, cnt_var;
extern uint32_t *all_var;
extern uint32_t *uniq_var;
extern uint64_t *uniq_parts;



extern unsigned char *s_ch;

extern uint32_t kl; //k-mer length

extern uint32_t suffixArrayCount;
extern uint32_t suffixArray1Count ;
extern uint32_t suffixArray1Max ;
extern uint32_t suffixArray2Count ;
extern uint32_t suffixArray2Max ;
extern uint32_t suffixArray3Count ;
extern uint32_t suffixArray3Max ;
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
struct variant
{
    uint32_t prePos;  //position 'one before variant is found'
    uint32_t alt;  //alt sign for SNP or delLen for DEL or idx in ins_sv_Arr for INS/SV
};

extern unsigned char * ins_sv_arr;
//for INS: 4 bytes for ins len, len bytes for string
//for SV: 4 bytes for delLen, 4 bytes for ins len, len bytes for string


struct singleDelFromRight
{
    uint32_t nextRefPos;
    uint32_t vtNo;
};

extern singleDelFromRight * delFromRight;


extern variant * vt;
extern uint32_t * prevSV, *prevDeletion;


extern unsigned char* reference;
extern uint32_t refSize;

unsigned char getRef(uint32_t n);
unsigned char getRefFinal(uint32_t n);
unsigned char encodeChar(char x);


uint32_t getCoveredVariants(uint32_t  prevPos,  uint32_t vtListLoc, uint32_t * firstVar, uint32_t gap, uint32_t * chCov, uint32_t * nextRefPos,  uint32_t * pos_vtToCheck, uint32_t  lenToCheck, uint32_t * remainedToInsert);
uint32_t getCoveredVariants(uint32_t  prevPos,  uint32_t * vtListLoc, uint32_t * firstVar, uint32_t gap, uint32_t * chCov, uint32_t * nextRefPos, uint32_t * currVt, uint32_t * pos_vtToCheck, uint32_t  lenToCheck, uint32_t * remainedToInsert);



uint64_t checkVTconfiguration1(uint32_t prePos,  uint32_t noVt,  uint32_t firstVar, uint32_t currVt, uint32_t lastVt, uint32_t  firstVarRead1, uint32_t refPos, uint32_t chCovered, unsigned char * read1, int32_t readLen, uint32_t delta,  int32_t preChToCheck,  uint64_t ** passed_result, uint64_t pos, uint64_t * max, uint32_t prevDelCount, uint32_t firstSV, uint32_t lastSV, uint32_t  pos_vtToCheck, unsigned char * wholeRead, uint32_t wholeReadLen,   int32_t leftPos, uint32_t rightPos, uint32_t error, uint32_t maxError, uint32_t nextRefPos, uint32_t gap, uint32_t remainedToInsert, uint32_t startPosFoundExact, uint32_t prePosExact, uint32_t prevDelVar);


uint64_t checkVTconfigurationExact(uint32_t prePos, uint32_t noVt,  uint32_t firstVar, uint32_t currVt,  uint32_t lastVt,  uint32_t firstVarRead1, uint32_t refPos, uint32_t chCovered, unsigned char * read1, int32_t readLen, uint32_t delta,  int32_t preChToCheck,  uint64_t ** passed_result, uint64_t pos, uint64_t * max, uint32_t prevDelCount, uint32_t firstSV, uint32_t lastSV,  uint32_t  pos_vtToCheck, uint32_t nextRefPos, uint32_t gap, uint32_t remainedToInsert, uint32_t prevDelVar);


uint64_t checkDelRegion(uint32_t currfFirstVar, uint32_t prevFirstVar, int32_t leftPos, int32_t gapFirstVar, int32_t error, int32_t maxError, unsigned char * wholeRead,  uint32_t  pos_vtToCheck, uint32_t vtCovered, int32_t  prev_firstSV, int32_t prev_lastSV, uint32_t  prev_prevDelCount, uint64_t ** passed_result, uint64_t  pos, uint64_t *max, uint32_t lastVt, uint32_t nextRefPos, uint32_t remainedToInsert, uint32_t chCovered, unsigned char * read1, int32_t readLen, uint32_t delta,  int32_t preChToCheck,  uint32_t wholeReadLen, uint32_t rightPos);


extern uint64_t **vtVector;

extern uint32_t *LUT, *LUT1, *LUT2, *LUT3, *VT_LUT;

extern int32_t encodeSuffix(const char * kmer);



extern uint64_t *resV;
extern uint32_t * vtList;
extern uint32_t max_vtToCheck, vtListSize;
extern uint32_t * vtToCheck;
extern uint64_t resSize;

extern uint32_t ileSame;



#endif
