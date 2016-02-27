/*
 This file is a part of MuGI software distributed under GNU GPL 2 licence.
 The homepage of the MuGI project is http://sun.aei.polsl.pl/mugi/
 
 Authors: Agnieszka Danek, Sebastian Deorowicz
 
 Version: 1.0
 Date   : 2015-Mar
 */

#ifndef __approximateSeachBetter_index__SuffixArray__
#define __approximateSeachBetter_index__SuffixArray__

#include <iostream>
#include "common.h"


extern uint32_t * suffixArray;;
void  searchString(unsigned char * read1, int32_t len, uint32_t encodedPattern,  int32_t delta, uint64_t ** passed_result, uint64_t * pos, uint64_t *max, unsigned char * wholeRead, uint32_t wholeReadLen,  int32_t leftPos, uint32_t rightPos, uint32_t maxError, int32_t lenToCheck, uint32_t rightBorder);


uint32_t getCoveredVariantsSA(uint32_t  prevPos,  uint32_t * firstVar,  uint32_t * chCov, uint32_t * nextRefPos, uint32_t * currVT, uint32_t * pos_vtToCheck, uint32_t  lenToCheck);



void searchStringExact(unsigned char * read1, int32_t len, uint32_t encodedPattern,  int32_t delta, uint64_t ** passed_result, uint64_t * pos, uint64_t *max, int32_t lenToCheck, uint32_t rightBorder);

#endif /* defined(__approximateSeachBetter_index__SuffixArray__) */
