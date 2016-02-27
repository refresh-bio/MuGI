/*
 This file is a part of MuGI software distributed under GNU GPL 2 licence.
 The homepage of the MuGI project is http://sun.aei.polsl.pl/mugi/
 
 Authors: Agnieszka Danek, Sebastian Deorowicz
 
 Version: 1.0
 Date   : 2015-Mar
 */

#ifndef __INDEX_v2__SuffixArray3__
#define __INDEX_v2__SuffixArray3__

#include <iostream>
#include "common.h"

struct suffix3
{
    uint32_t prevPos;
    uint32_t vtListSize;
    uint32_t *vtList; //variants found/not found list (bits)
    uint32_t gap; //before prevPos
    
};
extern suffix3 * suffixArray3;

int32_t compareReadSuffix(unsigned char* read, uint32_t readLen,  suffix3 suffix);

bool compareSuffix3(const suffix3& s1, const suffix3& s2);

void  searchSA3(unsigned char * read1, int32_t len, uint32_t encodedPattern, uint32_t delta, uint64_t  ** passed_result, uint64_t * pos, uint64_t *max, unsigned char * wholeRead, uint32_t wholeReadLen,  int32_t leftPos, uint32_t rightPos, uint32_t maxError, int32_t lenToCheck, uint32_t rightBorder);

void  searchSA3Exact(unsigned char * read1, int32_t len, uint32_t encodedPattern, uint32_t delta, uint64_t  ** passed_result, uint64_t * pos, uint64_t *max, int32_t lenToCheck, uint32_t rightBorder);


#endif /* defined(__INDEX_v2__SuffixArray3__) */
