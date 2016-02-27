/*
 This file is a part of MuGI software distributed under GNU GPL 2 licence.
 The homepage of the MuGI project is http://sun.aei.polsl.pl/mugi/
 
 Authors: Agnieszka Danek, Sebastian Deorowicz
 
 Version: 1.0
 Date   : 2015-Mar
 */

#ifndef __approximateSeachBetter_index__coding__
#define __approximateSeachBetter_index__coding__

#include <iostream>
#include "common.h"


#endif /* defined(__approximateSeachBetter_index__coding__) */


void StoreUInt(unsigned char *dest, uint32_t  data);

void LoadUInt(unsigned char *dest, uint32_t  &data);

unsigned char getRef(uint32_t n);

char decodeCode(unsigned char x);


char *  decode(uint32_t prevPos, uint32_t length, int32_t id,  uint32_t vtStart, int32_t gap);