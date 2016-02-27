/*
 This file is a part of MuGI software distributed under GNU GPL 2 licence.
 The homepage of the MuGI project is http://sun.aei.polsl.pl/mugi/
 
 Authors: Agnieszka Danek, Sebastian Deorowicz
 
 Version: 1.0
 Date   : 2015-Mar
 */

#ifndef __build_INDEX__coding__
#define __build_INDEX__coding__

#include <iostream>
#include "defs.h"


namespace Coding{
    
    
    uchar encodeChar(char x);
    char decodeCode(uchar x);
    int encodeSuffix(char * kmer);
    
    void StoreUInt(unsigned char *dest, uint32_t  data);
    void LoadUInt(unsigned char *dest, uint32_t  &data);
    
};
#endif /* defined(__build_INDEX__coding__) */
