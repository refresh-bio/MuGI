/*
 This file is a part of MuGI software distributed under GNU GPL 2 licence.
 The homepage of the MuGI project is http://sun.aei.polsl.pl/mugi/
 
 Authors: Agnieszka Danek, Sebastian Deorowicz
 
 Version: 1.0
 Date   : 2015-Mar
 */

#ifndef __build_INDEX__suffix1__
#define __build_INDEX__suffix1__

#include <iostream>
#include "defs.h"
#include "data.h"
#include "coding.h"


struct suffix1
{
    unsigned int prevPos;
    unsigned int vtList; //variants found/not found list (bits)
#ifdef CREATE_ALL_INDEXES
    unsigned int startWithVariant;
#endif
    
    
};



bool compareSuffix1(const suffix1& s1, const suffix1& s2);

char * getSuffix(suffix1 * suffix, unsigned int * nextPos);

char  getFirstChar(suffix1 * suffix);

#endif /* defined(__build_INDEX__suffix1__) */
