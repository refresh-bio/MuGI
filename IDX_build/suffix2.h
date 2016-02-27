/*
 This file is a part of MuGI software distributed under GNU GPL 2 licence.
 The homepage of the MuGI project is http://sun.aei.polsl.pl/mugi/
 
 Authors: Agnieszka Danek, Sebastian Deorowicz
 
 Version: 1.0
 Date   : 2015-Mar
 */

#ifndef __build_INDEX__suffix2__
#define __build_INDEX__suffix2__

#include <iostream>
#include "defs.h"
#include "data.h"
#include "coding.h"




struct suffix2
{
    unsigned int prevPos;
    unsigned int vtList; //variants found/not found list (bits)
    unsigned int gap; //before prevPos
  
};


bool compareSuffix2(const suffix2& s1, const suffix2& s2);
char * getSuffix(suffix2 * suffix, unsigned int * nextPos);
char  getFirstChar(suffix2 * suffix);

#endif /* defined(__build_INDEX__suffix2__) */
