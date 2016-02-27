/*
 This file is a part of MuGI software distributed under GNU GPL 2 licence.
 The homepage of the MuGI project is http://sun.aei.polsl.pl/mugi/
 
 Authors: Agnieszka Danek, Sebastian Deorowicz
 
 Version: 1.0
 Date   : 2015-Mar
 */

#ifndef __build_INDEX__suffix3__
#define __build_INDEX__suffix3__

#include <iostream>
#include "defs.h"
#include "data.h"
#include "coding.h"
#include "omp.h"

struct suffix3
{
    unsigned int prevPos;
    unsigned int vtListSize;
    unsigned int *vtList; //variants found/not found list (bits)
    unsigned int gap; 
    
};


bool compareSuffix3(const suffix3& s1, const suffix3& s2);
char * getSuffix(suffix3 * suffix, unsigned int * nextPos);

char  getFirstChar(suffix3 * suffix);




#endif /* defined(__build_INDEX__suffix3__) */
