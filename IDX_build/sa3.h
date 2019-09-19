/*
 This file is a part of MuGI software distributed under GNU GPL 2 licence.
 The homepage of the MuGI project is http://sun.aei.polsl.pl/mugi/
 
 Authors: Agnieszka Danek, Sebastian Deorowicz
 
 Version: 1.0
 Date   : 2015-Mar
 */

#ifndef __build_INDEX__sa3__
#define __build_INDEX__sa3__

#include <iostream>
#include <fstream>
#include <algorithm>
#include "tbb/parallel_sort.h"
#include "defs.h"
#include "suffix3.h"

#include "data.h"



class SA3{
private:
	///unsigned int n;
    unsigned int n_sparse;
    
    unsigned int suffixArray3Max;
    uint32_t filtered;
    
    unsigned int *LUT = NULL;
    suffix3 * suffixArray3 = NULL;
    
    
public:
    SA3(){n_sparse = 1; suffixArray3Count = 0; suffixArray3Max = 0; filtered = 0;}
    ~SA3(){if(suffixArray3){ free(suffixArray3); suffixArray3 = NULL;} if(LUT){ free(LUT); LUT = NULL;}}
    SA3(uint32_t sparse, uint32_t _max);

    unsigned int suffixArray3Count;
    
    void newSuffix3(unsigned int  prePos,  unsigned int vtListSize, unsigned int *newVtList, unsigned int gap,  unsigned int vtCovered, unsigned int firstVar, unsigned int lastFilt, unsigned int currSV);

    
    void sort();

    
    void  write_sa(FILE * file);
    void  write_sa_lut(FILE * file);
    
    
    void readData(FILE *in);
};

#endif /* defined(__build_INDEX__sa3__) */
