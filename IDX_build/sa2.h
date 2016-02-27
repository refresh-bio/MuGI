/*
 This file is a part of MuGI software distributed under GNU GPL 2 licence.
 The homepage of the MuGI project is http://sun.aei.polsl.pl/mugi/
 
 Authors: Agnieszka Danek, Sebastian Deorowicz
 
 Version: 1.0
 Date   : 2015-Mar
 */

#ifndef __build_INDEX__sa2__
#define __build_INDEX__sa2__

#include <iostream>
#include <fstream>
#include "tbb/parallel_sort.h"
#include "defs.h"
#include "suffix2.h"

class SA2{
private:
	///unsigned int n;
    unsigned int n_sparse;
    
    unsigned int suffixArray2Max;
    uint32_t filtered;
    
    unsigned int *LUT;
    suffix2 * suffixArray2;
    
    
public:
    SA2(){n_sparse = 1; suffixArray2Count = 0; suffixArray2Max = 0; filtered = 0;}
    ~SA2(){if(suffixArray2) free(suffixArray2);if(LUT) free(LUT);};
    SA2(uint32_t sparse, uint32_t _max);
	
    
    unsigned int suffixArray2Count;
    void newSuffix2(unsigned int  prePos,  unsigned int newVtList, unsigned int gap,  unsigned int vtCovered, unsigned int firstVar, unsigned int currSV);
    
    void  sort();
    void  write_sa(FILE * file);
    void  write_sa_lut(FILE * file);
    
    
    
    void readData(FILE *in);
};

#endif /* defined(__build_INDEX__sa2__) */
