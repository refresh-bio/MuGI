/*
 This file is a part of MuGI software distributed under GNU GPL 2 licence.
 The homepage of the MuGI project is http://sun.aei.polsl.pl/mugi/
 
 Authors: Agnieszka Danek, Sebastian Deorowicz
 
 Version: 1.0
 Date   : 2015-Mar
 */

#ifndef __build_INDEX__sa1__
#define __build_INDEX__sa1__

#include <iostream>
#include <fstream>
#include "tbb/parallel_sort.h"
#include "defs.h"
#include "variant_list.h"
#include "bit_vectors.h"
#include "fasta.h"
#include "suffix1.h"
#include "coding.h"





class SA1{
private:
	///unsigned int n;
    uint32_t n_sparse;
    uint32_t suffixArray1Max;
    uint32_t filtered;
    uint32_t kl;
    uint32_t no_var;
    uint32_t * LUT = NULL;
    suffix1 * suffixArray1 = NULL;
    
    BitVectors * bv = NULL;
    VariantList * vl = NULL;
    CFastaFile * fasta = NULL;
    
public:
    SA1(){n_sparse = 1; suffixArray1Count = 0; suffixArray1Max = 0; filtered = 0;}
    ~SA1(){if(suffixArray1){ free(suffixArray1); suffixArray1 = NULL;} if(LUT){ free(LUT); LUT = NULL;}}
    SA1(uint32_t sparse, uint32_t _max, uint32_t _kl, uint32_t _no_var,  BitVectors * _bv, VariantList * _vl, CFastaFile * _fasta);
    
    
    
    uint32_t suffixArray1Count;
    
    
    void newSuffix1(unsigned int  prePos,  unsigned int newVtList, unsigned int vtCovered, unsigned int firstVar, unsigned int currSV);
    
    void  sort();
    void  write_sa(FILE * file);
    void  write_sa_lut(FILE * file);
    void  write_sa(FILE * file, int32_t sparse);
    
    void readData(FILE *in);
    
    
    
};

#endif /* defined(__build_INDEX__sa1__) */
