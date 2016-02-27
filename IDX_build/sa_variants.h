/*
 This file is a part of MuGI software distributed under GNU GPL 2 licence.
 The homepage of the MuGI project is http://sun.aei.polsl.pl/mugi/
 
 Authors: Agnieszka Danek, Sebastian Deorowicz
 
 Version: 1.0
 Date   : 2015-Mar
 */

#ifndef __build_INDEX__sa_variants__
#define __build_INDEX__sa_variants__

#include <iostream>
#include "variant_list.h"
#include "bit_vectors.h"
#include "sa1.h"
#include "sa2.h"
#include "sa3.h"
#include "fasta.h"

class SAVariants{
private:
    uint32_t n_sparse;
    uint32_t kl;
    uint32_t no_var;
   // uint32_t data_size;
    
    BitVectors * bv;
    VariantList * vl;
    CFastaFile * fasta;
    
    SA1 * sa1;
    SA2 * sa2;
    SA3 * sa3;
    
  //  int encodeSuffix(char * kmer);
    
    void addSuffixes(unsigned int prePos, unsigned int vtListSize, unsigned int * vtList, unsigned int gap, unsigned int noVt, unsigned int currVt, unsigned int refPos, unsigned short chCovered, bool delRegVariatns, unsigned int skippedVt, unsigned int lastFiltFound, unsigned int meaningfulVTinDelReg, unsigned int currentSV);
    
public:
    SAVariants(){n_sparse = 1;}
    ~SAVariants(){};
   
    SAVariants(uint32 _n_sparse){n_sparse = _n_sparse;}
    int make_sa_variants(VariantList * _vl, BitVectors * _bv, CFastaFile * _fasta, uint32_t _kl);
    
    
    void  sort();
    
    void  write_sa1(FILE * file);
    void  write_sa1(FILE * file, int32_t sparse);
    void  write_sa1_lut(FILE * file);
    
    void  write_sa2(FILE * file);
    void  write_sa2_lut(FILE * file);
    
    void  write_sa3(FILE * file);
    void  write_sa3_lut(FILE * file);
    
    
    
    void readSA123andSort(FILE *file);
};



#endif /* defined(__build_INDEX__sa_variants__) */
