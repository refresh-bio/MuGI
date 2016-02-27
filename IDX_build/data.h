/*
 This file is a part of MuGI software distributed under GNU GPL 2 licence.
 The homepage of the MuGI project is http://sun.aei.polsl.pl/mugi/
 
 Authors: Agnieszka Danek, Sebastian Deorowicz
 
 Version: 1.0
 Date   : 2015-Mar
 */

#ifndef build_INDEX_data_h
#define build_INDEX_data_h

#include "fasta.h"
#include "variant_list.h"
#include "bit_vectors.h"

namespace Data
{
    extern CFastaFile * fasta;
    extern vector<CFastaSequence> data;
    extern VariantList * vl;
    extern BitVectors * bv;
    
    // parameters to build INDEX
    extern uint32_t kl;
    extern uint32_t n_sparse;
    extern uint32_t blocking;
    extern uint32_t ploidity;
}


#endif
