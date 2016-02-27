/*
 This file is a part of MuGI software distributed under GNU GPL 2 licence.
 The homepage of the MuGI project is http://sun.aei.polsl.pl/mugi/
 
 Authors: Agnieszka Danek, Sebastian Deorowicz
 
 Version: 1.0
 Date   : 2015-Mar
 */

#include "data.h"

namespace Data
{
    CFastaFile *fasta = new CFastaFile();
    vector<CFastaSequence> data;
    VariantList * vl;
    BitVectors * bv;
    
    uint32_t kl;
    uint32_t n_sparse;
    uint32_t blocking;
    uint32_t ploidity;
}