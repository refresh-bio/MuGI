/*
 This file is a part of MuGI software distributed under GNU GPL 2 licence.
 The homepage of the MuGI project is http://sun.aei.polsl.pl/mugi/
 
 Authors: Agnieszka Danek, Sebastian Deorowicz
 
 Version: 1.0
 Date   : 2015-Mar
 */


#ifndef __build_INDEX__bit_vectors__
#define __build_INDEX__bit_vectors__

#include <iostream>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <array>
#include <cstring>
#include <stdio.h>
#include <cstring>
#include <stdlib.h>
#include <vector>

#include "defs.h"
#include "variant_list.h"
using namespace std;

namespace std{
  
    template<typename T>
    struct hash<vector<T>>
    {
        typedef vector<T> argument_type;
        typedef size_t result_type;
        
        result_type operator()(const argument_type& a) const
        {
            hash<T> hasher;
            result_type h = 0;
            for (result_type i = 0; i < a.size(); ++i)
            {
                h = h * 31 + hasher(a[i]);
            }
            return h;
        }
    };
}




class BitVectors{
private:
    
    //metadata
    uint32_t no_gen;
    uint32_t ploidity;
    uint32_t no_var;
    
    
    
    //bit_vectors_per_genome
    uchar ** bit_gen;
    uint32_t * bit_pos;
    uint32_t no_bytes_per_gen;
    
    
    //variant_vectors
    uint64_t ** bit_var;
    uint32_t no_uint64_per_var;
    
    
    
    
    
    
public:
    
    BitVectors(){
        no_gen = 0;
        no_var = 0;
        ploidity = 1;
        bit_gen = NULL;
        bit_var = NULL;
        bit_pos = NULL;
        
	}
    
	~BitVectors() {
		if(bit_gen)
        {
            for (int32_t i=0; i < no_gen*ploidity; i++)
            {
                delete bit_gen[i];
            }
            delete bit_gen;
        }
        
        if(bit_var)
        {
            for(uint32_t i = 0; i < no_var; ++i)
            {
                delete bit_var[i];
            }
            delete bit_var;
        }
    }
    
    //genomes names (from VCF)
    char ** genomesNames;
    
    BitVectors(uint32 _ploidity, uint32 _no_gen, uint32 _no_var, string * _genomesNames);
    void writeBits(uint32_t gen_id, unsigned char newBits, int32 number);

    uint32_t getNoId(){return no_gen*ploidity;}
    
    
    // bit vectors per genome/individual set
    void writeBitBufToFile(uint32_t gen_id);
    void removeVariantsSetInDelRegions(VariantList * vl);
    
    
    // variant vectors
    void makeVarVec();
    void saveVarVec(FILE * fileName);
    
    
    // compressed
    int createCompressedCollection(uint32_t _blocking);
    void writeCompressedCollectionToFile(FILE * file);
    uint32_t * all_var;
    uint32_t * uniq_var;
    uint64_t * uniq_parts;
    // compressed collection
    uint32_t blocking;
    uint32_t parts;
    uint32_t cnt_uniq_var;
    uint32_t cnt_uniq_parts;
    
    
    void freeBitPos();
    
    
    /*******/
      void readData(FILE *in);
};



#endif /* defined(__build_INDEX__bit_vectors__) */
