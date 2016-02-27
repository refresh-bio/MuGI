/*
 This file is a part of MuGI software distributed under GNU GPL 2 licence.
 The homepage of the MuGI project is http://sun.aei.polsl.pl/mugi/
 
 Authors: Agnieszka Danek, Sebastian Deorowicz
 
 Version: 1.0
 Date   : 2015-Mar
 */

#ifndef __build_INDEX__variant_list__
#define __build_INDEX__variant_list__

#include <iostream>

#include <algorithm>
#include "defs.h"
#include "coding.h"

#define SNP 0
#define INS 1
#define DEL 2
#define SV  3
#define END  0xFF

using namespace std;

struct variantAuxiliary
{
    char * alt;
    uint32_t delLen;
};

struct variant
{
    uint32_t prePos;  //position 'one before variant is found'
    uint32_t alt;      //alt sign for SNP or delLen for DEL or idx in ins_sv_Arr for INS/SV
};

struct singleDelFromRight
{
    uint32_t nextRefPos;
    uint32_t vtNo;
};
bool compareDelFromRight(singleDelFromRight a, singleDelFromRight b);

class VariantList{
private:
    
    
    uchar * ins_sv_arr;
    
    
    // auxiliary structures for fast searching tool
    struct singleDelFromRight * delFromRight;
    uint32_t * prevDeletion;
    uint32_t * VT_LUT;


    
public:
    
    uint32_t no_var;
    uint32_t ref_size;
    
    uint32_t no_sv;
    uint32_t no_chars_inserted;
    uint32_t no_ins;
    uint32_t no_del;
    uint32_t max_del;
    uint32_t max_sv_del;
    
    variant * vt;
    variantAuxiliary * vtAux;
    uchar * vtType;
    uint32_t * listSV;
    uint32_t * prevSV;
    
    VariantList()
    {
        ref_size = 0;
        no_var = 0;
        no_sv = 0;
        no_ins = 0;
        no_del = 0;
        no_chars_inserted = 0;
        max_del = 0;
        max_sv_del = 0;
        vt = NULL;
        vtAux = NULL;
        prevSV = NULL;
        listSV = NULL;
        delFromRight = NULL;
        prevDeletion = NULL;
        VT_LUT = NULL;
        
    }
    
    ~VariantList()
    {
        if(vt)
            delete vt;
        if(vtAux)
            delete vtAux;
        if(vtType)
            delete vtType;
        if(prevSV)
            delete prevSV;
    }
    
    
    
    
    
    VariantList(uint32_t _no_var, uint32_t _ref_size);
    
    
    
    bool transform();
    
    bool write_compact_vl(FILE * file);
    
    uchar getType(uint32_t var_no);
    
    uint32_t getPos(uint32_t var_no);
    
    uint32_t getDelLen(uint32_t var_no);
    uint32_t getInsLen(uint32_t var_no);
    
    uint32_t getNoVar(){return no_var;}
    
    
    uint32_t get_no_ins(){return no_ins;}
    uint32_t get_no_sv(){return no_sv;}
    uint32_t get_maxSVdel(){return max_sv_del;}
    uint32_t get_maxDel(){return max_del;}
    uint32_t get_no_var(){return no_var;}
    
    void set_no_var(uint32_t _no_var){no_var = _no_var;}
    void set_ref_size(uint32_t _ref_size);
    
    
    uchar * getPtrToFirstInsChar(uint32_t var_no);

    /****/
     void readData(FILE *in);

    
};




#endif /* defined(__build_INDEX__variant_list__) */

