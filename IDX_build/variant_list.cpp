/*
 This file is a part of MuGI software distributed under GNU GPL 2 licence.
 The homepage of the MuGI project is http://sun.aei.polsl.pl/mugi/
 
 Authors: Agnieszka Danek, Sebastian Deorowicz
 
 Version: 1.0
 Date   : 2015-Mar
 */


#include "variant_list.h"



VariantList::VariantList(uint32_t _no_var, uint32_t _ref_size)
{
    ref_size = _ref_size;
    no_var = _no_var;
    no_sv = 0;
    no_ins = 0;
    no_del = 0;
    max_del = 0;
    max_sv_del = 0;
    no_chars_inserted = 0;
    
    vtAux = new variantAuxiliary[no_var];
    vt = new variant[no_var+1];
    vtType = new uchar[no_var+1];
    
    prevSV = new uint32_t[no_var];
    prevDeletion = new uint32_t[no_var];
    
    delFromRight = NULL;
    VT_LUT = NULL;
    
}

bool VariantList::transform()
{
    ins_sv_arr = new uchar[2*sizeof(uint32_t)*no_sv + sizeof(uint32_t)*no_ins + no_chars_inserted];
    uint32_t pos = 0;
    uint32_t tempLen;
    
    for(int32_t i=0; i < no_var; ++i)
    {
        switch (vtType[i])
        {
            case SNP:
            {
            
                vt[i].alt = Coding::encodeChar(vtAux[i].alt[0]);
                break;
            
            }
            case DEL:
            {
                
                vt[i].alt = vtAux[i].delLen;
                break;
                
            }
            case INS:
            {
                
                vt[i].alt = pos;
                
                
                //store insertion length
                tempLen = (uint32_t)strlen(vtAux[i].alt);
                //memcpy(ins_sv_arr+pos, &tempLen, sizeof(uint32_t));
                Coding::StoreUInt(ins_sv_arr+pos, tempLen);
                pos += sizeof(uint32);
                
                //store inserted chars
                for(int32_t n = 0; n < tempLen; n++)
                {
                    ins_sv_arr[pos + n]= Coding::encodeChar(vtAux[i].alt[n]);
                }
                pos += tempLen;
                
                break;
                
            }
            case SV:
            {
                
                vt[i].alt = pos;
                
                //store deletion length
                //memcpy(ins_sv_arr+pos, &vtAux[i].delLen, sizeof(uint32_t));
                Coding::StoreUInt(ins_sv_arr+pos, vtAux[i].delLen);
                
                pos += sizeof(uint32_t);

                
                //store insertion length
                if(vtAux[i].alt[0] == '-')
                    tempLen = 0;
                else
                    tempLen = (uint32_t)strlen(vtAux[i].alt);
                //memcpy(ins_sv_arr+pos, &tempLen, sizeof(uint32_t));
                Coding::StoreUInt(ins_sv_arr+pos, tempLen);
                pos += sizeof(uint32_t);
                
                //store inserted chars
                for(int32_t n = 0; n < tempLen; n++)
                {
                    ins_sv_arr[pos + n]= Coding::encodeChar(vtAux[i].alt[n]);
                }
                pos += tempLen;
                
                break;
                
            }
                
        }
    }
    
    vtType[no_var] = END;
    vt[no_var].prePos = ref_size - 1;
    vt[no_var].alt = 0xFFFFFFFF;
    
    listSV = new uint32_t [no_sv+1];
    uint32_t currSV = 0;
    for(uint32_t i=0; i < no_var; i++)
    {
        if(vtType[i] == SV)
            listSV[currSV++] = i;
    }
    
    if(no_sv == 0)
        listSV[0] = no_var;
    
    ///
    delFromRight = new struct singleDelFromRight[no_del+no_sv+1];
    unsigned int delLen;
    

    
    delFromRight[0].vtNo = no_var;
    delFromRight[0].nextRefPos = 0;
    
    int cnt = 1;
    
    for(int i=0; i < no_var; ++i)
    {
        if(vtType[i] == DEL)
        {
            delFromRight[cnt].vtNo = i;
            delFromRight[cnt++].nextRefPos = vt[i].prePos + vt[i].alt + 1;
            
        }
        else if(vtType[i] == SV)
        {
            delFromRight[cnt].vtNo = i;
            //LoadUint(ins_sv_arr+vt[i].alt, delLen);
            memcpy(&delLen, ins_sv_arr+vt[i].alt, sizeof(unsigned int));
            delFromRight[cnt++].nextRefPos = vt[i].prePos + delLen + 1;
        }
    }
    
    
    sort(delFromRight, delFromRight+no_del+no_sv+1, compareDelFromRight);
    
    
    
   // prevDeletion = (unsigned int * ) malloc((no_var)*sizeof(unsigned int));
    
    cnt = no_del+no_sv;
    
    for(int i = no_var-1; i >= 0; --i)
    {
        //    cout << i << " " << vt[i].prePos << " " << cnt << " " << delFromRight[cnt].nextRefPos << endl;
        
        while(vt[i].prePos + 1 < delFromRight[cnt].nextRefPos)
        {
            --cnt;
        }
        prevDeletion[i] = cnt;//delFromRight[cnt].vtNo;
        
    }
    
    
    unsigned int curr, lut_size =  ref_size/VT_LUT_STEP ; //vt[noVar-1].prePos/VT_LUT_STEP + 1;
    
    VT_LUT = new unsigned int[lut_size];
    
    unsigned int currPrePos = 0;
    for(int i=0; i < no_var && lut_size > 0; ++i)
    {
        if(vt[i].prePos/VT_LUT_STEP >= currPrePos)
        {
            curr = vt[i].prePos/VT_LUT_STEP;
            for(int j = currPrePos; j <= curr; j++)
                VT_LUT[j] = i;
            currPrePos = vt[i].prePos/VT_LUT_STEP + 1;
        }
    }
    while(currPrePos < lut_size)
        VT_LUT[currPrePos++] = no_var;
    

    delete vtAux;
    vtAux = NULL;
    return true;
}



bool VariantList::write_compact_vl(FILE * file)
{
    uint32_t temp;
    
    fwrite(&no_var, 1, sizeof(uint32_t), file);
    fwrite(&no_sv, 1, sizeof(uint32_t), file);
    fwrite(&no_ins, 1, sizeof(uint32_t), file);
    fwrite(&no_del, 1, sizeof(uint32_t), file);
    fwrite(&no_chars_inserted, 1, sizeof(uint32_t), file);
    
    
    fwrite(&max_sv_del, 1, sizeof(uint32_t), file);
    fwrite(&max_del, 1, sizeof(uint32_t), file);
    
    
    fwrite(vtType, sizeof(uchar), no_var+1, file);
    fwrite(vt, sizeof(variant), no_var+1, file);
    fwrite(ins_sv_arr, sizeof(uchar), (2*sizeof(uint32_t)*no_sv + sizeof(uint32_t)*no_ins + no_chars_inserted), file);
    
    
    
    fwrite(listSV, sizeof(uint32_t), no_sv, file);
    fwrite(prevSV, sizeof(uint32_t), no_var, file);
    fwrite(prevDeletion, sizeof(uint32_t), no_var, file);
    fwrite(delFromRight, sizeof(struct singleDelFromRight), (no_del+no_sv) + 1, file);
    temp = VT_LUT_STEP;
    fwrite(&temp, sizeof(uint32_t), 1, file);
    fwrite(VT_LUT, sizeof(uint32_t), ref_size/VT_LUT_STEP, file);

    return true;
}


uchar VariantList::getType(uint32_t var_no)
{
    
    return vtType[var_no];
}

uint32_t VariantList::getPos(uint32_t var_no)
{
    return vt[var_no].prePos;
}

uint32_t VariantList::getDelLen(uint32_t var_no)
{
    
    //   return vtAux[var_no].delLen;
    
    
    if(vtType[var_no] == SNP || vtType[var_no] == INS)
        return 0;
    else if(vtType[var_no] == DEL)
        return vt[var_no].alt;
    else //sv
    {
        uint32_t delLen;
        memcpy(&delLen, ins_sv_arr+vt[var_no].alt, sizeof(unsigned int));
        return delLen;
    }
}

uint32_t VariantList::getInsLen(uint32_t var_no)
{
    
    //   return vtAux[var_no].delLen;
    
    
    if(vtType[var_no] == SNP || vtType[var_no] == DEL)
        return 0;
    else if(vtType[var_no] == INS)
    {
        uint32_t insLen;
        memcpy(&insLen, ins_sv_arr+vt[var_no].alt, sizeof(unsigned int));
        return insLen;
    }
    else //sv
    {
        uint32_t insLen;
        memcpy(&insLen, ins_sv_arr+vt[var_no].alt+sizeof(unsigned int), sizeof(unsigned int));
        return insLen;
    }
}



uchar * VariantList::getPtrToFirstInsChar(uint32_t var_no)
{
    if(vtType[var_no] == INS)
    {
        return ins_sv_arr+vt[var_no].alt+sizeof(unsigned int);
      
    }
    else if(vtType[var_no] == SV)
    {
       return ins_sv_arr+vt[var_no].alt+2*sizeof(unsigned int);
    }
    else //SNP or DEL
        return NULL;
    
}




bool compareDelFromRight(singleDelFromRight a, singleDelFromRight b)
{
    if(a.nextRefPos < b.nextRefPos)
        return true;
    else if (a.nextRefPos > b.nextRefPos)
        return false;
    else
    {
        if(a.vtNo < b.vtNo)
            return true;
        else
            return false;
    }
}

void::VariantList::set_ref_size(uint32_t _ref_size){ref_size = _ref_size;}
/*********/

void VariantList::readData(FILE *in)
{
    cout << "vl "<< ref_size << endl;
    fflush(stdout);
    
    fread(&no_var, 1, sizeof(uint32_t), in);
    fread(&no_sv, 1, sizeof(uint32_t), in);
    fread(&no_ins, 1, sizeof(uint32_t), in);
    fread(&no_del, 1, sizeof(uint32_t), in);
    fread(&no_chars_inserted, 1, sizeof(uint32_t), in);
    
    cout << "no_var"<< endl;
    fflush(stdout);
    
    fread(&max_sv_del, 1, sizeof(uint32_t), in);
    fread(&max_del, 1, sizeof(uint32_t), in);
    
    cout << no_chars_inserted << " " <<  no_var <<
     " " << max_del << " " << no_del << " " << no_sv << endl;
    fflush(stdout);
    
    vtType = (uchar *)malloc((no_var+1)*sizeof(uchar)); //for sentinel
    vt = (variant*) malloc((no_var+1)*sizeof(variant)); //for sentinel
    
    
    uint32_t size_of_ins_sv_arr = 2*sizeof(uint32_t)*no_sv + sizeof(uint32_t)*no_ins + no_chars_inserted;
    ins_sv_arr = (unsigned char*) malloc((size_of_ins_sv_arr)*sizeof(unsigned char));
    
    
    fread(vtType, sizeof(uchar), no_var+1, in);
    fread(vt, sizeof(variant), no_var+1, in);
    fread(ins_sv_arr, sizeof(uchar), (2*sizeof(uint32_t)*no_sv + sizeof(uint32_t)*no_ins + no_chars_inserted), in);
    

    
    prevSV = (uint32_t * ) malloc((no_var)*sizeof(uint32_t));
    prevDeletion = (uint32_t * ) malloc((no_var)*sizeof(uint32_t));
    listSV = (uint32_t *) malloc (no_sv * sizeof(uint32_t));
    delFromRight = (struct singleDelFromRight * ) malloc (((no_del+no_sv) + 1) * sizeof(struct singleDelFromRight));
    
    
    fread(listSV, sizeof(uint32_t), no_sv, in);
    fread(prevSV, sizeof(uint32_t), no_var, in);
    fread(prevDeletion, sizeof(uint32_t), no_var, in);
    fread(delFromRight, sizeof(struct singleDelFromRight), (no_del+no_sv) + 1, in);
    
  
    
    
    uint32_t lut_step;
    fread(&lut_step, sizeof(uint32_t), 1, in);
    cout << lut_step<< endl;
    fflush(stdout);
    VT_LUT = (uint32_t * ) malloc ((uint32_t)(ref_size/VT_LUT_STEP) * sizeof(uint32_t));
    fread(VT_LUT, sizeof(uint32_t), ref_size/VT_LUT_STEP, in);
    
    
}


