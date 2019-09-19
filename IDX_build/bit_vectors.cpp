/*
 This file is a part of MuGI software distributed under GNU GPL 2 licence.
 The homepage of the MuGI project is http://sun.aei.polsl.pl/mugi/
 
 Authors: Agnieszka Danek, Sebastian Deorowicz
 
 Version: 1.0
 Date   : 2015-Mar
 */

#include "bit_vectors.h"



BitVectors::BitVectors(uint32 _ploidity, uint32 _no_gen, uint32 _no_var, string * _genomesNames)
{
    int32_t i;
    
    no_gen = _no_gen;
    ploidity = _ploidity;
    no_var = _no_var;
    
    no_bytes_per_gen = no_var/8 + ((no_var&7) > 0);

    
    
    
    genomesNames = new char*[no_gen*ploidity];
    
    bit_gen = new uchar*[no_gen*ploidity];
    bit_pos = new uint32_t[no_gen*ploidity];
    
    
    for (i=0; i < no_gen*ploidity; i++)
    {
        bit_gen[i] = new unsigned char[no_bytes_per_gen];
        for(int a = 0; a < no_bytes_per_gen; a++)
        {
            bit_gen[i][a] = 0;
        }
        genomesNames[i] = new char[_genomesNames[i].length()];
        strcpy(genomesNames[i], _genomesNames[i].c_str());
        bit_pos[i] = 0;
    }
    
    
    
    no_uint64_per_var = (no_gen*ploidity)/64 + (((no_gen*ploidity)%64) > 0);
    
    bit_var = new uint64_t*[no_var];
    for(uint32_t i = 0; i < no_var; ++i)
    {
        bit_var[i] = new uint64_t[no_uint64_per_var];
        for(int32_t j = 0; j < no_uint64_per_var; j++)
        {
            bit_var[i][j]=0LL;
        }
    }
    
    
}


void BitVectors::makeVarVec()
{
    uint64_t idx;
    uint64_t bitValue;
    
    for(uint32_t i = 0; i < no_var; ++i)
    {

        for(int32_t b = 0; b < no_gen*ploidity; b++)
        {
            idx = (unsigned long int)i/8UL; // * inBitV + (unsigned long int)b;
            bitValue = (bit_gen[b][idx]  & (1 << (7 - (i&7)))) != 0 ;
            
            bit_var[i][b/64] = bit_var[i][b/64] | (bitValue << (63-(b&63))) ;
            
        }
        
     
    }
    
    for(int32_t b = 0; b < no_gen*ploidity; b++)
    {
        delete bit_gen[b];
    }
    delete bit_gen;
    bit_gen = NULL;
    
}

void BitVectors::saveVarVec(FILE * file)
{
   
    uint32 no_id = (no_gen*ploidity);
    
    fwrite(&no_var, sizeof(uint32_t), 1,  file);
    fwrite(&no_id, sizeof(uint32_t), 1,  file);
    for(uint32_t i = 0; i < no_var; ++i)
    {
        for(int32_t b = 0; b < no_uint64_per_var; b++)
        {
            fwrite(&bit_var[i][b], sizeof(uint64_t), 1, file);
        }
    }

}


void  BitVectors::writeBits(uint32_t gen_id, unsigned char newBits, int32 number) //number:1-8
{
    int16 written;
    int32 leftToPut = number;
    
    
    while (leftToPut >0)
    {
        bit_gen[gen_id][bit_pos[gen_id]/8] = bit_gen[gen_id][bit_pos[gen_id]/8] | (newBits >> (bit_pos[gen_id]&7));
        
        
        written = ((8-(bit_pos[gen_id]&7)) <leftToPut? (8-(bit_pos[gen_id]&7)) : leftToPut);
        bit_pos[gen_id] = bit_pos[gen_id] + written;
        
        
        leftToPut = leftToPut - written;
        newBits = newBits << written;
    }
    
}



void BitVectors::writeBitBufToFile(uint32_t gen_id)
{
    string fileName;
    fileName.assign(genomesNames[gen_id]);
    fileName.append(".vt");
    FILE *fp = fopen(fileName.c_str(), "wb");
    
    
    
    
    
    fwrite(bit_gen[gen_id], 1, no_bytes_per_gen, fp);
    
    
    fclose(fp);
}



void BitVectors::removeVariantsSetInDelRegions(VariantList * vl)
{
    uint32_t delToPos = 0, currVt;
    uint64_t idx;
    bool bitValue;
    
    for(uint32_t b = 0; b < no_gen*ploidity; b++)
    {
       
        delToPos = 0;
        for (currVt = 0; currVt < no_var; currVt++)
        {
            
            
            idx = (unsigned long int)currVt/8UL; // * inBitV + (unsigned long int)b;
            bitValue = bit_gen[b][idx]  & (1 << (7 - (currVt&7)));
            
            if(bitValue)
            switch (vl->getType(currVt))
            {
                case SNP:
                    if(vl->getPos(currVt) < delToPos)
                        bit_gen[b][idx] = bit_gen[b][idx] & (~(1 << (7 - (currVt&7))));
                    break;
                case INS:
                    if(vl->getPos(currVt) < delToPos)
                        bit_gen[b][idx] = bit_gen[b][idx] & (~(1 << (7 - (currVt&7))));
                    break;
                case DEL:
                    if(vl->getPos(currVt)  < delToPos && vl->getPos(currVt)  + vl->getDelLen(currVt)  <= delToPos)
                        bit_gen[b][idx] = bit_gen[b][idx] & (~(1 << (7 - (currVt&7))));
                    else
                    {
                        delToPos = vl->getPos(currVt)  + vl->getDelLen(currVt) ;
                    }
                    break;
                case SV:
                    if(vl->getPos(currVt) < delToPos && vl->getPos(currVt) + vl->getDelLen(currVt) <= delToPos)
                        bit_gen[b][idx] = bit_gen[b][idx] & (~(1 << (7 - (currVt&7))));
                    else
                    {
                        delToPos = vl->getPos(currVt)  + vl->getDelLen(currVt) ;
                    }
                    break;
            }
            
            
        }
    }
    
}



int BitVectors::createCompressedCollection(uint32_t _blocking)
{
    makeVarVec();
    
    uint32_t max_uniq_parts = 4000000, max_uniq_var = 300000, max_var = 400000;
    uint32_t  cnt_var = 0, cnt_unique = 0, counter_0 = 0;
    uint32_t * temp_ui;
    uint64_t *temp_ull;
    
    blocking = _blocking;
    parts = (no_uint64_per_var + blocking - 1) / blocking;
    cnt_uniq_parts = 0, cnt_uniq_var = 0;
  
    all_var = (uint32_t *) malloc (max_var * sizeof(uint32_t));
    uniq_var = (uint32_t *) malloc (max_uniq_var * parts * sizeof(uint32_t));
    uniq_parts = (uint64_t *) malloc (max_uniq_parts * blocking * sizeof(uint64_t));
    
    int32 i, j;
    

    //typedef array<uint64_t, BLOCKING> part_t;
    typedef vector<unsigned long long> part_t;
    
    
    //typedef array<unsigned long long, NO_UINT64_PER_VARIANT> var_t;
    typedef vector<unsigned long long> var_t;
    
    unordered_map<part_t, uint32_t, hash<part_t>> data_part;
    
    //unordered_map<var_t, uint32_t, hash<var_t, 35>> data_var;
    
    unordered_map<var_t, unsigned int, hash<var_t>> data_var;

    vector<unsigned long long> arr_var;
    
    part_t arr_part, arr_part_0;
    for(j = 0; j < blocking; ++j)
		arr_part_0.push_back(0);
    
    for(cnt_var = 1; cnt_var <= no_var; cnt_var++)
	{
        
        arr_var.clear();
        arr_part.clear();
        
        
        
		for(i = 0; i < no_uint64_per_var; ++i)
		{
			
			//arr_var[i] = x;
            arr_var.vector::push_back(bit_var[cnt_var-1][i]);
            
        }
        
       
        if(data_var[arr_var] > 0)
        {
            //a = all_var[data_var[arr_var]-1];
            //all_var[cnt_var-1] = all_var[data_var[arr_var]-1];
            all_var[cnt_var-1] = all_var[data_var[arr_var]-1];
            if(cnt_var >= max_var)
            {
                max_var = max_var * 2;
                temp_ui = (unsigned int*) realloc(all_var, max_var*sizeof(unsigned int));
                if(temp_ui)
                    all_var = temp_ui;
                else
                    return 9;
            }
			continue;
		}

        
		data_var[arr_var] = cnt_var;
        
        all_var[cnt_var-1] = cnt_uniq_var;
        if(cnt_var >= max_var)
        {
            max_var = max_var * 2;
            temp_ui = (uint32_t*) realloc(all_var, max_var*sizeof(uint32_t));
            if(temp_ui)
                all_var = temp_ui;
            else
                return 9;
        }
        ++cnt_unique;
        
        for(i = 0; i < no_uint64_per_var; i += blocking)
		{
			arr_part = arr_part_0;
			for(j = 0; j < blocking && i+j < no_uint64_per_var; ++j)
				arr_part[j] = bit_var[cnt_var-1][i+j];
			
            if(arr_part == arr_part_0)
				counter_0++;
            
            if(data_part[arr_part] > 0)
            {
                uniq_var[cnt_uniq_var*parts+i/blocking] = (data_part[arr_part]-1);
                if(cnt_uniq_var + 1 >= max_uniq_var)
                {
                    max_uniq_var = max_uniq_var * 2;
                    temp_ui = (uint32_t*) realloc(uniq_var, max_uniq_var * parts *sizeof(uint32_t));
                    if(temp_ui)
                        uniq_var = temp_ui;
                    else
                        return 9;
                }
                
                continue;
            }
            else
            {
                
                cnt_uniq_parts++;
                
                //                uniq_var[cnt_uniq_var*parts+i/blocking] = arr_part
                //  for(j = 0; j < blocking && i+j < ull_per_var; ++j)
                //    uniq_partj[cnt_uniq_parts*3+j] = arr_part[j];
                
                memcpy(uniq_parts + (cnt_uniq_parts-1)*blocking, arr_part.data(),  j*sizeof(uint64_t));
                if(cnt_uniq_parts  >= max_uniq_parts)
                {
                    max_uniq_parts = max_uniq_parts * 2;
                    temp_ull = (uint64_t *) realloc(uniq_parts, max_uniq_parts * blocking * sizeof(uint64_t));
                    if(temp_ull)
                        uniq_parts = temp_ull;
                    else
                        return 9;
                }
                
                
                
                data_part[arr_part] = cnt_uniq_parts;
                uniq_var[cnt_uniq_var*parts+i/blocking] = cnt_uniq_parts-1;
                if(cnt_uniq_var + 1 >= max_uniq_var)
                {
                    max_uniq_var = max_uniq_var * 2;
                    temp_ui = (uint32_t*) realloc(uniq_var, max_uniq_var * parts *sizeof(uint32_t));
                    if(temp_ui)
                        uniq_var = temp_ui;
                    else
                        return 9;
                }
                
               
            }

		}
        cnt_uniq_var++;
        
    }
    
    for(uint32_t i = 0; i < no_var; ++i)
    {
        delete bit_var[i];
    }
    delete bit_var;
    bit_var = NULL;
    
    return 0;
}

void BitVectors::writeCompressedCollectionToFile(FILE * file)
{
    fwrite(&no_gen, sizeof(uint32_t), 1, file);
    
    fwrite(&ploidity, sizeof(uint32_t), 1, file);
    
    fwrite(&blocking, sizeof(uint32_t), 1, file);
    fwrite(&parts, sizeof(uint32_t), 1, file);
    fwrite(&no_var, sizeof(uint32_t), 1, file);
    fwrite(&cnt_uniq_var, sizeof(uint32_t), 1, file);
    fwrite(&cnt_uniq_parts, sizeof(uint32_t), 1, file);
    
    fwrite(all_var, sizeof(uint32_t), no_var, file);
    fwrite(uniq_var, sizeof(uint32_t), cnt_uniq_var * parts, file);
    fwrite(uniq_parts, sizeof(unsigned long long), cnt_uniq_parts * blocking, file);
    
}


void BitVectors::freeBitPos()
{
    delete bit_pos;
}


/******/
void BitVectors::readData(FILE *in)
{
    
    fread(&no_gen, sizeof(uint32_t), 1, in);
    fread(&ploidity, sizeof(uint32_t), 1, in);
    // inBitV = inBitV * ploidity;
    
    fread(&blocking, sizeof(uint32_t), 1, in);
    fread(&parts, sizeof(uint32_t), 1, in);
    fread(&no_var, sizeof(uint32_t), 1, in);
    fread(&cnt_uniq_var, sizeof(uint32_t), 1, in);
    fread(&cnt_uniq_parts, sizeof(uint32_t), 1, in);
    
    all_var = new uint32_t[no_var];
    uniq_var = new uint32_t[cnt_uniq_var * parts];
    uniq_parts= new uint64_t[cnt_uniq_parts * blocking];
    
    
    fread(all_var, sizeof(uint32_t), no_var, in);
    fread(uniq_var, sizeof(uint32_t), cnt_uniq_var * parts, in);
    fread(uniq_parts, sizeof(unsigned long long), cnt_uniq_parts * blocking, in);
    
}



