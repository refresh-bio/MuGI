/*
 This file is a part of MuGI software distributed under GNU GPL 2 licence.
 The homepage of the MuGI project is http://sun.aei.polsl.pl/mugi/
 
 Authors: Agnieszka Danek, Sebastian Deorowicz
 
 Version: 1.0
 Date   : 2015-Mar
 */

#include "sa2.h"

SA2::SA2(uint32_t sparse, uint32_t _max)
{
    n_sparse = sparse;
    suffixArray2Count = 0;
    suffixArray2Max = _max;
   
    
    suffixArray2 = (suffix2 *) malloc (suffixArray2Max * sizeof(suffix2));
    
    LUT = (uint32_t * ) malloc ((uint32_t)(LUT_MAX+1) * sizeof(uint32_t));
}



void SA2::newSuffix2(unsigned int  prePos,  unsigned int newVtList, unsigned int gap,  unsigned int vtCovered, unsigned int firstVar, unsigned int currSV)
{
    //  cout  << " " << prePos << " " << newVtList << " " <<  gap << " " <<  vtCovered << " " << firstVar  << endl;
    //  fflush(stdout);
   
    bool firstVarNotZero = firstVar > 0 ? true : false;
    
    bool includeVar = true;
     unsigned int back;
    
#ifdef FILTERING

    
    if(Data::vl->vtType[firstVar] != SV)  //if== SV, includeVar = true (because it may be SV+ins in del region
    {
        
        //meaningful SV
        bool inDelReg;
        int cSV = currSV, lastSV, firstSV;
        while(cSV >= 0 && Data::vl->vt[Data::vl->listSV[cSV]].prePos  >= prePos)
            cSV--;
        lastSV = cSV;
        if(cSV < 0)
        {
            firstSV = 0;
        }
        else
        {
            while (cSV > 0 && prePos - Data::vl->vt[Data::vl->listSV[cSV]].prePos <= Data::vl->get_maxSVdel())
                cSV--;
            while (cSV < Data::vl->get_no_sv() && Data::vl->vt[Data::vl->listSV[cSV]].prePos + Data::vl->getDelLen(Data::vl->listSV[cSV]) <= prePos)
                cSV++;
            firstSV = cSV;
        }
        unsigned int dels, prevDelCount = 0;
        
        //check how many meaningful DELs before
        back = 1;
        if(firstVarNotZero)
        {
            while (firstVar >= back && prePos - Data::vl->vt[firstVar - back].prePos <=  Data::vl->get_maxDel())
            {
                if(Data::vl->getType(firstVar - back) == DEL)// && vt[firstVar - back].prePos + vt[firstVar - back].delLen > prePos)
                {
                    prevDelCount++;
                }
                back++;
            }
            
            back--;
            while(back)
            {
                if(Data::vl->getType(firstVar - back) == DEL)
                {
                    if(Data::vl->getPos(firstVar - back) + Data::vl->getDelLen(firstVar - back) <= prePos)
                        prevDelCount--;
                    else
                        break;
                }
                back--;
            }
        }


        
        
        includeVar = false;
        int i;
        bool suffix, bitvector;
        unsigned long int temp;//, idx;
        for(int b=0; b < Data::bv->getNoId(); b++)
        {
            
            //check if first chars not deleted
            inDelReg = false;
            back = 1;
            dels = prevDelCount;
            while(dels)
            {
                if(Data::vl->getType(firstVar - back) == DEL)
                {
                    dels--;
                    
                    //idx = ((unsigned long int)(firstVar - back)/8UL)*(unsigned long int)bv->getNoGen() + (unsigned long int)b;
                    //inDelReg = inputBuffer_t[idx] & (1 << (7 - (firstVar - back)%8));
                    
                    inDelReg = Data::bv->uniq_parts[Data::bv->uniq_var[Data::bv->all_var[firstVar - back]*Data::bv->parts+b/64/Data::bv->blocking]*Data::bv->blocking+(b/64)%Data::bv->blocking] & (1LL << (63-(b&63)));
                    
                    
                    
                    if(inDelReg && Data::vl->vt[firstVar - back].prePos + Data::vl->getDelLen(firstVar - back) > prePos)
                    {
                        break;
                    }
                    else
                        inDelReg = false;
                }
                back++;
            }
            
            if(!inDelReg)
                for(int sv = lastSV; sv >= (int)firstSV; sv--)
                {
                    {
                       
                        
                        inDelReg = Data::bv->uniq_parts[Data::bv->uniq_var[Data::bv->all_var[Data::vl->listSV[sv]]*Data::bv->parts+b/64/Data::bv->blocking]*Data::bv->blocking+(b/64)%Data::bv->blocking] & (1LL << (63-(b&63)));
                        
                        if(inDelReg)
                        {
                            if(Data::vl->getPos(Data::vl->listSV[sv]) +  Data::vl->getDelLen(Data::vl->listSV[sv]) <= prePos)
                                inDelReg = false;
                            else
                                break;
                        }
                    }
                    
                }
            
            ///

            
            if(!inDelReg)
            {
                for(i=0; i < vtCovered; i++)
                {
                    suffix = newVtList & (1 << (NO_BITS_UINT - i - 1));
                    
                    temp = (unsigned long int)firstVar + (unsigned long int)i;
                 
                    bitvector = Data::bv->uniq_parts[Data::bv->uniq_var[Data::bv->all_var[temp]*Data::bv->parts + b/64/Data::bv->blocking]*Data::bv->blocking + (b/64)%Data::bv->blocking] & (1LL << (63-(b&63)));
                    
                    if(suffix != bitvector)
                        break;
                }
                if(i == vtCovered)
                {
                    includeVar = true;
                    break;
                }
            }
            
            
        }
    }
#endif

    char ch;
    if(includeVar)
    {
        
        
        suffixArray2[suffixArray2Count].prevPos = prePos;
        suffixArray2[suffixArray2Count].gap = gap;
        suffixArray2[suffixArray2Count].vtList = newVtList;
        
        ch = getFirstChar(&suffixArray2[suffixArray2Count]);
        if(ch !='N' && ch != 'n')
        {
            suffixArray2Count++;
            if(suffixArray2Count == suffixArray2Max)
            {
                suffixArray2 = (suffix2 *) realloc (suffixArray2, suffixArray2Max * 2 * sizeof(suffix2));
                suffixArray2Max = (unsigned int) (2 * suffixArray2Max);
            }
        }
        
    
    }
}


void SA2::sort()
{
    
    tbb::parallel_sort(suffixArray2, suffixArray2+suffixArray2Count, compareSuffix2);
    
}



void  SA2::write_sa(FILE * file)
{
    
    fwrite(&suffixArray2Count, 1, sizeof(uint32_t), file);
    
    
   // ofstream myfile;
  //  myfile.open ("check-sa2.txt");
    
    
    char *kmer;
    uint32_t nextPos, currWrittenSuffix = 0, i;
    int32_t encoded, currIdx = 0;
    for(i = 0; i < suffixArray2Count; i++)
    {
        kmer = getSuffix(&suffixArray2[i], &nextPos);
        
        encoded = Coding::encodeSuffix(kmer);
 
        //if(encoded != 0xFFFFFFFF)  //only ACGTacgt
        {
            // if(suffixArray2[i].prevPos % n_sparse == 0)
            {
                fwrite(&suffixArray2[i].prevPos, sizeof(unsigned int), 1, file);
                fwrite(&suffixArray2[i].vtList, sizeof(unsigned int), 1, file);
                fwrite(&suffixArray2[i].gap, sizeof(unsigned int), 1, file);
            
                
                //create LUT for SA2
                //encoded = encodeSuffix(kmer);
                if(encoded >= currIdx)
                {
                    for(int j = currIdx; j <=encoded; j++)
                        LUT[j] = currWrittenSuffix;
                    currIdx = encoded + 1;
                }
                currWrittenSuffix++;
            }
            //  currSuffix++;
            
        }
        free(kmer);
  
    }
    for(int j = currIdx; j <=LUT_MAX; j++)
        LUT[j] = currWrittenSuffix;
    
    //myfile.close();
}


void  SA2::write_sa_lut(FILE * file)
{
    fwrite(LUT, sizeof(uint32_t), LUT_MAX+1, file);
}


void SA2::readData(FILE *in)
{
    fread(&suffixArray2Count, 1, sizeof(uint32_t), in);
    suffixArray2 = (suffix2 *) malloc (suffixArray2Count * sizeof(suffix2));
    fread(suffixArray2, sizeof(unsigned int), suffixArray2Count*3, in);
    
    LUT = (uint32_t * ) malloc ((uint32_t)(LUT_MAX+1) * sizeof(uint32_t));
    fread(LUT, sizeof(uint32_t), LUT_MAX+1, in);
    
    
}



