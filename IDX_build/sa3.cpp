/*
 This file is a part of MuGI software distributed under GNU GPL 2 licence.
 The homepage of the MuGI project is http://sun.aei.polsl.pl/mugi/
 
 Authors: Agnieszka Danek, Sebastian Deorowicz
 
 Version: 1.0
 Date   : 2015-Mar
 */

#include "sa3.h"

SA3::SA3(uint32_t sparse, uint32_t _max)
{
    n_sparse = sparse;
    suffixArray3Count = 0;
    suffixArray3Max = _max;
    
    
    suffixArray3 = (suffix3 *) malloc (suffixArray3Max * sizeof(suffix3));
    LUT = (uint32_t * ) malloc ((uint32_t)(LUT_MAX+1) * sizeof(uint32_t));
}

void SA3::newSuffix3(unsigned int  prePos,  unsigned int vtListSize, unsigned int *newVtList, unsigned int gap,  unsigned int vtCovered, unsigned int firstVar, unsigned int lastFilt, unsigned int currSV)
{
    
    //    cout  << " " << prePos << " " << newVtList[0] << " " << gap << " " << vtCovered << " " << firstVar << " " << currSV << endl;
    //  fflush(stdout);
    
    bool firstVarNotZero = firstVar > 0 ? true : false;
    bool includeVar = true;
    
    unsigned long int refPos;
    unsigned int back;
#ifdef FILTERING
    
    bool inDelReg;
    int cSV = currSV, lastSV, firstSV;
    while(cSV >= 0 && Data::vl->vt[Data::vl->listSV[cSV]].prePos >= prePos)
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
        
        while (cSV < Data::vl->get_no_sv() && Data::vl->vt[Data::vl->listSV[cSV]].prePos +  Data::vl->getDelLen(Data::vl->listSV[cSV])  <= prePos)
            cSV++;
        firstSV = cSV;
    }
    
    unsigned int dels, prevDelCount = 0;
    //check how many meaningful DELs before
    back = 1;
    if(firstVarNotZero)
    {
        while (firstVar >= back &&  prePos - Data::vl->vt[firstVar - back].prePos <= Data::vl->get_maxDel())
        {
            if(Data::vl->getType(firstVar - back) == DEL)// && vt[firstVar - back].prePos + Data::vl->getDelLen(firstVar - back) > prePos)
            {
                prevDelCount++;
            }
            back++;
        }
        back--;
        while(back)
        {
            if(Data::vl->getType(firstVar - back)== DEL)
            {
                if(Data::vl->vt[firstVar - back].prePos + Data::vl->getDelLen(firstVar - back) <= prePos)
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
    unsigned long int temp,  skippedVt;//idx,
    
    
    
    
    for(int b = 0; b <  Data::bv->getNoId(); b++)
    {
        //check if first chars not deleted
        inDelReg = false;
        back = 1;
        dels = prevDelCount;
        while(dels)
        {
            if(Data::vl->getType(firstVar - back)== DEL)
            {
                dels--;
                //idx = ((unsigned long int)(firstVar - back)/8UL)*(unsigned long int) Data::bv->getNoGen() + (unsigned long int)b;
                //inDelReg = inputBuffer_t[idx] & (1 << (7 - (firstVar - back)%8));

                inDelReg = Data::bv->uniq_parts[Data::bv->uniq_var[Data::bv->all_var[firstVar - back]*Data::bv->parts+b/64/Data::bv->blocking]*Data::bv->blocking+(b/64)%3] & (1LL << (63-(b&63)));

                
                if(inDelReg && Data::vl->vt[firstVar - back].prePos + Data::vl->getDelLen(firstVar - back) > prePos)
                {
                    //cout << "del: " << prePos << " "<< back << " " << vt[firstVar - back].prePos << " " <<Data::vl->getDelLen(firstVar - back) << endl;
                    break;
                }
                else
                    inDelReg = false;
            }
            back++;
        }
        
        if(!inDelReg && Data::vl->no_sv > 0)
            for(int sv = lastSV; sv >= (int)firstSV; sv--)
            {
                //if(Data::vt[Data::vl->listSV[sv]].prePos +  Data::vl->getDelLen(Data::vl->listSV[sv]) > prePos)
                {
                    //idx = ((unsigned long int)Data::vl->listSV[sv]/8UL)*(unsigned long int) Data::bv->getNoGen() + (unsigned long int)b;
                    //inDelReg = inputBuffer_t[idx] & (1 << (7 - Data::vl->listSV[sv]%8));
                    
                    inDelReg = Data::bv->uniq_parts[Data::bv->uniq_var[Data::bv->all_var[Data::vl->listSV[sv]]*Data::bv->parts+b/64/Data::bv->blocking]*Data::bv->blocking+(b/64)%3] & (1LL << (63-(b&63)));
                    
                    if(inDelReg)
                    {
                        if(Data::vl->vt[Data::vl->listSV[sv]].prePos +  Data::vl->getDelLen(Data::vl->listSV[sv]) <= prePos)
                            inDelReg = false;
                        else
                            break;
                    }
                }
            }
        
        
        
        
        skippedVt = 0; refPos = 0;
        if(!inDelReg)
        {   for(i=0; i < vtCovered; i++)
            {
                suffix = newVtList[i/NO_BITS_UINT] & (1 << (NO_BITS_UINT - i%NO_BITS_UINT  - 1));
                
                temp = (uint32_t)firstVar + (uint32_t)i + (uint32_t)skippedVt;
                
                while ((Data::vl->getType(temp)== SNP && Data::vl->vt[temp].prePos + 1  < refPos) || (Data::vl->getType(temp)== INS && Data::vl->vt[temp].prePos + 1 + Data::vl->getDelLen(temp) < refPos) || ((Data::vl->getType(temp)== DEL || Data::vl->getType(temp)== SV) && Data::vl->vt[temp].prePos + 1 + Data::vl->getDelLen(temp) <= refPos))
                {
                    skippedVt++;
                    temp++;
                }
                if (suffix && Data::vl->vt[temp].prePos + 1 + Data::vl->getDelLen(temp) > refPos)
                    refPos = Data::vl->vt[temp].prePos + 1 +Data::vl->getDelLen(temp);
                
                
                
                //idx = ((unsigned long int)temp/8UL)*(unsigned long int) Data::bv->getNoGen() + (unsigned long int)b;
                //bitvector = inputBuffer_t[idx] & (1 << (7 - temp%8));
           
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
#endif
    
    /*
     CREATED NORMAL NOT  SPARSE
    if(prePos % n_sparse != 0)
        if(prePos != vt[firstVar].prePos || vt[firstVar].type == SNP)
            includeVar = false;
     */
    
    char ch;
    
    if(gap > 0 && Data::vl->vtType[firstVar] == SV)
    {
        //if(includeVar == false)
        //    cout << "*SA3 new suffix: "<< prePos << " " << gap << " " << newVtList[0] << endl;
        includeVar = true;
    }
    
    if(includeVar)
    {
        
        
        suffixArray3[suffixArray3Count].prevPos = prePos;
        suffixArray3[suffixArray3Count].gap = gap;
        suffixArray3[suffixArray3Count].vtListSize = vtListSize;
        suffixArray3[suffixArray3Count].vtList = (unsigned int *) malloc (vtListSize * sizeof(unsigned int));
        
        for(int i = 0; i < vtListSize; i++)
            suffixArray3[suffixArray3Count].vtList[i] = newVtList[i];
        
        
        ch = getFirstChar(&suffixArray3[suffixArray3Count]);
        if(ch !='N' && ch != 'n')
        {
            suffixArray3Count++;
            if(suffixArray3Count == suffixArray3Max)
            {
                suffixArray3 = (suffix3 *) realloc (suffixArray3, suffixArray3Max * 2 * sizeof(suffix3));
                suffixArray3Max = (unsigned int) (2 * suffixArray3Max);
            }
        }
        else
        {
            free(suffixArray3[suffixArray3Count].vtList);
        }
    }
    
    
}

void SA3::sort()
{
    

    tbb::parallel_sort(suffixArray3, suffixArray3+suffixArray3Count, compareSuffix3);
    
 
}

void  SA3::write_sa(FILE * file)
{
    
    fwrite(&suffixArray3Count, 1, sizeof(uint32_t), file);
    
   // ofstream myfile;
   // myfile.open ("check-sa3.txt");
    
    
    char *kmer;
    uint32_t nextPos, currWrittenSuffix = 0, i;
    int32_t encoded, currIdx = 0;
    for(i = 0; i < suffixArray3Count; i++)
    {
      
        kmer = getSuffix(&suffixArray3[i], &nextPos);
        //create LUT for SA3
        encoded = Coding::encodeSuffix(kmer);
        
      //  if(encoded != 0xFFFFFFFF)  //only ACGTacgt
        {
            //if(suffixArray3[i].prevPos % n_sparse == 0)
            {
                fwrite(&suffixArray3[i].prevPos, sizeof(unsigned int), 1, file);
                fwrite(&suffixArray3[i].gap, sizeof(unsigned int), 1, file);
                fwrite(&suffixArray3[i].vtListSize, sizeof(unsigned int), 1, file);
                for(int j = 0; j < suffixArray3[i].vtListSize; j++)
                    fwrite(&suffixArray3[i].vtList[j], sizeof(unsigned int), 1, file);
                
              //  myfile << currWrittenSuffix << " " << kmer << " " << suffixArray3[i].prevPos<< " " << suffixArray3[i].gap<< " " << suffixArray3[i].vtListSize<< " " << suffixArray3[i].vtList[0] <<  endl;
                
                //create LUT for SA3
                if(encoded >= currIdx)
                {
                    for(int j = currIdx; j <=encoded; j++)
                        LUT[j] = currWrittenSuffix;
                    currIdx = encoded + 1;
                }
                currWrittenSuffix++;
            }
            // currSuffix++;
        }
        free(kmer);
    }
    for(int j = currIdx; j <=LUT_MAX; j++)
        LUT[j] = currWrittenSuffix;
    
   // myfile.close();
}


void  SA3::write_sa_lut(FILE * file)
{
    fwrite(LUT, sizeof(uint32_t), LUT_MAX+1, file);
}



void SA3::readData(FILE *in)
{
    fread(&suffixArray3Count, 1, sizeof(uint32_t), in);
    suffixArray3 = (suffix3 *) malloc (suffixArray3Count * sizeof(suffix3));
    for(uint32_t i = 0; i < suffixArray3Count; i++)
    {
        
        
        fread(&suffixArray3[i].prevPos, sizeof(unsigned int), 1, in);
        fread(&suffixArray3[i].gap, sizeof(unsigned int), 1, in);
        fread(&suffixArray3[i].vtListSize, sizeof(unsigned int), 1, in);
        suffixArray3[i].vtList = (uint32_t * ) malloc (suffixArray3[i].vtListSize * sizeof(uint32_t));
        for(int j = 0; j < suffixArray3[i].vtListSize; j++)
            fread(&suffixArray3[i].vtList[j], sizeof(unsigned int), 1, in);
        
    }
    
    
    LUT = (uint32_t * ) malloc ((uint32_t)(LUT_MAX +1) * sizeof(uint32_t));
    fread(LUT, sizeof(uint32_t), LUT_MAX+1, in);
    
    
    
    
}



