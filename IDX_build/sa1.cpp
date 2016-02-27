/*
 This file is a part of MuGI software distributed under GNU GPL 2 licence.
 The homepage of the MuGI project is http://sun.aei.polsl.pl/mugi/
 
 Authors: Agnieszka Danek, Sebastian Deorowicz
 
 Version: 1.0
 Date   : 2015-Mar
 */

#include "sa1.h"



SA1::SA1(uint32_t sparse, uint32_t _max, uint32_t _kl, uint32_t _no_var,  BitVectors * _bv, VariantList * _vl, CFastaFile * _fasta)
{
    n_sparse = sparse;
    suffixArray1Count = 0;
    suffixArray1Max = _max;
    kl = _kl;
    no_var = _no_var;
    vl = _vl;
    bv = _bv;
    fasta = _fasta;
    filtered = 0;
    
    suffixArray1 = (suffix1 *) malloc (suffixArray1Max * sizeof(suffix1));
    
    LUT = (uint32_t * ) malloc ((uint32_t)(LUT_MAX+1) * sizeof(uint32_t));
}


void SA1::newSuffix1(unsigned int  prePos,  unsigned int newVtList, unsigned int vtCovered, unsigned int firstVar, unsigned int currSV)
{
    
    bool firstVarNotZero = firstVar > 0 ? true : false;
    bool includeVar = true;
    
    unsigned int back;
#ifdef FILTERING
    
    bool inDelReg;
    int cSV = currSV, lastSV, firstSV;
    while(cSV >= 0 && (vl->vt[vl->listSV[cSV]].prePos  >= prePos || prePos == 0xFFFFFFFF))
        cSV--;
    lastSV = cSV;
    if(cSV < 0)
    {
        firstSV = 0;
    }
    else
    {
        while (cSV > 0 && prePos - vl->vt[vl->listSV[cSV]].prePos <= vl->get_maxSVdel())
            cSV--;
        while (cSV < vl->get_no_sv() && vl->vt[vl->listSV[cSV]].prePos + vl->getDelLen(vl->listSV[cSV]) <= prePos)
            cSV++;
        firstSV = cSV;
    }
    unsigned int dels, prevDelCount = 0;
    
    //check how many meaningful DELs before
    back = 1;
    if(firstVarNotZero)
    {
        while (firstVar >= back && prePos - vl->vt[firstVar - back].prePos <=  vl->get_maxDel())
        {
            if(vl->getType(firstVar - back) == DEL)// && vt[firstVar - back].prePos + vt[firstVar - back].delLen > prePos)
            {
                prevDelCount++;
            }
            back++;
        }
        
        back--;
        while(back)
        {
            if(vl->getType(firstVar - back) == DEL)
            {
                if(vl->getPos(firstVar - back) + vl->getDelLen(firstVar - back) <= prePos)
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
    for(int b=0; b <  bv->getNoId(); b++)
    {
        
        //check if first chars not deleted
        inDelReg = false;
        back = 1;
        dels = prevDelCount;
        while(dels)
        {
            if(vl->getType(firstVar - back) == DEL)
            {
                dels--;
                
                //idx = ((unsigned long int)(firstVar - back)/8UL)*(unsigned long int)bv->getNoGen() + (unsigned long int)b;
                //inDelReg = inputBuffer_t[idx] & (1 << (7 - (firstVar - back)%8));
                
                inDelReg = bv->uniq_parts[bv->uniq_var[bv->all_var[firstVar - back]*bv->parts+b/64/bv->blocking]*bv->blocking+(b/64)%bv->blocking] & (1LL << (63-(b&63)));

               
                
                if(inDelReg && vl->vt[firstVar - back].prePos + vl->getDelLen(firstVar - back) > prePos)
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
                   // idx = ((unsigned long int)vl->listSV[sv]/8UL)*(unsigned long int)bv->getNoGen() + (unsigned long int)b;
                   // inDelReg = inputBuffer_t[idx] & (1 << (7 - vl->listSV[sv]%8));
                    
                    inDelReg = bv->uniq_parts[bv->uniq_var[bv->all_var[vl->listSV[sv]]*bv->parts+b/64/bv->blocking]*bv->blocking+(b/64)%bv->blocking] & (1LL << (63-(b&63)));
                    
                    if(inDelReg)
                    {
                        if(vl->getPos(vl->listSV[sv]) +  vl->getDelLen(vl->listSV[sv]) <= prePos)
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
                //idx = ((unsigned long int)temp/8UL)*(unsigned long int)bv->getNoGen() + (unsigned long int)b;
                //bitvector = inputBuffer_t[idx] & (1 << (7 - temp%8));
                
                bitvector = bv->uniq_parts[bv->uniq_var[bv->all_var[temp]*bv->parts + b/64/bv->blocking]*bv->blocking + (b/64)%bv->blocking] & (1LL << (63-(b&63)));
                
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
    
    
    char  ch;
    
    
    if(includeVar)
    {
     
        //cout << suffixArray1Count << " " <<prePos << " " << newVtList << endl;
        //if(suffixArray1Count == 88)
        //    cout << endl;
        suffixArray1[suffixArray1Count].prevPos = prePos;
        suffixArray1[suffixArray1Count].vtList = newVtList;
        
        
#ifdef CREATE_ALL_INDEXES
        if(Data::vl->vt[firstVar].prePos == prePos && Data::vl->getType(firstVar) != SNP)
            suffixArray1[suffixArray1Count].startWithVariant = 1;
        else
            suffixArray1[suffixArray1Count].startWithVariant = 0;
#endif
        
        
        ch = getFirstChar(&suffixArray1[suffixArray1Count]);
        if(ch !='N' && ch != 'n')
        {
            suffixArray1Count++;
            if(suffixArray1Count == suffixArray1Max)
            {
                suffixArray1 = (suffix1 *) realloc (suffixArray1, suffixArray1Max * 1.5 * sizeof(suffix1));
                suffixArray1Max = (unsigned int) (1.5 * suffixArray1Max);
            }
        }

    }
    else
        filtered++;
    
}


void SA1::sort()
{
  
    
    tbb::parallel_sort(suffixArray1, suffixArray1+suffixArray1Count, compareSuffix1);
}




void  SA1::write_sa(FILE * file)
{
    
    fwrite(&suffixArray1Count, 1, sizeof(uint32_t), file);
    
    //ofstream myfile;
    //myfile.open ("check-sa1.txt");
    
    
    char *kmer;
    uint32_t nextPos, currWrittenSuffix = 0;
    int64_t encoded, currIdx = 0;
    for(uint32_t i = 0; i < suffixArray1Count; i++)
    {
        kmer = getSuffix(&suffixArray1[i], &nextPos);
    
        encoded = Coding::encodeSuffix(kmer);
        //if(suffixArray1[i].prevPos >= 25597709 || suffixArray1[i].prevPos <= 25597712)
        //    cout << "";
       // if(encoded != 0xFFFFFFFF)  //only ACGTacgt
        {
            //if(suffixArray1[i].prevPos % n_sparse == 0)
            {
                fwrite(&suffixArray1[i].prevPos, sizeof(uint32_t), 1, file);
                fwrite(&suffixArray1[i].vtList, sizeof(uint32_t), 1, file);
                
             //  myfile << i << " " << kmer << " " << suffixArray1[i].prevPos<< " "<< suffixArray1[i].vtList <<  endl;
                
                //create LUT for SA1
               // encoded = Coding::encodeSuffix(kmer);
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
        //  fprintf(res, "%s\n", suffixArray1[i].kmer);
        //  cout << suffixArray1[i].kmer << " " << suffixArray1[i].nextPos << endl;
        // kmer = getSuffix(&suffixArray1[i], &nextPos);
        // cout << kmer << endl;
        // free(kmer);
    }
    for(int64 j = currIdx; j <=LUT_MAX; j++)
        LUT[j] = currWrittenSuffix;
    
  //  myfile.close();
}


void  SA1::write_sa_lut(FILE * file)
{
     fwrite(LUT, sizeof(uint32_t), LUT_MAX+1, file);
}



void  SA1::write_sa(FILE * file, int32_t sparse)
{
    
    uint32_t count = 0;
   // fwrite(&suffixArray1Count, 1, sizeof(uint32_t), file);
    
    //  ofstream myfile;
    //  myfile.open ("check-sa1.txt");
    
    for(uint32_t i = 0; i < suffixArray1Count; i++)
    {
#ifdef CREATE_ALL_INDEXES
        if(suffixArray1[i].prevPos % sparse == 0 || suffixArray1[i].startWithVariant == 1)
#else
        if(suffixArray1[i].prevPos % sparse == 0)
#endif
            count++;
    }
    fwrite(&count, 1, sizeof(uint32_t), file);
    cout << "sa1 sparse " << sparse << " count: " << count << endl;
    char *kmer;
    uint32_t nextPos, currWrittenSuffix = 0;
    int64_t encoded, currIdx = 0;
    for(uint32_t i = 0; i < suffixArray1Count; i++)
    {
#ifdef CREATE_ALL_INDEXES
        if(suffixArray1[i].prevPos % sparse == 0 || suffixArray1[i].startWithVariant == 1)
#else
        if(suffixArray1[i].prevPos % sparse == 0)
#endif
        {
            kmer = getSuffix(&suffixArray1[i], &nextPos);
        
            encoded = Coding::encodeSuffix(kmer);
             // if(encoded != 0xFFFFFFFF)  //only ACGTacgt
            {
                //if(suffixArray1[i].prevPos % n_sparse == 0)
                {
                    fwrite(&suffixArray1[i].prevPos, sizeof(uint32_t), 1, file);
                    fwrite(&suffixArray1[i].vtList, sizeof(uint32_t), 1, file);
                    
                    //     myfile << currWrittenSuffix << " " << kmer << " " << suffixArray1[i].prevPos<< " "<< suffixArray1[i].vtList <<  endl;
                    
                    //create LUT for SA1
                    if(encoded >= currIdx)
                    {
                        for(int j = currIdx; j <=encoded; j++)
                            LUT[j] = currWrittenSuffix;
                        currIdx = encoded + 1;
                    }
                    currWrittenSuffix++;
                }
                
                
            }
            free(kmer);
          
        }
    }
    for(int64 j = currIdx; j <=LUT_MAX; j++)
        LUT[j] = currWrittenSuffix;
    cout << "sa1 sparse " << sparse << " currWrittenSuffix: " << currWrittenSuffix << endl;
    
    //   myfile.close();
}


void SA1::readData(FILE *in)
{
    fread(&suffixArray1Count, 1, sizeof(uint32_t), in);
      uint64_t temp_size = (uint64_t)suffixArray1Count*2*sizeof(uint32_t) ;
    
    suffixArray1 = (suffix1 *) malloc (temp_size);
    
    if(!suffixArray1)
        cout << strerror(errno) << endl;
    
    
    //to work on all mac os -> bug in fread (some mac os), which cannot read more than 2GB at once
    int i = 0;
    size_t  a ;
  
    
    while (temp_size > 2147483647)
    {
          a = fread(((char *)suffixArray1)+i*2147483647, 1, 2147483647, in);
        temp_size = temp_size - 2147483647;
        i++;
	}
     a = fread(((char *)suffixArray1)+i*2147483647, 1, temp_size, in);
    
    
 
    LUT = (uint32_t * ) malloc ((uint32_t)(LUT_MAX+1) * sizeof(uint32_t));
    fread(LUT, sizeof(uint32_t), LUT_MAX+1, in);
    
  }

