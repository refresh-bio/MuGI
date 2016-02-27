/*
 This file is a part of MuGI software distributed under GNU GPL 2 licence.
 The homepage of the MuGI project is http://sun.aei.polsl.pl/mugi/
 
 Authors: Agnieszka Danek, Sebastian Deorowicz
 
 Version: 1.0
 Date   : 2015-Mar
 */


#include "SuffixArray1.h"


int32_t compareReadSuffix(unsigned char* read, uint32_t lenToCheck, suffix1 suffix)
{
    uint32_t insLen, delLen;
    int32_t pos = 0;
   
    uint32_t refPos = suffix.prevPos+1;
    
    ///BINARY SEARCH FOR FIRST VT
    uint32_t vtNo = 0,m, l, r; // l = 0, r = noVar-1, m;
    

    uint32_t temp = suffix.prevPos/VT_LUT_STEP;

    
    
    if(refSize >= VT_LUT_STEP)
    {
        l = VT_LUT[temp++];
        r = VT_LUT[temp];
    }
    else
    {
        l = 0;
        r = noVar;
    }
    
    
    while (l < r) {
        m = (r + l) / 2;
        if ((vt[m].prePos <  suffix.prevPos || vt[m].prePos == 0xFFFFFFFF) && suffix.prevPos != 0xFFFFFFFF)            l = m + 1;
        else
            r = m;
    }
    vtNo = l;
    
    
    uint32_t vtCovered = 0, len;
    
    
    
    while (pos < lenToCheck && refPos < refSize)
    {
       
        if( vt[vtNo].prePos+1 > refPos)
        {
            
            len = vt[vtNo].prePos - refPos + 1 < lenToCheck - pos ? vt[vtNo].prePos - refPos + 1 : lenToCheck - pos ;
            
        
            for(int32_t i = 0; i < len; i++)
                if(read[pos] == getRefFinal(refPos))
                {
                    pos++;
                    //chCovered++;
                    refPos++;
                }
                else
                {
                    if(refPos == refSize)
                        return 1;
                    return read[pos] > getRefFinal(refPos) ? 1 : -1;
                }
               
            
            
            
          
        }
        else
        {
            if(suffix.vtList & (1 << (NO_BITS_UINT - vtCovered - 1))) //variant found/introduced
            {
                switch (vtType[vtNo])
                {
                    case SNP:
                        //introduce variant
                        if(read[pos] == (vt[vtNo].alt))
                        {
                            //chCovered++;
                            refPos++;
                            vtCovered++;
                            vtNo++;
                            pos++;
                        }
                        else
                            return read[pos] > (vt[vtNo].alt) ? 1 : -1;
                        

                        break;
                    case DEL:
                        //introduce variant
                        refPos = refPos + vt[vtNo].alt;
                        vtCovered++;
                        vtNo++;
                        
                        break;
                    case INS:
                        //introduce variant
                        LoadUInt(ins_sv_arr+vt[vtNo].alt, insLen);
                        
                        for (int32_t c = 0; c < insLen && pos < lenToCheck; c++)
                        {                            if(read[pos] == (ins_sv_arr[vt[vtNo].alt+sizeof(uint32_t)+c]))
                            {
                                pos++;
                            }
                            else
                                return read[pos] > (ins_sv_arr[vt[vtNo].alt+sizeof(uint32_t)+c]) ? 1 : -1;

                            
                        }
                        
                        vtCovered++;
                        vtNo++;
                        
                        break;
                    case SV:
                        
                        LoadUInt(ins_sv_arr+vt[vtNo].alt, delLen);
                        LoadUInt(ins_sv_arr+vt[vtNo].alt+sizeof(uint32_t), insLen);
                       
                        
                        //introduce variant
                        if(insLen == 0)
                        {
                            refPos = refPos + delLen;
                            vtCovered++;
                            vtNo++;
                        }
                        else
                        {
                            for (int32_t c = 0; c < insLen && pos < lenToCheck; c++)
                            {
                                if(read[pos] == (ins_sv_arr[vt[vtNo].alt+(sizeof(uint32_t)<<1)+c]))
                                {
                                    pos++;
                                }
                                else
                                    return read[pos] > (ins_sv_arr[vt[vtNo].alt+(sizeof(uint32_t)<<1)+c]) ? 1 : -1;
                            }
                            refPos = refPos + delLen;
                            vtCovered++;
                            vtNo++;
                        }
                        break;
                }
                
            }
            else //variant not found/not introduced
            {
                vtCovered++;
                vtNo++;
            }
        }
    }
    
    
    
    return 0;
}



void  searchSA1Exact(unsigned char * read1, int32_t len, uint32_t encodedPattern, uint32_t delta, uint64_t  ** passed_result, uint64_t * pos, uint64_t *max, int32_t lenToCheck, uint32_t rightBorder)
{

    uint64_t l, r, m;
    bool firstIndAtPos;
    
    uint64_t * result = *passed_result;

    uint64_t * tempResult;

    
    
    
    uint32_t delLen;
        
    l = LUT1[encodedPattern];
    r = LUT1[rightBorder];
    
   
    uint64_t ll, rr;
    
    while (l < r) {
        m =(r + l) / 2UL;
        
        if (compareReadSuffix(read1+delta, lenToCheck, suffixArray1[m]) > 0)
            l = m + 1;
        else
            r = m;
    }
    ll = l;
    
    r = LUT1[rightBorder];
    while (l < r) {
        m = (r + l) / 2UL;
        if (compareReadSuffix(read1+delta, lenToCheck, suffixArray1[m]) < 0)
            r = m;
        else
            l = m + 1;
    }
    rr = l;

    
    if(rr>ll)
    {
        

        uint32_t  pos_vtToCheck;
        
        uint64_t  j;
        uint32_t vtCovered,  firstVar, chCovered, nextRefPos, remainedToInsert;
       
        int32_t lastSV = 0, firstSV = 0;
        uint32_t  prevDelCount = 0;
       
        
        for(j = ll; j < rr; j++)
        {
            remainedToInsert=0;
            
            firstIndAtPos = true;
            
            vtCovered = getCoveredVariants(suffixArray1[j].prevPos, suffixArray1[j].vtList , &firstVar, 0, &chCovered, &nextRefPos, &pos_vtToCheck, lenToCheck, &remainedToInsert);
            if(!vtCovered)
                continue;
            
            
            

            /** Check if not in deleted region - find meaningful SVs and DELs **/
            bool firstVarNotZero = firstVar > 0 ? true : false;
            uint32_t back;
            
            int32_t cSV = prevSV[firstVar];//noSV - 1;
            while(cSV >= 0 && vt[listSV[cSV]].prePos >= suffixArray1[j].prevPos)
                cSV--;
            lastSV = cSV;
            if(cSV < 0)
            {
                firstSV = 0;
            }
            else
            {
                while (cSV > 0 && suffixArray1[j].prevPos - vt[listSV[cSV]].prePos <= maxSVdel)
                    cSV--;
                
                /*
                LoadUInt(ins_sv_arr+vt[listSV[cSV]].alt, delLen);
                while (cSV < noSV && vt[listSV[cSV]].prePos +  delLen <= suffixArray1[j].prevPos)
                {
                    cSV++;
                    LoadUInt(ins_sv_arr+vt[listSV[cSV]].alt, delLen);
                }*/
                
                if(cSV < noSV)
                {
                    LoadUInt(ins_sv_arr+vt[listSV[cSV]].alt, delLen);
                    while (vt[listSV[cSV]].prePos +  delLen  <= suffixArray1[j].prevPos)
                    {
                        cSV++;
                        if(cSV < noSV)
                            LoadUInt(ins_sv_arr+vt[listSV[cSV]].alt, delLen);
                        else
                            break;
                    }
                }
                
                
                
                firstSV = cSV;
            }
            prevDelCount = 0;
            //check how many meaningful DELs before
            back = 1;
            if(firstVarNotZero)
            {
                while (firstVar >= back && suffixArray1[j].prevPos - vt[firstVar - back].prePos <= maxDel)
                {
                    if(vtType[firstVar - back] == DEL)
                    {
                        prevDelCount++;
                    }
                    back++;
                }
                
                back--;
                while(back)
                {
                    if(vtType[firstVar - back] == DEL)
                    {
                        if(vt[firstVar - back].prePos + vt[firstVar - back].alt <= suffixArray1[j].prevPos)
                            prevDelCount--;
                        else
                            break;
                    }
                    back--;
                }
            }
                
            
            vtList[0]  = 0;
            for(int32_t v=0; v < vtCovered; ++v)
                vtList[0] = vtList[0] | (1 << (NO_BITS_UINT - v - 1));
            vtList[0] = vtList[0] & suffixArray1[j].vtList;
            
            for (int32_t i = 1; i < vtListSize; i++)
                vtList[i]  = 0;
            
            uint32_t lastVt = remainedToInsert > 0 ? firstVar+vtCovered-1 : firstVar+vtCovered;
            
            int32_t preVDelNo = prevDeletion[firstVar];
            while (delFromRight[preVDelNo].nextRefPos > suffixArray1[j].prevPos + 1) {
                --preVDelNo;
            }
            
            *pos = checkVTconfigurationExact(suffixArray1[j].prevPos,  vtCovered, firstVar, firstVar+vtCovered, lastVt, firstVar, nextRefPos, chCovered, read1, len, delta, delta, &result, *pos, max, prevDelCount, firstSV, lastSV, pos_vtToCheck,  nextRefPos, 0, remainedToInsert,  preVDelNo);
            if (*pos  == *max)
            {
                *max = *max * 2;
                

                tempResult = ( uint64_t *) realloc (result, *max * sizeof( uint64_t));

                if (NULL == tempResult)
                {
                    printf("SA: Error in realloc, the message is - %s\n max = %llu\n", strerror(errno), *max);
                    fflush(stdout);
                    //free(result);
                    exit(8);
                }
                else
                {
                    result = tempResult;
                }
            }
            
        }
    }
    
    
    result[*pos] = 0LL;

    *passed_result=result;
    
}


void searchSA1(unsigned char * read1, int32_t len, uint32_t encodedPattern, uint32_t delta, uint64_t ** passed_result, uint64_t * pos, uint64_t *max, unsigned char * wholeRead, uint32_t wholeReadLen,  int32_t leftPos, uint32_t rightPos, uint32_t maxError, int32_t lenToCheck, uint32_t rightBorder)
{
  
    uint64_t l, r, m;
    bool firstIndAtPos;
    

    uint64_t * result = *passed_result;


    uint64_t * tempResult;

    
   
    
    uint32_t delLen;
   
    
    l = LUT1[encodedPattern];
    r = LUT1[rightBorder];
    
  
    /** test **/
    uint64_t ll, rr;
    
    while (l < r) {
        m =(r + l) / 2UL;
        
        if (compareReadSuffix(read1+delta, lenToCheck, suffixArray1[m]) > 0)
            l = m + 1;
        else
            r = m;
    }
    ll = l;
    
    r = LUT1[rightBorder];
    
    while (l < r) {
        m = (r + l) / 2UL;
        if (compareReadSuffix(read1+delta, lenToCheck, suffixArray1[m]) < 0)
            r = m;
        else
            l = m + 1;
    }
    rr = l;
    /***/
    
    
    if(rr>ll)
    {
       
        uint32_t  pos_vtToCheck;
        
        uint64_t  j;
        uint32_t vtCovered,  firstVar, chCovered, nextRefPos, remainedToInsert;
      
        int32_t lastSV = 0, firstSV = 0;
        uint32_t  prevDelCount = 0, auxiliary;
        uint32_t * nn;//p = (uint32_t *)malloc(1*sizeof(uint32_t));
        
        for(j = ll; j < rr; j++)
        {
            remainedToInsert=0;
            
            firstIndAtPos = true;
            
            vtCovered = getCoveredVariants(suffixArray1[j].prevPos, suffixArray1[j].vtList , &firstVar, 0, &chCovered, &nextRefPos, &pos_vtToCheck,  lenToCheck, &remainedToInsert);
            if(!vtCovered)
                continue;

            /** Check if not in deleted region - find meaningful SVs and DELs **/
            bool firstVarNotZero = firstVar > 0 ? true : false;
            uint32_t back;
            
            int32_t cSV = prevSV[firstVar];//noSV - 1;
            while(cSV >= 0 && vt[listSV[cSV]].prePos >= suffixArray1[j].prevPos)
                cSV--;
            lastSV = cSV;
            if(cSV < 0)
            {
                firstSV = 0;
            }
            else
            {
                while (cSV > 0 && suffixArray1[j].prevPos - vt[listSV[cSV]].prePos <= maxSVdel)
                    cSV--;
            
 /*               LoadUInt(ins_sv_arr+vt[listSV[cSV]].alt, delLen);
                while (cSV < noSV && vt[listSV[cSV]].prePos +  delLen <= suffixArray1[j].prevPos)
                {
                    cSV++;
                    LoadUInt(ins_sv_arr+vt[listSV[cSV]].alt, delLen);
                }
   */
                
                if(cSV < noSV)
                {
                    LoadUInt(ins_sv_arr+vt[listSV[cSV]].alt, delLen);
                    while (vt[listSV[cSV]].prePos +  delLen  <= suffixArray1[j].prevPos)
                    {
                        cSV++;
                        if(cSV < noSV)
                            LoadUInt(ins_sv_arr+vt[listSV[cSV]].alt, delLen);
                        else
                            break;
                    }
                }
                
                
                firstSV = cSV;
            }
            prevDelCount = 0;
            //check how many meaningful DELs before
            back = 1;
            if(firstVarNotZero)
            {
                while (firstVar >= back && suffixArray1[j].prevPos - vt[firstVar - back].prePos <= maxDel)
                {
                    if(vtType[firstVar - back] == DEL)
                    {
                        prevDelCount++;
                    }
                    back++;
                }
                
                back--;
                while(back)
                {
                    if(vtType[firstVar - back] == DEL)
                    {
                        if(vt[firstVar - back].prePos + vt[firstVar - back].alt <= suffixArray1[j].prevPos)
                            prevDelCount--;
                        else
                            break;
                    }
                    back--;
                }
            }
            
            
            vtList[0]  = 0;
            for(int32_t v=0; v < vtCovered; ++v)
                vtList[0] = vtList[0] | (1 << (NO_BITS_UINT - v - 1));
            vtList[0] = vtList[0] & suffixArray1[j].vtList;
            
            for (int32_t i = 1; i < vtListSize; i++)
                vtList[i]  = 0;
            
            uint32_t lastVt = remainedToInsert > 0 ? firstVar+vtCovered-1 : firstVar+vtCovered;
            
            int32_t preVDelNo = prevDeletion[firstVar];
            while (delFromRight[preVDelNo].nextRefPos > suffixArray1[j].prevPos + 1) {
                --preVDelNo;
            }
            
            nn = &auxiliary;
            *pos = checkVTconfiguration1(suffixArray1[j].prevPos, vtCovered, firstVar, firstVar+vtCovered, lastVt, firstVar, nextRefPos, chCovered, read1, len, delta, delta, &result, *pos, max, prevDelCount, firstSV, lastSV, pos_vtToCheck,  wholeRead, wholeReadLen, leftPos, rightPos, 0, maxError, nextRefPos, 0, remainedToInsert, leftPos+1, suffixArray1[j].prevPos, preVDelNo);
            if (*pos  == *max)
            {
                *max = *max * 2;
                

                tempResult = ( uint64_t *) realloc (result, *max * sizeof( uint64_t));

                if (NULL == tempResult)
                {
                    printf("SA: Error in realloc, the message is - %s\n max = %llu\n", strerror(errno), *max);
                    fflush(stdout);
                    //free(result);
                    exit(8);
                }
                else
                {
                    result = tempResult;
                }
            }
            
            
            
            
            
            
        }
    }
    
    
    result[*pos] = 0LL;

    *passed_result=result;
    
    
}



