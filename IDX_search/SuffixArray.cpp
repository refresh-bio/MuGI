/*
 This file is a part of MuGI software distributed under GNU GPL 2 licence.
 The homepage of the MuGI project is http://sun.aei.polsl.pl/mugi/
 
 Authors: Agnieszka Danek, Sebastian Deorowicz
 
 Version: 1.0
 Date   : 2015-Mar
 */

#include "SuffixArray.h"

void searchStringExact(unsigned char * read1, int32_t len, uint32_t encodedPattern,  int32_t delta, uint64_t ** passed_result, uint64_t * pos, uint64_t *max, int32_t lenToCheck, uint32_t rightBorder)
{
    
    
    uint64_t * result = *passed_result;
    
    
    uint64_t l, r, m, s;
    uint64_t * tempResult;
    

    
    

    l = LUT[encodedPattern];
    r = LUT[rightBorder];
    int32_t p, idx;
    int32_t greater;
    uint32_t c, delLen;
    
    
    while (l < r)
    {
        m = ((l+r) / 2UL);
        
        idx = 0; greater = 0;
        for(p = delta; idx < lenToCheck; p++, idx++)
        {
            c = getRefFinal(suffixArray[m]+idx);
            if(read1[p] == c)
                ;
            else if(read1[p] > c  || c == 0xF)
            {
                greater = 1;
                l = m + 1;
                break;
            }
            else
            {
                break;
            }
        }
        if(!greater)
            r = m;
    }
    s = l;
    
    r = LUT[rightBorder];
    
    while (l < r) {
        m = (r + l) / 2UL;
        idx = 0; greater = 1;
        for(p = delta; idx < lenToCheck; p++, idx++)
        {
            c = getRefFinal(suffixArray[m]+idx);
            if(read1[p] == c)
                ;
            else if(read1[p] < c || c == 0xF)
            {
                greater = 0;
                r = m;
                break;
            }
            else
            {
                break;
            }
        }
        if(greater)
            l = m + 1;
        
    }
    r = l;
    

    if(r > s)
    {
        uint32_t coveredVariants, firstVar, chCov, currVt;
        
        
        
       
        uint32_t nextRefPos;
        
        
        
        uint32_t  pos_vtToCheck;
        
        
        
        
        
        int32_t lastSV = 0, firstSV = 0;
        uint32_t  prevDelCount = 0;
        
        for(uint64_t  j = s; j < r; j++)
        {
            
            for(int u = 0; u < vtListSize; u++)
                vtList[u]=0;

            
            
            pos_vtToCheck = 0;
            coveredVariants = getCoveredVariantsSA(suffixArray[j]-1,  &firstVar,  &chCov, &nextRefPos ,&currVt,  &pos_vtToCheck,  lenToCheck);
            
            
            
            

            /** Check if not in deleted region - find meaningful SVs and DELs **/
            bool firstVarNotZero = firstVar > 0 ? true : false;
            uint32_t back;
            
            int32_t cSV = prevSV[firstVar];
            while(cSV >= 0 && vt[listSV[cSV]].prePos >= suffixArray[j]-1)
                cSV--;
            lastSV = cSV;
            if(cSV < 0)
            {
                firstSV = 0;
            }
            else
            {
                while (cSV > 0 && suffixArray[j]-1- vt[listSV[cSV]].prePos <= maxSVdel)
                    cSV--;
                
                
                /*LoadUInt(ins_sv_arr+vt[listSV[cSV]].alt, delLen);
                
                while (cSV < noSV && vt[listSV[cSV]].prePos +  delLen < suffixArray[j])
                {
                    cSV++;
                    LoadUInt(ins_sv_arr+vt[listSV[cSV]].alt, delLen);
                }*/
                
                if(cSV < noSV)
                {
                    LoadUInt(ins_sv_arr+vt[listSV[cSV]].alt, delLen);
                    while (vt[listSV[cSV]].prePos +  delLen  < suffixArray[j])
                    {
                        cSV++;
                        if(cSV < noSV)
                            LoadUInt(ins_sv_arr+vt[listSV[cSV]].alt, delLen);
                        else
                            break;
                    }
                }
                
                
                firstSV = cSV < noSV? cSV : noSV+1;
            }
            
            prevDelCount = 0;
            
            //check how many meaningful DELs before
            back = 1;
            if(firstVarNotZero)
            {
                while (firstVar >= back && suffixArray[j]-1 - vt[firstVar - back].prePos <= maxDel)
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
                        if(vt[firstVar - back].prePos + vt[firstVar - back].alt < suffixArray[j])
                            prevDelCount--;
                        else
                            break;
                    }
                    back--;
                }
            }
            /***/

            
            int32_t preVDelNo = prevDeletion[firstVar];
            while (delFromRight[preVDelNo].nextRefPos > suffixArray[j]-1 + 1) {
                --preVDelNo;
            }
            
         
            *pos = checkVTconfigurationExact(suffixArray[j]-1,   coveredVariants, firstVar, currVt, currVt,  firstVar, nextRefPos, chCov, read1, len, delta, delta, &result, *pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck,  nextRefPos, 0, 0, preVDelNo);
            
            
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

void  searchString(unsigned char * read1, int32_t len, uint32_t encodedPattern,  int32_t delta, uint64_t ** passed_result, uint64_t * pos, uint64_t *max, unsigned char * wholeRead, uint32_t wholeReadLen,  int32_t leftPos, uint32_t rightPos, uint32_t maxError, int32_t lenToCheck, uint32_t rightBorder)
{
    
    
    uint64_t * result = *passed_result;
    
    
    uint64_t l, r, m, s;
    uint64_t * tempResult;

    
  
    l = LUT[encodedPattern];
    r = LUT[rightBorder];
    int32_t p, idx;
    int32_t greater;
    uint32_t c, delLen;
    
    
    while (l < r)
    {
        m = ((l+r) / 2UL);
        
        idx = 0; greater = 0;
        for(p = delta; idx < lenToCheck; p++, idx++)
        {
            c = getRefFinal(suffixArray[m]+idx);
            if(read1[p] == c)
                ;
            else if(read1[p] > c || c == 0xF)
            {
                greater = 1;
                l = m + 1;
                break;
            }
            else
            {
                break;
            }
        }
        if(!greater)
            r = m;
    }
    s = l;
    
    r = LUT[rightBorder];
    
    while (l < r) {
        m = (r + l) / 2UL;
        idx = 0; greater = 1;
        for(p = delta; idx < lenToCheck; p++, idx++)
        {
            c = getRefFinal(suffixArray[m]+idx);
            if(read1[p] == c)
                ;
            else if(read1[p] < c || c == 0xF)
            {
                greater = 0;
                r = m;
                break;
            }
            else
            {
                break;
            }
        }
        if(greater)
            l = m + 1;
        
    }
    r = l;
    
    if(r > s)
    {
        uint32_t coveredVariants, firstVar, chCov, currVt;
        
        
        
        uint32_t nextRefPos;
        
        uint32_t  pos_vtToCheck;
        
        
        
        
        
        int32_t lastSV = 0, firstSV = 0;
        uint32_t  prevDelCount = 0;
        
        for(uint64_t  j = s; j < r ; j++)
        {
            for(int u = 0; u < vtListSize; u++)
                vtList[u]=0;
            
            
            pos_vtToCheck = 0;
            coveredVariants = getCoveredVariantsSA(suffixArray[j]-1,  &firstVar,  &chCov, &nextRefPos ,&currVt,  &pos_vtToCheck, lenToCheck);
            
            
            

            /** Check if not in deleted region - find meaningful SVs and DELs **/
            bool firstVarNotZero = firstVar > 0 ? true : false;
            uint32_t back;
            
           
            int32_t cSV = prevSV[firstVar];
            while(cSV >= 0 && vt[listSV[cSV]].prePos >= suffixArray[j]-1)
                cSV--;
            lastSV = cSV;
            if(cSV < 0)
            {
                firstSV = 0;
            }
            else
            {
                while (cSV > 0 && suffixArray[j]-1- vt[listSV[cSV]].prePos <= maxSVdel)
                    cSV--;
                
                
               /* LoadUInt(ins_sv_arr+vt[listSV[cSV]].alt, delLen);
                
                while (cSV < noSV && vt[listSV[cSV]].prePos +  delLen < suffixArray[j])
                {
                    cSV++;
                    LoadUInt(ins_sv_arr+vt[listSV[cSV]].alt, delLen);
                }
                */
                
                if(cSV < noSV)
                {
                    LoadUInt(ins_sv_arr+vt[listSV[cSV]].alt, delLen);
                    while (vt[listSV[cSV]].prePos +  delLen  < suffixArray[j])
                    {
                        cSV++;
                        if(cSV < noSV)
                            LoadUInt(ins_sv_arr+vt[listSV[cSV]].alt, delLen);
                        else
                            break;
                    }
                }
                
                firstSV = cSV < noSV? cSV : noSV+1;
            }
            
            prevDelCount = 0;
            
            //check how many meaningful DELs before
            back = 1;
            if(firstVarNotZero)
            {
                while (firstVar >= back && suffixArray[j]-1 - vt[firstVar - back].prePos <= maxDel)
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
                        if(vt[firstVar - back].prePos + vt[firstVar - back].alt < suffixArray[j])
                            prevDelCount--;
                        else
                            break;
                    }
                    back--;
                }
            }

            
            int32_t preVDelNo = prevDeletion[firstVar];
            while (delFromRight[preVDelNo].nextRefPos > suffixArray[j]-1 + 1) {
                --preVDelNo;
            }
            
            *pos = checkVTconfiguration1(suffixArray[j]-1,   coveredVariants, firstVar, currVt, currVt,  firstVar, nextRefPos, chCov, read1, len, delta, delta, &result, *pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck, wholeRead, wholeReadLen, leftPos, rightPos, 0, maxError, nextRefPos, 0, 0, leftPos+1, suffixArray[j]-1,  preVDelNo);
            
            
            
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




uint32_t getCoveredVariantsSA(uint32_t  prevPos,  uint32_t * firstVar,  uint32_t * chCov, uint32_t * nextRefPos, uint32_t * currVT,  uint32_t * pos_vtToCheck,  uint32_t  lenToCheck)
{
    
    uint32_t vtCovered = 0, temp;
    uint32_t *temp_vtToCheck;
    uint32_t refPos = prevPos+1;
    ///BINARY SEARCH FOR FIRST VT
    uint32_t vtNo = 0,  m, l, r;
    
    
    temp = prevPos/VT_LUT_STEP;
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
    
    while (l < r)
    {
        m = ((l+r) / 2);
        if ((vt[m].prePos < prevPos || vt[m].prePos == 0xFFFFFFFF) && prevPos != 0xFFFFFFFF)
        {
            l = m + 1;
        }
        else
        {
            r = m;
        }
    }
    vtNo = l;
    *firstVar = vtNo;
    
    
    uint32_t nextPos =  lenToCheck < kl ? refPos +  lenToCheck :  refPos + kl;
    
    while ( vt[vtNo].prePos+1 < nextPos &&  vt[vtNo].prePos+1 < refSize)
    {
        (vtToCheck)[(*pos_vtToCheck)++] = vtNo;
        if(*pos_vtToCheck == max_vtToCheck)
        {
            max_vtToCheck = (max_vtToCheck) * 2;
            
            temp_vtToCheck = ( uint32_t *) realloc (vtToCheck, (max_vtToCheck) * sizeof( uint32_t));
            if (NULL == temp_vtToCheck)
            {
                printf("SA: Error in realloc, the message is - %s\n max = %d\n", strerror(errno), max_vtToCheck);
                fflush(stdout);
                //free(result);
                exit(8);
            }
            else
            {
                vtToCheck = temp_vtToCheck;
            }
            
        }
        vtCovered++;
        vtNo++;
    }
    
    if( nextPos <= refSize )
    {
        *nextRefPos = nextPos;
        *chCov = lenToCheck;
    }
    else
    {
        *nextRefPos = refSize;
        *chCov = refSize - refPos;
    }
    
   
    
    
    *currVT = vtNo;
    
    return vtCovered;
}