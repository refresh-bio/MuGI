/*
 This file is a part of MuGI software distributed under GNU GPL 2 licence.
 The homepage of the MuGI project is http://sun.aei.polsl.pl/mugi/
 
 Authors: Agnieszka Danek, Sebastian Deorowicz
 
 Version: 1.0
 Date   : 2015-Mar
 */

#include "SuffixArray2.h"





int32_t compareReadSuffix(unsigned char* read, uint32_t lenToCheck, suffix2 suffix)
{
    
    int32_t pos = 0;
    uint32_t insLen, delLen;
   
    bool firstVT = true;
   
    uint32_t refPos = suffix.prevPos+1;
    
    ///BINARY SEARCH FOR FIRST VT
    uint32_t vtNo = 0, l, r, m;
 

    
    if(refSize >= VT_LUT_STEP)
    {
        l = VT_LUT[suffix.prevPos/VT_LUT_STEP];
        r = VT_LUT[suffix.prevPos/VT_LUT_STEP+1];
    }
    else
    {
        l = 0;
        r = noVar;
    }
    
    
    while (l < r) {
        m = (r + l) / 2;
        if ((vt[m].prePos <  suffix.prevPos || vt[m].prePos == 0xFFFFFFFF) && suffix.prevPos != 0xFFFFFFFF)
            l = m + 1;
        else
            r = m;
    }
      vtNo = l;
    
    
    uint32_t vtCovered = 0;

    
    
    while (pos < lenToCheck)
    {
        while( vt[vtNo].prePos+1 < refPos)
        {
            cout << vtNo << " "<< vt[vtNo].prePos <<" "<< refPos<< " l=" << l << " r=" << r << " m = " << m << " " << suffix.prevPos << " "<< suffix.gap << " " << suffix.vtList <<  endl;
            cout << "-shouldn't end up here with suffix2 (no var in del)\n";
            vtNo++;
            
        }
        
        if( vt[vtNo].prePos+1 > refPos)  //no variant at currPos //vtNo >= noVar ||
        {
            
            if(read[pos] == getRefFinal(refPos))
            {
                pos++;
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
                        
                        if(firstVT)
                        {
                            for (int32_t c = suffix.gap; c < insLen && pos < lenToCheck; c++)
                            {
                                if(read[pos] == (ins_sv_arr[vt[vtNo].alt+sizeof(uint32_t)+c]))
                                {
                                    pos++;
                                }
                                else
                                    return read[pos] > (ins_sv_arr[vt[vtNo].alt+sizeof(uint32_t)+c]) ? 1 : -1;

                               
                            }
                        }
                        else
                        {
                            for (int32_t c = 0; c < insLen && pos < lenToCheck; c++)
                            {
                                
                                if(read[pos] == (ins_sv_arr[vt[vtNo].alt+sizeof(uint32_t)+c]))
                                {
                                    pos++;
                                }
                                else
                                    return read[pos] > (ins_sv_arr[vt[vtNo].alt+sizeof(uint32_t)+c]) ? 1 : -1;
                                
                            }
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
                            if(firstVT)
                            {;
                                for (int32_t c = suffix.gap; c < insLen && pos < lenToCheck; c++)
                                {
                                    if(read[pos] == (ins_sv_arr[vt[vtNo].alt+(sizeof(uint32_t)<<1)+c]))
                                    {
                                        pos++;
                                    }
                                    else
                                        return read[pos] > (ins_sv_arr[vt[vtNo].alt+(sizeof(uint32_t)<<1)+c]) ? 1 : -1;
                                }
                                
                             
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
            firstVT = false;
        }
    }
    
    
    return 0;
}





void searchSA2Exact(unsigned char * read1, int32_t len, uint32_t encodedPattern, uint32_t delta, uint64_t ** passed_result, uint64_t * pos, uint64_t *max, int32_t lenToCheck, uint32_t rightBorder)
{
    
    
    
    uint64_t l, r, m, s;

    uint64_t * result = *passed_result;

    uint64_t * tempResult;

    
    
    l = LUT2[encodedPattern];
    r = LUT2[rightBorder];
    
    while (l < r)
    {
        m = ((l+r) / 2UL);
        if (compareReadSuffix(read1+delta, lenToCheck, suffixArray2[m]) > 0)
        {
            l = m + 1;
        }
        else
        {
            r = m;
        }
    }
    s = l;
    
    r = LUT2[rightBorder];
    while (l < r) {
        m = (r + l) / 2UL;
        if (compareReadSuffix(read1+delta, lenToCheck, suffixArray2[m]) < 0)
        {
            r = m;
        }
        else
        {
            l = m + 1;
        }
    }
    r = l;
    
    
    if (r>s)
    {

        
        uint64_t  j;
        uint32_t vtCovered, firstVar, chCovered, nextRefPos, remainedToInsert;
        
        uint32_t pos_vtToCheck;
      
        int32_t lastSV = 0, firstSV = 0;
        uint32_t  prevDelCount = 0;
        
        uint32_t  delLen;
        
        for(j = s; j < r; j++)
        {
            remainedToInsert = 0;
            vtCovered = getCoveredVariants(suffixArray2[j].prevPos, suffixArray2[j].vtList , &firstVar, suffixArray2[j].gap, &chCovered, &nextRefPos, &pos_vtToCheck, lenToCheck, &remainedToInsert);
            if(!vtCovered)
                continue;
            
            /** Check if not in deleted region - find meaningful SVs and DELs **/
            bool firstVarNotZero = firstVar > 0 ? true : false;
            uint32_t back;
            
            int32_t cSV = prevSV[firstVar];//noSV - 1;
            while(cSV >= 0 && vt[listSV[cSV]].prePos >= suffixArray2[j].prevPos)
                cSV--;
            lastSV = cSV;
            if(cSV < 0)
            {
                firstSV = 0;
            }
            else
            {
                while (cSV > 0 && suffixArray2[j].prevPos - vt[listSV[cSV]].prePos <= maxSVdel)
                    cSV--;
                
                
                if(cSV < noSV)
                {
                    LoadUInt(ins_sv_arr+vt[listSV[cSV]].alt, delLen);
                    while (vt[listSV[cSV]].prePos +  delLen  <= suffixArray2[j].prevPos)
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
                while (firstVar >= back && suffixArray2[j].prevPos - vt[firstVar - back].prePos <= maxDel)
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
                        if(vt[firstVar - back].prePos + vt[firstVar - back].alt <= suffixArray2[j].prevPos)
                            prevDelCount--;
                        else
                            break;
                    }
                    back--;
                }
            }
            /************************************/
            
            
            vtList[0]  = 0;
            for(int32_t v=0; v < vtCovered; ++v)
                vtList[0] = vtList[0] | (1 << (NO_BITS_UINT - v - 1));
            vtList[0] = vtList[0] & suffixArray2[j].vtList;
            
            for (int32_t i = 1; i < vtListSize; i++)
                vtList[i]  = 0;
            
            // int32_t c = 1;
            int32_t d = delta;
            int32_t gap = suffixArray2[j].gap;
            
            
                
            
            uint32_t lastVt = remainedToInsert > 0 ? firstVar+vtCovered-1 : firstVar+vtCovered;
            
            int32_t preVDelNo = prevDeletion[firstVar];
            while (delFromRight[preVDelNo].nextRefPos > suffixArray2[j].prevPos + 1) {
                --preVDelNo;
            }
            
            *pos = checkVTconfigurationExact(suffixArray2[j].prevPos, vtCovered, firstVar, firstVar+vtCovered, lastVt,firstVar,  nextRefPos, chCovered, read1, len, delta, d, &result, *pos, max, prevDelCount, firstSV, lastSV, pos_vtToCheck, nextRefPos, gap, remainedToInsert,  preVDelNo);
            
            if (*pos == *max)
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
            
            //check if current is a SV with insertion in some other del region
            if(vtType[firstVar] == SV) // so there is ins (gap > 0 for SA2) that may be present in other deletion region
            {
                
                
                uint32_t dels, delLenFirst;
                LoadUInt(ins_sv_arr+vt[firstVar].alt, delLenFirst);
                dels = prevDelCount;
                back = 1;
                while(dels)
                {
                    if(vtType[firstVar - back] == DEL)
                    {
                        dels--;
                        //if it is indeed in del region of other DEL
                        if(vt[firstVar - back].prePos + vt[firstVar - back].alt > vt[firstVar].prePos \
                           && vt[firstVar - back].prePos + vt[firstVar - back].alt < vt[firstVar].prePos + delLenFirst)
                        {
                            
                            *pos = checkDelRegion(firstVar - back, firstVar, -1, gap, 0, 0, 0, pos_vtToCheck, vtCovered, firstSV, lastSV, prevDelCount, &result, *pos, max, lastVt, nextRefPos, remainedToInsert, chCovered, read1, len, delta, d, 0, 0);
                            
                            
                            
                            
                        }
                    }
                    back++;
                }
                
                
                if(firstSV < noSV)
                    for(int32_t sv = lastSV; sv >= (int32_t)firstSV; sv--)
                    {
                        LoadUInt(ins_sv_arr+vt[listSV[sv]].alt, delLen);
                        
                        //if it is indeed in del region of other DEL
                        if(vt[listSV[sv]].prePos +  delLen > vt[firstVar].prePos \
                           && vt[listSV[sv]].prePos +  delLen < vt[firstVar].prePos + delLenFirst)
                        {
                            *pos = checkDelRegion(listSV[sv], firstVar, -1, gap, 0, 0, 0, pos_vtToCheck, vtCovered, firstSV, lastSV, prevDelCount, &result, *pos, max, lastVt, nextRefPos, remainedToInsert, chCovered, read1, len, delta, d, 0, 0);
                        }
                    }
                }
            
        }
    }
    
    result[*pos] = 0LL;

    *passed_result=result;
    
    
}






void  searchSA2(unsigned char * read1, int32_t len, uint32_t encodedPattern, uint32_t delta, uint64_t ** passed_result, uint64_t * pos, uint64_t *max, unsigned char * wholeRead, uint32_t wholeReadLen,  int32_t leftPos, uint32_t rightPos, uint32_t maxError, int32_t lenToCheck, uint32_t rightBorder)
{
    
    
    
    uint64_t l, r, m, s;
    

    uint64_t * result = *passed_result;

    uint64_t * tempResult;

    
    
    
    l = LUT2[encodedPattern];
    r = LUT2[rightBorder];
    
    while (l < r)
    {
        m = ((l+r) / 2UL);
        if (compareReadSuffix(read1+delta, lenToCheck, suffixArray2[m]) > 0)
        {
            l = m + 1;
        }
        else
        {
            r = m;
        }
    }
    s = l;
    
    r = LUT2[rightBorder];
    while (l < r) {
        m = (r + l) / 2UL;
        if (compareReadSuffix(read1+delta, lenToCheck, suffixArray2[m]) < 0)
        
        {
            r = m;
        }
        else
        {
            l = m + 1;
        }
    }
    r = l;
    
    
    if (r>s)
    {
    
        uint64_t  j;
        uint32_t vtCovered, firstVar, chCovered, nextRefPos, remainedToInsert;
        
        uint32_t pos_vtToCheck;
     
        int32_t lastSV = 0, firstSV = 0;
        uint32_t  prevDelCount = 0;
        
        uint32_t  delLen;
        
        for(j = s; j < r; j++)
        {
            remainedToInsert = 0;
            vtCovered = getCoveredVariants(suffixArray2[j].prevPos, suffixArray2[j].vtList , &firstVar, suffixArray2[j].gap, &chCovered, &nextRefPos, &pos_vtToCheck, lenToCheck, &remainedToInsert);
            if(!vtCovered)
                continue;
            
            uint32_t lastVt = remainedToInsert > 0 ? firstVar+vtCovered-1 : firstVar+vtCovered;
            
            
            vtList[0]  = 0;
            for(int32_t v=0; v < vtCovered; ++v)
                vtList[0] = vtList[0] | (1 << (NO_BITS_UINT - v - 1));
            vtList[0] = vtList[0] & suffixArray2[j].vtList;
            
            for (int32_t i = 1; i < vtListSize; i++)
                vtList[i]  = 0;
            
            // int32_t c = 1;
            int32_t d = delta;

            
            
            
            /** Check if not in deleted region - find meaningful SVs and DELs **/
            bool firstVarNotZero = firstVar > 0 ? true : false;
            uint32_t back;
            
            int32_t cSV = prevSV[firstVar];//noSV - 1;
            while(cSV >= 0 && vt[listSV[cSV]].prePos >= suffixArray2[j].prevPos)
                cSV--;
            lastSV = cSV;
            if(cSV < 0)
            {
                firstSV = 0;
            }
            else
            {
                while (cSV > 0 && suffixArray2[j].prevPos - vt[listSV[cSV]].prePos <= maxSVdel)
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
                    while (vt[listSV[cSV]].prePos +  delLen  <= suffixArray2[j].prevPos)
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
                while (firstVar >= back && suffixArray2[j].prevPos - vt[firstVar - back].prePos <= maxDel)
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
                        if(vt[firstVar - back].prePos + vt[firstVar - back].alt <= suffixArray2[j].prevPos)
                            prevDelCount--;
                        else
                            break;
                    }
                    back--;
                }
            }
/************************************/
            
            
        
            

            
            
            
            
            
            
            
            
            /*//always delta==0 for SA2
             int32_t gap = suffixArray2[j].gap - 1;
             if(suffixArray2[j].gap > 0 && delta > 0)
            {
                if(vtType[firstVar] == INS)
                    tempLen = vt[firstVar].alt+sizeof(uint32_t);
                else //SV
                    tempLen = vt[firstVar].alt+(sizeof(uint32_t)<<1);
                
                
                d = delta - 1;
                
                while (d >= 0 && gap >= 0)
                {
                    if ((ins_sv_arr[tempLen+gap]) == read1[d])
                    {
                        d--; gap--;
                    }
                    else
                    {
                        firstIndAtPos = false;
                        break;
                    }
                }
                d++;
            }
            gap++;*/
            
            
        
                
            // check variant configuration, do not let any previous deletion (even del regions)
            int32_t preVDelNo = prevDeletion[firstVar];
            while (delFromRight[preVDelNo].nextRefPos > suffixArray2[j].prevPos + 1) {
                --preVDelNo;
            }
            *pos = checkVTconfiguration1(suffixArray2[j].prevPos, vtCovered, firstVar, firstVar+vtCovered, lastVt,firstVar,  nextRefPos, chCovered, read1, len, delta, d, &result, *pos, max, prevDelCount, firstSV, lastSV,  pos_vtToCheck, wholeRead, wholeReadLen, leftPos, rightPos, 0, maxError, nextRefPos, suffixArray2[j].gap, remainedToInsert, leftPos+1, suffixArray2[j].prevPos, preVDelNo);
            
            if (*pos == *max)
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
            
            
            
            
             //check if current is a SV with insertion in some other del region
             if(vtType[firstVar] == SV) // so there is ins (gap > 0 for SA2) that may be present in other deletion region
             {
                 
                 
                 int32_t gap = suffixArray2[j].gap;
                 uint32_t dels, delLenFirst;
                 LoadUInt(ins_sv_arr+vt[firstVar].alt, delLenFirst);
                 dels = prevDelCount;
                 back = 1;
                 while(dels)
                 {
                     if(vtType[firstVar - back] == DEL)
                     {
                         dels--;
                          //if it is indeed in del region of other DEL
                         if(vt[firstVar - back].prePos + vt[firstVar - back].alt > vt[firstVar].prePos \
                         && vt[firstVar - back].prePos + vt[firstVar - back].alt < vt[firstVar].prePos + delLenFirst)
                         {
                             
                             *pos = checkDelRegion(firstVar - back, firstVar, leftPos, gap, 0, maxError, wholeRead, pos_vtToCheck, vtCovered, firstSV, lastSV, prevDelCount, &result, *pos, max, lastVt, nextRefPos, remainedToInsert, chCovered, read1, len, delta, d, wholeReadLen, rightPos);
                             
                             

                         
                         }
                     }
                     back++;
                 }
                

                 if(firstSV < noSV)
                 for(int32_t sv = lastSV; sv >= (int32_t)firstSV; sv--)
                 {
                     LoadUInt(ins_sv_arr+vt[listSV[sv]].alt, delLen);
                     
                      //if it is indeed in del region of other DEL
                     if(vt[listSV[sv]].prePos +  delLen > vt[firstVar].prePos \
                     && vt[listSV[sv]].prePos +  delLen < vt[firstVar].prePos + delLenFirst)
                     {
                        *pos = checkDelRegion(listSV[sv], firstVar, leftPos, gap, 0, maxError, wholeRead, pos_vtToCheck, vtCovered, firstSV, lastSV, prevDelCount, &result, *pos, max, lastVt, nextRefPos, remainedToInsert, chCovered, read1, len, delta, d, wholeReadLen, rightPos);
                     }
                 }
             }

        }
    }
    
    result[*pos] = 0LL;

    *passed_result=result;
    
}


