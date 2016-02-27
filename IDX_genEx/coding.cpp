/*
 This file is a part of MuGI software distributed under GNU GPL 2 licence.
 The homepage of the MuGI project is http://sun.aei.polsl.pl/mugi/
 
 Authors: Agnieszka Danek, Sebastian Deorowicz
 
 Version: 1.0
 Date   : 2015-Mar
 */

#include "coding.h"


void StoreUInt(unsigned char *dest, uint32_t  data)
{
    *dest++ = (data      ) & 0xFF;
    *dest++ = (data >>  8) & 0xFF;
    *dest++ = (data >> 16) & 0xFF;
    *dest++ = (data >> 24) & 0xFF;
    
}

void LoadUInt(unsigned char *dest, uint32_t  &data)
{
    
  // ver1
    data = 0;
    data = *dest++;
    data = data + ((*dest++) << 8);
    data = data + ((*dest++) << 16);
    data = data + ((*dest  ) << 24);
   
    
    //ver2
    //A_memcpy(&data, dest, sizeof(uint32_t));
    
    
    //ver3
    //memcpy(&data, dest, sizeof(uint32_t));

}




unsigned char getRef(uint32_t n)
{
    unsigned char c;
    //if(n%2)
    if(n&1)
        c = reference[n/2] & 0xF;
    else
        c = (reference[n/2] >> 4) & 0x0F ;
    
    
    switch(c)
    {
        case 0:
            c = 'A';
            break;
        case 1:
            c = 'C';
            break;
        case 2:
            c = 'G';
            break;
        case 8:
            c = 'T';
            break;
        case 5:
            c = 'N';
            break;
        case 3:
            c = 'K';
            break;
        case 4:
            c = 'M';
            break;
        case 6:
            c = 'R';
            break;
        case 7:
            c = 'S';
            break;
        case 9:
            c = 'U';
            break;
        case 10:
            c = 'W';
            break;
        case 11:
            c = 'Y';
            break;
        default:
          //  printf(" Unknown nucleotide: %d\n", c);
            c = 'N';
            break;
            
    }
    return c;
    
}

char decodeCode(unsigned char x)
{
    char c;
    switch(x)
    {
        case 0:
            c = 'A';
            break;
        case 1:
            c = 'C';
            break;
        case 2:
            c = 'G';
            break;
        case 8:
            c = 'T';
            break;
        case 5:
            c = 'N';
            break;
        case 3:
            c = 'K';
            break;
        case 4:
            c = 'M';
            break;
        case 6:
            c = 'R';
            break;
        case 7:
            c = 'S';
            break;
        case 9:
            c = 'U';
            break;
        case 10:
            c = 'W';
            break;
        case 11:
            c = 'Y';
            break;
        default:
            printf("decodeCode() Unknown nucleotide: %d\n", x);
            c = 0XF;
            break;
            
    }
    
    return c;
    
}

char *  decode(uint32_t prevPos, uint32_t length, int32_t id,  uint32_t vtStart, int32_t gap)
{
    uint32_t insLen, delLen, tempLen;
    
    char * result = (char *) malloc ((length+1) * sizeof( char));
    char * tempResult;
    int32_t pos = 0;
    int32_t max = length+1;
    bool firstVT = true;
    uint32_t chCovered = 0;
    uint32_t refPos = prevPos+1;
    int32_t back, prevIns;
    
    uint32_t vtNo;
    ///BINARY SEARCH FOR FIRST VT
    if(vtStart != 0xFFFFFFFF)
    {
        vtNo = vtStart;
    }
    else
    {
        uint32_t l = 0, r = noVar-1, m;
        vtNo = 0;
        while (l < r)
        {
            m = ((l+r) / 2);
            if (vt[m].prePos < prevPos)
            {
                l = m + 1;
            }
            else
            {
                r = m;
            }
        }
        vtNo = l;
        vtStart = vtNo;
    }
    
    bool vtPresent;
    uint32_t newRefPos;
    
    /** Check if not in deleted region - find meaningful SVs and DELs **/
    bool firstVarNotZero = vtStart > 0 ? true : false;
    //uint32_t back;
    
    int32_t cSV = prevSV[vtStart];//noSV - 1;
    int32_t lastSV, firstSV;
    while(cSV >= 0 && vt[listSV[cSV]].prePos >= prevPos)
        cSV--;
    lastSV = cSV;
    if(cSV < 0)
    {
        firstSV = 0;
    }
    else
    {
        
        while (cSV > 0 && prevPos - vt[listSV[cSV]].prePos <= maxSVdel)
            cSV--;
        
        
        
        if(cSV < noSV)
        {
            LoadUInt(ins_sv_arr+vt[listSV[cSV]].alt, delLen);
            while (vt[listSV[cSV]].prePos +  delLen  <= prevPos)
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
    uint32_t dels, prevDelCount = 0;
    //check how many meaningful DELs before
    back = 1;
    if(firstVarNotZero)
    {
        while (vtStart >= back && prevPos - vt[vtStart - back].prePos <= maxDel)
        {
            if(vtType[vtStart - back] == DEL)// && vt[firstVar - back].prePos + vt[firstVar - back].delLen > prePos)
            {
                prevDelCount++;
            }
            back++;
        }
        
        back--;
        while(back)
        {
            if(vtType[vtStart - back] == DEL)
            {
                if(vt[vtStart - back].prePos + vt[vtStart - back].alt <= prevPos)
                    prevDelCount--;
                else
                    break;
            }
            back--;
        }
    }
    /***/
    
    //check if first chars not deleted
    bool inDelReg = false;
    back = 1;
    dels = prevDelCount;
    while(dels)
    {
        if(vtType[vtStart - back] == DEL)
        {
            dels--;
            inDelReg =  uniq_parts[uniq_var[all_var[vtStart - back]*parts+id/64/blocking]*blocking+(id/64)%blocking] & (1LL << (63-id%64));
            
            if(inDelReg && vt[vtStart - back].prePos + vt[vtStart - back].alt > prevPos)
            {
                break;
            }
            else
                inDelReg = false;
        }
        back++;
    }
    if(!inDelReg)
        for(int32_t sv = lastSV; sv >= (int32_t)firstSV; sv--)
        {
            {

                inDelReg =  uniq_parts[uniq_var[all_var[listSV[sv]]*parts+id/64/blocking]*blocking+(id/64)%blocking] & (1LL << (63-id%64));
                if(inDelReg)
                {
                    LoadUInt(ins_sv_arr+vt[listSV[sv]].alt, delLen);
                    //A_memcpy(&delLen, ins_sv_arr+vt[listSV[sv]].alt, sizeof(uint32_t));
                    
                    if(vt[listSV[sv]].prePos + delLen <= prevPos)
                        inDelReg = false;
                    break;
                }
            }
            
        }
    
    ///
    if(!inDelReg)
    {
        
        
        while (chCovered < length)
        {
            
            
            while( vt[vtNo].prePos+1 < refPos && chCovered < length && vtNo < noVar)
            {
                
               vtPresent = uniq_parts[uniq_var[all_var[vtNo]*parts+id/64/blocking]*blocking+(id/64)%blocking] & (1LL << (63-id%64));
              
                if(vtPresent)
                {
                    if(vtType[vtNo] == DEL || vtType[vtNo] == SV)
                    {
                        if(vtType[vtNo] == DEL)
                        {
                            delLen = vt[vtNo].alt;
                            insLen = 0;
                        }
                        else
                        {
                            LoadUInt(ins_sv_arr+vt[vtNo].alt, delLen);
                            //A_memcpy(&delLen, ins_sv_arr+vt[vtNo].alt, sizeof(uint32_t));
                            LoadUInt(ins_sv_arr+vt[vtNo].alt+sizeof(uint32_t), insLen);
                            //A_memcpy(&insLen, ins_sv_arr+vt[vtNo].alt+sizeof(uint32_t), sizeof(uint32_t));
                        }
                        
                        newRefPos = vt[vtNo].prePos+1+delLen > refPos ? vt[vtNo].prePos+1+delLen : refPos;
                        if(newRefPos > refPos)
                        {
                            refPos = newRefPos;
                            
                           
                            if(insLen != 0)
                            {
                                
                                back = 1;
                                prevIns = 0;
                                if(vtType[vtNo-back] == SV)
                                {
                                    LoadUInt(ins_sv_arr+vt[vtNo-back].alt+sizeof(uint32_t), tempLen);
                                }
                                while(vt[vtNo].prePos == vt[vtNo-back].prePos && (vtType[vtNo-back] == INS || (vtType[vtNo-back] == SV && tempLen != 0)))
                                {
                                    prevIns =  insLen > prevIns?  insLen : prevIns;
                                    back++;
                                }
                                if(prevIns == 0 || insLen > prevIns)
                                {
                                    if(vtType[vtNo] == INS)
                                        tempLen = vt[vtNo].alt+sizeof(uint32_t);
                                    else //SV
                                        tempLen = vt[vtNo].alt+(sizeof(uint32_t)<<1);
                                    
                                    if(firstVT)
                                    {
                                        for (int32_t c = gap; c < insLen && pos < length; c++)
                                        {
                                            result[pos++] = decodeCode(ins_sv_arr[tempLen+c]);
                                            if (pos == max)
                                            {
                                                max = max * 2;
                                                tempResult = ( char *) realloc (result, max * sizeof( char));
                                                if (NULL == tempResult)
                                                {
                                                    printf("Error in realloc, the message is - %s\n", strerror(errno));
                                                    exit(8);
                                                }
                                                else
                                                {
                                                    result = tempResult;
                                                }
                                            }
                                        }
                                        chCovered = chCovered +  insLen-gap;
                                    }
                                    else
                                    {
                                        for (int32_t c = 0; c < insLen && pos < length; c++)
                                        {
                                            result[pos++] = decodeCode(ins_sv_arr[tempLen+c]);
                                            if (pos == max)
                                            {
                                                max = max * 2;
                                                tempResult = ( char *) realloc (result, max * sizeof( char));
                                                if (NULL == tempResult)
                                                {
                                                    printf("Error in realloc, the message is - %s\n", strerror(errno));
                                                    exit(8);
                                                }
                                                else
                                                {
                                                    result = tempResult;
                                                }
                                            }
                                        }
                                        chCovered = chCovered + insLen;
                                    }
                                }
                            }
                        }
                        
                    }                }
                
                vtNo++;
            }
            
            if(chCovered < length )
            {
                if(vtNo >= noVar || vt[vtNo].prePos+1 > refPos)  //no variant at currPos
                {
                    result[pos++] = getRef(refPos++);
                    if (pos  == max)
                    {
                        max = max * 2;
                        
                        tempResult = ( char *) realloc (result, max * sizeof( char));
                        if (NULL == tempResult)
                        {
                            printf("Error in realloc, the message is - %s\n", strerror(errno));
                            exit(8);
                        }
                        else
                        {
                            result = tempResult;
                        }
                    }
                    
                    chCovered++;
                }
                else
                {
                    
                    vtPresent = uniq_parts[uniq_var[all_var[vtNo]*parts+id/64/blocking]*blocking+(id/64)%blocking] & (1LL << (63-id%64));
                    
                    if(vtPresent) //variant found/introduced
                    {
                        switch (vtType[vtNo])
                        {
                            case SNP:
                                //introduce variant
                                
                                result[pos++] = decodeCode(vt[vtNo].alt);
                                if (pos  == max)
                                {
                                    max = max * 2;
                                    tempResult = ( char *) realloc (result, max * sizeof( char));
                                    if (NULL == tempResult)
                                    {
                                        printf("Error in realloc, the message is - %s\n", strerror(errno));
                                        exit(8);
                                    }
                                    else
                                    {
                                        result = tempResult;
                                    }
                                }
                                
                                chCovered++;
                                refPos++;
                                vtNo++;
                                
                                break;
                            case DEL:
                                //introduce variant
                                refPos = refPos + vt[vtNo].alt;
                                vtNo++;
                                
                                break;
                            case INS:
                                //introduce variant
                                LoadUInt(ins_sv_arr+vt[vtNo].alt, insLen);
                                
                                
                                if(firstVT)
                                {
                                    for (int32_t c = gap; c < insLen && pos < length; c++)
                                    {
                                        result[pos++] = decodeCode(ins_sv_arr[vt[vtNo].alt+sizeof(uint32_t)+c]);
                                        if (pos == max)
                                        {
                                            max = max * 2;
                                            tempResult = ( char *) realloc (result, max * sizeof( char));
                                            if (NULL == tempResult)
                                            {
                                                printf("Error in realloc, the message is - %s\n", strerror(errno));
                                                exit(8);
                                            }
                                            else
                                            {
                                                result = tempResult;
                                            }
                                        }
                                    }
                                    chCovered = chCovered +  insLen-gap;
                                    
                                }
                                else
                                {
                                    for (int32_t c = 0; c < insLen && pos < length; c++)
                                    {
                                        result[pos++] = decodeCode(ins_sv_arr[vt[vtNo].alt+sizeof(uint32_t)+c]);
                                        if (pos  == max)
                                        {
                                            max = max * 2;
                                            
                                            tempResult = ( char *) realloc (result, max * sizeof( char));
                                            if (NULL == tempResult)
                                            {
                                                printf("Error in realloc, the message is - %s\n", strerror(errno));
                                                exit(8);
                                            }
                                            else
                                            {
                                                result = tempResult;
                                            }
                                        }
                                    }
                                    chCovered = chCovered + insLen;
                                }
                                
                                vtNo++;
                                
                                break;
                            case SV:
                                LoadUInt(ins_sv_arr+vt[vtNo].alt, delLen);
                               
                                LoadUInt(ins_sv_arr+vt[vtNo].alt+sizeof(uint32_t), insLen);
                                
                                
                                //introduce variant
                                if(insLen == 0)
                                {
                                    refPos = refPos + delLen;
                                    vtNo++;
                                }
                                else
                                {
                                    if(firstVT)
                                    {
                                        for (int32_t c = gap; c < insLen && pos < length; c++)
                                        {
                                            result[pos++] = decodeCode(ins_sv_arr[vt[vtNo].alt+(sizeof(uint32_t)<<1)+c]);
                                            if (pos == max)
                                            {
                                                max = max * 2;
                                                
                                                tempResult = ( char *) realloc (result, max * sizeof( char));
                                                if (NULL == tempResult)
                                                {
                                                    printf("Error in realloc, the message is - %s\n", strerror(errno));
                                                    exit(8);
                                                }
                                                else
                                                {
                                                    result = tempResult;
                                                }
                                                
                                            }
                                        }
                                        
                                        chCovered = chCovered +  insLen-gap;
                                        firstVT = false;
                                    }
                                    else
                                    {
                                        for (int32_t c = 0; c < insLen&& pos < length; c++)
                                        {
                                            result[pos++] =  decodeCode(ins_sv_arr[vt[vtNo].alt+(sizeof(uint32_t)<<1)+c]);
                                            if (pos  == max)
                                            {
                                                max = max * 2;
                                                
                                                tempResult = ( char *) realloc (result, max * sizeof( char));
                                                if (NULL == tempResult)
                                                {
                                                    printf("Error in realloc, the message is - %s\n", strerror(errno));
                                                    exit(8);
                                                }
                                                else
                                                {
                                                    result = tempResult;
                                                }
                                            }
                                        }
                                        chCovered = chCovered + insLen;
                                        
                                    }
                                    refPos = refPos + delLen;
                                    vtNo++;
                                }
                                break;
                        }
                        
                    }
                    else //variant not found/not introduced
                    {
                        vtNo++;
                    }
                    firstVT = false;
                }
            }
            
        }
        
        result[pos++] = '\0';
        
        
        return result;
    }
    else
    {
        free(result);
        return NULL;
    }
}
