/*
 This file is a part of MuGI software distributed under GNU GPL 2 licence.
 The homepage of the MuGI project is http://sun.aei.polsl.pl/mugi/
 
 Authors: Agnieszka Danek, Sebastian Deorowicz
 
 Version: 1.0
 Date   : 2015-Mar
 */

#include "sa_variants.h"


int SAVariants::make_sa_variants(VariantList * _vl, BitVectors * _bv, CFastaFile * _fasta, uint32_t _kl)
{
    kl = _kl;
    vl = _vl;
    bv = _bv;
    fasta = _fasta;

    no_var = vl->getNoVar();
    
    uint32_t logTreshold = no_var/10;
    
    
    sa1 = new SA1(n_sparse, no_var*kl, kl, no_var, bv, vl, fasta);
    
    sa2 = new SA2(n_sparse, vl->get_no_ins() * 10 + 10);
    sa3 = new SA3(n_sparse, vl->get_no_sv() * 10 + 10);
    
    uint32_t currSV = 0, lastPrePos = 0, currPrePos;
    
    //for first variant and zero position
    uint32_t i = 0;
    while(vl->getPos(i) == 0xFFFFFFFF) //wariants at first pos, so prePos = -1, uint: 0xFFFFFFFF)
    {
        if(i > logTreshold)
        {
            cout  << (((int)(((float)(logTreshold)/(float)(no_var))*100))+1)/10*10 << "% done" << endl;
            logTreshold = logTreshold + no_var/10;
        }

        
        if(vl->getType(i) == SV && i != vl->listSV[0])
            currSV++;
        
        unsigned int vtListInit[1];
        vtListInit[0] = 0;
        currPrePos = 0xFFFFFFFF;
        addSuffixes(currPrePos, 1U, vtListInit, 0U, 0U, i, currPrePos+1, 0U, false, 0U, 0, 0, currSV);
        
        i++;
    }
    if(i == 0)
    {
        
        if(i > logTreshold)
        {
            cout  << (((int)(((float)(logTreshold)/(float)(no_var))*100))+1)/10*10 << "% done" << endl;
            logTreshold = logTreshold + no_var/10;
        }
        
        
        if(vl->getType(i) == SV && i != vl->listSV[0])
            currSV++;
        
        currPrePos = vl->getPos(i) > (kl-1) ? vl->getPos(i) - (kl-1) : 0;
        for(; currPrePos <= vl->getPos(i); currPrePos++)
        {
            if(currPrePos > lastPrePos || currPrePos == 0)
            {
                lastPrePos = currPrePos;
                
                unsigned int vtListInit[1];
                vtListInit[0] = 0;
                addSuffixes(currPrePos, 1U, vtListInit, 0U, 0U, i, currPrePos+1, 0U, false, 0U, 0, 0, currSV);
                
            }
        }
        
    }
    
    for (i = 1; i < no_var ; i++)
    {
        
        if(i > logTreshold)
        {
            cout  << (((int)(((float)(logTreshold)/(float)(no_var))*100))+1)/10*10 << "% done" << endl;
            logTreshold = logTreshold + no_var/10;
        }
        
        
        if(vl->getType(i) == SV && i != vl->listSV[0])
            currSV++;
        
        currPrePos = vl->getPos(i) > (kl-1) ? vl->getPos(i) - (kl-1) : 0;
        for(; currPrePos <= vl->getPos(i); currPrePos++)
        {
            if(currPrePos > lastPrePos)
            {
                lastPrePos = currPrePos;
             
                unsigned int vtListInit[1];
                vtListInit[0] = 0;
                addSuffixes(currPrePos, 1U, vtListInit, 0U, 0U, i, currPrePos+1, 0U, false, 0U, 0, 0, currSV);
                
            }
        }
        
    }

   
    sort();
    
    
    return 0;
}


void SAVariants::addSuffixes(unsigned int prePos, unsigned int vtListSize, unsigned int * vtList, unsigned int gap, unsigned int noVt, unsigned int currVt, unsigned int refPos, unsigned short chCovered, bool delRegVariatns, unsigned int skippedVt, unsigned int lastFiltFound, unsigned int meaningfulVTinDelReg, unsigned int currentSV)
{
    
  
    
    
    unsigned int * newVtList;
  
    newVtList = (unsigned int *) malloc (vtListSize * sizeof(unsigned int));
    if (NULL == newVtList) {
        printf("Error in malloc, the message is - %s\n", strerror(errno));
        exit(8);
    }
    for (int i = 0; i < vtListSize; i++)
    {
        newVtList[i] = vtList[i];
    }
    
    if(noVt == vtListSize * NO_BITS_UINT)
    {
        vtListSize++;
        unsigned int * tempList = (unsigned int *) realloc(newVtList, vtListSize * sizeof(unsigned int));
        if (NULL ==tempList)
        {
            printf("Error in realloc, the message is - %s\n", strerror(errno));
            exit(8);
        }
        else
        {
            newVtList = tempList;
        }
        newVtList[vtListSize-1] = 0;
    }
    
    unsigned int last = 0;//, prevIns;
    //int i;
    //bool suffix, bitvector;
    //unsigned long int temp, idx, skippedVtNow, refPosTemp;
    //unsigned long tNo;
    unsigned int lastFiltFoundTemp;
    unsigned int newRefPos;
    
    
    while (chCovered < kl)
    {
        
        
        //in case of previous del, some succiding variants may not have impact on suffix...
        while(currVt < no_var && vl->getPos(currVt)+1 < refPos && chCovered < kl)
        {
            delRegVariatns = true;
            
            if(vl->getType(currVt) == DEL || vl->getType(currVt) == SV)
            {
                newRefPos = vl->getPos(currVt)+1+vl->getDelLen(currVt) > refPos ? vl->getPos(currVt)+1+vl->getDelLen(currVt) : refPos;
                
                if(newRefPos > refPos)
                {
                    
                    last = newVtList[noVt/NO_BITS_UINT];
                    newVtList[noVt/NO_BITS_UINT] = newVtList[noVt/NO_BITS_UINT] | (1U << (NO_BITS_UINT -  noVt%NO_BITS_UINT -  1U));
                    
                    if (vl->getType(currVt) == SV)
                        currentSV++;
                    lastFiltFoundTemp = lastFiltFound;
                    
                    if(vl->getInsLen(currVt) == 0)
                    {
                        addSuffixes(prePos, vtListSize, newVtList, gap, noVt+1, currVt+1, newRefPos, chCovered, delRegVariatns, skippedVt,lastFiltFoundTemp, meaningfulVTinDelReg+1, currentSV);
                    }
                    else
                    {
                        ///
                        if(chCovered + vl->getInsLen(currVt) >= kl)
                        {
                            sa3->newSuffix3(prePos, vtListSize, newVtList, gap, noVt+1, currVt-noVt-skippedVt, lastFiltFound, currentSV);
                            
                        }
                        else
                            addSuffixes(prePos, vtListSize, newVtList, gap, noVt+1, currVt+1, refPos+vl->getDelLen(currVt), chCovered+vl->getInsLen(currVt), delRegVariatns, skippedVt, lastFiltFound, meaningfulVTinDelReg, currentSV);
                        ///
                    }
                    
                    newVtList[noVt/NO_BITS_UINT] = last;
                    noVt++;
                    
                }
                else
                    skippedVt++;
                
            }
            else
            {
                skippedVt++;
            }
            
            currVt++;
            
        }
        
        if( chCovered < kl)
        {
            if(currVt < no_var && vl->getPos(currVt)+1 == refPos)
            {
                switch (vl->getType(currVt))
                {
                    case SNP:
                        //introduce variant
                        //
                        last = newVtList[noVt/NO_BITS_UINT];
                        newVtList[noVt/NO_BITS_UINT] = newVtList[noVt/NO_BITS_UINT] | (1U << (NO_BITS_UINT -  noVt%NO_BITS_UINT -  1U));
                        addSuffixes(prePos, vtListSize, newVtList, gap, noVt+1, currVt+1, refPos+1, chCovered+1, delRegVariatns, skippedVt, lastFiltFound, meaningfulVTinDelReg, currentSV);
                        
                        break;
                    case DEL:
                        //introduce variant
                        // if(chCovered > 0)
                    {
                        last = newVtList[noVt/NO_BITS_UINT];
                        newVtList[noVt/NO_BITS_UINT] = newVtList[noVt/NO_BITS_UINT] | (1U << (NO_BITS_UINT -  noVt%NO_BITS_UINT -  1U));
                        addSuffixes(prePos, vtListSize, newVtList, gap, noVt+1, currVt+1, refPos+vl->getDelLen(currVt), chCovered, delRegVariatns, skippedVt, lastFiltFound, meaningfulVTinDelReg, currentSV);
                    }
                        break;
                    case INS:
                        //introduce variant
                        //
                        last = newVtList[noVt/NO_BITS_UINT];
                        newVtList[noVt/NO_BITS_UINT] = newVtList[noVt/NO_BITS_UINT] | (1U << (NO_BITS_UINT -  noVt%NO_BITS_UINT -  1));
                        if(noVt > 0 || chCovered > 0)
                        {
                            if(chCovered + vl->getInsLen(currVt) >= kl)
                            {
                                if(vtListSize == 1 && delRegVariatns == false)
                                {
                                    if (gap == 0)
                                    {
                                        if(prePos % n_sparse == 0 || (vl->getType(currVt-noVt) != SNP && prePos == vl->getPos(currVt-noVt)))
                                            sa1->newSuffix1(prePos, newVtList[0], noVt+1, currVt-noVt, currentSV);
                                    }
                                    else
                                        sa2->newSuffix2(prePos, newVtList[0], gap, noVt+1, currVt-noVt, currentSV);
                                }
                                else
                                {
                                    sa3->newSuffix3(prePos, vtListSize, newVtList, gap, noVt+1, currVt-noVt-skippedVt, lastFiltFound, currentSV);
                                    // if(delRegVariatns && k > 8)
                        
                                }
                            }
                            else
                                addSuffixes(prePos, vtListSize, newVtList, gap, noVt+1, currVt+1, refPos, chCovered + vl->getInsLen(currVt), delRegVariatns, skippedVt, lastFiltFound, meaningfulVTinDelReg, currentSV);
                        }
                        else if (gap == 0)//suffix start with insertion // //it's first chCovered, so no prev variants (vtListSize==1)
                        {
                            //no gap, start suffix at firts pos of alt
                            if (vl->getInsLen(currVt) >= kl)
                            {
                                if(prePos % n_sparse == 0 || (vl->getType(currVt-noVt) != SNP && prePos == vl->getPos(currVt-noVt)))
                                    sa1->newSuffix1(prePos, newVtList[0], noVt+1, currVt-noVt, currentSV);
                            }
                            else
                            {
                                addSuffixes(prePos, vtListSize, newVtList, gap, noVt+1, currVt+1, refPos, chCovered+vl->getInsLen(currVt), delRegVariatns, skippedVt, lastFiltFound, meaningfulVTinDelReg, currentSV);
                            }
                            
                            //some gap between prevPos and suffix start
                            unsigned int i;
                            for(i = 1; i <= vl->getInsLen(currVt); i++)
                            {
                                //unsigned long a = strlen(&vt[currVt].alt[i]);
                                if (vl->getInsLen(currVt)-i >= kl)
                                {
                                    sa2->newSuffix2(prePos, newVtList[0], i, noVt+1, currVt-noVt, currentSV);
                                }
                                else
                                {
                                    addSuffixes(prePos, vtListSize, newVtList, i, noVt+1, currVt+1, refPos, chCovered+vl->getInsLen(currVt)-i, delRegVariatns, skippedVt, lastFiltFound, meaningfulVTinDelReg, currentSV);
                                }
                            }
                        }
                        else
                        {
                            cout << "sth wrong with vtList/gap\n";
                            exit(9);
                        }
                        break;
                    case SV:
                        //introduce variant
                        //no insertion -> treat like simple deletion
                        if(vl->getInsLen(currVt) == 0)
                        {
                            //  if(chCovered > 0)
                            {
                                last = newVtList[noVt/NO_BITS_UINT];
                                newVtList[noVt/NO_BITS_UINT] = newVtList[noVt/NO_BITS_UINT] | (1 << (NO_BITS_UINT -  noVt%NO_BITS_UINT -  1));
                                addSuffixes(prePos, vtListSize, newVtList, gap, noVt+1, currVt+1, refPos+vl->getDelLen(currVt), chCovered, delRegVariatns, skippedVt, lastFiltFound, meaningfulVTinDelReg, currentSV);
                            }
                        }
                        else    // deletion+insertion
                        {
                            last = newVtList[noVt/NO_BITS_UINT];
                            newVtList[noVt/NO_BITS_UINT] = newVtList[noVt/NO_BITS_UINT] | (1 << (NO_BITS_UINT -  noVt%NO_BITS_UINT -  1));
                            if(noVt > 0 || chCovered > 0)
                            {
                                if(chCovered + vl->getInsLen(currVt) >= kl)
                                {
                                    if(vtListSize == 1 && delRegVariatns == false)
                                    {
                                        if (gap == 0)
                                        {
                                            if(prePos % n_sparse == 0 || (vl->getType(currVt-noVt) != SNP && prePos == vl->getPos(currVt-noVt)))
                                                sa1->newSuffix1(prePos, newVtList[0], noVt+1, currVt-noVt, currentSV);
                                        }
                                        else
                                            sa2->newSuffix2(prePos, newVtList[0], gap, noVt+1, currVt-noVt, currentSV);
                                    }
                                    else
                                    {
                                        sa3->newSuffix3(prePos, vtListSize, newVtList, gap, noVt+1, currVt-noVt-skippedVt, lastFiltFound, currentSV);
                            
                                    }
                                }
                                else
                                    addSuffixes(prePos, vtListSize, newVtList, gap, noVt+1, currVt+1, refPos+vl->getDelLen(currVt), chCovered+vl->getInsLen(currVt), delRegVariatns, skippedVt, lastFiltFound, meaningfulVTinDelReg, currentSV);
                            }
                            else if(gap == 0)//suffix start with insertion
                            {
                                //no gap, start suffix at firts pos of alt
                                if (vl->getInsLen(currVt) >= kl)
                                {
                                    if(prePos % n_sparse == 0 || (vl->getType(currVt-noVt) != SNP && prePos == vl->getPos(currVt-noVt)))
                                        sa1->newSuffix1(prePos, newVtList[0], noVt+1, currVt-noVt, currentSV);
                                }
                                else
                                {
                                    addSuffixes(prePos, vtListSize, newVtList, gap, noVt+1, currVt+1, refPos+vl->getDelLen(currVt), chCovered+vl->getInsLen(currVt), delRegVariatns, skippedVt, lastFiltFound, meaningfulVTinDelReg, currentSV);
                                }
                                
                                //some gap between prevPos and suffix start
                                unsigned int i;
                                for(i = 1; i <= vl->getInsLen(currVt); i++)
                                {
                                    //unsigned long a = strlen(&vt[currVt].alt[i]);
                                    //if (strlen(&vt[currVt].alt[i]) >= kl)
                                    if (vl->getInsLen(currVt)-i >= kl)
                                    {
                                        sa2->newSuffix2(prePos, newVtList[0], i, noVt+1, currVt-noVt, currentSV);
                                    }
                                    else
                                    {
                                        addSuffixes(prePos, vtListSize, newVtList, i, noVt+1, currVt+1, refPos+vl->getDelLen(currVt), chCovered+vl->getInsLen(currVt)-i, delRegVariatns, skippedVt, lastFiltFound, meaningfulVTinDelReg, currentSV);
                                    }
                                }
                            }
                            else
                            {
                                cout << "sth wrong with vtList/gap\n";
                                exit(9);
                            }
                        }
                        
                        break;
                }
                
                //don't introduce variant
                newVtList[noVt/NO_BITS_UINT] = last;
                currVt++;
                noVt++;
                
                if(noVt == vtListSize * NO_BITS_UINT)
                {
                    
                    vtListSize++;
                    unsigned int * tempList = (unsigned int *) realloc(newVtList, vtListSize * sizeof(unsigned int));
                    if (NULL ==tempList) {
                        printf("Error in realloc, the message is - %s\n", strerror(errno));
                        exit(8);
                    }
                    else
                    {
                        newVtList = tempList;
                    }
                    newVtList[vtListSize-1]=0;
                }
                
            }
            else
            {
                //	cout << "chCovered++ " << chCovered << endl;
                if (refPos != fasta->data_size)
                {
                    refPos++;
                    chCovered++;
                }
                else
                    break;
                //cout << " cc = " << chCovered << endl;
                //fflush(stdout);
            }
        }
    }
    
    if(vtListSize == 1 && delRegVariatns == false)
    {
        if(newVtList[0] > 0)
        {
            if (gap == 0)
            {
                if(prePos % n_sparse == 0 || (vl->getType(currVt-noVt) != SNP && prePos == vl->getPos(currVt-noVt)))
                    sa1->newSuffix1(prePos, newVtList[0], noVt, currVt-noVt, currentSV);
            }
            else
                sa2->newSuffix2(prePos, newVtList[0], gap, noVt, currVt-noVt, currentSV);
        }
        
        free(newVtList);
    }
    else
    {
        bool newVt = false;
        for (int j = 0; j < vtListSize; j++)
            if(newVtList[j] > 0)
                newVt = true;
        if(newVt)
        {
            sa3->newSuffix3(prePos, vtListSize, newVtList, gap, noVt, currVt-noVt-skippedVt, lastFiltFound, currentSV);
            
        }
        
        free(newVtList);
    }
}
void  SAVariants::sort()
{
    

    
    cout << "Sorting " << sa3->suffixArray3Count << " suffixes in SA3..." <<  endl;
    fflush(stdout);
    sa3->sort();
    
    cout << "Sorting " << sa2->suffixArray2Count << " suffixes in SA2..." <<  endl;
    fflush(stdout);
    sa2->sort();
    
    cout << "Sorting " << sa1->suffixArray1Count << " suffixes in SA1..." <<  endl;
    fflush(stdout);
    sa1->sort();
}


void  SAVariants::write_sa1(FILE * file)
{
    cout << "Writing sa1.. " <<endl;
    sa1->write_sa(file);
}

void  SAVariants::write_sa1(FILE * file, int32_t sparse)
{
    cout << "Writing sa1.. " <<endl;
    sa1->write_sa(file, sparse);
}


void  SAVariants::write_sa1_lut(FILE * file)
{
    //cout << "Writing sa1 LUT.. " <<endl;
    sa1->write_sa_lut(file);
}


void  SAVariants::write_sa2(FILE * file)
{
    cout << "Writing sa2.. " <<endl;
    sa2->write_sa(file);
}


void  SAVariants::write_sa2_lut(FILE * file)
{
    //cout << "Writing sa2 LUT.. " <<endl;
    sa2->write_sa_lut(file);
}


void  SAVariants::write_sa3(FILE * file)
{
    cout << "Writing sa3.. " <<endl;
    sa3->write_sa(file);
}


void  SAVariants::write_sa3_lut(FILE * file)
{
    //cout << "Writing sa3 LUT.. " <<endl;
    sa3->write_sa_lut(file);
}



/*******/
void SAVariants::readSA123andSort(FILE *in)
{
    
    cout << "sizeof(suffix1)  " << sizeof(suffix1) << endl;
    sa1 = new SA1();
    
    cout << "sizeof(suffix2)  " << sizeof(suffix2) << endl;
    sa2 = new SA2();
    
    cout << "sizeof(suffix3)  " << sizeof(suffix3) << endl;
    sa3 = new SA3();
    
    sa1->readData(in);
    sa2->readData(in);
    sa3->readData(in);
    
    
    sort();
}
