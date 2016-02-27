/*
 This file is a part of MuGI software distributed under GNU GPL 2 licence.
 The homepage of the MuGI project is http://sun.aei.polsl.pl/mugi/
 
 Authors: Agnieszka Danek, Sebastian Deorowicz
 
 Version: 1.0
 Date   : 2015-Mar
 */

#include "sa.h"

SA_index::SA_index()
{
    suffix_count = 0;
    n_sparse = 1;
    LUT = (unsigned int * ) malloc ((unsigned int)(LUT_MAX+1) * sizeof(unsigned int));
    SA = NULL;
    
    
}


SA_index::SA_index(unsigned int sparse)
{
    suffix_count = 0;
    n_sparse = sparse;
    LUT = (unsigned int * ) malloc ((unsigned int)(LUT_MAX+1) * sizeof(unsigned int));
    SA = NULL;
    
}


SA_index::~SA_index()
{
    
   // delete [] s_ch;
    if(LUT) free(LUT);
    if(SA)
        delete [] SA;
}







int SA_index::make_sa(CFastaFile * fasta)
{
    if(fasta->data_size > UINT32_MAX)
    {
		fprintf(stderr, "\nSize of input sequence not supported!\n");
		return 0;
	}
    unsigned int n = (unsigned int)fasta->data_size ;
	if(n == 0) {
		fprintf(stderr, "\nEmpty string, nothing to sort!");
		return 0;
	}
	n++; // append the virtual sentinel
//	fprintf(stderr, "Allocating input and output space: %llu bytes = %.2lf MB", (uint64) 4*n, (double)4*n/1024/1024);
//	s_ch=new unsigned char[n];
	unsigned int * SA_temp = new unsigned int[n];
	//if(s_ch==NULL || SA==NULL) {
    if( SA_temp==NULL) {
        //    delete [] s_ch;
        delete [] SA_temp;
		fprintf(stderr, "\nInsufficient memory, exit!");
		return 0;
	}
    
	
    /// create suffix array for reference
	SA_IS(fasta->raw_data, SA_temp, n, 255, sizeof(char), 0);
    
    
    /****/
    
    SA = new unsigned int[n/n_sparse+1];
    char  seq[LUT_ENCODED_CHARS+1];
    seq[LUT_ENCODED_CHARS]='\0';
    
    int32_t encoded , currIdx = 0;
    //unsigned int *LUT;
    
    //From i = 1, because fasta->raw_data + SA_temp[0] is '\0'
	for(unsigned int i = 1; i < n; i++)
	{
        strncpy(seq, (char*)(fasta->raw_data + SA_temp[i]), LUT_ENCODED_CHARS);
       
        if(seq[0] !='N' && seq[0] != 'n')
        //if(encoded != 0xFFFFFFFF)  //only ACGTacgt
        {
            if((SA_temp[i]-1) % n_sparse == 0)
            {
                encoded = Coding::encodeSuffix(seq);
                
                SA[suffix_count] = SA_temp[i];
                //create LUT for SA
                //if(encoded != 0xFFFFFFFF)  //only ACGTacgt
                    if(encoded >= currIdx)
                    {
                        for(int32 j = currIdx; j <=encoded; j++)
                            LUT[j] = suffix_count;
                        currIdx = encoded + 1;
                    }
                
                suffix_count++;
            }
        }
	}
    for(int j = currIdx; j <=LUT_MAX; j++)
        LUT[j] = suffix_count;
    
    delete [] SA_temp;
    

    
    return 0;
    
}







unsigned int SA_index::get_count()
{
    return suffix_count;
}

void SA_index::write_sa(FILE *file)
{
    
    fwrite(&suffix_count,  sizeof(uint32_t), 1, file);
    fwrite(SA, sizeof(uint32_t), suffix_count, file);
}

void SA_index::write_sa(FILE *file, int32_t sparse)
{
    uint32_t count = 0, curr_count = 0;
    for(uint32_t i = 0; i < suffix_count; i++)
    {
        if((SA[i]-1) % sparse == 0)
            count++;
    }
    fwrite(&count, sizeof(uint32_t), 1,  file);

    
    char  seq[LUT_ENCODED_CHARS+1];
    seq[LUT_ENCODED_CHARS]='\0';
    int32_t encoded , currIdx = 0;
    for(uint32_t i = 0; i < suffix_count; i++)
    {
        if((SA[i]-1) % sparse == 0)
        {
            strncpy(seq, (char*)(Data::fasta->raw_data + SA[i]), LUT_ENCODED_CHARS);
       
            encoded = Coding::encodeSuffix(seq);
            
            fwrite(&SA[i], sizeof(uint32_t), 1, file);
            if(encoded >= currIdx)
            {
                for(int32 j = currIdx; j <=encoded; j++)
                    LUT[j] = curr_count;
                currIdx = encoded + 1;
            }
            
            curr_count++;
            
        }
                
                
    }
    
    
    for(int32 j = currIdx; j <=LUT_MAX; j++)
        LUT[j] = curr_count;

    
    
}

void  SA_index::write_sa_lut(FILE *file)
{
    fwrite(LUT, sizeof(unsigned int), LUT_MAX+1, file);
}


/*********/

void SA_index::readData(FILE *in)
{
    fread(&suffix_count, 1, sizeof(uint32_t), in);
    SA = (uint32_t *) malloc (suffix_count * sizeof(uint32_t));
    fread(SA, sizeof(uint32_t), suffix_count, in);
    LUT = (uint32_t * ) malloc ((uint32_t)(LUT_MAX+1) * sizeof(uint32_t));
    fread(LUT, sizeof(unsigned int), LUT_MAX+1, in);
    

}


void SA_index::calcLUT(CFastaFile * fasta)
{
     int32_t encoded , currIdx = 0, count = 0;
    char  seq[LUT_ENCODED_CHARS+1];
    seq[LUT_ENCODED_CHARS]='\0';
    for(uint32_t i = 0; i < suffix_count; i++)
	{
        for(int j = 0; j < LUT_ENCODED_CHARS; j++)
            seq[j] = fasta->getRefFromCompressed(SA[i]+j);
        encoded = Coding::encodeSuffix(seq);
        
        //if(seq[0] !='N' && seq[0] != 'n')
            //if(encoded != 0xFFFFFFFF)  //only ACGTacgt
        {
          //  if((SA_temp[i]-1) % n_sparse == 0)
            {
            //    SA[suffix_count] = SA_temp[i];
                //create LUT for SA
                //if(encoded != 0xFFFFFFFF)  //only ACGTacgt
                if(encoded >= currIdx)
                {
                    for(int32 j = currIdx; j <=encoded; j++)
                        LUT[j] = count;
                    currIdx = encoded + 1;
                }
                
                count++;
            }
        }
	}
    for(int j = currIdx; j <= LUT_MAX; j++)
        LUT[j] = count;
    
    
    
    for(int i = 0; i <= LUT_MAX; i++)
        cout << i << " " << LUT[i] << endl;

}