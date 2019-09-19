/*
 This file is a part of MuGI software distributed under GNU GPL 2 licence.
 The homepage of the MuGI project is http://sun.aei.polsl.pl/mugi/
 
 Authors: Agnieszka Danek, Sebastian Deorowicz
 
 Version: 1.0
 Date   : 2015-Mar
 */

#include <iostream>
#include "defs.h"
#include "data.h"
#include "sa.h"
#include "variant_list.h"
#include "bit_vectors.h"
#include "vcf.h"
#include "sa_variants.h"

#include "fasta.h"

#include "tbb/task_scheduler_init.h"

#define LUT_ENCODED_CHARS 12
#define LUT_MAX  ((1U<<(LUT_ENCODED_CHARS*2)))//16777216


uint32_t n_sparse, ploidity, kl, blocking;


#define BLOCKING 3
#define DEFAULT_PLOIDITY 2

using namespace std;
//input: REF + VCF
int main(int argc, const char * argv[])
{

 //
    
  
    tbb::task_scheduler_init TBBinit;
    //create index based on VCF and REF file
#ifndef CREATE_ALL_INDEXES
    blocking = BLOCKING ;
    // Read basic parameters
    if(!argv[1] || !argv[2] ||!argv[3] || !argv[4])
    {
        cout << "Wrong number of parameters." << endl;
        cout << "Required parameters:\n[fasta] [vcf] [k_len] [sparsity] {[ploidity]}" << endl;
        cout << "where:\n[fasta] path to file with the reference sequence," << endl;
        cout << "[vcf] path to file with the VCF file" << endl;
        cout << "[k_len] length of created k-mers" << endl;
        cout << "[sparsity] sparsity for SA and SA1 suffix arrays (1-16 recommended)" << endl;
        cout << "[ploidity] optional, 1 or 2 (default 2)." << endl;
        exit(1);
    }
    kl = atoi(argv[3]);
    n_sparse = atoi(argv[4]);
    
    if(argv[5])
    {
        ploidity = atoi(argv[5]);
        if(ploidity != 1 && ploidity != 2)
        {
            cout << "Ploidity " << ploidity << " is not supported." << endl;
            exit(1);
        }

    }
    else
        ploidity = DEFAULT_PLOIDITY;
    
    
    
    Data::kl = kl;
    Data::ploidity = ploidity;
    Data::n_sparse = n_sparse;
    Data::blocking = blocking;
    
    
    std::stringstream out_name;
    out_name << "index-" << kl << "-" << n_sparse;
    
    
    // Open file to write all index data
    FILE * out = fopen(out_name.str().c_str() , "wb");
    if (out == NULL) {
        printf("Cannot open output file %s\n", out_name.str().c_str());
        printf("The message is - %s\n", strerror(errno));
        exit(1);
    }
    // Write index parameters
    fwrite(&Data::kl, 1, sizeof(uint32_t), out);
    fwrite(&Data::n_sparse, 1, sizeof(uint32_t), out);
    uint32_t temp = LUT_ENCODED_CHARS;
    fwrite(&temp, 1, sizeof(uint32_t), out);
    
    
    //Open reference sequence
    if(!Data::fasta->Open(argv[1]))
    {
        cout << "File " << argv[1] << " does not exists or is not in FASTA format\n";
        exit(1);
    }
    cout << "Reference sequence: " << argv[1] << endl;
    //Read reference sequence
    Data::fasta->Read(Data::data);
    // Compress reference sequence
    cout << "Compressing reference sequence..." << endl;
    Data::fasta->compressSequence();
    // Write compressed sequence to file
    Data::fasta->writeCompressedSequenceToFile(out);
    cout << "Reference sequence compressed and written to the output file." << endl;
     fflush(out);
    
    cout << "VCF file: " << argv[2] << endl;
    //Open VCF file and read it's metadata
    VCF * vcf = new VCF();
    if(!vcf->Open(argv[2]))
    {
        cout << "File " << argv[2] << " does not exists or is not in VCF format\n";
        exit(1);
    }
    cout << "There are " << vcf->get_no_var() << " variant calls for " << vcf->get_no_gen();
    if(ploidity == 2)
        cout << " diploid individuals." << endl;
    else
        cout << " haploid individuals." << endl;
    
    // Create Container for bit vectors and for variant list
    Data::bv = new BitVectors(ploidity, vcf->get_no_gen() , vcf->get_no_var(), vcf->get_gen_names());
    Data::vl = new VariantList(vcf->get_no_var(), Data::fasta->data_size);
    
    
    // Create a dictionary file with data about all chromosomes (name, start pos, size) and
    // all individuals IDs (diploid individuals divided totwo chromosome sets named: ID-1 and ID-2)
    cout << "Creating a dictionary file." << endl;
    ofstream dict;
    dict.open("index-dict");
    dict << Data::fasta->sequences.size() << endl;
    for(int32 i = 0; i < Data::fasta->sequences.size(); ++i)
    {
        dict << Data::fasta->sequences[i].seq_name << endl;
        dict << Data::fasta->sequences[i].start_pos << "\t" << Data::fasta->sequences[i].size << endl;
    }
    dict << Data::bv->getNoId() << endl;
    for(int32_t i = 0; i < Data::bv->getNoId(); ++i)
    {
        dict << Data::bv->genomesNames[i] << endl;
    }
    dict.close();
    
    
    // Process VCF file to create variant list file and bit vectors
    cout << "Processing VCF data... " << endl;
    if(!vcf->ProcessVCFData(Data::vl, Data::bv))
        exit(1);
    
    Data::vl->set_no_var(vcf->get_no_var());
    // Close vcf file
    vcf->Close();
    
    cout << "Compressing variant list... " << endl;
    // Transform variant list to more compact representation
    Data::vl->transform();
    //Write compact variant list to the output file
    Data::vl->write_compact_vl(out);
    cout << "Variant list compressed and written to the output file." << endl;
     fflush(out);
    
    // Remove info about meaningless variants found in every individual form bit vectors representing them
    cout << "Filtering info about found variants... " << endl;
    Data::bv->removeVariantsSetInDelRegions(Data::vl);
    // compress bit variant collection
    cout << "Compressing data about found variants... " << endl;
    Data::bv->createCompressedCollection(blocking);
    cout << "Writing compressed collection to the output file.. " << endl;
    Data::bv->writeCompressedCollectionToFile(out);
    cout << "Data about found variants compressed and written to the output file." << endl;
    fflush(out);
   
    cout << "Creating suffix array (sparsity = " << n_sparse << ") for the reference... " << endl;
    // Create suffix array for the reference sequence
    SA_index * sa = new SA_index(n_sparse);
    sa->make_sa(Data::fasta);
    
    // Write LUT and suffix array to the output file
    sa->write_sa(out);
    sa->write_sa_lut(out);
    delete sa;
    cout << "Suffix array created and written to the output file." << endl;
    fflush(out);

    
    cout << "Creating suffix arrays (sparsity = " << Data::n_sparse << ") for all subsequences of length " << Data::kl << " with variants... " << endl;
    // Create suffix arrays for sequences with variants
    SAVariants * saVar = new SAVariants(Data::n_sparse);
    saVar->make_sa_variants(Data::vl, Data::bv, Data::fasta, Data::kl);
    
    // Write LUTs and suffix arrays to the output file
    cout << "Writing to file..." <<endl;
    saVar->write_sa1(out);
    saVar->write_sa1_lut(out);
    
    //cout << "size without sa2, sa3 " << ftell(out) << endl;
    
    saVar->write_sa2(out);
    saVar->write_sa2_lut(out);
    
    saVar->write_sa3(out);
    saVar->write_sa3_lut(out);
    cout << "Suffix arrays created and written to the output file." << endl;
    
    fclose(out);
    
    
    
    
#else
    
    if(!argv[1] || !argv[2] ||!argv[3] || !argv[4])
    {
        cout << "Wrong number of parameters." << endl;
        cout << "Required parameters:\n[fasta] [vcf] [k_len] [max_sparsity] {[ploidity]}" << endl;
        cout << "where:\n[fasta] path to file with the reference sequence," << endl;
        cout << "[vcf] path to file with the VCF file" << endl;
        cout << "[k_len] length of created k-mers" << endl;
        cout << "[max_sparsity] max sparsity for SA and SA1 suffix arrays (max_sparsity indexes will be created, with sparsities from 1 to max_sparsity)" << endl;
        cout << "[ploidity] optional, 1 or 2 (default 2)." << endl;
        exit(1);
    }
    int max_sparsity = atoi(argv[4]);
    
    if(argv[5])
    {
        ploidity = atoi(argv[5]);
        if(ploidity != 1 && ploidity != 2)
        {
            cout << "Ploidity " << ploidity << " is not supported." << endl;
            exit(1);
        }
        
    }
    else
        ploidity = DEFAULT_PLOIDITY;
    
    
    blocking = BLOCKING ;
    // Read basic parameters
    kl = atoi(argv[3]);
    n_sparse = 1;
    
    if(argv[5])
    {
        ploidity = atoi(argv[5]);
        if(ploidity != 1 && ploidity != 2)
        {
            cout << "Ploidity " << ploidity << " is not supported." << endl;
            exit(1);
        }
        
    }
    else
        ploidity = DEFAULT_PLOIDITY;
    
    
    
    Data::kl = kl;
    Data::ploidity = ploidity;
    Data::n_sparse = n_sparse;
    Data::blocking = blocking;
    
    
    std::stringstream out_name[max_sparsity];
    FILE * out[max_sparsity];
    
    for(int i = 1; i <= max_sparsity; i++)
    {
        out_name[i-1] << "index-" << kl << "-" << i;
    
        // Open file to write all index data
            out[i-1] = fopen(out_name[i-1].str().c_str() , "wb");
        
        // Write index parameters
        fwrite(&Data::kl, 1, sizeof(uint32_t), out[i-1]);
        fwrite(&i, 1, sizeof(uint32_t), out[i-1]);
        uint32_t temp = LUT_ENCODED_CHARS;
        fwrite(&temp, 1, sizeof(uint32_t), out[i-1]);
    }
    
    //Open reference sequence
    if(!Data::fasta->Open(argv[1]))
    {
        cout << "File " << argv[1] << " does not exists or is not in FASTA format\n";
        exit(1);
    }
    cout << "Reference sequence: " << argv[1] << endl;
    //Read reference sequence
    Data::fasta->Read(Data::data);
    // Compress reference sequence
    cout << "Compressing reference sequence..." << endl;
    Data::fasta->compressSequence();
    // Write compressed sequence to file
    for(int i = 0; i < max_sparsity; i++)
    {
        Data::fasta->writeCompressedSequenceToFile(out[i]);
        fflush(out[i]);
    }
    cout << "Reference sequence compressed and written to the output file." << endl;
    
    
    cout << "VCF file: " << argv[2] << endl;
    //Open VCF file and read it's metadata
    VCF * vcf = new VCF();
    if(!vcf->Open(argv[2]))
    {
        cout << "File " << argv[2] << " does not exists or is not in VCF format\n";
        exit(1);
    }
    cout << "There are " << vcf->get_no_var() << " variant calls for " << vcf->get_no_gen();
    if(ploidity == 2)
        cout << " diploid individuals." << endl;
    else
        cout << " haploid individuals." << endl;
    
    // Create Container for bit vectors and for variant list
    Data::bv = new BitVectors(ploidity, vcf->get_no_gen() , vcf->get_no_var(), vcf->get_gen_names());
    Data::vl = new VariantList(vcf->get_no_var(), Data::fasta->data_size);
    
    
    // Create a dictionary file with data about all chromosomes (name, start pos, size) and
    // all individuals IDs (diploid individuals divided totwo chromosome sets named: ID-1 and ID-2)
    cout << "Creating a dictionary file." << endl;
    ofstream dict;
    dict.open("index-dict");
    dict << Data::fasta->sequences.size() << endl;
    for(int32 i = 0; i < Data::fasta->sequences.size(); ++i)
    {
        dict << Data::fasta->sequences[i].seq_name << endl;
        dict << Data::fasta->sequences[i].start_pos << "\t" << Data::fasta->sequences[i].size << endl;
    }
    dict << Data::bv->getNoId() << endl;
    for(int32_t i = 0; i < Data::bv->getNoId(); ++i)
    {
        dict << Data::bv->genomesNames[i] << endl;
    }
    dict.close();
    
    
    // Process VCF file to create variant list file and bit vectors
    cout << "Processing VCF data... " << endl;
    if(!vcf->ProcessVCFData(Data::vl, Data::bv))
        exit(1);
    // Close vcf file
    vcf->Close();
    
    cout << "Compressing variant list... " << endl;
    // Transform variant list to more compact representation
    Data::vl->transform();
    
    
    //Write compact variant list to the output file
    for(int i = 0; i < max_sparsity; i++)
    {
        Data::vl->write_compact_vl(out[i]);
        fflush(out[i]);
    }
    cout << "Variant list compressed and written to the output file." << endl;
    
    
    // Remove info about meaningless variants found in every individual form bit vectors representing them
    cout << "Filtering info about found variants... " << endl;
    Data::bv->removeVariantsSetInDelRegions(Data::vl);
    // compress bit variant collection
    cout << "Compressing data about found variants... " << endl;
    Data::bv->createCompressedCollection(blocking);
    cout << "Writing compressed collection to the output file.. " << endl;
    
    for(int i = 0; i < max_sparsity; i++)
    {
        Data::bv->writeCompressedCollectionToFile(out[i]);
        fflush(out[i]);
    }
    cout << "Data about found variants compressed and written to the output file." << endl;
    
    cout << "Creating suffix array (sparsity = " << n_sparse << ") for the reference... " << endl;
    // Create suffix array for the reference sequence
    SA_index * sa = new SA_index(n_sparse);
    sa->make_sa(Data::fasta);
    
    // Write LUT and suffix array to the output file
    sa->write_sa(out[0]);
    sa->write_sa_lut(out[0]);
    fflush(out[0]);
    for(uint32_t i = 1; i < max_sparsity; i++)
    {
        sa->write_sa(out[i], i+1);
        sa->write_sa_lut(out[i]);
        fflush(out[i]);
    }
    delete sa;
    cout << "Suffix array created and written to the output file." << endl;
    
    
    
    cout << "Creating suffix arrays (sparsity = " << Data::n_sparse << ") for all subsequences of length " << Data::kl << " with variants... " << endl;
    // Create suffix arrays for sequences with variants
    SAVariants * saVar = new SAVariants(Data::n_sparse);
    saVar->make_sa_variants(Data::vl, Data::bv, Data::fasta, Data::kl);
    
    // Write LUTs and suffix arrays to the output file
    cout << "Writing to file..." <<endl;
    saVar->write_sa1(out[0]);
    saVar->write_sa1_lut(out[0]);
    
    saVar->write_sa2(out[0]);
    saVar->write_sa2_lut(out[0]);
    
    saVar->write_sa3(out[0]);
    saVar->write_sa3_lut(out[0]);
    for(int32_t i = 1; i < max_sparsity; i++)
    {
        saVar->write_sa1(out[i], i+1);
        saVar->write_sa1_lut(out[i]);
        
        saVar->write_sa2(out[i]);
        saVar->write_sa2_lut(out[i]);
        
        saVar->write_sa3(out[i]);
        saVar->write_sa3_lut(out[i]);
    }
    cout << "Suffix arrays created and written to the output file." << endl;
    
    for(int32_t i = 0; i < max_sparsity; i++)
        fclose(out[i]);
    
    

    
    
#endif
     /**/


    delete vcf;
    delete saVar;
    delete Data::bv;
    delete Data::vl;
      
    
    
    return 0;
}

