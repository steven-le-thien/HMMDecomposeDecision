#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "msa.h"
#include "tree.h"
#include "utilities.h"
#include "options.h"
#include "tools.h"

#define DEBUG

// Helper function that determines how many structures are completely allocated (as opposed to aborted by malloc failure) and free the allocation
void clean_up(int allocated, msa_t * msa, msa_t * msa1, msa_t * msa2, option_t * options){
    switch(allocated){
        case 4:
            destroy_msa(msa2);
        case 3:
            destroy_msa(msa1);
        case 2:
            destroy_msa(msa);
        case 1:
            destroy_options(options);
        default:
            return;
    }
}

#define ALLOCATED_INFO allocated, &msa, &msa1, &msa2, &options

// Main function
int main(int argc, char ** argv){
    option_t options;           //commnd line options
    msa_t msa, msa1, msa2;      // msa struct for the single HMM and 2 msa structs for the double HMM
    int allocated;              // heap allocation counter (to prevent mem leak)
    int left_root, right_root;  // endpoints of centroid decomposition
    float L, L1, L2;            // bit score for the single HMM and 2 HMMs for the double HMM 

    allocated =  0;

    // Read options
    if(init_options(&options)                   != SUCCESS)         PRINT_AND_EXIT("init option failed in main",                GENERAL_ERROR, ALLOCATED_INFO);
    if(read_cmd_arg(argc, argv, &options)       != SUCCESS)         PRINT_AND_EXIT("read command line args failed in main",     GENERAL_ERROR, ALLOCATED_INFO);
    allocated = 1;
    if(options.input_index                      == NULL_OPTION)     PRINT_AND_EXIT("must have valid input name",                GENERAL_ERROR, ALLOCATED_INFO);
    else{
        single_model_build_option.input_name            = options.input_name;
        fasttree_options.input_name                     = options.input_name;
        single_model_search_option.input_sequences_name = options.input_name;
    }


    // Allocate MSA and HMM structs on the stack
    if(parse_input(&msa, options.input_name)    != SUCCESS)         PRINT_AND_EXIT("init msa failed in main",                   GENERAL_ERROR, ALLOCATED_INFO);
    allocated = 2;

    // Call HMMbuild job for first model
    if(hmmbuild_job(&single_model_build_option) != SUCCESS)         PRINT_AND_EXIT("first model HMMbuild failed in main",       GENERAL_ERROR, ALLOCATED_INFO);


    // Call FastTree to build ML tree for the second model
    if(fasttree_job(&fasttree_options)          != SUCCESS)         PRINT_AND_EXIT("fast tree failed in main",                  GENERAL_ERROR, ALLOCATED_INFO);


    // Parse result of FastTree
    if(read_newick(fasttree_options.output_name)!= SUCCESS)         PRINT_AND_EXIT("read newick failed in main",                GENERAL_ERROR, ALLOCATED_INFO);

    for(int i = 0; i < 200; i++){
        printf("%d %s\n", parent_map[i], name_map[i]);
    }
    // Do centroid decomposition on the second model
    if(centroid_decomposition(&left_root, &right_root)
                                                != SUCCESS)         PRINT_AND_EXIT("centroid decomposition failed in main",     GENERAL_ERROR, ALLOCATED_INFO);
    
    if(make_smaller_msa(&msa, &msa1)            != SUCCESS)         PRINT_AND_EXIT("make small msa 1 failed in main",           GENERAL_ERROR, ALLOCATED_INFO);
    allocated = 3;
    if(make_smaller_msa(&msa, &msa2)            != SUCCESS)         PRINT_AND_EXIT("make small msa 2 failed in main",           GENERAL_ERROR, ALLOCATED_INFO);
    allocated = 4;
    if(retrieve_msa_from_root(left_root, right_root, &msa1, &msa2, &msa)
                                                != SUCCESS)         PRINT_AND_EXIT("retrieve_msa_from_root failed in main",     GENERAL_ERROR, ALLOCATED_INFO);

    // Write MSA to a file in FASTA format
    if(write_msa(&msa1, DEFAULT_DOUBLE_FIRST_MSA_NAME)
                                                != SUCCESS)         PRINT_AND_EXIT("write msa 1 failed in main",                GENERAL_ERROR, ALLOCATED_INFO);
    if(write_msa(&msa2, DEFAULT_DOUBLE_SECOND_MSA_NAME)
                                                != SUCCESS)         PRINT_AND_EXIT("write msa 2 failed in main",                GENERAL_ERROR, ALLOCATED_INFO);

    // Call HMMbuild for the second model
    if(hmmbuild_job(&double_model_first_build_option) 
                                                != SUCCESS)         PRINT_AND_EXIT("double model 1st hmmbuild failed in main",  GENERAL_ERROR, ALLOCATED_INFO);
    if(hmmbuild_job(&double_model_second_build_option)
                                                != SUCCESS)         PRINT_AND_EXIT("double model 2nd hmmbuild failed in main",  GENERAL_ERROR, ALLOCATED_INFO);

    // Compute Viterbi for first model
    if(hmmsearch_job(&single_model_search_option) 
                                                != SUCCESS)         PRINT_AND_EXIT("single model hmmsearch failed in main",     GENERAL_ERROR, ALLOCATED_INFO);

    L = compute_likelihood(DEFAULT_SINGLE_SEARCH_NAME, msa.num_seq);

    // Compute Viterbi for second model
    if(hmmsearch_job(&double_model_fist_search_option) 
                                                != SUCCESS)         PRINT_AND_EXIT("double model first hmmsearch failed in main",GENERAL_ERROR, ALLOCATED_INFO);
    if(hmmsearch_job(&double_model_second_search_option) 
                                                != SUCCESS)         PRINT_AND_EXIT("double model second hmmsearch failed in main",GENERAL_ERROR, ALLOCATED_INFO);
    L1 = compute_likelihood(DEFAULT_DOUBLE_FIRST_SEARCH_NAME,   msa1.num_seq);
    L2 = compute_likelihood(DEFAULT_DOUBLE_SECOND_SEARCH_NAME,  msa2.num_seq);

    // Perform statistical test 
    printf("%f %f %f\n", L, L1, L2);

    // Clean up

    return 0;
}