#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "msa.h"
#include "tree.h"
#include "utilities.h"
#include "options.h"
#include "tools.h"

#define DEBUG

int main(int argc, char ** argv){
	// Allocate commnd line options
	option_t options;
	if(init_options(&options) 					!= SUCCESS) 		PRINT_AND_RETURN("init option failed in main", 				GENERAL_ERROR);
	if(read_cmd_arg(argc, argv, &options) 		!= SUCCESS) 		PRINT_AND_RETURN("read command line args failed in main", 	GENERAL_ERROR);

	if(options.input_index 						== NULL_OPTION) 	PRINT_AND_RETURN("must have valid input name",			 	GENERAL_ERROR);
	else{
		single_model_build_option.input_name 	= options.input_name;
		fasttree_options.input_name 	= options.input_name;
		single_model_search_option.input_sequences_name = options.input_name;
	}

	// Allocate MSA and HMM structs on the stack
	msa_t msa;
	if(parse_input(&msa, options.input_name) 	!= SUCCESS) 		PRINT_AND_RETURN("init msa failed in main",				 	GENERAL_ERROR);

	// Call HMMbuild job for first model
	if(hmmbuild_job(&single_model_build_option) != SUCCESS) 		PRINT_AND_RETURN("first model HMMbuild failed in main",	 	GENERAL_ERROR);

	// Call FastTree to build ML tree for the second model
	if(fasttree_job(&fasttree_options) 			!= SUCCESS) 		PRINT_AND_RETURN("fast tree failed in main",	 			GENERAL_ERROR);

	// Parse result of FastTree
	if(read_newick(fasttree_options.output_name)!= SUCCESS) 		PRINT_AND_RETURN("read newick failed in main",	 			GENERAL_ERROR);

	// Do centroid decomposition on the second model
	int left_root, right_root;
	if(centroid_decomposition(&left_root, &right_root)
												!= SUCCESS) 		PRINT_AND_RETURN("centroid decomposition failed in main",	GENERAL_ERROR);

	// Allocate space for the 2 sub msa
	msa_t msa1, msa2;
	if(make_smaller_msa(&msa, &msa1) 			!= SUCCESS) 		PRINT_AND_RETURN("make small msa 1 failed in main",	 		GENERAL_ERROR);
	if(make_smaller_msa(&msa, &msa2) 			!= SUCCESS) 		PRINT_AND_RETURN("make small msa 2 failed in main",	 		GENERAL_ERROR);
	if(retrieve_msa_from_root(left_root, right_root, &msa1, &msa2, &msa, name_map, parent_map)
												!= SUCCESS)			PRINT_AND_RETURN("retrieve_msa_from_root failed in main",	GENERAL_ERROR);
	printf("testing msa is parsed correctly\n");
	printf("num seq is %d, N is %d\n", msa1.num_seq, msa1.N);
	for(int i = 0; i < msa1.num_seq; i++){
		printf("%s ", msa1.name[i]);
	}
	printf("\n");
	for(int i = 0; i < msa1.num_seq; i++){
		for(int j = 0; j < msa1.N; j++){
			printf("%c", msa1.msa[i][j]);
		}
		printf("\n");
	}

	

	// Write MSA to a file in FASTA format
	if(write_msa(&msa1, DEFAULT_DOUBLE_FIRST_MSA_NAME)
												!= SUCCESS)			PRINT_AND_RETURN("write msa 1 failed in main",	 			GENERAL_ERROR);
	if(write_msa(&msa2, DEFAULT_DOUBLE_SECOND_MSA_NAME)
												!= SUCCESS)			PRINT_AND_RETURN("write msa 2 failed in main",	 			GENERAL_ERROR);

	// Call HMMbuild for the second model
	if(hmmbuild_job(&double_model_first_build_option) != SUCCESS) 	PRINT_AND_RETURN("double model 1st hmmbuild failed in main",GENERAL_ERROR);
	if(hmmbuild_job(&double_model_second_build_option)!= SUCCESS) 	PRINT_AND_RETURN("double model 2nd hmmbuild failed in main",GENERAL_ERROR);

	// Compute Viterbi for first model
	if(hmmsearch_job(&single_model_search_option) 
												!= SUCCESS)			PRINT_AND_RETURN("single model hmmsearch failed in main",	GENERAL_ERROR);

	float L = compute_likelihood(DEFAULT_SINGLE_SEARCH_NAME, &msa);

	// Compute Viterbi for second model
	if(hmmsearch_job(&double_model_fist_search_option) 
												!= SUCCESS)			PRINT_AND_RETURN("double model first hmmsearch failed in main",	GENERAL_ERROR);
	if(hmmsearch_job(&double_model_second_search_option) 
												!= SUCCESS)			PRINT_AND_RETURN("double model second hmmsearch failed in main",GENERAL_ERROR);
	float L1 = compute_likelihood(DEFAULT_DOUBLE_FIRST_SEARCH_NAME, &msa1);
	float L2 = compute_likelihood(DEFAULT_DOUBLE_SECOND_SEARCH_NAME, &msa2);

	// Perform statistical test 
	printf("%f %f %f\n", L, L1, L2);
	return 0;
}