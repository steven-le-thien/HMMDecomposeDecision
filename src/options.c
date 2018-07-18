// File in HMMDecompositionDecision, created by Thien Le in July 2018

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "options.h"
#include "utilities.h"


// Constants
int DEFAULT_NUM_OPTIONS = 3;

char DEFAULT_SINGLE_HMM_NAME             []  = "defaultjob.single_hmm";
char DEFAULT_SYMFRAC                     []  = "--symfrac=0.0";
char DEFAULT_HMM_MOLECULE                []  = "--dna";

char DEFAULT_TREE_OUTPUT                 []  = "defaultjob.fasttree.out";
char DEFAULT_TREE_MODEL                  []  = "-gtr";
char DEFAULT_TREE_MOLECULE               []  = "-nt";
char DEFAULT_SUPPORT                     []  = "-nosupport";

char DEFAULT_DOUBLE_FIRST_HMM_NAME       []  = "defaultjob.double_first_hmm";
char DEFAULT_DOUBLE_SECOND_HMM_NAME      []  = "defaultjob.double_second_hmm";
char DEFAULT_DOUBLE_FIRST_MSA_NAME       []  = "defaultjob.double_first_msa";
char DEFAULT_DOUBLE_SECOND_MSA_NAME      []  = "defaultjob.double_second_msa";

char DEFAULT_NOALI                       []  = "--noali";
char DEFAULT_HEURISTICS_FILTER           []  = "--max";
char DEFAULT_E_VAL                       []  = "-E Infinity";

char DEFAULT_SINGLE_SEARCH_NAME          []  = "defaultjob.single_search_out";
char DEFAULT_DOUBLE_FIRST_SEARCH_NAME    []  = "defaultjob.double_first_search_out";
char DEFAULT_DOUBLE_SECOND_SEARCH_NAME   []  = "defaultjob.double_second_search_out";

char DEFAULT_SINGLE_SEARCH_FLAG          []  = "--tblout defaultjob.single_search_out";
char DEFAULT_DOUBLE_FIRST_SEARCH_FLAG    []  = "--tblout defaultjob.double_first_search_out";
char DEFAULT_DOUBLE_SECOND_SEARCH_FLAG   []  = "--tblout defaultjob.double_second_search_out";

char DEFAULT_HMMBUILD_OUT_SINGLE_NAME    []  = "defaultjob.hmmbuild_single.stdout";

char DEFAULT_HMMBUILD_OUT_SINGLE         []  = "defaultjob.hmmbuild_single.stdout";
char DEFAULT_HMMBUILD_OUT_FIRST_DOUBLE   []  = "defaultjob.hmmbuild_first_double.stdout";
char DEFAULT_HMMBUILD_OUT_SECOND_DOUBLE  []  = "defaultjob.hmmbuild_second_double.stdout";

char DEFAULT_HMMSEARCH_OUT_SINGLE         []  = "defaultjob.hmmsearch_single.stdout";
char DEFAULT_HMMSEARCH_OUT_FIRST_DOUBLE   []  = "defaultjob.hmmsearch_first_double.stdout";
char DEFAULT_HMMSEARCH_OUT_SECOND_DOUBLE  []  = "defaultjob.hmmsearch_second_double.stdout";

// Fields                                                               input_name                      symfrac             output_name                     molecule_name                       
hmmbuild_option_t        single_model_build_option           = { NULL,                           DEFAULT_SYMFRAC,    DEFAULT_SINGLE_HMM_NAME,        DEFAULT_HMM_MOLECULE,   DEFAULT_HMMBUILD_OUT_SINGLE};
hmmbuild_option_t        double_model_first_build_option     = { DEFAULT_DOUBLE_FIRST_MSA_NAME,  DEFAULT_SYMFRAC,    DEFAULT_DOUBLE_FIRST_HMM_NAME,  DEFAULT_HMM_MOLECULE,   DEFAULT_HMMBUILD_OUT_FIRST_DOUBLE};
hmmbuild_option_t        double_model_second_build_option    = { DEFAULT_DOUBLE_SECOND_MSA_NAME, DEFAULT_SYMFRAC,    DEFAULT_DOUBLE_SECOND_HMM_NAME, DEFAULT_HMM_MOLECULE,   DEFAULT_HMMBUILD_OUT_SECOND_DOUBLE};

// Fields                                                           input_sequences_name            input_hmm_name                  output_name                         no_ali_option       e_value_threshold       heuristics_filtering_threshold      
hmmsearch_options_t  single_model_search_option          = { NULL,                           DEFAULT_SINGLE_HMM_NAME,        DEFAULT_SINGLE_SEARCH_FLAG,         DEFAULT_NOALI,      DEFAULT_E_VAL,          DEFAULT_HEURISTICS_FILTER,              DEFAULT_HMMSEARCH_OUT_SINGLE}; 
hmmsearch_options_t  double_model_fist_search_option     = { DEFAULT_DOUBLE_FIRST_MSA_NAME,  DEFAULT_DOUBLE_FIRST_HMM_NAME,  DEFAULT_DOUBLE_FIRST_SEARCH_FLAG,   DEFAULT_NOALI,      DEFAULT_E_VAL,          DEFAULT_HEURISTICS_FILTER,              DEFAULT_HMMSEARCH_OUT_FIRST_DOUBLE}; 
hmmsearch_options_t  double_model_second_search_option   = { DEFAULT_DOUBLE_SECOND_MSA_NAME, DEFAULT_DOUBLE_SECOND_HMM_NAME, DEFAULT_DOUBLE_SECOND_SEARCH_FLAG,  DEFAULT_NOALI,      DEFAULT_E_VAL,          DEFAULT_HEURISTICS_FILTER,              DEFAULT_HMMSEARCH_OUT_SECOND_DOUBLE}; 

// Fields                                       input_name      output_name             model_name              molecule_name               support
fasttree_options_t fasttree_options  = { NULL,           DEFAULT_TREE_OUTPUT,    DEFAULT_TREE_MODEL,     DEFAULT_TREE_MOLECULE,      DEFAULT_SUPPORT};

/* Function to check the flag tag and assign its content to the appropriate field in option structure
 * This function also stores the index in argv to the content
 * Input    the flag, the content and pointer to an option struct as well as index to argv
 * Output   0 on success, ERROR otherwise
 * Effect   set fields in option struct
 */
int find_arg_index(char * flag, char * content, option_t * options, int i){
    if(strcmp(flag, "-i") == 0){//reading input name
        options->input_index = i;
        options->input_name = malloc(strlen(content));  

        if(!options->input_name)        
            PRINT_AND_RETURN("malloc failure for input name in find_arg_index",     MALLOC_ERROR);
        else
            strcpy(options->input_name, content);

    } else if(strcmp(flag, "-o") == 0){
        options->output_index = i;
        options->output_name = malloc(strlen(content));

        if(!options->output_name)       
            PRINT_AND_RETURN("malloc failure for output name in find_arg_index",    MALLOC_ERROR);
        else
            strcpy(options->output_name, content);

    } else if(strcmp(flag, "--symfrac") == 0){
        options->symfrac_index = i;
        options->symfrac = malloc(strlen(content));

        if(!options->symfrac)           
            PRINT_AND_RETURN("malloc failure for symfrac in find_arg_index",        MALLOC_ERROR);
        else
            strcpy(options->symfrac, content);
    } else PRINT_AND_RETURN("unrecognized argument", GENERAL_ERROR); 

    return 0;
}

/* Wrapper function to check for command length and loop through command to parse options
 * Input    argc, argv from command line and an option struct pointer
 * Output   0 on success, ERROR otherwise
 * Effect   set fields in option struct
 */
int read_cmd_arg(int argc, char ** argv, option_t * options){
    if(!options) 
        PRINT_AND_RETURN("options is null in read_cmd_arg", GENERAL_ERROR);

    // Make sure argc is in the correct range
    if(!argc % 2 || !(1 < argc && argc < 2 * options->num_options + 2)) 
        PRINT_AND_RETURN("incorrect number of argument", GENERAL_ERROR);

    // Loop through command
    for(int i = 0; i < argc; i++)
        if(i != 0 && i % 2) //currently reading a flag
            if(find_arg_index(argv[i], argv[i + 1], options, i) != SUCCESS)
                PRINT_AND_RETURN("failed reading of the arguments", GENERAL_ERROR);

    return 0;
}

/* Constructor for option struct
 * Input    pointer to option struct
 * Output   0 on success
 * Effect   set fields in option struct to default value
 */
int init_options(option_t * options){
    options->num_options = DEFAULT_NUM_OPTIONS;

    options->input_index = -1;
    options->output_index = -1;
    options->symfrac_index = -1;

    options->input_name = NULL;
    options->output_name = NULL;
    options->symfrac = NULL;

    return 0;
}

/* Destructor for option struct
 * Input    pointer to option struct
 * Output   0 on success, ERROR otherwise
 * Effect   set fields in option struct
 */
void destroy_options(option_t * options){
    // If the struct passed in is not allocated then nothing happens
    if(!options) return;
    if(options->input_name)     free(options->input_name);
    if(options->output_name)    free(options->output_name);
    if(options->symfrac)        free(options->symfrac);

    options->input_index = options->output_index = options->symfrac_index = 0;
}