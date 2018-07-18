// File in HMMDecompositionDecision, created by Thien Le in July 2018

#ifndef OPTION_H
#define OPTION_H

#define NULL_OPTION             -1

#include <stdlib.h>

extern char DEFAULT_SINGLE_HMM_NAME             [];
extern char DEFAULT_SYMFRAC                     [];
extern char DEFAULT_HMM_MOLECULE                [];

extern char DEFAULT_TREE_OUTPUT                 [];
extern char DEFAULT_TREE_MODEL                  [];
extern char DEFAULT_TREE_MOLECULE               [];
extern char DEFAULT_SUPPORT                     [];

extern char DEFAULT_DOUBLE_FIRST_HMM_NAME       [];
extern char DEFAULT_DOUBLE_SECOND_HMM_NAME      [];
extern char DEFAULT_DOUBLE_FIRST_MSA_NAME       [];
extern char DEFAULT_DOUBLE_SECOND_MSA_NAME      [];

extern char DEFAULT_NOALI                       [];
extern char DEFAULT_HEURISTICS_FILTER           [];
extern char DEFAULT_E_VAL                       [];

extern char DEFAULT_SINGLE_SEARCH_NAME          [];
extern char DEFAULT_DOUBLE_FIRST_SEARCH_NAME    [];
extern char DEFAULT_DOUBLE_SECOND_SEARCH_NAME   [];

extern char DEFAULT_SINGLE_SEARCH_FLAG          [];
extern char DEFAULT_DOUBLE_FIRST_SEARCH_FLAG    [];
extern char DEFAULT_DOUBLE_SECOND_SEARCH_FLAG   [];

extern char DEFAULT_HMMBUILD_OUT_SINGLE_NAME    [];

extern char DEFAULT_HMMBUILD_OUT_SINGLE         [];
extern char DEFAULT_HMMBUILD_OUT_FIRST_DOUBLE   [];
extern char DEFAULT_HMMBUILD_OUT_SECOND_DOUBLE  [];

extern char DEFAULT_HMMSEARCH_OUT_SINGLE         [];
extern char DEFAULT_HMMSEARCH_OUT_FIRST_DOUBLE   [];
extern char DEFAULT_HMMSEARCH_OUT_SECOND_DOUBLE  [];


extern int DEFAULT_NUM_OPTIONS;

typedef struct options{
    int num_options;

    int input_index;
   char * input_name;

    int output_index;
   char * output_name;

    int symfrac_index;
   char * symfrac;
} option_t;

typedef struct hmm_options{
   char * input_name;
   char * symfrac;
   char * output_name;
   char * molecule_name;
   char * stdout;
} hmmbuild_option_t;

typedef struct hmm_options_2{
   char * input_sequences_name;
   char * input_hmm_name;
   char * output_name; 
   char * no_ali_option;
   char * e_value_threshold;
   char * heuristics_filtering_threshold;
   char * stdout;
} hmmsearch_options_t;

typedef struct tree_options{
   char * input_name;
   char * output_name;
   char * model_name;
   char * molecule_name;
   char * support;
} fasttree_options_t;

// Fields                                   
extern hmmbuild_option_t        single_model_build_option;
extern hmmbuild_option_t        double_model_first_build_option;
extern hmmbuild_option_t        double_model_second_build_option;

// Fields                                  
extern hmmsearch_options_t  single_model_search_option;
extern hmmsearch_options_t  double_model_fist_search_option;
extern hmmsearch_options_t  double_model_second_search_option;

// Fields                                 
extern fasttree_options_t fasttree_options;

// Functions
extern int read_cmd_arg(int argc,char ** argv, option_t * options);
extern int init_options(option_t * options);
extern void destroy_options(option_t * options);

#endif