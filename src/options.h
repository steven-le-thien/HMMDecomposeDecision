// File in HMMDecompositionDecision, created by Thien Le in July 2018

#ifndef OPTION_H
#define OPTION_H

#define NULL_OPTION             -1

#include <stdlib.h>

static char DEFAULT_SINGLE_HMM_NAME             []  = "defaultjob.single_hmm";
static char DEFAULT_SYMFRAC                     []  = "--symfrac=0.0";
static char DEFAULT_HMM_MOLECULE                []  = "--dna";

static char DEFAULT_TREE_OUTPUT                 []  = "defaultjob.fasttree.out";
static char DEFAULT_TREE_MODEL                  []  = "-gtr";
static char DEFAULT_TREE_MOLECULE               []  = "-nt";
static char DEFAULT_SUPPORT                     []  = "-nosupport";

static char DEFAULT_DOUBLE_FIRST_HMM_NAME       []  = "defaultjob.double_first_hmm";
static char DEFAULT_DOUBLE_SECOND_HMM_NAME      []  = "defaultjob.double_second_hmm";
static char DEFAULT_DOUBLE_FIRST_MSA_NAME       []  = "defaultjob.double_first_msa";
static char DEFAULT_DOUBLE_SECOND_MSA_NAME      []  = "defaultjob.double_second_msa";

static char DEFAULT_NOALI                       []  = "--noali";
static char DEFAULT_HEURISTICS_FILTER           []  = "--max";
static char DEFAULT_E_VAL                       []  = "-E Infinity";

static char DEFAULT_SINGLE_SEARCH_NAME          []  = "defaultjob.single_search_out";
static char DEFAULT_DOUBLE_FIRST_SEARCH_NAME    []  = "defaultjob.double_first_search_out";
static char DEFAULT_DOUBLE_SECOND_SEARCH_NAME   []  = "defaultjob.double_second_search_out";

static char DEFAULT_SINGLE_SEARCH_FLAG          []  = "--tblout defaultjob.single_search_out";
static char DEFAULT_DOUBLE_FIRST_SEARCH_FLAG    []  = "--tblout defaultjob.double_first_search_out";
static char DEFAULT_DOUBLE_SECOND_SEARCH_FLAG   []  = "--tblout defaultjob.double_second_search_out";

static char DEFAULT_HMMBUILD_OUT_SINGLE_NAME    []  = "defaultjob.hmmbuild_single.stdout";

static char DEFAULT_HMMBUILD_OUT_SINGLE         []  = "> defaultjob.hmmbuild_single.stdout";
static char DEFAULT_HMMBUILD_OUT_FIRST_DOUBLE   []  = "> defaultjob.hmmbuild_first_double.stdout";
static char DEFAULT_HMMBUILD_OUT_SECOND_DOUBLE  []  = "> defaultjob.hmmbuild_second_double.stdout";

static char DEFAULT_HMMSEARCH_OUT_SINGLE         []  = "> defaultjob.hmmsearch_single.stdout";
static char DEFAULT_HMMSEARCH_OUT_FIRST_DOUBLE   []  = "> defaultjob.hmmsearch_first_double.stdout";
static char DEFAULT_HMMSEARCH_OUT_SECOND_DOUBLE  []  = "> defaultjob.hmmseacrh_second_double.stdout";


static const char DEFAULT_NUM_OPTIONS = 3;

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

// Fields                                                               input_name                      symfrac             output_name                     molecule_name                       
static hmmbuild_option_t        single_model_build_option           = { NULL,                           DEFAULT_SYMFRAC,    DEFAULT_SINGLE_HMM_NAME,        DEFAULT_HMM_MOLECULE,   DEFAULT_HMMBUILD_OUT_SINGLE};
static hmmbuild_option_t        double_model_first_build_option     = { DEFAULT_DOUBLE_FIRST_MSA_NAME,  DEFAULT_SYMFRAC,    DEFAULT_DOUBLE_FIRST_HMM_NAME,  DEFAULT_HMM_MOLECULE,   DEFAULT_HMMBUILD_OUT_FIRST_DOUBLE};
static hmmbuild_option_t        double_model_second_build_option    = { DEFAULT_DOUBLE_SECOND_MSA_NAME, DEFAULT_SYMFRAC,    DEFAULT_DOUBLE_SECOND_HMM_NAME, DEFAULT_HMM_MOLECULE,   DEFAULT_HMMBUILD_OUT_SECOND_DOUBLE};

// Fields                                                           input_sequences_name            input_hmm_name                  output_name                         no_ali_option       e_value_threshold       heuristics_filtering_threshold      
static hmmsearch_options_t  single_model_search_option          = { NULL,                           DEFAULT_SINGLE_HMM_NAME,        DEFAULT_SINGLE_SEARCH_FLAG,         DEFAULT_NOALI,      DEFAULT_E_VAL,          DEFAULT_HEURISTICS_FILTER,              DEFAULT_HMMSEARCH_OUT_SINGLE}; 
static hmmsearch_options_t  double_model_fist_search_option     = { DEFAULT_DOUBLE_FIRST_MSA_NAME,  DEFAULT_DOUBLE_FIRST_HMM_NAME,  DEFAULT_DOUBLE_FIRST_SEARCH_FLAG,   DEFAULT_NOALI,      DEFAULT_E_VAL,          DEFAULT_HEURISTICS_FILTER,              DEFAULT_HMMSEARCH_OUT_FIRST_DOUBLE}; 
static hmmsearch_options_t  double_model_second_search_option   = { DEFAULT_DOUBLE_SECOND_MSA_NAME, DEFAULT_DOUBLE_SECOND_HMM_NAME, DEFAULT_DOUBLE_SECOND_SEARCH_FLAG,  DEFAULT_NOALI,      DEFAULT_E_VAL,          DEFAULT_HEURISTICS_FILTER,              DEFAULT_HMMSEARCH_OUT_SECOND_DOUBLE}; 

// Fields                                       input_name      output_name             model_name              molecule_name               support
static fasttree_options_t fasttree_options  = { NULL,           DEFAULT_TREE_OUTPUT,    DEFAULT_TREE_MODEL,     DEFAULT_TREE_MOLECULE,      DEFAULT_SUPPORT};

// Functions
extern int read_cmd_arg(int argc, char ** argv, option_t * options);
extern int init_options(option_t * options);
extern void destroy_options(option_t * options);

#endif