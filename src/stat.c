// File in HMMDecompositionDecision, created by Thien Le in July 2018

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "stat.h"
#include "utilities.h"

#define CONST_E 2.7182818284590452353602874713

int bic(msa_t * single, msa_t * first_double, msa_t* second_double, float * L, float * L1, float * L2, char * hmmbuild_stdout, int * best_model){
    float prior_first;      //prior value for the first hmm
    float prior_second;     //prior value for the second hmm
    float first_model_log_odd;
    float second_model_log_odd;
    int i;
    FILE * f;
    char buffer[GENERAL_BUFFER_SIZE];       //placeholder for bit score (must be longer than the longest width of hmmsearch output)


    prior_first = 1.0 * first_double->num_seq / single->num_seq;
    prior_second = 1.0 - prior_first;

    first_model_log_odd = 0.0;
    second_model_log_odd = 0.0; 

    for(i = 0; i < single->num_seq; i++){
        first_model_log_odd += L[i];
    }

    for(i = 0; i < first_double->num_seq; i++){
        second_model_log_odd += L1[i] + log2f(prior_first);
    }

    for(i = 0; i < second_double->num_seq; i++){
        second_model_log_odd += L2[i] + log2f(prior_second);
    }

    f = fopen(hmmbuild_stdout, "r");
    if(!f) PRINT_AND_RETURN("open file error in bic calculation", OPEN_ERROR);

    for(i = 0; i < 13; i++){
        fgets(buffer, sizeof(buffer), f);
    }
    for(int i = 0; i < 4; i++){
        fscanf(f, "%s", buffer);
    }
    int len = atoi(buffer);

    fclose(f);
    // Copied from hmmer.h in HMMER, 
    /* Some notes:
 *   0. The model might be either in counts or probability form.
 *   1. t[0] is special: t[0][TMM,TMI,TMD] are the begin->M_1,I_0,D_1 entry probabilities,
 *      t[0][TIM,TII] are the I_0 transitions, and delete state 0 doesn't
 *      exist. Therefore D[0] transitions and mat[0] emissions are unused.
 *      To simplify some normalization code, we adopt a convention that these are set
 *      to valid probability distributions: 1.0 for t[0][TDM] and mat[0][0],
 *      and 0 for the rest.
 *   2. t[M] is also special: TMD and TDD are 0 because there is no next delete state;
 *      TDM is therefore 1.0 by definition. TMM and TDM are interpreted as the
 *      M->E and D->E end transitions. t[M][TDM] must be 1.0, therefore.
 */
    // Thus, intial state has 3 free transition parameters, for each of the next M - 1 states, we have 4 more free transition parameters and the end state does not have any free transition  parameter
    // Thus, total number of free transition parameter is 4 * M - 1

    // Since there are 4 nucleotides, there are 3 free emission parameters for each match state or insertion state (which we actually don't have)
    // Thus, total number of free parameter is 7 * M - 1
    // Thus the difference in number of free parameter is 7 * M (the double model has twice as many paramters, plus the `prior' parameter)
    int delta_k = 7 * len;

    float bic_double_over_single = log2f(single->num_seq) / log2f(CONST_E) * delta_k - 2.0 / log2f(CONST_E) * (second_model_log_odd - first_model_log_odd); 
    printf("Delta BIC is %f, log odd of double model is %f, log odd of single model is %f\n", bic_double_over_single, second_model_log_odd, first_model_log_odd);

    if(bic_double_over_single > 0) *best_model = 2;
    else *best_model = 1;
    return 0;
}