// File in HMMDecompositionDecision, created by Thien Le in July 2018

#ifndef STAT_H
#define STAT_H

#include "msa.h"

extern int bic(msa_t * single, msa_t * first_double, msa_t* second_double, float * L, float * L1, float * L2, char * hmmbuild_stdout, int * best_model);
extern int aic(msa_t * single, msa_t * first_double, msa_t* second_double, float * L, float * L1, float * L2, char * hmmbuild_stdout, int * best_model);

#endif