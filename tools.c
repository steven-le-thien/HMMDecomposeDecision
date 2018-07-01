#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "tools.h"
#include "utilities.h"

void add_command(char * current_command, char * new_command){
	str_add_space(current_command);
	strcat(current_command, new_command);
}

int hmmsearch_job(hmmsearch_options_t * hmm_options){
	char command[CMD_BUFFER_SIZE];

	strclr(command);
	strcat(command, "hmmsearch");
	add_command(command, hmm_options->no_ali_option);
	add_command(command, hmm_options->e_value_threshold);
	add_command(command, hmm_options->heuristics_filtering_threshold);
	add_command(command, hmm_options->output_name);
	add_command(command, hmm_options->input_hmm_name);
	add_command(command, hmm_options->input_sequences_name);

	// Call hmmbuild
	if(system(command) != SUCCESS)			PRINT_AND_RETURN("error in calling hmmsearch job\n", GENERAL_ERROR);

	return 0;
}

int hmmbuild_job(hmmbuild_option_t * hmm_options){
	char command[CMD_BUFFER_SIZE];

	// Build the command string
	strclr(command);
	strcat(command, "hmmbuild");
	add_command(command, hmm_options->molecule_name);
	add_command(command, hmm_options->symfrac);
	add_command(command, hmm_options->output_name);
	add_command(command, hmm_options->input_name);

	// Call hmmbuild
	if(system(command) != SUCCESS)			PRINT_AND_RETURN("error in calling hmmbuild job\n", GENERAL_ERROR);

	return 0;
}

int fasttree_job(fasttree_options_t * fasttree_options){
	char command[CMD_BUFFER_SIZE];

	// Build the command string
	strclr(command);
	strcat(command, "./FastTree");
	add_command(command, fasttree_options->model_name);
	add_command(command, fasttree_options->molecule_name);
	add_command(command, fasttree_options->support);
	add_command(command, fasttree_options->input_name);
	add_command(command, ">");
	add_command(command, fasttree_options->output_name);

	// Call FastTree
	if(system(command) != SUCCESS) 			PRINT_AND_RETURN("error in calling fasttree job\n", GENERAL_ERROR);
	return 0;
}

