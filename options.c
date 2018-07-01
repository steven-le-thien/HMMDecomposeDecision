#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "options.h"
#include "utilities.h"

int find_arg_index(char * flag, char * content, option_t * options, int i){
	if(strcmp(flag, "-i") == 0){//reading input name
		options->input_index = i;
		options->input_name = malloc(strlen(content));
		strcpy(options->input_name, content);
	} else if(strcmp(flag, "-o") == 0){
		options->output_index = i;
		options->input_name = malloc(strlen(content));
		strcpy(options->output_name, content);
	} else if(strcmp(flag, "--symfrac") == 0){
		options->symfrac_index = i;
		options->input_name = malloc(strlen(content));
		strcpy(options->symfrac, content);
	} else PRINT_AND_RETURN("unrecognized argument", GENERAL_ERROR); 

	return 0;
}

int read_cmd_arg(int argc, char ** argv, option_t * options){
	// Make sure argc is in the correct range
	if(!argc % 2 || !(1 < argc && argc < 2 * options->num_options + 2)) 
		PRINT_AND_RETURN("incorrect number of argument", GENERAL_ERROR);
	for(int i = 0; i < argc; i++)
		if(i != 0 && i % 2) //currently reading a flag
			if(find_arg_index(argv[i], argv[i + 1], options, i) != SUCCESS)
				PRINT_AND_RETURN("failed reading of the arguments", GENERAL_ERROR);

	return 0;
}


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

void destroy_options(option_t * options){
	if(options->input_name) 	free(options->input_name);
	if(options->output_name)	free(options->output_name);
	if(options->symfrac)		free(options->symfrac);

	options->input_index = options->output_index = options->symfrac_index = 0;
}