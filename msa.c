#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "msa.h"
#include "utilities.h"
#include "tree.h"

/* Initialize fields in the MSA to predetermined values
 * Input: 	msa 		pointer to the MSA srtuct
 * 			N 			length of each sequence
 * 			num_seq 	number of sequences in the MSA
 * 			msa_core 	array of sequences 
 * 			msa_name 	array of name for each sequence
 * Output: 	0 on success, ERROR otherwise
 * Effect 	set fields in the input msa
 */
int init_msa(msa_t* msa, int N, int num_seq, char** msa_core, char** msa_name){
	// Safety check
	if(!msa) 		PRINT_AND_RETURN("msa is NULL from init_msa", 		GENERAL_ERROR);
	if(!msa_core)	PRINT_AND_RETURN("msa_core is NULL from init_msa", 	GENERAL_ERROR);
	if(!msa_name) 	PRINT_AND_RETURN("msa_name is NULL from init_msa", 	GENERAL_ERROR);

	// Set the fields
	msa->N = N;
	msa->num_seq = num_seq;
	msa->msa = msa_core;
	msa->name = msa_name;
	return 0;
}

/* Wrapper to malloc name holder and copy name into that place
 * Input: 	the name and the location to be copied to
 * Output: 	0 on success, ERROR otherwise	
 * Effect: 	calls malloc
 */
int setup_name(char ** name_holder, char* name){
	if(!name) 	PRINT_AND_RETURN("name is NULL when trying to parse", 	GENERAL_ERROR);

	name_holder[0] = (char *) malloc(strlen(name));
	if(!name_holder[0]) PRINT_AND_RETURN("malloc failed for name holder", MALLOC_ERROR);
	strcpy(name_holder[0], name);
	return 0;
}

/* Wrapper to malloc sequence holder and copy sequence into that place
 * Input: 	the sequence and the location to be copied to
 * Output: 	0 on success, ERROR otherwise	
 * Effect: 	calls malloc
 */
int setup_sequence(char ** sequence_holder, char* sequence){
	// The syntax is exactly the same as setup_name so we will just call that function instead
	return setup_name(sequence_holder, sequence);
}

/* Parse FASTA file into msa_t struct
 * Input: 	pointer to the msa and name of the FASTA file
 * Output: 	0 on sucess, GENERAL_ERROR otherwise
 * Effect: 	may print onto outstream if error occurs, set fields in the msa struct, 
 *			open instream, assuming filename is in the same directory as the binary file
 * 			calls malloc
 */
int parse_input(msa_t * msa_ptr, char * filename){
	// May switch to dynamically allocation when dealing with longer sequences (e.g. concatenations)
	char 	line	[MAX_SEQUENCE_LENGTH];
	char 	seq 	[MAX_SEQUENCE_LENGTH];
	strclr(seq);
	strclr(line);

	char** msa_core = malloc(MAX_NUM_SEQUENCE * sizeof(char*));
	char** msa_name = malloc(MAX_NUM_SEQUENCE * sizeof(char*)); 
	if(!msa_core) 	PRINT_AND_RETURN("malloc for msa_core failed in parse_input",	MALLOC_ERROR);
	if(!msa_name) 	PRINT_AND_RETURN("malloc for msa_name failed in parse_input", 	MALLOC_ERROR);

	// Redirecting stdin
#ifdef DEBUG
	printf("Trying to redirect standard input to %s\n", filename);
#endif
	FILE * f = freopen(filename, "r", stdin);
	if(!f) PRINT_AND_RETURN("fail to open input file",	OPEN_ERROR);

	// FASTA parsing section
	int sequence_counter = 0;

	// Ignoring comments
	while(1){
		if(scanf("%s", line) < 0) PRINT_AND_RETURN("input file contains no sequence",	GENERAL_ERROR);
		if(str_start_with(line, '>')) break; 
	}

	// Set up name for the first sequence
	setup_name(&msa_name[sequence_counter], &line[1]);

	// Read the stdin line by line until EOF signal
	while(scanf("%s", line) >= 0){
		switch(line[0]){
			case ';': break;
			case '>': // Finished previous seuqence
				setup_sequence(&msa_core[sequence_counter], seq);
				sequence_counter++;
				setup_name(&msa_name[sequence_counter], &line[1]);
				strclr(seq);
				break;
			default:
				strcat(seq, line);
		}
	}

	// Final iteration for the last sequence
	setup_sequence(&msa_core[sequence_counter], seq);
	sequence_counter++;

	// Close file
	fclose(f);

	// Fill in the fields
	return(init_msa(msa_ptr, strlen(seq), sequence_counter, msa_core, msa_name));
}


/* Extending an original msa to a blank msa with same meta
 * Input: 	pointer to the two msa
 * Output: 	0 on sucess, ERROR otherwise
 * Effect: 	calls malloc, set fields in the new msa
 */
int make_smaller_msa(msa_t * original, msa_t * new){
	new->N 			= original->N;
	new->num_seq 	= 0;

	new->msa 		= (char **) malloc(original->num_seq * sizeof(char *));
	new->name 		= (char **) malloc(original->num_seq * sizeof(char *));
	if(!new->msa) 	PRINT_AND_RETURN("malloc for new->msa failed in make_smaller_msa",		MALLOC_ERROR);
	if(!new->name) 	PRINT_AND_RETURN("malloc for new->name failed in make_smaller_msa", 	MALLOC_ERROR);

	for(int i = 0; i < original->num_seq; i++){
		new->msa[i] 		= (char *) malloc(new->N);
		new->name[i]	 	= (char *) malloc(MAX_NAME_LENGTH);
		if(!new->msa[i]) 	PRINT_AND_RETURN("malloc for new msa sequence failed in make_smaller_msa",	MALLOC_ERROR);
		if(!new->msa[i]) 	PRINT_AND_RETURN("malloc for new msa name failed in make_smaller_msa", 		MALLOC_ERROR);

		new->msa[i][0] 		= 0;
		new->name[i][0] 	= 0;
	}
	return 0;
}

/* Brute force algorithm to find sequence based on name
 * Input: 	msa structure and the name of the sequence
 * Output: 	pointer to the sequence on success, NULL otherwise
 * Effect: 	none
 */
char * find_sequence_by_name(msa_t * msa, char * name){
	if(!msa) 		PRINT_AND_RETURN("msa is NULL in find_sequence_by_name",	NULL);
	if(!name) 		PRINT_AND_RETURN("name is NULL in find_sequence_by_name",	NULL);
	for(int i = 0; i < msa->num_seq; i++)
		if(strcmp(name, msa->name[i]) == 0)
			return msa->msa[i];
	return NULL;
}

/* Brute force algorithm to add a sequence from an msa struct to another msa struct
 * Input: 	both msa structures and the name of the sequence
 * Output: 	0 on success, ERROR otherwise
 * Effect: 	set fields in the new msa
 */
int add_to_msa_from_msa(msa_t * original, msa_t * new, char * name){
	if(!original) 	PRINT_AND_RETURN("original is NULL in add_to_msa_from_msa",	GENERAL_ERROR);
	if(!new) 		PRINT_AND_RETURN("new is NULL in add_to_msa_from_msa",		GENERAL_ERROR);
	if(!name) 		PRINT_AND_RETURN("name is NULL in add_to_msa_from_msa",		GENERAL_ERROR);

	char* sequence = find_sequence_by_name(original, name);
	if(!sequence) 	PRINT_AND_RETURN("cannot find sequence with said name in add_to_msa_from_msa",		GENERAL_ERROR);

	strcpy(new->msa[new->num_seq], sequence);
	strcpy(new->name[new->num_seq], name);
	new->num_seq++;
	return 0;
}

/* Open a file and write msa struct in FASTA format
 * Input: 	msa structure and the filename
 * Output: 	0 on success, ERROR otherwise
 * Effect: 	set fields in the new msa
 */
int write_msa(msa_t * msa, char * filename){
	if(!msa) 		PRINT_AND_RETURN("msa is NULL in write_msa",				GENERAL_ERROR);

	// Try to open file
	FILE *f = fopen(filename, "w");
	if (!f) 		PRINT_AND_RETURN("cannot open file to write in write_msa", 	GENERAL_ERROR);

	// Write to file
	for(int i = 0; i < msa->num_seq; i++){
		fprintf(f, ">%s\n", msa->name[i]);
		fprintf(f, "%s\n", msa->msa[i]);
	}

	// Close file
	fclose(f);
	return 0;
}


void dfs_msa(int cur_node, int par_node, int centroid, msa_t * all_msa, msa_t * small_msa){
	if(is_leave(cur_node)) add_to_msa_from_msa(all_msa, small_msa, name_map[cur_node]);
	for(int i = 0; i < maxN; i++){
		if(!adj_mat[cur_node][i] || i == par_node || i == cur_node || i == centroid) continue;

		dfs_msa(i, cur_node, centroid, all_msa, small_msa);
	}
}

int retrieve_msa_from_root(int centroid1, int centroid2, msa_t * msa1, msa_t * msa2, msa_t * all_msa){
	dfs_msa(centroid1, parent_map[centroid1], centroid2, all_msa, msa1);
	dfs_msa(centroid2, parent_map[centroid2], centroid1, all_msa, msa2);
	return 0;
}

float compute_likelihood(char * filename, msa_t * msa){
	// Open file
	FILE *f = fopen(filename, "r");
	if (!f) 								PRINT_AND_RETURN("cannot open file to write in hmmsearch_output", 	GENERAL_ERROR);

	// Get to the 4th line of the hmmsearch output file 
	char buffer[1000];
	fgets(buffer, sizeof(buffer), stdin);
	fgets(buffer, sizeof(buffer), stdin);

	float result = 0.0;
	for(int i = 0; i < msa->num_seq; i++){
		fgets(buffer, sizeof(buffer), stdin);
		scanf("%s", buffer);			
		scanf("%s", buffer);			
		scanf("%s", buffer);			
		scanf("%s", buffer);
		scanf("%s", buffer);
		result += atof(buffer);
	}

	fclose(f);

	return result;
}

