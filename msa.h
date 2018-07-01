#ifndef MSA_H
#define MSA_H

// Constants
const static int MAX_SEQUENCE_LENGTH 	= (int) 1e6;
const static int MAX_NUM_SEQUENCE 		= (int) 1e6;
const static int MAX_NAME_LENGTH 		= (int) 1e3;

// Structure for the multiple sequence alignment
typedef struct msa {
	int 	N; 			// size of one sequence
	int 	num_seq; 	// number of sequences
	char** 	msa;
	char** 	name;
} msa_t;

// Public functions. Details are in definition

// Constructor & destructor
extern int init_msa(msa_t* msa, int N, int num_seq, char** msa_core,  char** msa_name);
extern void destroy_msa(msa_t * msa);

// IO functions
extern int parse_input(msa_t * msa, char * filename);
extern int write_msa(msa_t * msa, char * filename);

// Extending to a new MSA
extern int add_to_msa_from_msa(msa_t * original, msa_t * new, char * name);
extern int make_smaller_msa(msa_t * original, msa_t * new);

// With tree 
extern void dfs_msa(int cur_node, int par_node, int centroid, msa_t * all_msa, msa_t * small_msa);
extern int retrieve_msa_from_root(int centroid1, int centroid2, msa_t * msa1, msa_t * msa2, msa_t * all_msa);


extern float compute_likelihood(char * filename, msa_t * msa);

#endif
