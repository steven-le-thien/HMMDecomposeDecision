#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "msa.h"
#include "utilities.h"
#include "tree.h"

// Private function templates
int setup_name(char ** name_holder, char* name);
int setup_sequence(char ** sequence_holder, char* sequence);
char * find_sequence_by_name(msa_t * msa, char * name);
int add_to_msa_from_msa(msa_t * original, msa_t * new, char * name);
void dfs_msa(int cur_node, int par_node, int centroid, msa_t * all_msa, msa_t * small_msa);

/* Initialize fields in the MSA to predetermined values
 * Input:   msa         pointer to the MSA srtuct
 *          N           length of each sequence
 *          num_seq     number of sequences in the MSA
 *          msa_core    array of sequences 
 *          msa_name    array of name for each sequence
 * Output:  0 on success, ERROR otherwise
 * Effect   set fields in the input msa
 */
int init_msa(msa_t* msa, int N, int num_seq, char** msa_core, char** msa_name){
    // Safety check
    if(!msa)        PRINT_AND_RETURN("msa is NULL from init_msa",       GENERAL_ERROR);
    if(!msa_core)   PRINT_AND_RETURN("msa_core is NULL from init_msa",  GENERAL_ERROR);
    if(!msa_name)   PRINT_AND_RETURN("msa_name is NULL from init_msa",  GENERAL_ERROR);

    // Set the fields
    msa->N = N;
    msa->num_seq = num_seq;
    msa->msa = msa_core;
    msa->name = msa_name;
    return 0;
}

/* A destructor for the struct msa that loops through all mallocated fields and set other fields to 0
 * This does not deallocate the struct itself if the struct was from heap memory
 * Input:       msa struct
 * Output:      nothing
 * Effect:      freeing mallocated blocks
 */
void destroy_msa(msa_t * msa){
    int i; //loop variable

    // Safe cheking: if the structure was never allocated then does nothing
    if(!msa) return;

    // Free the names and the sequences containing in the msa
    for (i = 0; i < msa->num_seq; i++){
        free(msa->msa[i]);
        free(msa->name[i]);
    }

    // Set other fields to 0
    free(msa->msa);
    free(msa->name);
    msa->N = msa->num_seq = 0;
}


/* Parse FASTA file into msa_t struct. This function assumes a maximum length of a sequence and number of sequences currently
 * Input:   pointer to the msa and name of the FASTA file
 * Output:  0 on sucess, GENERAL_ERROR otherwise
 * Effect:  may print onto outstream if error occurs, set fields in the msa struct, 
 *          open instream, assuming filename is in the same directory as the binary file
 *          calls malloc
 */
int parse_input(msa_t * msa_ptr, char * filename){
    // Placeholder for contents of the msa
    char** msa_core;
    char** msa_name;

    int i; //loop variable

    // Input FASTA file
    FILE * f;

    // File reader variables
    int sequence_counter;
    char line[MAX_SEQUENCE_LENGTH];
    char seq[MAX_SEQUENCE_LENGTH];

    // Safe checking
    if(!msa_ptr)    PRINT_AND_RETURN("msa_ptr is NULL in parse input",          GENERAL_ERROR);
    if(!filename)   PRINT_AND_RETURN("filename is NULL in parse input",         GENERAL_ERROR);

    // Allocate space
    msa_core = malloc(MAX_NUM_SEQUENCE * sizeof(char*));
    msa_name = malloc(MAX_NUM_SEQUENCE * sizeof(char*)); 
    if(!msa_core)   PRINT_AND_RETURN("malloc for msa_core failed in parse_input",   MALLOC_ERROR);
    if(!msa_name){
        free(msa_core);
        PRINT_AND_RETURN("malloc for msa_name failed in parse_input",   MALLOC_ERROR);
    }

    // Clear strings    
    sequence_counter = 0;
    strclr(seq);
    strclr(line);

    // Redirecting stdin
    f = freopen(filename, "r", stdin);
    if(!f) PRINT_AND_RETURN("fail to open input file",  OPEN_ERROR);    

    // Ignoring comments
    while(1){
        if(scanf("%s", line) < 0) PRINT_AND_RETURN("input file contains no sequence",   GENERAL_ERROR);
        if(str_start_with(line, '>')) break; 
    }

    // Set up name for the first sequence
    if(setup_name(&msa_name[sequence_counter], &line[1]) != SUCCESS){ //roll back and deallocate
        free(msa_core);
        free(msa_name);
    }

    // Read the stdin line by line until EOF signal
    while(scanf("%s", line) >= 0){
        switch(*line){ //read the first character of the word
            case ';': break;
            case '>': // Finished previous seuqence
                if(setup_sequence(&msa_core[sequence_counter], seq) != SUCCESS){ // roll back and deallocate
                    for(i = 0; i < sequence_counter; i++){
                        if(i != sequence_counter - 1) free(msa_core[i]);
                        free(msa_name[i]);
                    }
                    free(msa_core);
                    free(msa_name);
                    PRINT_AND_RETURN("sequence allocation faiiled in no sequence", MALLOC_ERROR);
                }

                sequence_counter++;
                if(setup_name(&msa_name[sequence_counter], &line[1]) != SUCCESS){ // roll back and deallocate
                    for(i = 0; i < sequence_counter; i++){
                        free(msa_core[i]);
                        free(msa_name[i]);
                    }
                    free(msa_core);
                    free(msa_name);
                    PRINT_AND_RETURN("name allocation faiiled in no sequence", MALLOC_ERROR);
                }
                strclr(seq);
                break;
            default:
                strcat(seq, line);
        }
    }

    // Final iteration for the last sequence
    if(setup_sequence(&msa_core[sequence_counter], seq) != SUCCESS){ // roll back and deallocate
        for(i = 0; i < sequence_counter; i++){
            if(i != sequence_counter - 1) free(msa_core[i]);
            free(msa_name[i]);
        }
        free(msa_core);
        free(msa_name);
        PRINT_AND_RETURN("sequence allocation faiiled in no sequence", MALLOC_ERROR);
    }
    sequence_counter++;

    // Close file
    fclose(f);

    // Fill in the fields
    return(init_msa(msa_ptr, strlen(seq), sequence_counter, msa_core, msa_name));
}


/* Open a file and write msa struct in FASTA format
 * Input:   msa structure and the filename
 * Output:  0 on success, ERROR otherwise
 * Effect:  set fields in the new msa
 */
int write_msa(msa_t * msa, char * filename){
    FILE * f; //file to the opened
    int i; //loop variable

    if(!msa)        PRINT_AND_RETURN("msa is NULL in write_msa",                GENERAL_ERROR);
    if(!filename)   PRINT_AND_RETURN("filename is NULL is write_msa",           GENERAL_ERROR);

    // Try to open file
    f = fopen(filename, "w");
    if (!f)         PRINT_AND_RETURN("cannot open file to write in write_msa",  GENERAL_ERROR);

    // Write to file
    for(i = 0; i < msa->num_seq; i++){
        fprintf(f, ">%s\n", msa->name[i]);
        fprintf(f, "%s\n", msa->msa[i]);
    }

    // Close file
    fclose(f);
    return 0;
}

/* Extending an original msa to a blank msa with same meta
 * Input:   pointer to the two msa
 * Output:  0 on sucess, ERROR otherwise
 * Effect:  calls malloc, set fields in the new msa
 */
int make_smaller_msa(msa_t * original, msa_t * new){
    int i, j; //loop variable

    // Safe-checking 
    if(!original)           PRINT_AND_RETURN("original is NULL in make_smaller_msa",        GENERAL_ERROR);
    if(!new)                PRINT_AND_RETURN("new is NULL in make_smaller_msa",             GENERAL_ERROR);

    new->N          = original->N;
    new->num_seq    = 0;

    new->msa        = (char **) malloc(original->num_seq * sizeof(char *));
    new->name       = (char **) malloc(original->num_seq * sizeof(char *));
    if(!new->msa)   PRINT_AND_RETURN("malloc for new->msa failed in make_smaller_msa",      MALLOC_ERROR);
    if(!new->name)  {
        free(new->msa);
        PRINT_AND_RETURN("malloc for new->name failed in make_smaller_msa",     MALLOC_ERROR);
    }

    for(i = 0; i < original->num_seq; i++){
        new->msa[i]         = (char *) malloc(new->N);
        new->name[i]        = (char *) malloc(MAX_NAME_LENGTH);
        if(!new->msa[i]){ //roll back and free all allocated memories 
            for(j = 0; j < i; j++){
                free(new->msa[j]);
                free(new->name[j]);
            }
            free(new->msa);
            free(new->name);
            PRINT_AND_RETURN("malloc for new msa sequence failed in make_smaller_msa",  MALLOC_ERROR);
        }
        if(!new->name[i]){//roll back and free all allocated memories 
            for(j = 0; j < i; j++){
                free(new->msa[j]);
                free(new->name[j]);
            }
            free(new->msa[i]);
            free(new->msa);
            free(new->name);
            PRINT_AND_RETURN("malloc for new msa name failed in make_smaller_msa",      MALLOC_ERROR);
        }
        strclr(new->msa[i]);
        strclr(new->name[i]);
    }
    return 0;
}

/* Given a tree initialized with write newick, its centroid edge (from centroid_decomposition) and 2 empty msa structure, 
 * Fill the 2 msa with leave sequences obtained from breaking the tree into 2 subtrees at the centroid edge
 * Input:       2 centroids (in any order), 2 empty msa with meta set up (in any order) and the original msa
 * Output:      0 on success, ERROR otherwise 
 * Effect:      setting fields in the 2 empty msa and increases its sequence count
 */
int retrieve_msa_from_root(int centroid1, int centroid2, msa_t * msa1, msa_t * msa2, msa_t * all_msa){
    // Safe checking 
    if(!msa1)           PRINT_AND_RETURN("msa1 is NULL in retrieve_msa_from_root",      GENERAL_ERROR);
    if(!msa2)           PRINT_AND_RETURN("msa2 is NULL in retrieve_msa_from_root",      GENERAL_ERROR);
    if(!all_msa)        PRINT_AND_RETURN("all_msa is NULL in retrieve_msa_from_root",   GENERAL_ERROR);

    if(!(0 <= centroid1 && centroid1 < maxN)) PRINT_AND_RETURN("centroid1 is out of range in retrieve_msa_from_root",   GENERAL_ERROR);
    if(!(0 <= centroid2 && centroid2 < maxN)) PRINT_AND_RETURN("centroid2 is out of range in retrieve_msa_from_root",   GENERAL_ERROR);

    dfs_msa(centroid1, parent_map[centroid1], centroid2, all_msa, msa1);
    dfs_msa(centroid2, parent_map[centroid2], centroid1, all_msa, msa2);
    return 0;
}

/* Given a file created by hmmsearch, computing the joint likelihood of observing all the sequences in that file by adding up the bitscores
 * Input:       name of the file created by hmmsearch and expected number of sequences in the file
 * Output:      0 on success, ERROR otherwise
 * Effect:      open a file and read its content
 */
float compute_likelihood(char * filename, int num_seq){
    FILE * f;                               //file to be opened
    char buffer[GENERAL_BUFFER_SIZE];       //placeholder for bit score (must be longer than the longest width of hmmsearch output)
    int i;                                  //loop variable
    float result;                           //storing the resulting likelihood

    // Open file
    f = fopen(filename, "r");
    if (!f)                                 PRINT_AND_RETURN("cannot open file to write in hmmsearch_output",   GENERAL_ERROR);

    // Get to the 3rd line of the hmmsearch output file 
    fgets(buffer, sizeof(buffer), f);
    fgets(buffer, sizeof(buffer), f);

    result = 0.0;
    for(i = 0; i < num_seq; i++){
        // Skip this line
        fgets(buffer, sizeof(buffer), f);

        // The bit score is the 5th word on the current line
        scanf("%s", buffer);            
        scanf("%s", buffer);            
        scanf("%s", buffer);            
        scanf("%s", buffer);
        scanf("%s", buffer);
        result += atof(buffer); //TODO: add error checking for this part
    }

    // Close the file
    fclose(f);

    return result;
}


//INTERNAL FUNCTIONS IMPLEMENTATIONS

/* Wrapper to malloc name holder and copy name into that place
 * Input:   the name and the location to be copied to
 * Output:  0 on success, ERROR otherwise   
 * Effect:  calls malloc
 */
int setup_name(char ** name_holder, char* name){
    if(!name)   PRINT_AND_RETURN("name is NULL when trying to parse",   GENERAL_ERROR);

    *name_holder = (char *) malloc(strlen(name));
    if(!*name_holder) PRINT_AND_RETURN("malloc failed for name holder", MALLOC_ERROR);
    strcpy(*name_holder, name);
    return 0;
}

/* Wrapper to malloc sequence holder and copy sequence into that place
 * Input:   the sequence and the location to be copied to
 * Output:  0 on success, ERROR otherwise   
 * Effect:  calls malloc
 */
int setup_sequence(char ** sequence_holder, char* sequence){
    // The syntax is exactly the same as setup_name so we will just call that function instead
    return setup_name(sequence_holder, sequence);
}


/* Brute force algorithm to find sequence based on name
 * Input:   msa structure and the name of the sequence
 * Output:  pointer to the sequence on success, NULL otherwise
 * Effect:  none
 */
char * find_sequence_by_name(msa_t * msa, char * name){
    int i; //loop variable

    if(!msa)        PRINT_AND_RETURN("msa is NULL in find_sequence_by_name",    NULL);
    if(!name)       PRINT_AND_RETURN("name is NULL in find_sequence_by_name",   NULL);
    for(i = 0; i < msa->num_seq; i++)
        if(strcmp(name, msa->name[i]) == 0)
            return msa->msa[i];
    return NULL;
}

/* Brute force algorithm to add a sequence from an msa struct to another msa struct
 * Input:   both msa structures and the name of the sequence
 * Output:  0 on success, ERROR otherwise
 * Effect:  set fields in the new msa
 */
int add_to_msa_from_msa(msa_t * original, msa_t * new, char * name){
    char * sequence; //hold the pointer to the sequene to be copied

    if(!original)   PRINT_AND_RETURN("original is NULL in add_to_msa_from_msa", GENERAL_ERROR);
    if(!new)        PRINT_AND_RETURN("new is NULL in add_to_msa_from_msa",      GENERAL_ERROR);
    if(!name)       PRINT_AND_RETURN("name is NULL in add_to_msa_from_msa",     GENERAL_ERROR);

    sequence = find_sequence_by_name(original, name);
    if(!sequence)   PRINT_AND_RETURN("cannot find sequence with said name in add_to_msa_from_msa",      GENERAL_ERROR);

    strcpy(new->msa[new->num_seq], sequence);
    strcpy(new->name[new->num_seq], name);
    new->num_seq++;
    return 0;
}

/* DFS for traversing the tree from one centroid edge avoiding the other, adding the msa at each leave traversed to a provided msa.
 * Input:       cur_node    numeric node in the tree (not names) currently at
 *              par_node    numeric node of the parents last visited (instead of an array of visited nodes)
 *              all_msa     the original msa where the sequences come from
 *              small_msa   the new msa to which the sequences are added 
 * Output:  none
 * Effect:  write to the small_msa, incrementing its node count
 */
void dfs_msa(int cur_node, int par_node, int centroid, msa_t * all_msa, msa_t * small_msa){
    if(is_leave(cur_node)) add_to_msa_from_msa(all_msa, small_msa, name_map[cur_node]);
    for(int i = 0; i < maxN; i++){
        if(!adj_mat[cur_node][i] || i == par_node || i == cur_node || i == centroid) continue;

        dfs_msa(i, cur_node, centroid, all_msa, small_msa);
    }
}




