#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tree.h"
#include "utilities.h"

// Private functions
int size_dfs(int node, int parent);
int centroid_search(int node, int parent);
int check(int node_a);
int make_adjacent(int node_a, int node_b);
int make_parent(int node_p, int node_c);
int save_name(int cur_node, char * name);
int init();

/* Brute force algorithm to test whether a node is a leaf (has degree exactly 1 in a binary tree)
 * Input:   node to be tested
 * Output:  1 if the node is a leave, 0 otherwise
 * Effect:  none
 */ 
int is_leave(int cur_node){
    int seen_one_neighbor = 0;
    for(int i = 0; i < maxN; i++){
        if(adj_mat[cur_node][i]){
            if(seen_one_neighbor) return 0;
            else seen_one_neighbor = 1;
        }
    }
    return seen_one_neighbor;
}

/* Centroid decomposition of a tree into 2 subtrees by deleting the centroid edge. This requires the fields are already initialized by read_newick
 * Input:   pointer to write the 2 endpoints of the centroid edge to
 * Output:  0 on success, ERROR otherwise
 * Effect:  none
 */ 
int centroid_decomposition(int * left_subtree_root, int * right_subtree_root){
    if(!left_subtree_root)              PRINT_AND_RETURN("left_subtree_root is NULL in cetroid decomposition",  GENERAL_ERROR);
    if(!right_subtree_root)             PRINT_AND_RETURN("right_subtree_root is NULL in centroid decomposition",GENERAL_ERROR);

    // Assuming input is already read by read_newick and is by definition a tree
    size_dfs(0, -1);
    *left_subtree_root      = centroid_search(0, -1);
    *right_subtree_root     = parent_map[*left_subtree_root];
    return 0;
}

#define READ_NAME_STATE 0
#define OTHER_STATE     1

#define incr_level(node_p, node_c, max_node) do{node_p = node_c; incr_node(node_c, max_node);} while(0)
#define decr_level(node_p, node_c, max_node) do{node_c = node_p; node_p = parent_map[node_c];} while(0)
#define incr_node(node, max_node)               node = ++max_node

int read_newick(char * filename){
    init();
    FILE * f = freopen(filename, "r", stdin);
    if(!f) PRINT_AND_RETURN("fail to open read_newick file",    OPEN_ERROR);

    // FSM for reading newick format
    int cur_state = READ_NAME_STATE;
    int cur_node = 0;
    int max_node = 0;
    int parent_node = -1;
    char cur_char = 0;
    char cur_name[maxNameSize];
    strclr(cur_name);

    while(scanf("%c", &cur_char) == 1){
        printf("cur char is %c cur node is %d cur par is %d cur state is %d\n", cur_char, cur_node, parent_node, cur_state);
        switch(cur_char){
            case '(': //start of a new level and start of a new node
                incr_level(parent_node, cur_node, max_node);
                make_parent(parent_node, cur_node);
                cur_state = READ_NAME_STATE;
                break;
            case ',': // start of new node end of old node
                save_name(cur_node, cur_name);
                incr_node(cur_node, max_node);
                make_parent(parent_node, cur_node);
                cur_state = READ_NAME_STATE;
                break;
            case ')': // end of level
                save_name(cur_node, cur_name);
                decr_level(parent_node, cur_node, max_node);
                break;
            case ':': 
                cur_state = OTHER_STATE;
                break;
            default:
                // printf("%c %d %s\n", cur_char, cur_state, cur_name);
                if(cur_state == READ_NAME_STATE){
                    char buf[2];
                    buf[0] = cur_char;
                    buf[1] = 0;
                    strcat(cur_name, buf);
                }
        }
    }

    fclose(f);
    return 0;
}

// INTERNAL FUNCTION IMPLEMENTATIONS

/* Helper function to check whether a node is in range
 * Input:   node to check
 * Output:  1 if the node is in range, 0 otherwise
 * Effect:  none
 */ 
int check(int node_a){
    if(0 <= node_a && node_a < maxN) return 1;
    else return 0;
}

/* Helper function to make a node adjacent to another in the adjacency matrix
 * Input:   2 nodes
 * Output:  0 on success, ERROR otherwise
 * Effect:  none
 */ 
int make_adjacent(int node_a, int node_b){
    if(!check(node_a) || !check(node_b))  PRINT_AND_RETURN("node out of range in make_adjacent", GENERAL_ERROR);

    adj_mat[node_a][node_b] = 1;
    adj_mat[node_b][node_a] = 1;

    return 0;
}

/* Helper function to make a node a parent of another 
 * Input:   2 nodes (order matters, the parent node comes first)
 * Output:  0 on success, ERROR otherwise
 * Effect:  none
 */ 
int make_parent(int node_p, int node_c){
    // printf("make paret %d is paret f %d\n", node_p, node_c);
    if(make_adjacent(node_p, node_c) != SUCCESS)    PRINT_AND_RETURN("node out of range in make_parent", GENERAL_ERROR);

    parent_map[node_c] = node_p;
    return 0;
}

/* Helper function to save name for a numeric node
 * Input:   node number and its name 
 * Output:  0 on success
 * Effect:  none
 */ 
int save_name(int cur_node, char * name){
    printf("save_name %d %s\n", cur_node, name);
    if(strempty(name)){
        char buf[1000];
        snprintf(buf, 1000, "%d", cur_node);
        strcpy(name_map[cur_node], buf);
    } else {
        strcpy(name_map[cur_node], name);
    }
    strclr(name); 

    return 0;
}

/* Helper function to intialize global variables to some starting state
 * Input:   noner
 * Output:  0 on success
 * Effect:  none
 */ 
int init(){
    for(int i = 0; i < maxN; i++){
        subtree_size[i] = 0;
        parent_map[i] = -1;
        for(int j = 0; j < maxNameSize; j++){
            name_map[i][j] = 0;
        }
        for(int j = 0; j < maxN; j++){
            adj_mat[i][j] = 0;
        }
    }
    return 0;
}

/* DFS routine to compute the size of the leaf set of the subtree rooted at some particular node
 * Input:   a node and it's parent (recursive state)
 * Output:  0 on success
 * Effect:  none
 */ 
int size_dfs(int node, int parent){ 
        // printf("%d\n", node);

    if(is_leave(node)) subtree_size[node] = 1;

    for(int i = 0; i < maxN; i++){
        if(!adj_mat[node][i] || i == parent || i == node) continue;
        size_dfs(i, node);
        subtree_size[node] += subtree_size[i];
    }
    return 0;
}

/* Greedy routine to find the centroid node of a tree. This requires subtree_size to have been initialized
 * Input:   a node and it's parent (recursive state)
 * Output:  the centroid node of the tree
 * Effect:  none
 */ 
int centroid_search(int node, int parent){ 
    int is_centroid = 1;
    int heaviest_child = 0;
    for(int i = 0; i < maxN; i++){
        if(!adj_mat[node][i] || i == parent || i == node) continue;
        if(subtree_size[i] > subtree_size[0] / 2) 
            is_centroid = 0;

        if(heaviest_child == 0 || subtree_size[i] > subtree_size[heaviest_child]) 
            heaviest_child = i;
    }
    if(is_centroid && subtree_size[node] <= subtree_size[0] / 2)
        return node;
    return centroid_search(heaviest_child, node);
}

