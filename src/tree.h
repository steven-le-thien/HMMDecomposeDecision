// File in HMMDecompositionDecision, created by Thien Le in July 2018

#ifndef TREE_H
#define TREE_H

#define maxN 			10000 // these are number of nodes, not just number of leaves
#define maxNameSize 	1000

// Size of subtree rooted at a particular node
int subtree_size[maxN]; 
int parent_map[maxN];
char name_map[maxN][maxNameSize];

int adj_mat[maxN][maxN];

int debug_counter;

extern int centroid_decomposition(int * left_subtree_root, int * right_subtree_root);
extern int read_newick(char * filename);
extern int is_leave(int cur_node);

#endif