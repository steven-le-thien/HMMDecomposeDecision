#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tree.h"
#include "utilities.h"

extern int centroid_decomposition(int * left_subtree_root, int * right_subtree_root);
extern int read_newick(char * filename);
extern int is_leave(int cur_node);

int is_leave(int cur_node){
	int seen_one_neighbor = 0;
	for(int i = 0; i < maxN; i++){
		if(adj_mat[cur_node][i]){
			if(seen_one_neighbor) return -1;
			else seen_one_neighbor = 1;
		}
	}
	return 0;
}

int check(int node_a){
	if(0 <= node_a && node_a < maxN) return 1;
	else return 0;
}

int make_adjacent(int node_a, int node_b){
	if(!check(node_a) || !check(node_b))  PRINT_AND_RETURN("node out of range in make_adjacent", GENERAL_ERROR);

	adj_mat[node_a][node_b] = 1;
	adj_mat[node_b][node_a] = 1;

	return 0;
}

int make_parent(int node_p, int node_c){
	if(!make_adjacent(node_p, node_c))	PRINT_AND_RETURN("node out of range in make_parent", GENERAL_ERROR);

	parent_map[node_c] = node_p;
	return 0;
}

#define READ_NAME_STATE 0
#define OTHER_STATE 	1

int read_newick(char * filename){
	FILE * f = freopen(filename, "r", stdin);
	if(!f) PRINT_AND_RETURN("fail to open read_newick file",	OPEN_ERROR);

	// FSM for reading newick format
	int cur_state = READ_NAME_STATE;
	int cur_node = 0;
	int parent_node = -1;
	char cur_char;
	char cur_name[maxNameSize];
	strclr(cur_name);
	while(scanf("%c", cur_char)){
		switch(cur_char){
			case '(': //start of a new level and start of a new node
				cur_state = READ_NAME_STATE;
				cur_node++;
		}
	}

	fclose(f);
}