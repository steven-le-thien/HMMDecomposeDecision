#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tree.h"
#include "utilities.h"

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

int size_dfs(int node, int parent){ 
	if(is_leave(node)) subtree_size[node] = 1;

	for(int i = 0; i < maxN; i++){
		if(!adj_mat[node][i] || i == parent || i == node) continue;
		size_dfs(i, node);
		subtree_size[node] += subtree_size[i];
	}
	return 0;
}

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

int centroid_decomposition(int * left_subtree_root, int * right_subtree_root){
	// Assuming input is already read by read_newick and is by definition a tree
	size_dfs(0, -1);
	left_subtree_root[0] = centroid_search(0, -1);
	right_subtree_root[0] = parent_map[left_subtree_root[0]];
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
	if(make_adjacent(node_p, node_c) != SUCCESS)	PRINT_AND_RETURN("node out of range in make_parent", GENERAL_ERROR);

	parent_map[node_c] = node_p;
	return 0;
}

int save_name(int cur_node, char * name){
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

void incr_level(int * parent_node, int * cur_node){
	parent_node[0] = cur_node[0];
}

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

#define READ_NAME_STATE 0
#define OTHER_STATE 	1

#define incr_level(node_p, node_c) do{node_p = node_c; cur_node++;} while(0)
#define decr_level(node_p, node_c) do{node_c = node_p; node_p = parent_map[node_c];} while(0)
#define incr_node(node) node++

int read_newick(char * filename){
	init();
	FILE * f = freopen(filename, "r", stdin);
	if(!f) PRINT_AND_RETURN("fail to open read_newick file",	OPEN_ERROR);

	// FSM for reading newick format
	int cur_state = READ_NAME_STATE;
	int cur_node = 0;
	int parent_node = -1;
	char cur_char = 0;
	char cur_name[maxNameSize];
	strclr(cur_name);

	while(scanf("%c", &cur_char)){
		switch(cur_char){
			case '(': //start of a new level and start of a new node
				incr_level(parent_node, cur_node);
				make_parent(cur_node, parent_node);
				cur_state = READ_NAME_STATE;
				break;
			case ',': // start of new node end of old node
				save_name(cur_node, cur_name);
				incr_node(cur_node);
				make_parent(cur_node, parent_node);
				cur_state = READ_NAME_STATE;
				break;
			case ')': // end of level
				save_name(cur_node, cur_name);
				decr_level(parent_node, cur_node);
				break;
			case ':': 
				cur_state = OTHER_STATE;
				break;
			default:
				if(cur_node == READ_NAME_STATE){
					char buf[1];
					buf[0] = cur_char;
					strcat(cur_name, buf);
				}
		}
	}

	fclose(f);
	return 0;
}