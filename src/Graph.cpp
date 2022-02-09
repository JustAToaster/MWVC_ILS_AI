#include <fstream>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <algorithm>

#include "Graph.h"

using namespace std;

Graph::Graph(){
	this->n = 0;
}

Graph::~Graph(){
	delete w;
	delete degree_over_weight;
	delete edges;
}

void Graph::initialize_graph_nodes(int n){
	this->n = n;
	w = new int[n];
	fill_n(w, n, 0);
}

void Graph::add_node_weight(int index, int weight){
	w[index] = weight;
}

void Graph::add_edges(vector< pair<int, int> > vec_edges){
	this->num_edges = vec_edges.size();
	this->edges = new pair<int,int>[num_edges];
	copy(vec_edges.begin(), vec_edges.end(), this->edges);
}

int Graph::get_weight(int i){
	/*if(i < 0 || i > n){
		cerr << "Out of bounds" << endl;
		exit(1);
	}*/
	return w[i];
}

int Graph::min_weight(){
	return *min_element(w, w + n);
}

void Graph::print_edges(){
	pair<int, int> edge;
	cout << "Numero di archi: " << num_edges << endl;
	for(int i = 0; i < num_edges; ++i){
		edge = edges[i];
		cout << "(" << edge.first << "," << edge.second << ")|";
	}
	cout << endl;
}