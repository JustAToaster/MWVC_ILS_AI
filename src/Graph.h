#include <fstream>
#include <cstdlib>
#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <list>

using namespace std;

class Graph{
	private:
		int n;	//Numero di nodi

		int num_edges;	//Numero di archi

		pair<int, int> *edges; //Array di archi

		list<int> *adj_lists; //Liste di adiacenza

		int *w;	//Funzione peso dei nodi

		int *degrees; //Degree dei nodi

		float *degree_over_weight; //Rapporto grado/peso dei nodi
		
	public:
		Graph();
		
		~Graph();
		
		void initialize_graph_nodes(int n);
		
		void add_node_weight(int index, int weight);
		
		void add_edges(vector< pair<int, int> > vec_edges);

		void build_adj_lists(vector< pair<int, int> > vec_edges);
		
		int get_weight(int i);
		
		int min_weight();
		
		void print_edges();
		
		int total_weight(bool* solution);

		void compute_degrees();
				
		void compute_degree_weight_ratio();

		//Seguono le funzioni specifiche al problema

		vector< pair<int, int> > compute_uncovered_edges(bool* solution);

		void remove_covered_edges(int added_node, vector< pair<int, int> >& uncovered_edges);
		
		bool valid_solution(bool* solution);

		bool valid_solution_node(int removed_node, bool* solution);

		bool valid_solution_two_nodes(int removed_node1, int removed_node2, bool* solution);

		void greedy_heuristic(bool* solution);

		void greedy_heuristic_prob(bool* solution);

		void greedy_heuristic_stochastic(bool* solution);

		void greedy_heuristic_queue(bool* solution);
		
		void greedy_heuristic_queue_prob(bool* solution);

		int random_neighbor(int node, bool* solution, bool neighbor_present);

};