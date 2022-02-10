#include <fstream>
#include <cstdlib>
#include <iostream>
#include <string>
#include <cmath>
#include <vector>

using namespace std;

class Graph{
	private:
		int n;	//Numero di nodi

		int num_edges;	//Numero di archi

		pair<int, int> *edges; //Array di archi

		int *w;	//Funzione peso dei nodi

		int *degrees; //Degree dei nodi

		float *degree_over_weight; //Rapporto grado/peso dei nodi
		
	public:
		Graph();
		
		~Graph();
		
		void initialize_graph_nodes(int n);
		
		void add_node_weight(int index, int weight);
		
		void add_edges(vector< pair<int, int> > vec_edges);
		
		int get_weight(int i);
		
		int min_weight();
		
		void print_edges();
		
		int total_weight(bool* solution);

		void compute_degrees();
				
		void compute_degree_weight_ratio();

		bool valid_solution(bool* solution);

		void greedy_heuristic(bool* solution);

		void greedy_heuristic_stochastic(bool* solution);

		void fix_invalid_solution(bool* solution);

		void fix_invalid_solution_stochastic(bool* solution);
};