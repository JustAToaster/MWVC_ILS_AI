#include <fstream>
#include <cstdlib>
#include <iostream>
#include <string>
#include <cmath>
#include <queue>

#include "Graph.h"
#include "utilities.h"

#define MAX_ITER_WITHOUT_IMPROVEMENTS 256
#define MAX_OBJECTIVE_FUNCTION_EVAL 2e4
#define MAX_LS_ITER 512

using namespace std;

float random_01();

int random_index(int size);

void perturbation(bool* perturbed_solution, int num_elems_changed, int n);

bool acceptance_criteria(int curr_best_weight, int prev_weight, int curr_weight, int eps, int iter_without_improvements);

int epsilon_scheduling(int eps, int curr_weight, int curr_best_weight, int min_weight, int min_eps);

int perturbation_scheduling(int num_elems_changed, float eps, float min_eps, float max_eps, int pert_min_changes, int pert_max_changes);

void print_solution(bool* solution, int n);

bool local_search(bool* curr_solution, int& curr_weight, Graph& graph, int n, int iter, int& num_obj_func_eval, bool VERBOSE_FLAG);

void separate_nodes(bool* curr_solution, int* present_nodes, int* external_nodes, int& n_present, int& n_external, int n);

void separate_nodes_vector(bool* curr_solution, vector<int>& present_nodes, vector<int>& external_nodes, int n);

bool local_search_stochastic(bool* curr_solution, int& curr_weight, Graph& graph, int n, int num_swaps, int& actual_swaps, int iter, int& num_obj_func_eval, bool VERBOSE_FLAG);

bool local_search_stochastic_removals(bool* curr_solution, int& curr_weight, Graph& graph, int n, int num_removals, int num_swaps, int& actual_swaps, int iter, int& num_obj_func_eval, bool VERBOSE_FLAG);

bool local_search_stochastic_vector(bool* curr_solution, int& curr_weight, Graph& graph, int n, int num_swaps, int iter, int& num_obj_func_eval, bool VERBOSE_FLAG);

void build_node_queue(priority_queue< pair<float, int> >& queue, int* weights, int n, vector< pair<int, int> > uncovered_edges);

void build_node_queue(priority_queue< pair<float, int> >& queue, int* weights, int n, vector< pair<int, int> > uncovered_edges, float& min_ratio, float& max_ratio);

void build_node_vec(vector< pair<float, int> >& node_vec, int* weights, int n, vector< pair<int, int> > uncovered_edges);