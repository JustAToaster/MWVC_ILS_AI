#include <fstream>
#include <cstdlib>
#include <iostream>
#include <string>
#include <cmath>

#include "Graph.h"

#define MAX_ITER_WITHOUT_IMPROVEMENTS 8192
#define MAX_OBJECTIVE_FUNCTION_EVAL 2e4
#define MAX_LS_ITER 512

using namespace std;

void perturbation(bool* perturbed_solution, int num_elems_changed, int n);

bool acceptance_criteria(int curr_best_weight, int prev_weight, int curr_weight, int eps, int iter_without_improvements);

int epsilon_scheduling(int eps, int curr_weight, int curr_best_weight, int min_weight, int iter);

int perturbation_scheduling(int num_elems_changed, int pert_min_changes, int pert_tries);

void print_solution(bool* solution, int n);

bool local_search(bool* curr_solution, int& curr_weight, Graph& graph, int n, int iter, int& num_obj_func_eval, bool VERBOSE_FLAG);

void separate_nodes(bool* curr_solution, int* present_nodes, int* external_nodes, int& n_present, int& n_external, int n);

bool local_search_stochastic(bool* curr_solution, int& curr_weight, Graph& graph, int n, int num_swaps, int iter, int& num_obj_func_eval, bool VERBOSE_FLAG);