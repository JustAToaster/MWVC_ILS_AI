#include <fstream>
#include <cstdlib>
#include <iostream>
#include <string>
#include <cmath>
#include <utility>
#include <time.h>

#include "optimization.h"

using namespace std;

//Funzione obiettivo, la somma dei pesi dei nodi di una soluzione
int Graph::total_weight(bool* solution){
	int total_weight = 0;
	for(int i = 0; i < n; ++i){
		total_weight += (solution[i] ? w[i] : 0);
	}
	return total_weight;
}

//Metodo per verificare la validita' di una soluzione rispetto al vincolo del problema
bool Graph::valid_solution(bool* solution){
	pair<int, int> edge;
	for(int i = 0; i < num_edges; ++i){
		edge = edges[i];
		if (!solution[edge.first] && !solution[edge.second]) return false;
	}
	return true;
}


void Graph::compute_degree_weight_ratio(){
	degree_over_weight = new float[n]();
	int degrees[n];
	fill_n(degrees, n, 0);
	pair<int, int> edge;
	for(int i = 0; i < num_edges; ++i){
		edge = edges[i];
		degrees[edge.first]++;
	}
	for(int i = 0; i < n; ++i){
		degree_over_weight[i] = (float)degrees[i]/w[i];
	}
}

//Per ogni arco non coperto, inserisci l'estremo con rapporto degree/peso piu' alto
void Graph::greedy_solution(bool* solution){
	compute_degree_weight_ratio();
	pair<int, int> edge;
	for(int i = 0; i < num_edges; ++i){
		edge = edges[i];
		if (!solution[edge.first] && !solution[edge.second]){
			if (degree_over_weight[edge.first] > degree_over_weight[edge.second]) solution[edge.first] = 1;
			else solution[edge.second] = 1;
		}
	}
}

//Usa l'euristica greedy per sistemare una soluzione invalida (ottenuta ad esempio da una perturbazione)
void Graph::fix_invalid_solution(bool* solution){
	pair<int, int> edge;
	for(int i = 0; i < num_edges; ++i){
		edge = edges[i];
		if (!solution[edge.first] && !solution[edge.second]){
			if (degree_over_weight[edge.first] > degree_over_weight[edge.second]) solution[edge.first] = 1;
			else solution[edge.second] = 1;
		}
	}
}

//Per ogni arco non coperto, il nodo tra i due da aggiungere alla soluzione viene scelto con probabilita' proporzionale al rapporto degree/peso
//Sembra dare risultati peggiori
void Graph::fix_invalid_solution_stochastic(bool* solution){
	pair<int, int> edge;
	float first_node_prob;
	for(int i = 0; i < num_edges; ++i){
		edge = edges[i];
		if (!solution[edge.first] && !solution[edge.second]){
			first_node_prob = degree_over_weight[edge.first]/(degree_over_weight[edge.first] + degree_over_weight[edge.second]);
			if (rand()/RAND_MAX <= first_node_prob) solution[edge.first] = 1;
			else solution[edge.second] = 1;
		}
	}
}

void print_solution(bool* solution, int n){
	for(int i = 0; i < n; ++i){
		if(solution[i]) cout << i << "|";
	}
	cout << endl;
}

//Bit flip casuali su una soluzione. Potrebbe portare a soluzioni invalide, quindi va sistemata la soluzione o ripetuta la perturbazione.
void perturbation(bool* perturbed_solution, int num_elems_changed, int n){
	int rand_index;
	for(int i = 0; i < num_elems_changed; ++i){
		rand_index = rand() % n;
		perturbed_solution[rand_index] = !perturbed_solution[rand_index];
	}
}

//Criterio di accettazione di un nuovo ottimo locale per ILS
bool acceptance_criteria(int curr_best_weight, int prev_weight, int curr_weight, int eps, int iter_without_improvements){
	// Se la soluzione risulta migliore della precedente o dell'ottimo attuale, accettala
	if(curr_weight < curr_best_weight || curr_weight < prev_weight) return true;
	// Se risulta peggiore di un epsilon e non stiamo deviando per troppo tempo, accettala
	if(curr_weight < curr_best_weight + eps && iter_without_improvements < MAX_ITER_WITHOUT_IMPROVEMENTS) return true;
	return false;
}

//Scheduling della variabile epsilon, che definisce la deviazione massima concessa rispetto al migliore attuale
int epsilon_scheduling(int eps, int curr_weight, int curr_best_weight, int min_weight, int iter){
	int new_eps;
	if (curr_weight < curr_best_weight) new_eps = eps - 2*min_weight;
	else new_eps = eps;
	if (new_eps > min_weight) return new_eps;
	else return min_weight;
}

//Scheduling del numero di bit da flippare in una perturbazione, deve diminuire col tempo per non generare soluzioni troppo lontane
int perturbation_scheduling(int num_elems_changed, int pert_min_changes, int iter){
	int new_num_elems_changed = num_elems_changed;
	if(!(iter & 31)) new_num_elems_changed--;
	if(new_num_elems_changed <= pert_min_changes) new_num_elems_changed = pert_min_changes;
	return new_num_elems_changed;
}

//Ricerca locale con sostituzioni o rimozioni di nodi dalla soluzione. Ritorna true se il numero di valutazioni della funzione obiettivo supera il massimo
bool local_search(bool* curr_solution, int& curr_weight, Graph& graph, int n, int iter, int& num_obj_func_eval, bool VERBOSE_FLAG){

	bool ls_stuck = false;
	curr_weight = graph.total_weight(curr_solution);
	num_obj_func_eval++;
	if(num_obj_func_eval >= MAX_OBJECTIVE_FUNCTION_EVAL) return true;
	int prev_weight = curr_weight;

	bool curr_local_best_solution[n];
	copy(curr_solution, curr_solution+n, curr_local_best_solution);
	int curr_local_best_weight = curr_weight;

	bool neighbor[n];
	int neighbor_weight = curr_weight;

	for(int ls_iter = 1; ls_iter <= MAX_LS_ITER && !ls_stuck; ++ls_iter){
		if(VERBOSE_FLAG) cout << "Local search numero " << iter << ", iterazione " << ls_iter << ", esploro il vicinato..." << endl;

		//Esplorazione vicini
		for (int i = 0; i < n; ++i){
			for (int j = i; j < n; ++j){
				//Vicini rimuovendo un elemento dalla soluzione attuale
				if (i == j && curr_solution[i]){
					copy(curr_solution, curr_solution+n, neighbor);
					neighbor[i] = 0;
				}
				//Se gli indici sono sullo stesso bit e il valore e' 0 oppure i due bit non sono uguali, non si puo' fare rimozione o swap
				else if ((i == j && !curr_solution[i]) || (curr_solution[i] == curr_solution[j])){
					continue;
				}
				//Un bit e' 0 e l'altro 1, effettua uno swap dei nodi
				else{
					copy(curr_solution, curr_solution+n, neighbor);
					neighbor[i] = !neighbor[i];
					neighbor[j] = !neighbor[j];
				}
				if (graph.valid_solution(neighbor)){
					if(neighbor[i]) neighbor_weight = curr_weight + graph.get_weight(i) - graph.get_weight(j);
					else{
						if(neighbor[j]) neighbor_weight = curr_weight - graph.get_weight(i) + graph.get_weight(j);
						else neighbor_weight = curr_weight - graph.get_weight(i);
					}
					num_obj_func_eval++;
					if (neighbor_weight < curr_local_best_weight){
						copy(neighbor, neighbor+n, curr_local_best_solution);
						curr_local_best_weight = neighbor_weight;
					}
					if(num_obj_func_eval >= MAX_OBJECTIVE_FUNCTION_EVAL){
						copy(curr_local_best_solution, curr_local_best_solution+n, curr_solution);
						curr_weight = curr_local_best_weight;
						return true;
					}
				}
			}
		}

		//Aggiornamento soluzione attuale con il migliore vicino
		copy(curr_local_best_solution, curr_local_best_solution+n, curr_solution);
		curr_weight = curr_local_best_weight;

		//cout << "Soluzione attuale: " << endl;
		//print_solution(curr_solution, n);
		if(VERBOSE_FLAG) cout << "Peso soluzione attuale: " << curr_weight << endl;

		if(curr_weight == prev_weight){
			if (VERBOSE_FLAG) cout << "Local search numero " << iter << ", iterazione " << ls_iter << ". Bloccato in un ottimo locale. Uscita..." << endl;
			//print_solution(curr_local_best_solution, n);
			ls_stuck = true;
		}

		prev_weight = curr_weight;

	}

	curr_weight = curr_local_best_weight;

	return false;
}