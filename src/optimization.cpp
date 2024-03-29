#include <fstream>
#include <cstdlib>
#include <iostream>
#include <string>
#include <cmath>
#include <utility>
#include <time.h>
#include <vector>
#include <algorithm>
#include <numeric>
#include <queue>
#include <unordered_set>
#include <iterator>
#include <float.h>

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
	for(int first_node = 0; first_node < n; ++first_node){
		for (const int & second_node : adj_lists[first_node]){
			if(!solution[first_node] && !solution[second_node]) return false;
		}
	}
	return true;
}

//Se viene rimosso un solo nodo dalla soluzione, verifica che tutti gli archi della sua lista di adiacenza siano ancora coperti
bool Graph::valid_solution_node(int removed_node, bool* solution){
	for (const int & second_node : adj_lists[removed_node]){
		if(!solution[second_node]) return false;
	}
	return true;
}

//Se vengono rimossi due nodi dalla soluzione, verifica che tutti gli archi delle loro liste di adiacenza siano ancora coperti
bool Graph::valid_solution_two_nodes(int removed_node1, int removed_node2, bool* solution){
	for (const int & second_node : adj_lists[removed_node1]){
		if(!solution[second_node]) return false;
	}
	for (const int & second_node : adj_lists[removed_node2]){
		if(!solution[second_node]) return false;
	}
	return true;
}

//Salva i gradi dei nodi
void Graph::compute_degrees(){
	degrees = new int[n];
	fill_n(degrees, n, 0);
	for(int first_node = 0; first_node < n; ++first_node){
		degrees[first_node] = adj_lists[first_node].size();
	}
}

//Dividi i gradi per i pesi, per ottenere una stima di importanza dei nodi
void Graph::compute_degree_weight_ratio(){
	degree_over_weight = new float[n]();
	compute_degrees();
	for(int i = 0; i < n; ++i){
		degree_over_weight[i] = (float)degrees[i]/w[i];
	}
}

//Per ogni arco non coperto, inserisci l'estremo con grado piu' alto
//Si può usare pure per sistemare una soluzione non valida (ottenuta ad esempio da una perturbazione)
void Graph::greedy_heuristic(bool* solution){
	for(int first_node = 0; first_node < n; ++first_node){
		for (const int & second_node : adj_lists[first_node]){
			if(!solution[first_node] && !solution[second_node]){
				if (degree_over_weight[first_node] > degree_over_weight[second_node]) solution[first_node] = 1;
				else solution[second_node] = 1;
			}
		}
	}
}

//Per ogni arco non coperto, il nodo tra i due da aggiungere alla soluzione viene scelto con probabilita' proporzionale al rapporto degree/peso
void Graph::greedy_heuristic_prob(bool* solution){
	float first_node_prob;
	for(int first_node = 0; first_node < n; ++first_node){
		for (const int & second_node : adj_lists[first_node]){
			if(!solution[first_node] && !solution[second_node]){
				first_node_prob = degree_over_weight[first_node]/(degree_over_weight[first_node] + degree_over_weight[second_node]);
				if (random_01() <= first_node_prob) solution[first_node] = 1;
				else solution[second_node] = 1;
			}
		}
	}
}

//Calcola tutti gli archi non coperti da una soluzione
vector< pair<int, int> > Graph::compute_uncovered_edges(bool* solution){
	vector< pair<int, int> > uncovered_edges;
	for(int first_node = 0; first_node < n; ++first_node){
		for (const int & second_node : adj_lists[first_node]){
			if(!solution[first_node] && !solution[second_node]) uncovered_edges.push_back(make_pair(first_node, second_node));
		}
	}
	return uncovered_edges;
}

//Dopo un aggiornamento della soluzione aggiungendo un nodo, rimuovi gli archi che ora copre dal vettore di archi non coperti
void Graph::remove_covered_edges(int added_node, vector< pair<int, int> >& uncovered_edges){
	pair<int, int> edge;
	for (const int & adj_node : adj_lists[added_node]){
		uncovered_edges.erase(remove(uncovered_edges.begin(), uncovered_edges.end(), make_pair(added_node, adj_node)), uncovered_edges.end());
		uncovered_edges.erase(remove(uncovered_edges.begin(), uncovered_edges.end(), make_pair(adj_node, added_node)), uncovered_edges.end());
	}
}

//Presi gli archi non coperti, continua a pescarne uno a caso a caso e aggiungere il nodo ottimale fino a quando non è una soluzione valida
void Graph::greedy_heuristic_stochastic(bool* solution){
	pair<int, int> rand_edge;
	vector< pair<int, int> > uncovered_edges = compute_uncovered_edges(solution);
	int added_node;
	//int iter = 0;
	float first_node_prob;
	while(uncovered_edges.size() > 0){
		//iter++;
		rand_edge = uncovered_edges.at(random_index(uncovered_edges.size()));
		first_node_prob = degree_over_weight[rand_edge.first]/(degree_over_weight[rand_edge.first] + degree_over_weight[rand_edge.second]);
		if (random_01() <= first_node_prob) added_node = rand_edge.first;
		else added_node = rand_edge.second;
		solution[added_node] = 1;
		remove_covered_edges(added_node, uncovered_edges);
	}
}

//Costruisci una coda di priorità basata sul rapporto grado/peso nel sottografo degli archi non coperti
void build_node_queue(priority_queue< pair<float, int> >& queue, int* weights, int n, vector< pair<int, int> > uncovered_edges){
	queue = priority_queue< pair<float, int> >();
	unordered_set<int> nodes;
	int degrees[n] = {};
	for(pair<int, int> edge: uncovered_edges){
		degrees[edge.first]++;
		degrees[edge.second]++;
		nodes.insert(edge.first);
		nodes.insert(edge.second);
	}
	for(const int& node: nodes){
		queue.push(make_pair((float)degrees[node]/(float)weights[node], node));
	}
}

//Costruisci la coda e salva il minimo e il massimo valore
void build_node_queue(priority_queue< pair<float, int> >& queue, int* weights, int n, vector< pair<int, int> > uncovered_edges, float& min_ratio, float& max_ratio){
	queue = priority_queue< pair<float, int> >();
	unordered_set<int> nodes;
	int degrees[n] = {};
	float node_ratio;
	for(pair<int, int> edge: uncovered_edges){
		degrees[edge.first]++;
		degrees[edge.second]++;
		nodes.insert(edge.first);
		nodes.insert(edge.second);
	}
	for(const int& node: nodes){
		node_ratio = (float)degrees[node]/(float)weights[node];
		if (node_ratio < min_ratio) min_ratio = node_ratio;
		if (node_ratio > max_ratio) max_ratio = node_ratio;
		queue.push(make_pair(node_ratio, node));
	}
}

//Inserisci ad ogni istante il nodo ottimale che massimizza il rapporto grado/peso nel sottografo degli archi non coperti
void Graph::greedy_heuristic_queue(bool* solution){
	vector< pair<int, int> > uncovered_edges = compute_uncovered_edges(solution);
	priority_queue< pair<float, int> > node_queue;
	build_node_queue(node_queue, w, n, uncovered_edges);
	int node_index;
	while(uncovered_edges.size() > 0){
		node_index = node_queue.top().second;
		node_queue.pop();
		solution[node_index] = 1;
		remove_covered_edges(node_index, uncovered_edges);
		build_node_queue(node_queue, w, n, uncovered_edges);
	}
}

//Versione randomizzata del miglior algoritmo greedy, per dare più possibilità a nodi apparentemente meno ottimali
void Graph::greedy_heuristic_queue_prob(bool* solution){
	pair<int, int> rand_edge;
	vector< pair<int, int> > uncovered_edges = compute_uncovered_edges(solution);
	priority_queue< pair<float, int> > node_queue;
	float min_ratio = FLT_MAX;
	float max_ratio = FLT_MIN;
	build_node_queue(node_queue, w, n, uncovered_edges, min_ratio, max_ratio);
	int node_index;
	int best_node_index = node_queue.top().second;
	float node_ratio;
	float node_importance = 1.0f;
	while(uncovered_edges.size() > 0){
		node_ratio = node_queue.top().first;
		node_index = node_queue.top().second;
		node_queue.pop();
		if(max_ratio > min_ratio) node_importance = (node_ratio - min_ratio)/(max_ratio - min_ratio);
		else node_importance = 1.0f;
		//Usa una funzione per trasformare l'importanza del nodo in probabilità di essere inserito
		if (random_01() <= probability_function(node_importance, -0.5f, 's')){
			solution[node_index] = 1;
			remove_covered_edges(node_index, uncovered_edges);
			build_node_queue(node_queue, w, n, uncovered_edges, min_ratio, max_ratio);
			best_node_index = node_queue.top().second;
		}
		//Gestisci caso limite in cui nessun nodo viene scelto
		if(node_queue.empty()){
			solution[best_node_index] = 1;
			remove_covered_edges(best_node_index, uncovered_edges);
			build_node_queue(node_queue, w, n, uncovered_edges, min_ratio, max_ratio);
			best_node_index = node_queue.top().second;
		}
	}
}

void build_node_vec(vector< pair<float, int> >& node_vec, int* weights, int n, vector< pair<int, int> > uncovered_edges){
	node_vec = vector< pair<float, int> >();
	unordered_set<int> nodes;
	int degrees[n] = {};
	float node_ratio;
	for(pair<int, int> edge: uncovered_edges){
		degrees[edge.first]++;
		degrees[edge.second]++;
		nodes.insert(edge.first);
		nodes.insert(edge.second);
	}
	for(const int& node: nodes){
		node_ratio = (float)degrees[node]/(float)weights[node];
		node_vec.push_back(make_pair(node_ratio, node));
	}
	sort(node_vec.begin(), node_vec.end(), [](const pair<float, int>& x, const pair<float, int>& y) { return x.first > y.first; } );
}

//Versione dell'algoritmo greedy randomizzato con vettore e prendendo un indice normale casuale
void Graph::greedy_heuristic_vec_prob(bool* solution){
	pair<int, int> rand_edge;
	vector< pair<int, int> > uncovered_edges = compute_uncovered_edges(solution);
	vector< pair<float, int> > node_vec;
	build_node_vec(node_vec, w, n, uncovered_edges);
	int index_std = 1;
	int vec_index, node_index;
	while(uncovered_edges.size() > 0 && node_vec.size() > 0){
		vec_index = random_normal_index(0, index_std);
		if(vec_index >= (int)node_vec.size()) vec_index = node_vec.size() - 1;
		node_index = node_vec.at(vec_index).second;
		solution[node_index] = 1;
		remove_covered_edges(node_index, uncovered_edges);
		build_node_vec(node_vec, w, n, uncovered_edges);
	}
}

//Bit flip casuali su una soluzione. Potrebbe portare a soluzioni invalide, quindi va sistemata la soluzione o ripetuta la perturbazione.
void perturbation(bool* perturbed_solution, int num_elems_changed, int n){
	int rand_index;
	for(int i = 0; i < num_elems_changed; ++i){
		rand_index = random_index(n);
		perturbed_solution[rand_index] = !perturbed_solution[rand_index];
	}
}

//Criterio di accettazione di un nuovo ottimo locale per ILS
bool acceptance_criteria(int curr_best_weight, int prev_weight, int curr_weight, int eps, int iter_without_improvements){
	// Se la soluzione risulta migliore della precedente, accettala a prescindere
	if(curr_weight < prev_weight) return true;
	// Se risulta peggiore di un epsilon del migliore attuale e non stiamo deviando per troppo tempo, accettala
	if(curr_weight < curr_best_weight + eps && iter_without_improvements < MAX_ITER_WITHOUT_IMPROVEMENTS) return true;
	return false;
}

//Scheduling della variabile epsilon, che definisce la deviazione massima concessa rispetto al migliore attuale
int epsilon_scheduling(int eps, int curr_weight, int curr_best_weight, int min_weight, int min_eps){
	int new_eps;
	//Decrementi costanti o in base al miglioramento non sembrano avere chiare differenze
	//if (curr_weight < curr_best_weight) new_eps = eps - 2*min_weight;
	if (curr_weight < curr_best_weight) new_eps = eps - ((curr_best_weight - curr_weight) << 1);
	else new_eps = eps;
	if (new_eps > min_eps) return new_eps;
	else return min_eps;
}

//Scheduling del numero di bit da flippare in una perturbazione, deve diminuire col tempo per non generare soluzioni troppo lontane
int perturbation_scheduling(int num_elems_changed, float eps, float min_eps, float max_eps, int pert_min_changes, int pert_max_changes){
	float improvement = (float)(1.0f - (eps - min_eps)/(max_eps - min_eps));
	int new_num_elems_changed = pert_max_changes - (int)(improvement*(pert_max_changes - pert_min_changes));
	return new_num_elems_changed;
}

//Local search con sostituzioni o rimozioni di nodi dalla soluzione. Ritorna true se il numero di valutazioni della funzione obiettivo supera il massimo
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

//Separa nodi presenti e non presenti in due array di interi con gli indici dei nodi
void separate_nodes_array(bool* curr_solution, int* present_nodes, int* external_nodes, int& n_present, int& n_external, int n){
	fill_n(present_nodes, n_present, 0);
	fill_n(external_nodes, n_external, 0);
	n_present = 0;
	n_external = 0;
	for(int i = 0; i < n; ++i){
		if(curr_solution[i]){
			present_nodes[n_present] = i;
			n_present++;
		}
		else{
			external_nodes[n_external] = i;
			n_external++;
		}
	}
}

//Prendi un vicino casuale dalla lista di adiacenza
int Graph::random_neighbor(int node, bool* solution, bool neighbor_present){
	list<int> node_list = adj_lists[node];
	int list_length = node_list.size();
	int curr_index = -1;
	int rand_index;
	int curr_neighbor = 0;
	const int MAX_TRIES = 8;
	for(int tries = 0; tries < MAX_TRIES; ++tries){
		rand_index = random_index(list_length);
		for (list<int>::iterator it = node_list.begin(); it != node_list.end() && curr_index != rand_index; ++it){
			curr_neighbor = *it;
			curr_index++;
		}
		if(solution[curr_neighbor] == neighbor_present){
			return curr_neighbor;
		}
	}
	//Se non riesci a trovare un vicino buono, ritorna il nodo stesso da reinserire o rimuovere
	return node;
}

//Local search campionando i vicini ottenuti con lo scambio
bool local_search_stochastic(bool* curr_solution, int& curr_weight, Graph& graph, int n, int num_swaps, int& actual_swaps, int iter, int& num_obj_func_eval, bool VERBOSE_FLAG){

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

	int rand_node, rand_remove, rand_insert;

	bool node_removed = false;

	for(int ls_iter = 1; ls_iter <= MAX_LS_ITER && !ls_stuck; ++ls_iter){
		if(VERBOSE_FLAG) cout << "Local search numero " << iter << ", iterazione " << ls_iter << ", esploro il vicinato..." << endl;

		node_removed = false;

		//Esplorazione vicini rimuovendo un elemento dalla soluzione attuale
		for (int i = 0; i < n; ++i){
			if(curr_solution[i]){
				copy(curr_solution, curr_solution+n, neighbor);
				neighbor[i] = 0;
				//if (graph.valid_solution(neighbor)){
				if (graph.valid_solution_node(i, neighbor)){
					neighbor_weight = curr_weight - graph.get_weight(i);
					num_obj_func_eval++;
					if (neighbor_weight < curr_local_best_weight){
						node_removed = true;
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
		
		//Esplorazione vicini scambiando un elemento della soluzione attuale con uno non presente che migliora il peso
		for(int i = 0; i < num_swaps && !node_removed; ++i){
			rand_node = random_index(n);
			if(curr_solution[rand_node]){
				rand_remove = rand_node;
				rand_insert = graph.random_neighbor(rand_remove, curr_solution, 0);
			}
			else{
				rand_insert = rand_node;
				rand_remove = graph.random_neighbor(rand_insert, curr_solution, 1);
			}
			if(graph.get_weight(rand_insert) >= graph.get_weight(rand_remove)) continue;
			copy(curr_solution, curr_solution+n, neighbor);
			neighbor[rand_remove] = 0;
			neighbor[rand_insert] = 1;
			//if (graph.valid_solution(neighbor)){
			if (graph.valid_solution_node(rand_remove, neighbor)){
				neighbor_weight = curr_weight - graph.get_weight(rand_remove) + graph.get_weight(rand_insert);
				num_obj_func_eval++;
				if (neighbor_weight < curr_local_best_weight){
					actual_swaps++;
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

//Local search campionando i vicini ottenuti con una rimozione e con due scambi casuali
bool local_search_stochastic_removals(bool* curr_solution, int& curr_weight, Graph& graph, int n, int num_removals, int num_swaps, int& actual_swaps, int iter, int& num_obj_func_eval, bool VERBOSE_FLAG){

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

	int node_index = 0, rand_s1_remove = 0, rand_s1_insert = 0, rand_s2_remove = 0, rand_s2_insert = 0;

	bool node_removed = false;

	for(int ls_iter = 1; ls_iter <= MAX_LS_ITER && !ls_stuck; ++ls_iter){
		if(VERBOSE_FLAG) cout << "Local search numero " << iter << ", iterazione " << ls_iter << ", esploro il vicinato..." << endl;

		node_removed = false;

		//Esplorazione vicini rimuovendo un elemento dalla soluzione attuale
		for (int i = 0; i < num_removals; ++i){
			node_index = random_index(n);
			if(curr_solution[node_index]){
				copy(curr_solution, curr_solution+n, neighbor);
				neighbor[node_index] = 0;
				//if (graph.valid_solution(neighbor)){
				if (graph.valid_solution_node(node_index, neighbor)){
					neighbor_weight = curr_weight - graph.get_weight(node_index);
					num_obj_func_eval++;
					if (neighbor_weight < curr_local_best_weight){
						node_removed = true;
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
		
		//Esplorazione vicini scambiando un elemento della soluzione attuale con uno non presente che migliora il peso
		for(int i = 0; i < num_swaps && !node_removed; ++i){
			copy(curr_solution, curr_solution+n, neighbor);
			node_index = random_index(n);
			if(curr_solution[node_index]){
				rand_s1_remove = node_index;
				rand_s1_insert = graph.random_neighbor(rand_s1_remove, curr_solution, 0);
			}
			else{
				rand_s1_insert = node_index;
				rand_s1_remove = graph.random_neighbor(rand_s1_insert, curr_solution, 1);
			}
			node_index = random_index(n);
			if(curr_solution[node_index]){
				rand_s2_remove = node_index;
				rand_s1_insert = graph.random_neighbor(rand_s1_remove, curr_solution, 0);
			}
			else{
				rand_s2_insert = node_index;
				rand_s2_remove = graph.random_neighbor(rand_s2_insert, curr_solution, 1);
			}
			rand_s2_remove = random_index(n);
			//Se si è pescato un solo nodo da rimuovere, effettua un singolo swap con s1
			if (rand_s1_remove == rand_s2_remove){
				neighbor_weight = curr_weight - graph.get_weight(rand_s1_remove) + graph.get_weight(rand_s1_insert);
				neighbor[rand_s1_remove] = 0;
				neighbor[rand_s1_insert] = 0;
			}
			else{
				//Inizialmente sfida la fortuna rimuovendo due nodi diversi e inserendone uno solo
				neighbor_weight = curr_weight - graph.get_weight(rand_s1_remove) - graph.get_weight(rand_s2_remove) + graph.get_weight(rand_s1_insert);
				neighbor[rand_s1_remove] = 0;
				neighbor[rand_s2_remove] = 0;
				neighbor[rand_s1_insert] = 1;
				if(!graph.valid_solution_two_nodes(rand_s1_remove, rand_s2_remove, neighbor) && rand_s1_insert != rand_s2_insert){
					neighbor_weight = curr_weight - graph.get_weight(rand_s1_remove) - graph.get_weight(rand_s2_remove) + graph.get_weight(rand_s1_insert) + graph.get_weight(rand_s2_insert);
					neighbor[rand_s1_remove] = 0;
					neighbor[rand_s2_remove] = 0;
					neighbor[rand_s1_insert] = 1;
					neighbor[rand_s2_insert] = 1;
				}
			}
			//if(graph.get_weight(rand_v2) >= graph.get_weight(rand_v1)) continue;
			//if (graph.valid_solution(neighbor)){
			if (graph.valid_solution_two_nodes(rand_s1_remove, rand_s2_remove, neighbor)){
				num_obj_func_eval++;
				if (neighbor_weight < curr_local_best_weight){
					actual_swaps++;
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

void separate_nodes_vector(bool* curr_solution, vector<int>& present_nodes, vector<int>& external_nodes, int n){
	for(int i = 0; i < n; ++i){
		if(curr_solution[i]) present_nodes.push_back(i);
		else external_nodes.push_back(i);
	}
}

//Versione con vector per aggiungere e rimuovere elementi facilmente, realtà non sembra più efficiente (alla fine le iterazioni delle singole local search sono poche)
bool local_search_stochastic_vector(bool* curr_solution, int& curr_weight, Graph& graph, int n, int num_swaps, int iter, int& num_obj_func_eval, bool VERBOSE_FLAG){
	bool ls_stuck = false;
	curr_weight = graph.total_weight(curr_solution);
	num_obj_func_eval++;
	if(num_obj_func_eval >= MAX_OBJECTIVE_FUNCTION_EVAL) return true;
	int prev_weight = curr_weight;

	vector<int> present_nodes;
	vector<int> external_nodes;

	separate_nodes_vector(curr_solution, present_nodes, external_nodes, n);

	bool curr_local_best_solution[n];
	copy(curr_solution, curr_solution+n, curr_local_best_solution);
	int curr_local_best_weight = curr_weight;

	bool neighbor[n];
	int neighbor_weight = curr_weight;

	int rand_v1, rand_v2;

	int remove_index, insert_index;

	int node_to_remove = -1;
	int node_to_insert = -1;
	int curr_node;

	bool node_removed = false;

	if(VERBOSE_FLAG){
		cout << "Local search numero " << iter << ", soluzione attuale: " << endl;
		print_solution(curr_solution, n);
		cout << "Peso: " << curr_weight << endl;
	}

	for(int ls_iter = 1; ls_iter <= MAX_LS_ITER && !ls_stuck; ++ls_iter){
		if(VERBOSE_FLAG) cout << "Local search numero " << iter << ", iterazione " << ls_iter << ", esploro il vicinato..." << endl;
		node_removed = false;

		//Esplorazione vicini rimuovendo un elemento dalla soluzione attuale
		for (int i = 0; i < (int)present_nodes.size(); ++i){
			copy(curr_solution, curr_solution+n, neighbor);
			curr_node = present_nodes.at(i);
			neighbor[curr_node] = 0;
			//cout << "Provo a rimuovere " << curr_node << endl;
			if (graph.valid_solution(neighbor)){
				neighbor_weight = curr_weight - graph.get_weight(curr_node);
				num_obj_func_eval++;
				if (neighbor_weight < curr_local_best_weight){
					node_to_remove = i;
					node_removed = true;
					copy(neighbor, neighbor+n, curr_local_best_solution);
					curr_local_best_weight = neighbor_weight;
				}
				if(num_obj_func_eval >= MAX_OBJECTIVE_FUNCTION_EVAL){
					copy(curr_local_best_solution, curr_local_best_solution+n, curr_solution);
					curr_weight = curr_local_best_weight;
					return true;
				}
			}
			//cout << "Fine iterazione " << i << endl;
		}
		//int total_swaps = 0;
		//Esplorazione vicini scambiando un elemento della soluzione attuale con uno non presente che migliora il peso
		for(int i = 0; i < num_swaps && !node_removed && present_nodes.size() > 0 && external_nodes.size() > 0; ++i){
			remove_index = random_index(present_nodes.size());
			insert_index = random_index(external_nodes.size());
			rand_v1 = present_nodes.at(remove_index);
			rand_v2 = external_nodes.at(insert_index);
			if(graph.get_weight(rand_v2) >= graph.get_weight(rand_v1)) continue;
			//total_swaps++;
			copy(curr_solution, curr_solution+n, neighbor);
			neighbor[rand_v1] = 0;
			neighbor[rand_v2] = 1;
			if (graph.valid_solution(neighbor)){
				neighbor_weight = curr_weight - graph.get_weight(rand_v1) + graph.get_weight(rand_v2);
				num_obj_func_eval++;
				if (neighbor_weight < curr_local_best_weight){
					node_to_remove = remove_index;
					node_to_insert = insert_index;
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
		//cout << "Local search " << iter << ", iterazione " << ls_iter << ", numero scambi " << total_swaps << endl;
		//Sistemazione dei vector di nodi
		//Se il vicino ottimale è dato da una rimozione, togli il nodo da present_nodes e mettilo in external_nodes
		if(node_to_insert == -1 && node_to_remove >= 0){
			external_nodes.push_back(present_nodes.at(node_to_remove));
			present_nodes.erase(present_nodes.begin() + node_to_remove);
		}
		//Se il vicino ottimale è dato da uno scambio, salva il nodo da rimuovere per poi inserirlo in external_nodes
		else if(node_to_insert >= 0 && node_to_remove >= 0){
			curr_node = present_nodes.at(node_to_remove);
			present_nodes.erase(present_nodes.begin() + node_to_remove);
			present_nodes.push_back(external_nodes.at(node_to_insert));
			external_nodes.erase(external_nodes.begin() + node_to_insert);
			external_nodes.push_back(curr_node);
		}

		//Aggiornamento soluzione attuale con il migliore vicino
		copy(curr_local_best_solution, curr_local_best_solution+n, curr_solution);
		curr_weight = curr_local_best_weight;

		//cout << "Soluzione attuale: " << endl;
		//print_solution(curr_solution, n);
		//if(VERBOSE_FLAG) cout << "Peso soluzione attuale: " << curr_weight << endl;

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