#include <fstream>
#include <cstdlib>
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include <climits>
#include <vector>
#include <utility>
#include <time.h>
#include <algorithm>

#include "optimization.h"

using namespace std;

#define MAX_ITER 16384
//#define MAX_ITER 2

int main(int argc, char *argv[]){
	if(argc == 1){
		cerr << "Specificare file di input!" << endl;
		cerr << "Utilizzo: " << argv[0] << " input [VERBOSE_FLAG]" << endl;
		cerr << "Esempio: ./" << argv[0] << " ../wvcp-instances/vc_200_750_02.txt 0" << endl;
		exit(1);
	}
	string input = argv[1];
	bool VERBOSE_FLAG = true;
	if(argc >= 3) VERBOSE_FLAG = atoi(argv[2]);
	ifstream read;
	ofstream write;
	ofstream write_benchmark;
	read.open(input);
	if(read.fail())
	{
		cout << "Apertura input fallita!\n";
		exit(1);
	}
	//write.open(input.substr(0, input.length()-4)+"_output.txt");
	write.open("output.txt");
	if(write.fail())
	{
		cout << "Non sono riuscito a scrivere l'output.txt!'\n";
		exit(1);
	}

	write_benchmark.open("./benchmark.csv", ios_base::app);
	if(write_benchmark.fail())
	{
		cout << "Non sono riuscito a scrivere il benchmark!'\n";
		exit(1);
	}

	srand(640);

	int n;
	int weight;
	bool curr_edge;
	Graph graph;
	vector< pair<int, int> > edges;
	
	// Lettura input
	read >> n;
	cout << "N: " << n << endl;
	graph.initialize_graph_nodes(n);
	for(int i = 0; i < n; ++i){
		read >> weight;
		graph.add_node_weight(i, weight);
	}
	for(int i = 0; i < n; ++i){
		for(int j = 0; j < n; ++j){
			read >> curr_edge;
			//if(j >= i && curr_edge){
			if(curr_edge){
				edges.push_back(make_pair(i, j));
			}
			//cout << curr_edge << "\t";
		}
		//cout << endl;
	}

	//Scrivi i valori per il benchmark
	//write_benchmark << "num_nodes,num_edges,min_eps,max_eps,pert_min_changes,pert_max_changes,num_swaps,iter,weight," << endl;
	//write_benchmark << "iter,weight,obj_fun_eval,eps,num_elems_changed,actual_swaps" << endl;

	//graph.add_edges(edges);
	graph.build_adj_lists(edges);
	graph.compute_degree_weight_ratio();
	//graph.print_edges();

	bool* curr_solution = new bool[n];

	//Soluzione banale, tutti i nodi. In alcuni casi potrebbe essere un miglior punto di partenza di quella greedy.
	//fill_n(curr_solution, n, 1);
	//graph.greedy_heuristic_queue_prob(curr_solution);
	//graph.compute_degree_weight_ratio();

	//Soluzione greedy
	fill_n(curr_solution, n, 0);
	graph.greedy_heuristic_queue(curr_solution);

	bool* curr_best_solution = new bool[n];
	copy(curr_solution, curr_solution+n, curr_best_solution);

	int curr_weight = graph.total_weight(curr_solution);
	print_solution(curr_solution, n);
	cout << "Peso iniziale: " << curr_weight << endl;
	
	int prev_ils_weight = curr_weight;
	int curr_best_weight = curr_weight;

	int iter_without_improvements = 0;

	// Il minimo peso di un nodo nel grafo. Usato per lo scheduling di epsilon.
	int min_weight = graph.min_weight();

	//Numero di valutazioni della funzione obiettivo
	int num_obj_func_eval = 0;

	//Vero se il numero di valutazioni della funzione obiettivo supera il massimo
	bool max_eval_reached = false;
	
	int max_eps = curr_weight >> 3;
	int min_eps = curr_weight >> 5;
	if (min_eps < 2*min_weight) min_eps = 2*min_weight;
	if (max_eps < min_eps) max_eps = 2*min_eps;

	int eps = max_eps;

	int pert_min_changes = n >> 8;
	if(pert_min_changes == 0) pert_min_changes = 1;

	int pert_max_changes = n >> 4;
	if(pert_max_changes <= pert_min_changes) pert_max_changes = pert_min_changes;

	int num_elems_changed = pert_max_changes;

	int num_swaps = n >> 3;
	if (num_swaps < 8) num_swaps = 8;

	cout << "Numero sostituzioni per iterazione della local search: " << num_swaps << endl;
	cout << "Valore iniziale di epsilon: " << eps << ", valore minimo: " << min_eps << endl;
	cout << "Numero di cambiamenti iniziali in una perturbazione: " << num_elems_changed << endl;
	cout << "Numero di cambiamenti minimi in una perturbazione: " << pert_min_changes << endl << endl;

	int final_iterations = MAX_ITER;


	clock_t start_opt, end_opt;
	start_opt = clock();

	int actual_swaps = 0;

	//int objfun_optimal = 2e4;

	//Ciclo di iterated local search
	for(int iter = 1; iter <= MAX_ITER; ++iter){
		
		//write_benchmark << iter << "," << curr_weight << "," << num_obj_func_eval << "," << eps << "," << num_elems_changed << "," << actual_swaps << endl;
		//actual_swaps = 0;
		//Calcola l'ottimo locale a partire dalla soluzione attuale
		
		//Local search esaustiva: troppo costosa
		//max_eval_reached = local_search(curr_solution, curr_weight, graph, n, iter, num_obj_func_eval, VERBOSE_FLAG);
		
		//Local search campionando gli scambi
		max_eval_reached = local_search_stochastic(curr_solution, curr_weight, graph, n, num_swaps, actual_swaps, iter, num_obj_func_eval, VERBOSE_FLAG);
		
		//Local search campionando scambi e rimozioni. Necessaria se il numero di nodi è eccessivamente grande o si parte dalla soluzione banale, altrimenti è meglio quella che campiona solo gli scambi
		//max_eval_reached = local_search_stochastic_removals(curr_solution, curr_weight, graph, n, 200, num_swaps, actual_swaps, iter, num_obj_func_eval, VERBOSE_FLAG);
		
		if(max_eval_reached){
			cout << "Esaurito il numero massimo di valutazioni della funzione obiettivo. Uscita..." << endl;
			if(curr_weight < curr_best_weight){
				curr_best_weight = curr_weight;
				//objfun_optimal = num_obj_func_eval;
				copy(curr_solution, curr_solution+n, curr_best_solution);
			}
			final_iterations = iter;
			break;
		}

		// Terminata la local search, verifico il criterio di accettazione
		eps = epsilon_scheduling(eps, curr_weight, curr_best_weight, min_weight, min_eps);

		if(VERBOSE_FLAG) cout << "Iterazione " << iter << ", valore di epsilon: " << eps  << ", numero di bit flip con una perturbazione: " << num_elems_changed << endl;
		if(acceptance_criteria(curr_best_weight, prev_ils_weight, curr_weight, eps, iter_without_improvements)){
			if(curr_weight < curr_best_weight){
				curr_best_weight = curr_weight;
				copy(curr_solution, curr_solution+n, curr_best_solution);
				iter_without_improvements = 0;
				//objfun_optimal = num_obj_func_eval;
			}
			else{
				iter_without_improvements++;
			}
			prev_ils_weight = curr_weight;
			perturbation(curr_solution, num_elems_changed, n);
			num_elems_changed = perturbation_scheduling(num_elems_changed, eps, min_eps, max_eps, pert_min_changes, pert_max_changes);
			if (!graph.valid_solution(curr_solution)) graph.greedy_heuristic_queue_prob(curr_solution);
		}
		else{
			if(VERBOSE_FLAG) cout << "Soluzione non accettata dal criterio di accettazione. Ritorno all'ottimo precedente." << endl;
			iter_without_improvements++;
			/*
			if(iter_without_improvements >= MAX_ITER_WITHOUT_IMPROVEMENTS){
				cout << "Esaurito il numero di tentativi per migliorare la soluzione. Uscita da ILS..." << endl;
				break;
			}*/
			copy(curr_best_solution, curr_best_solution+n, curr_solution);
			curr_weight = curr_best_weight;
		}
	}

	end_opt = clock();
	double runtime_ms = (end_opt - start_opt)*1.0e3/CLOCKS_PER_SEC;
	cout << "Terminato ciclo di ottimizzazione con " << final_iterations << " iterazioni di ILS in " << runtime_ms << " ms." << endl;

	cout << "Valore finale di epsilon: " << eps << ", numero cambiamenti in una perturbazione: " << num_elems_changed << endl;
	cout << "Scambi effettuati: " << actual_swaps << endl;
	cout << "Numero valutazioni della funzione obiettivo: " << num_obj_func_eval << endl;
	cout << "Soluzione migliore trovata: " << endl;
	print_solution(curr_best_solution, n);
	cout << "Peso: " << curr_best_weight << endl;

	cout << "La soluzione è valida? " << graph.valid_solution(curr_best_solution) << endl;
	cout << "Peso calcolato: " << graph.total_weight(curr_best_solution) << endl;

	//Scrivi la soluzione su file
	for(int i = 0; i < n; ++i){
		if(curr_best_solution[i]) write << i << " ";
	}
	write << endl << curr_best_weight << endl;

	//write_benchmark << n << "," << edges.size()/2 << "," << num_obj_func_eval << "," << actual_swaps << endl;
	//write_benchmark << n << "," << edges.size()/2 << "," << objfun_optimal << endl;

	delete curr_solution;
	delete curr_best_solution;

	read.close();
	write.close();
	write_benchmark.close();
	return 0;
}
