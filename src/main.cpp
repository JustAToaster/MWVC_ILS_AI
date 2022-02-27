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

//Massimo numero di iterazioni di ILS (tenendo comunque conto del numero di valutazioni massime della funzione obiettivo)
#define MAX_ITER 32768

int main(int argc, char *argv[]){
	//Processa l'input da riga di comando
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
	string graph_name = input.substr(input.find_last_of("/\\") + 1);
	graph_name = graph_name.substr(0, input.length()-4);
	graph_name = graph_name.substr(0, graph_name.length()-4);
	write.open(graph_name + "_output.txt");
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
	
	//Fixed seed
	//srand(640);

	//Processa il grafo in input
	int n;
	int weight;
	bool curr_edge;
	Graph graph;
	vector< pair<int, int> > edges;
	
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
			//Salva gli archi in entrambe le direzioni per la validazione ottimizzata, altrimenti aggiungere j >= i
			if(curr_edge){
				edges.push_back(make_pair(i, j));
			}
			//cout << curr_edge << "\t";
		}
		//cout << endl;
	}

	//Scrivi i valori per il benchmark
	//write_benchmark << "num_nodes,num_edges,min_eps,max_eps,pert_min_changes,pert_max_changes,num_swaps,iter,weight" << endl;
	//write_benchmark << "iter,weight,obj_fun_eval,eps,num_elems_changed,actual_swaps" << endl;

	graph.build_adj_lists(edges);
	graph.compute_degree_weight_ratio();
	//graph.print_edges();

	bool* curr_solution = new bool[n];

	//Soluzione banale, tutti i nodi. In alcuni casi potrebbe essere un miglior punto di partenza di quella greedy.
	//fill_n(curr_solution, n, 1);

	//Soluzione greedy
	fill_n(curr_solution, n, 0);
	graph.greedy_heuristic_queue(curr_solution);

	//Soluzione migliore attuale, da stampare e scrivere alla fine
	bool* curr_best_solution = new bool[n];
	copy(curr_solution, curr_solution+n, curr_best_solution);
	int curr_best_weight = graph.total_weight(curr_solution);

	//Peso della soluzione attuale processata dalle local search
	int curr_weight = curr_best_weight;
	
	cout << "Soluzione iniziale: " << endl;
	print_solution(curr_solution, n);
	cout << "Peso iniziale: " << curr_weight << endl;
	
	//Peso della soluzione dell'iterazione di ILS precedente
	int prev_ils_weight = curr_weight;

	//Numero di iterazioni di ILS da quando non è stata migliorata la funzione obiettivo
	int iter_without_improvements = 0;

	// Il minimo peso di un nodo nel grafo. Usato per lo scheduling di epsilon.
	int min_weight = graph.min_weight();

	//Contatore valutazioni della funzione obiettivo
	int num_obj_func_eval = 0;

	//Numero di valutazioni per raggiungere l'ultima soluzione ottimale
	int objfun_optimal = 2e4;

	//Vero se il numero di valutazioni della funzione obiettivo richieste supera il massimo
	bool max_eval_reached = false;
	
	//Valori di partenza e finale di epsilon, il parametro per lo scheduling del distacco concesso per il criterio di accettazione
	int max_eps = curr_weight >> 3;
	int min_eps = curr_weight >> 5;
	if (min_eps < 2*min_weight) min_eps = 2*min_weight;
	if (max_eps < min_eps) max_eps = 2*min_eps;
	int eps = max_eps;

	//Valori di partenza e finale del numero di bit flip in una perturbazione
	int pert_max_changes = n >> 4;
	int pert_min_changes = n >> 8;
	if(pert_min_changes == 0) pert_min_changes = 1;
	if(pert_max_changes <= pert_min_changes) pert_max_changes = pert_min_changes;
	int num_elems_changed = pert_max_changes;

	//Numero di sostituzioni massime in una iterazione di local search
	int num_swaps = n >> 3;
	if (num_swaps < 8) num_swaps = 8;

	//Numero di rimozioni massime se si effettua il campionamento di queste ultime
	//int num_removals = num_swaps << 1;

	//Contatore per il numero di sostituzioni totali effettuate
	int actual_swaps = 0;

	cout << "Numero sostituzioni per iterazione della local search: " << num_swaps << endl;
	cout << "Valore iniziale di epsilon: " << eps << ", valore minimo: " << min_eps << endl;
	cout << "Numero di cambiamenti iniziali in una perturbazione: " << num_elems_changed << endl;
	cout << "Numero di cambiamenti minimi in una perturbazione: " << pert_min_changes << endl << endl;

	//Numero di iterazioni finali di ILS
	int final_iterations = MAX_ITER;

	//Per salvare il runtime
	clock_t start_opt, end_opt;
	start_opt = clock();

	//Ciclo di iterated local search
	for(int iter = 1; iter <= MAX_ITER; ++iter){
		
		//write_benchmark << iter << "," << curr_weight << "," << num_obj_func_eval << "," << eps << "," << num_elems_changed << "," << actual_swaps << endl;

		//Calcola l'ottimo locale a partire dalla soluzione attuale
		
		//Local search esaustiva: troppo costosa
		//max_eval_reached = local_search(curr_solution, curr_weight, graph, n, iter, num_obj_func_eval, VERBOSE_FLAG);
		
		//Local search campionando gli scambi
		max_eval_reached = local_search_stochastic(curr_solution, curr_weight, graph, n, num_swaps, actual_swaps, iter, num_obj_func_eval, VERBOSE_FLAG);
		
		//Local search campionando due scambi totalmente casuali e rimozioni. Necessario campionare rimozioni se il numero di nodi è eccessivamente grande o si parte dalla soluzione banale
		//max_eval_reached = local_search_stochastic_removals(curr_solution, curr_weight, graph, n, num_removals, num_swaps, actual_swaps, iter, num_obj_func_eval, VERBOSE_FLAG);
		
		if(max_eval_reached){
			cout << "Esaurito il numero massimo di valutazioni della funzione obiettivo. Uscita..." << endl;
			if(curr_weight < curr_best_weight){
				curr_best_weight = curr_weight;
				objfun_optimal = num_obj_func_eval;
				copy(curr_solution, curr_solution+n, curr_best_solution);
			}
			final_iterations = iter;
			break;
		}

		// Terminata la local search, aggiorno epsilon per poi verificare il criterio di accettazione
		eps = epsilon_scheduling(eps, curr_weight, curr_best_weight, min_weight, min_eps);

		if(VERBOSE_FLAG) cout << "Iterazione " << iter << ", valore di epsilon: " << eps  << ", numero di bit flip con una perturbazione: " << num_elems_changed << endl;
		if(acceptance_criteria(curr_best_weight, prev_ils_weight, curr_weight, eps, iter_without_improvements)){
			//Se c'è un miglioramento, aggiorna la soluzione attuale e azzera il numero di iterazioni senza miglioramento
			if(curr_weight < curr_best_weight){
				curr_best_weight = curr_weight;
				objfun_optimal = num_obj_func_eval;
				copy(curr_solution, curr_solution+n, curr_best_solution);
				iter_without_improvements = 0;
			}
			else{
				iter_without_improvements++;
			}
			prev_ils_weight = curr_weight;
			perturbation(curr_solution, num_elems_changed, n);
			num_elems_changed = perturbation_scheduling(num_elems_changed, eps, min_eps, max_eps, pert_min_changes, pert_max_changes);
			if (!graph.valid_solution(curr_solution)) graph.greedy_heuristic_vec_prob(curr_solution);
		}
		else{
			if(VERBOSE_FLAG) cout << "Soluzione non accettata dal criterio di accettazione. Ritorno all'ottimo precedente." << endl;
			iter_without_improvements++;
			//Se non si dà un massimo di valutazioni della funzione obiettivo, si potrebbe decidere di uscire da ILS dopo un paio di tentativi di miglioramento
			/*
			if(iter_without_improvements >= TOTAL_MAX_ITER_WITHOUT_IMPROVEMENTS){
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
	cout << "Numero valutazioni della funzione obiettivo: " << objfun_optimal << ", totali: " << num_obj_func_eval << endl;
	cout << "Soluzione migliore trovata: " << endl;
	print_solution(curr_best_solution, n);
	cout << "Peso: " << curr_best_weight << endl;

	cout << "La soluzione è valida? " << graph.valid_solution(curr_best_solution) << endl;

	//Scrivi la soluzione su file
	for(int i = 0; i < n; ++i){
		if(curr_best_solution[i]) write << i << " ";
	}
	write << endl << curr_best_weight << endl;

	write_benchmark << graph_name << "," << n << "," << edges.size()/2 << "," << curr_best_weight << "," << objfun_optimal << "," << actual_swaps << "," << runtime_ms << endl;

	delete curr_solution;
	delete curr_best_solution;

	read.close();
	write.close();
	write_benchmark.close();
	return 0;
}
