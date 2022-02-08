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

#define MAX_ITER 4096

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
	read.open(input);
	if(read.fail())
	{
		cout<< "Apertura input fallita!\n";
		exit(1);
	}
	write.open(input.substr(0, input.length()-4)+"_output.txt");
	if(write.fail())
	{
		cout << "Non sono riuscito a scrivere l'output.txt!'\n";
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
			if(j >= i && curr_edge){
				edges.push_back(make_pair(i, j));
			}
			//cout << curr_edge << "\t";
		}
		//cout << endl;
	}

	graph.add_edges(edges);
	//graph.print_edges();

	bool* curr_solution = new bool[n];

	//Soluzione banale, tutti i nodi. In alcuni casi potrebbe essere un miglior punto di partenza di quella greedy.
	//fill_n(curr_solution, n, 1);

	//Soluzione greedy
	fill_n(curr_solution, n, 0);
	graph.greedy_solution(curr_solution);

	bool* curr_best_solution = new bool[n];
	fill_n(curr_best_solution, n, 1);

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
	
	int eps = curr_weight >> 4;
	if (eps < 0) eps = 0;

	int pert_min_changes = n >> 8;
	if(pert_min_changes == 0) pert_min_changes = 1;

	int num_elems_changed = n >> 6;
	if(num_elems_changed == 0) num_elems_changed = pert_min_changes;

	cout << "Valore iniziale di epsilon: " << eps << ", peso minimo: " << min_weight << endl;
	cout << "Numero di cambiamenti iniziali in una perturbazione: " << num_elems_changed << endl;
	cout << "Numero di cambiamenti minimi in una perturbazione: " << pert_min_changes << endl << endl;

	clock_t start_opt, end_opt;
	start_opt = clock();

	//Ciclo di iterated local search
	for(int iter = 1; iter <= MAX_ITER; ++iter){

		if(VERBOSE_FLAG){
			cout << "Local search numero " << iter << ", soluzione attuale: " << endl;
			print_solution(curr_solution, n);
			cout << "Peso: " << curr_weight << endl;
		}

		//Calcola l'ottimo locale a partire dalla soluzione attuale
		max_eval_reached = local_search(curr_solution, curr_weight, graph, n, iter, num_obj_func_eval, VERBOSE_FLAG);

		if(max_eval_reached){
			cout << "Esaurito il numero massimo di valutazioni della funzione obiettivo. Uscita..." << endl;
			break;
		}

		// Terminata la local search, verifico il criterio di accettazione
		eps = epsilon_scheduling(eps, curr_weight, curr_best_weight, min_weight, iter);

		if (VERBOSE_FLAG) cout << "Iterazione " << iter << ", valore di epsilon: " << eps  << ", numero di bit flip con una perturbazione: " << num_elems_changed << endl;
		if(acceptance_criteria(curr_best_weight, prev_ils_weight, curr_weight, eps, iter_without_improvements)){
			if(curr_weight < curr_best_weight){
				curr_best_weight = curr_weight;
				copy(curr_solution, curr_solution+n, curr_best_solution);
				iter_without_improvements = 0;
			}
			else{
				iter_without_improvements++;
			}
			prev_ils_weight = curr_weight;
			perturbation(curr_solution, num_elems_changed, n);
			num_elems_changed = perturbation_scheduling(num_elems_changed, pert_min_changes, iter);
			if (!graph.valid_solution(curr_solution)) graph.fix_invalid_solution(curr_solution);
		}
		else{
			if(VERBOSE_FLAG) cout << "Soluzione non accettata dal criterio di accettazione. Ritorno all'ottimo precedente." << endl;
			iter_without_improvements++;
			if(iter_without_improvements >= MAX_ITER_WITHOUT_IMPROVEMENTS){
				cout << "Esaurito il numero di tentativi per migliorare la soluzione. Uscita da ILS..." << endl;
				break;
			}
			copy(curr_best_solution, curr_best_solution+n, curr_solution);
			curr_weight = curr_best_weight;
		}
		
	}

	end_opt = clock();
	double runtime_ms = (end_opt - start_opt)*1.0e3/CLOCKS_PER_SEC;
	cout << "Terminato ciclo di ottimizzazione in " << runtime_ms << " ms." << endl;

	cout << "Numero valutazioni della funzione obiettivo: " << num_obj_func_eval << endl;
	cout << "Soluzione migliore trovata: " << endl;
	print_solution(curr_best_solution, n);
	cout << "Peso: " << curr_best_weight << endl;

	//Scrivi la soluzione su file
	for(int i = 0; i < n; ++i){
		if(curr_best_solution[i]) write << i << " ";
	}
	write << endl << curr_best_weight << endl;

	delete curr_solution;
	delete curr_best_solution;

	read.close();
	write.close();
	return 0;
}
