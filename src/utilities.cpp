#include <fstream>
#include <cstdlib>
#include <iostream>
#include <string>
#include <cmath>
#include <utility>
#include <time.h>
#include<float.h>
#include<random>

using namespace std;

float random_01(){
    random_device rd;
    default_random_engine re(rd());
    uniform_real_distribution<float> uniform_distribution(0.0f, 1.0f);
    return uniform_distribution(re);
}

int random_index(int size){
    random_device rd;
    default_random_engine re(rd());
    uniform_int_distribution<int> uniform_distribution(0, size-1);
    return uniform_distribution(re);
}

int random_normal_index(int mean, int std){
    random_device rd;
    default_random_engine re(rd());
    normal_distribution<float> normal_distribution((float)mean, (float)std);
    return (int)abs(normal_distribution(re));
}

/*
float random_01(){
	return (float)rand()/RAND_MAX;
}

int random_index(int size){
	return rand() % size;
}
*/
float probability_function(float x, float bias, char function){
	//Sigmoide
	if (function == 's') return 1/(1+exp(-x+bias));
	//Tangente iperbolica
	if (function == 't') return tanh(x+bias);
	//Esponenziale semplice
	if (function == 'e') return 1-exp(-x+bias);
	return x;
}

void print_solution(bool* solution, int n){
	for(int i = 0; i < n; ++i){
		if(solution[i]) cout << i << "|";
	}
	cout << endl;
}

void print_array(int* array, int n){
	for(int i = 0; i < n; ++i){
		cout << array[i] << "|";
	}
	cout << endl;
}