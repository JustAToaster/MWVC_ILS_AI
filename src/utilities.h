#include <fstream>
#include <cstdlib>
#include <iostream>
#include <string>
#include <cmath>

using namespace std;

float random_01();

int random_index(int size);

int random_normal_index(int mean, int std);

void print_solution(bool* solution, int n);

void print_array(int* array, int n);

float probability_function(float x, float threshold, char function);