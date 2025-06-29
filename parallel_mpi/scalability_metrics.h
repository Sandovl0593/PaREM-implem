#pragma once

#include <vector>
#include <string>
#include <chrono>
#include "optimized_parallel_dfa.h"

using namespace std;

// Función para generar string de prueba
string generateTestInput(int size);

// Función para crear DFA de prueba
OptimizedDFA createTestDFA();

// Función para medir tiempo secuencial
double measureSequentialTime(const OptimizedDFA& dfa, const string& input); 