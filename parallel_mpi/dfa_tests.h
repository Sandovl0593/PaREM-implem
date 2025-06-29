#ifndef DFA_TESTS_H
#define DFA_TESTS_H

#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <random>
#include "optimized_parallel_dfa.h"
using namespace std;

// Estructura para almacenar configuración de test
struct DFATest {
    string name;
    int numStates;
    int initialState;
    vector<int> finalStates;
    string alphabet;
    string input;
    bool expectedResult;
    string description;
};

// Estructura para almacenar resultados de tests
struct TestResult {
    string testName;
    bool passed;
    bool expectedResult;
    bool actualResult;
    long long executionTime; // en microsegundos
    string description;
};

// Declaraciones de funciones
vector<TestResult> runSequentialTestsWithResults(const OptimizedDFA& dfa, const vector<DFATest>& tests);
void printTestSummary(const vector<TestResult>& results);
void runSequentialTests(const OptimizedDFA& dfa, const vector<DFATest>& tests);
vector<DFATest> createBinaryTests();
vector<DFATest> createComplexTests();
vector<DFATest> createStressTests();
vector<DFATest> createABCTests();
vector<DFATest> createMultiFinalTests();
void runDetailedAnalysis(const OptimizedDFA& dfa, const string& input, const string& testName);

#endif // DFA_TESTS_H 