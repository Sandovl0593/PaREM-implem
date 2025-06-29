#pragma once
#include <iostream>
#include <mpi.h>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
using namespace std;

// Clase DFA optimizada
class OptimizedDFA {
private:
    int numStates;
    int initialState;
    vector<int> acceptingStates;
    string alphabet;
    map<pair<int, char>, int> transitions;
    
public:
    OptimizedDFA(int states, int initial, vector<int> accepting, string alpha);
    void addTransition(int from, char symbol, int to);
    bool run(const string& input) const;
    int getNextState(int currentState, char symbol) const;
    bool isAccepting(int state) const;
    
    // Métodos públicos para compatibilidad con parallel_dfa.cpp
    int getInitialState() const { return initialState; }
    int getNumStates() const { return numStates; }
    string getAlphabet() const { return alphabet; }
    int nextState(int currentState, char symbol) const { return getNextState(currentState, symbol); }
    vector<int> getAcceptingStates() const { return acceptingStates; }
    
    // Método para tracking detallado
    pair<bool, vector<int>> run_with_tracking(const string& input) const {
        vector<int> stateSequence;
        int currentState = initialState;
        stateSequence.push_back(currentState);
        
        for (char symbol : input) {
            currentState = getNextState(currentState, symbol);
            if (currentState == -1) {
                return {false, stateSequence};
            }
            stateSequence.push_back(currentState);
        }
        
        return {isAccepting(currentState), stateSequence};
    }
};

// Funciones de comunicación optimizada
void setupOptimizedBuffers(int bufferSize);
void cleanupOptimizedBuffers();

// Procesamiento con diferentes tipos de comunicación
void processWithBlockingCommunication(const string& input, int numProcs, int myRank);
void processWithAsyncCommunication(const string& input, int numProcs, int myRank);

// Estructura de resultado paralelo optimizado
struct OptimizedParallelDFAResult {
    bool accepted;
    double computation_time;
    double communication_time;
    int final_state;
    vector<int> process_states;
    
    OptimizedParallelDFAResult() : accepted(false), computation_time(0), communication_time(0), final_state(-1) {}
};

// Función principal de procesamiento paralelo optimizado
OptimizedParallelDFAResult processOptimizedDFAParallel(const OptimizedDFA& dfa, const string& input, 
                                                      int numProcs, int myRank, bool useAsyncComm = false);

// Funciones de testing de comunicación
void testPointToPointCommunication(int numProcs, int myRank);
void testCollectiveCommunication(int numProcs, int myRank);
void testDerivedDataTypes(int numProcs, int myRank);
void testOptimizedBuffers(int numProcs, int myRank); 