#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <omp.h>
using namespace std;

// Clase DFA optimizada
class DFA {
private:
    int numStates;
    int initialState;
    vector<int> acceptingStates;
    string alphabet;
    map<pair<int, char>, int> transitions;
    
public:
    DFA(int states, int initial, vector<int> accepting, string alpha);
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

// Estructura de resultado paralelo
struct DFAResult {
    bool accepted;
    double computationTime;
    int finalState;
    vector<int> processStates;
    int numMatches;

    DFAResult() : accepted(false), computationTime(0), numMatches(0), finalState(-1) {}
};

struct DFAOptimizedResult {
    bool accepted;
    double preComputeJumpTime;
    double computationTime;
    int finalState;
    vector<int> processStates;
    int numMatches;

    DFAOptimizedResult() : accepted(false), computationTime(0), preComputeJumpTime(0),
                        numMatches(0), finalState(-1) {}
};

struct BlockResult {
    std::vector<int> route;
    int matches;
    int startState;
    int endState;
};

struct BlockInfo {
    int endState;
    int matches;
};

// FUNCIONES DE OPTIMIZACION CON OMP
DFAResult processOMPParallel(const DFA& dfa, const std::string &T, int num_threads);

std::vector<std::vector<BlockInfo>> precomputeJump(const DFA& afd, const std::string &T, int num_threads, int chunk);
std::vector<int> computeStartStates(const std::vector<std::vector<BlockInfo>> &jumpTable);

DFAOptimizedResult runOptimizedOMPParallel(const DFA& afd, const std::string &T, int num_threads);
