#include "dfa.h"

// ================= IMPLEMENTACIÓN DEL DFA =================
DFA::DFA(int states, int initial, vector<int> accepting, string alpha) 
    : numStates(states), initialState(initial), acceptingStates(accepting), alphabet(alpha) {}

void DFA::addTransition(int from, char symbol, int to) {
    transitions[{from, symbol}] = to;
}

bool DFA::run(const string& input) const {
    int currentState = initialState;
    
    for (char symbol : input) {
        currentState = getNextState(currentState, symbol);
        if (currentState == -1) return false; // Estado inválido
    }
    
    return isAccepting(currentState);
}

int DFA::getNextState(int currentState, char symbol) const {
    auto it = transitions.find({currentState, symbol});
    return (it != transitions.end()) ? it->second : -1;
}

bool DFA::isAccepting(int state) const {
    return find(acceptingStates.begin(), acceptingStates.end(), state) != acceptingStates.end();
}
