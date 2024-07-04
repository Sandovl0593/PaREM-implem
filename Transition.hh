#ifndef TRANSITION_HH
#define TRANSITION_HH

#include <iostream>
#include <vector>
#include "State.hh"
#include "AFN.hh"
using namespace std;

class Transition {
private:
    State* initialState;
    State* finalState;
    string transitionSymbol;

public:
    Transition(string transitionSymbol) {
        this->transitionSymbol = transitionSymbol;
        this->initialState = new State(stateCount);
        this->finalState = new State(stateCount);
    }

    Transition(string transitionSymbol, State* initialState, State* finalState) {
        this->transitionSymbol = transitionSymbol;
        this->initialState = initialState;
        this->finalState = finalState;
    }

    State* getInitialState() {
        return this->initialState;
    }

    State* getFinalState() {
        return this->finalState;
    }

    string toString() {
        return initialState->toString() + " - " + transitionSymbol + " - " + finalState->toString();
    }

    string getTransitionSymbol() { return this->transitionSymbol; }

    ~Transition() {
        if (initialState != nullptr) {
          delete initialState;
          initialState = nullptr;
        } 
        if (finalState != nullptr) {
          delete finalState;
          finalState = nullptr;
        }
    }
};

#endif