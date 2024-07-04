#ifndef TRANSITION_HH
#define TRANSITION_HH

#include <iostream>
#include <vector>
#include "State.hh"
using namespace std;

class Transition {
public:
    State* initialState;
    State* finalState;
    string transitionSymbol;

    Transition(string transitionSymbol, int& currentState) {
        this->transitionSymbol = transitionSymbol;
        this->initialState = new State(currentState++);
        this->finalState = new State(currentState++);
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
          cout << "Deleted initialState from Transition" << endl;
          initialState = nullptr;
        } 
        if (finalState != nullptr) {
          delete finalState;
          cout << "Deleted finalState from Transition" << endl;
          finalState = nullptr;
        }
    }
};

#endif