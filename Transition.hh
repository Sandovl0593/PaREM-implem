#ifndef TRANSITION_HH
#define TRANSITION_HH

#include <iostream>
#include <vector>
using namespace std;

class State {
public:
  int stateId;

  State(int stateId) : stateId(stateId) {
  }

  string toString() {
    return to_string(this->stateId);
  }

  int getStateId() {
    return this->stateId;
  }

  ~State() {
  }
};

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
    }
};

#endif