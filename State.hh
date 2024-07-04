#ifndef STATE_HH
#define STATE_HH
#include <iostream>
using namespace std;

class State {
public:
  vector<State*> previousStates;
  vector<State*> nextStates;
  int stateId;

  State(int stateId) : stateId(stateId) {
  }

  State(int stateId, vector<State*> previousState, vector<State*> nextState)
    : stateId(stateId), previousStates(previousState), nextStates(nextState) {
  }

  State(int stateId, bool dfa) : stateId(stateId) {}

  void addPreviousState(State* previousState) {
    this->previousStates.push_back(previousState);
  }

  void addNextState(State* nextState) {
    this->nextStates.push_back(nextState);
  }

  vector<State*> getPreviousStates() {
    return this->previousStates;
  }

  vector<State*> getNextStates() {
    return this->nextStates;
  }

  string toString() {
    return to_string(this->stateId);
  }

  int getStateId() {
    return this->stateId;
  }

  ~State() {
    for (auto state : previousStates) {
      if (state != nullptr) {
        delete state;
        cout << "delete prev state from State" << endl;
        state = nullptr;
      } 
    }

    for (auto state : nextStates) {
      if (state != nullptr) {
        delete state;
        cout << "delete next state from State" << endl;
        state = nullptr;
      } 
    }
  }
};

#endif