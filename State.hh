#ifndef STATE_HH
#define STATE_HH
#include "AFN.hh"
using namespace std;

class State {
  vector<State*> previousStates;
  vector<State*> nextStates;
  int stateId;
  bool isInitial;
  bool isFinal;

public:
  State(int stateId) : stateId(stateId), isInitial(false), isFinal(false) {
    // AFN::stateCount++; // sets a unique ID which comes from AFN's total stateCount
  }

  State(int stateId, vector<State*> previousState, vector<State*> nextState)
    : stateId(stateId), previousStates(previousState), nextStates(nextState), isInitial(false), isFinal(false) {
    // AFN::stateCount++;
  }

  State(int stateId, bool dfa) : stateId(stateId), isInitial(false), isFinal(false) {}

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

  void setInitial(bool isInitial) {
    this->isInitial = isInitial;
  }

  bool getInitial() {
    return this->isInitial;
  }

  void setFinal(bool isFinal) {
    this->isFinal = isFinal;
  }

  bool getFinal() {
    return this->isFinal;
  }

  ~State() {
    for (auto state : previousStates) {
      if (state != nullptr) {
        delete state;
        state = nullptr;
      } 
    }

    for (auto state : nextStates) {
      if (state != nullptr) {
        delete state;
        state = nullptr;
      } 
    }
  }
};

#endif