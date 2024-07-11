#ifndef AFD_HH
#define AFD_HH
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <list>
#include "AFN.hh"
using namespace std;

class AFD {
public:
  vector<char> symbolList;
  vector<Transition*> transitionsList;
  vector<State*> finalStates;
  int initialState;
  vector<int> states;

  AFD(AFN* afn) {
    // Init AFD
    symbolList = afn->getSymbolList();
    initialState = 0;

    // Table of states and mapping of subsets of states of the AFN to states of the AFD
    map<vector<int>, State*> subsetToState;
    list<vector<int>> pendingSubsets;

    // Initial state of the AFD is the ε-closure of the initial state of the AFN
    vector<int> initialSubset = epsilonClosure({afn->getInitialState()}, afn);
    State* initialDFAState = new State(initialState);
    subsetToState[initialSubset] = initialDFAState;
    states.push_back(initialState);
    pendingSubsets.push_back(initialSubset);

    // Process all subsets of states of the AFN
    while (!pendingSubsets.empty()) {
      vector<int> currentSubset = pendingSubsets.front();
      pendingSubsets.pop_front();
      State* currentDFAState = subsetToState[currentSubset];

      for (char symbol : symbolList) {
        vector<int> newSubset = move(currentSubset, symbol, afn);
        newSubset = epsilonClosure(newSubset, afn);

        if (newSubset.empty()) continue;

        if (subsetToState.find(newSubset) == subsetToState.end()) {
          // New state in AFD
          int newStateId = states.size();
          State* newState = new State(newStateId);
          subsetToState[newSubset] = newState;
          states.push_back(newStateId);
          pendingSubsets.push_back(newSubset);
        }

        State* nextDFAState = subsetToState[newSubset];
        Transition* newTransition = new Transition(string(1, symbol), currentDFAState, nextDFAState);
        transitionsList.push_back(newTransition);
        // currentDFAState->nextStates.push_back(nextDFAState);
        // nextDFAState->previousStates.push_back(currentDFAState);
      }
    }

    // Determinates the final states of the AFD
    for (const auto& subsetStatePair : subsetToState) {
      const vector<int>& subset = subsetStatePair.first;
      State* dfaState = subsetStatePair.second;
      for (State* afnFinalState : afn->getFinalStates()) {
        if (find(subset.begin(), subset.end(), afnFinalState->stateId) != subset.end()) {
          finalStates.push_back(dfaState);
          break;
        }
      }
    }
  }

  vector<char> getSymbolList() const {
    return symbolList;
  }

  vector<Transition*> getTransitionsList() const {
    return transitionsList;
  }

  vector<State*> getFinalStates() const {
    return finalStates;
  }

  vector<int> getStates() const {
    return states;
  }

  int getInitialState() const {
    return initialState;
  }

private:
  vector<int> epsilonClosure(const vector<int>& states, AFN* afn) {
    vector<int> closure = states;
    list<int> toProcess(states.begin(), states.end());

    while (!toProcess.empty()) {
      int state = toProcess.front();
      toProcess.pop_front();
      for (Transition* transition : afn->getTransitionsList()) {
        if (transition->initialState->stateId == state && transition->transitionSymbol == "ε") {
          int nextStateId = transition->finalState->stateId;
          if (find(closure.begin(), closure.end(), nextStateId) == closure.end()) {
            closure.push_back(nextStateId);
            toProcess.push_back(nextStateId);
          }
        }
      }
    }
    return closure;
  }

  vector<int> move(const vector<int>& states, char symbol, AFN* afn) {
    vector<int> result;
    for (int state : states)
      for (Transition* transition : afn->getTransitionsList())
        if (transition->initialState->stateId == state && transition->transitionSymbol == string(1, symbol))
          if (find(result.begin(), result.end(), transition->finalState->stateId) == result.end()) 
            result.push_back(transition->finalState->stateId);
    return result;
  }
};

#endif