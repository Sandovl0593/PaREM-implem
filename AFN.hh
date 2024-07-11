#ifndef AFN_HH
#define AFN_HH

#include <iostream>
#include <vector>
#include <stack>
#include <algorithm>
#include <string>
#include <unordered_map>
#include "PostFix.hh"
#include "Simplifier.hh"
#include "Transition.hh"
using namespace std;

const char EPSILON = '\x040';

class AFN {
  PostFix* postFix; // to handle infix to postfix
  string regExp; // the regExp in infix
  string postFixRegExp; // the regExp in postfix

  vector<char> symbolList; // AFN's symbol list
  vector<Transition*> transitionsList; // AFN's transitions list
  vector<State*> finalStates; // AFN's acceptation states list
  int initialState; // AFN's initial state
  vector<int> states; // AFN's states

  int stateCount = 0; // id for States

  /* To save state reference */
  State* saveFinal;
  stack<State*> stackInitial;
  stack<State*> stackFinal;

  /* Abbreviations and juxtaposition concatenation */
  ExpressionSimplifier* expressionSimplifier = new ExpressionSimplifier();

  // Checks for symbols in postFixRegExp and adds them to symbolList
  void computeSymbolList() {
    for (char ch : postFixRegExp) {
      if (!PostFix::precedenceMap.count(ch)) {
        if (find(symbolList.begin(), symbolList.end(), ch) == symbolList.end()) {
          symbolList.push_back(ch);
          sort(symbolList.begin(), symbolList.end());
        }
      }
    }
  }

  // Adds AFN's states to stateList
  void computeStateList() {
    for (auto transition : transitionsList) {
      if (find(states.begin(), states.end(), transition->getInitialState()->getStateId()) == states.end()) {
        states.push_back(transition->getInitialState()->getStateId());
      }
      if (find(states.begin(), states.end(), transition->getFinalState()->getStateId()) == states.end()) {
        states.push_back(transition->getFinalState()->getStateId());
      }
    }
  }

  // Adds AFN's initial state to list
  void computeInitialState() {
    this->initialState = stackInitial.top()->getStateId();
    stackInitial.pop();
  }

  // The real deal; parses the expression to an AFN.
  void regExpToAFN() {
    for (size_t i = 0; i < postFixRegExp.length(); ++i) {
      char ch = postFixRegExp[i];
      if (find(symbolList.begin(), symbolList.end(), ch) != symbolList.end()) {
        Transition* tr1 = new Transition(string(1, ch), stateCount);
        transitionsList.push_back(tr1);

        State* initialState = tr1->getInitialState();
        State* finalState = tr1->getFinalState();

        stackInitial.push(initialState);
        stackFinal.push(finalState);
      } else if (ch == '|') {
        State* lowerInitial = stackInitial.top();
        stackInitial.pop();
        State* lowerFinal = stackFinal.top();
        stackFinal.pop();
        State* upperInitial = stackInitial.top();
        stackInitial.pop();
        State* upperFinal = stackFinal.top();
        stackFinal.pop();

        unify(upperInitial, upperFinal, lowerInitial, lowerFinal);
      } else if (ch == '*') {
        State* initialState = stackInitial.top();
        stackInitial.pop();
        State* finalState = stackFinal.top();
        stackFinal.pop();

        kleene(initialState, finalState);
      } else if (ch == '.') {
        saveFinal = stackFinal.top();
        stackFinal.pop();
        State* finalState = stackFinal.top();
        stackFinal.pop();
        State* initialState = stackInitial.top();
        stackInitial.pop();

        concatenate(finalState, initialState);
      }

      if (i == postFixRegExp.length() - 1) {
        finalStates.push_back(stackFinal.top());
        stackFinal.pop();
        auto it = find(symbolList.begin(), symbolList.end(), EPSILON);
        if (it != symbolList.end()) {
          symbolList.erase(it);
        }
      }
    }
  }

  // Union | in Thompson's algorithm
  void unify(State* upperInitialState, State* upperFinalState, State* lowerInitialState, State* lowerFinalState) {
    State* in = new State(stateCount++);
    State* out = new State(stateCount++);

    Transition* tr1 = new Transition("ε", in, upperInitialState);
    Transition* tr2 = new Transition("ε", in, lowerInitialState);
    Transition* tr3 = new Transition("ε", upperFinalState, out);
    Transition* tr4 = new Transition("ε", lowerFinalState, out);

    transitionsList.push_back(tr1);
    transitionsList.push_back(tr2);
    transitionsList.push_back(tr3);
    transitionsList.push_back(tr4);

    stackInitial.push(in);
    stackFinal.push(out);
  }

  // Concatenation in Thompson's algorithm
  void concatenate(State* initialState, State* finalState) {
    Transition* tr1 = new Transition("ε", initialState, finalState);
    transitionsList.push_back(tr1);

    stackFinal.push(saveFinal);
  }

  // Kleene star in Thompson's algorithm
  void kleene(State* initialState, State* finalState) {
    State* in = new State(stateCount++);
    State* out = new State(stateCount++);

    Transition* tr1 = new Transition("ε", finalState, initialState); // upper kleene part
    Transition* tr2 = new Transition("ε", in, out); // lower kleene part
    Transition* tr3 = new Transition("ε", in, initialState); // initial to initial of expression
    Transition* tr4 = new Transition("ε", finalState, out); // final of expression to out

    transitionsList.push_back(tr1);
    transitionsList.push_back(tr2);
    transitionsList.push_back(tr3);
    transitionsList.push_back(tr4);

    stackInitial.push(in);
    stackFinal.push(out);
  }

public:

  // static int stateCount; // id for States

  AFN(const string& regExp) : regExp(regExp) {
    expressionSimplifier = new ExpressionSimplifier(regExp);
    // stateCount = 0;
    postFixRegExp = PostFix::infixToPostfix(expressionSimplifier->getRegExp());
    computeSymbolList();
    regExpToAFN();
    computeStateList();
    computeInitialState();
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

  string getPostFixRegExp() const {
    return postFixRegExp;
  }

  ~AFN() {
    
  }
};

// int stateCount = 0; // id for States

#endif