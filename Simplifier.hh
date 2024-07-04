#ifndef SIMPLIFIER_HH
#define SIMPLIFIER_HH

#include <iostream>
#include <string>
using namespace std;

class ExpressionSimplifier {
    string regExp; // the regexp
  public:
    ExpressionSimplifier() {
      this->regExp = "";
    }
    ExpressionSimplifier(string regExp) {
      this->regExp = regExp;
      handleKleeneSum();
      handleLua();
    }

    string getRegExp() {
      return this->regExp;
    }

    // This function handles the '?' operator
    void handleLua() {
      for (int i = 0; i < regExp.size(); i++) {
        if (regExp[i] == '?') {
          if (regExp[i-1] != ')') {
            // ? is after a symbol, in which case ? afects only that symbol
            string symbol = string(1, regExp[i-1]);
            string subExpression = "(" + symbol + "|ε)";

            string left = regExp.substr(0, i-1);
            string right = regExp.substr(i+1);
            regExp = left + subExpression + right;
          } 
          else {
            // ? is after a ')', in which case ? afects a whole expresion inside '(' and ')'
            int j = i-1;
            while (j > 0) {
              if (regExp[j] == '(') {
                string symbolSequence = regExp.substr(j+1, i-j-2);
                string symbolSequenceWithBrackets = regExp.substr(j, i-j+1);
                string subExp = "(" + symbolSequenceWithBrackets + "|ε)";

                string left = regExp.substr(0, j);
                string right = regExp.substr(i+1);
                regExp = left + subExp + right;
                break;
              }
              j--;
            }
          }
        }
      }
    }

    // This function handles the '+' operator
    void handleKleeneSum() {
      for (int i = 0; i < regExp.length(); i++) {
        if (regExp[i] == '+') {
          if (i > 0 && regExp[i-1] != ')') {
            // + is after a symbol, in which case + afects only that symbol
            string symbol = string(1, regExp[i-1]);
            string subExpression = symbol + symbol + "*";

            string left = regExp.substr(0, i-1);
            string right = regExp.substr(i+1);
            regExp = left + subExpression + right;
          } 
          else {
            // + is after a ')', in which case + afects a whole expresion inside '(' and ')'
            int bracketCounter = 0;
            int j = i-1;
            while (j >= 0) {
              if (j != (i-1) && regExp[j] == ')') {
                bracketCounter++;
              }
              if (regExp[j] == '(') {
                if (bracketCounter == 0) {
                  bracketCounter--;
                } else {
                  string symbolSequence = regExp.substr(j+1, i-j-2);
                  string symbolSequenceWithBrackets = regExp.substr(j, i-j+1);
                  string subExp = symbolSequenceWithBrackets + symbolSequenceWithBrackets + "*";

                  string left = regExp.substr(0, j);
                  string right = regExp.substr(i+1);
                  regExp = left + subExp + right;
                }
              }
              j--;
            }
          }
        }
      }
    }

    ~ExpressionSimplifier() {}
};

#endif