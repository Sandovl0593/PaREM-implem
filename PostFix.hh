#ifndef POSTFIX_HH
#define POSTFIX_HH

#include <iostream>
#include <string>
#include <stack>
#include <unordered_map>
#include <vector>
#include <algorithm>
using namespace std;

class PostFix {
public:
  static const unordered_map<char, int> precedenceMap;

  /**
   * Get character precedence.
   * @param c character
   * @return corresponding precedence
   */
  static int getPrecedence(char c) {
    auto it = precedenceMap.find(c);
    return it == precedenceMap.end() ? 6 : it->second;
  }

  // Transform regular expression by inserting a '.' as explicit concatenation operator.
  static string formatRegEx(const string &regex) {
    string res;
    vector<char> allOperators = {'|', '?', '+', '*', '^'};
    vector<char> binaryOperators = {'^', '|'};

    for (size_t i = 0; i < regex.length(); ++i) {
      char c1 = regex[i];

      if (i + 1 < regex.length()) {
        char c2 = regex[i + 1];

        res += c1;

        if (c1 != '(' && c2 != ')' && find(allOperators.begin(), allOperators.end(), c2) == allOperators.end() && find(binaryOperators.begin(), binaryOperators.end(), c1) == binaryOperators.end()) {
          res += '.';
        }
      }
    }
    res += regex.back();

    return res;
  }

  /**
   * Convert regular expression from infix to postfix notation using Shunting-yard algorithm.
   * @param regex infix notation
   * @return postfix notation
   */
  static string infixToPostfix(const string &regex) {
    string postfix;
    stack<char> stack;

    string formattedRegEx = formatRegEx(regex);

    for (char c : formattedRegEx) {
      switch (c) {
        case '(':
          stack.push(c);
          break;

        case ')':
          while (stack.top() != '(') {
            postfix += stack.top();
            stack.pop();
          }
          stack.pop();
          break;

        default:
          while (!stack.empty()) {
            char peekedChar = stack.top();

            int peekedCharPrecedence = getPrecedence(peekedChar);
            int currentCharPrecedence = getPrecedence(c);

            if (peekedCharPrecedence >= currentCharPrecedence) {
              postfix += stack.top();
              stack.pop();
            } else {
              break;
            }
          }
          stack.push(c);
          break;
      }
    }

    while (!stack.empty()) {
      postfix += stack.top();
      stack.pop();
    }

    return postfix;
  }

  ~PostFix() {};
};

const unordered_map<char, int> PostFix::precedenceMap = {
  {'(', 1},
  {'|', 2},
  {'.', 3}, // explicit concatenation operator
  {'?', 4},
  {'*', 4},
  {'+', 4},
  {'^', 5}
};

#endif