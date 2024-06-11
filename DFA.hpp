#ifndef DFA_HPP
#define DFA_HPP

#include <iostream>
#include <vector>
#include <fstream>
#include <unordered_map>

using namespace std;

class DFA {
private:
    int n;
    vector<unordered_map<char, int>> transitions;
    vector<int> states;
    vector<char> alphabet;
    vector<int> F;
    int q0;

public:
    DFA();
    DFA(string filename);
    DFA(int n, vector<int> finalStates, unordered_map<int, unordered_map<char, int>> transitions);
    void readDFA(string filename);
    bool run(string input);
    void printDFA();
    vector<string> split(string str, char delimiter);
    vector<int> splitToInt(string str, char delimiter);
};

#endif // DFA_HPP