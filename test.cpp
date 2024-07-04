#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "AFN.hh"
#include "PostFix.hh"

using namespace std;

void writeFile(string regexp) {
    AFN afn(regexp);
    ofstream writer("AFNresult.txt");
    writer << "REGULAR EXPRESSION: " << regexp << endl;
    writer << "REGULAR EXPRESSION IN POSTFIX: " << afn.getPostFixRegExp() << endl;
    writer << "SYMBOL LIST: ";
    for (char ch : afn.getSymbolList()) {
        writer << ch << " ";
    } writer << endl;

    writer << "TRANSITIONS LIST: ";
    for (Transition* transition : afn.getTransitionsList()) {
        writer << transition->toString() << " ";
    } writer << endl;

    writer << "FINAL STATE: ";
    for (State* state : afn.getFinalStates()) {
        writer << state->toString() << " ";
    } writer << endl;

    writer << "STATES: ";
    for (string state : afn.getStates()) {
        writer << state << " ";
    } writer << endl;

    writer << "INITIAL STATE: ";
    for (string state : afn.getInitialState()) {
        writer << state << " ";
    } writer << endl;

    writer.close();
}


int main() {
    string regexp;
    cout << "Enter a regular expression: ";
    getline(cin, regexp);
    writeFile(regexp);

    return 0;
}
