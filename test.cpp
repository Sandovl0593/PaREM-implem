#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "AFN.hh"
#include "AFD.hh"

using namespace std;

void writeFile(string regexp) {
    AFN* afn = new AFN(regexp);
    ofstream writeAFN("AFNresult.txt");
    writeAFN << "REGULAR EXPRESSION: " << regexp << endl;
    writeAFN << "REGULAR EXPRESSION IN POSTFIX: " << afn->getPostFixRegExp() << endl;
    writeAFN << "SYMBOL LIST: ";
    for (char ch : afn->getSymbolList()) {
        writeAFN << ch << ", ";
    } writeAFN << endl;

    writeAFN << "TRANSITIONS LIST: ";
    for (Transition* transition : afn->getTransitionsList()) {
        writeAFN << transition->toString() << ", ";
    } writeAFN << endl;

    writeAFN << "FINAL STATE: ";
    for (State* state : afn->getFinalStates()) {
        writeAFN << state->toString() << ", ";
    } writeAFN << endl;

    writeAFN << "STATES: ";
    for (int state : afn->getStates()) {
        writeAFN << state << ", ";
    } writeAFN << endl;

    writeAFN << "INITIAL STATE: " << afn->getInitialState() << endl;

    writeAFN.close();


    AFD* afd = new AFD(afn);
    ofstream writeAFD("AFDresult.txt");
    writeAFD << "REGULAR EXPRESSION: " << regexp << endl;
    writeAFD << "SYMBOL LIST: ";
    for (char ch : afd->getSymbolList()) {
        writeAFD << ch << ", ";
    } writeAFD << endl;

    writeAFD << "TRANSITIONS LIST: ";
    for (Transition* transition : afd->getTransitionsList()) {
        writeAFD << transition->toString() << ", ";
    } writeAFD << endl;

    writeAFD << "FINAL STATE: ";
    for (State* state : afd->getFinalStates()) {
        writeAFD << state->toString() << ", ";
    } writeAFD << endl;

    writeAFD << "STATES: ";
    for (int state : afd->getStates()) {
        writeAFD << state << ", ";
    } writeAFD << endl;

    writeAFD << "INITIAL STATE: " << afd->getInitialState() << endl;

    delete afd;
    delete afn;
}


int main() {
    string regexp;
    cout << "Enter a regular expression: ";
    getline(cin, regexp);
    writeFile(regexp);
    cout << "Done!" << endl;

    return 0;
}
