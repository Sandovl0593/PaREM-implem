// PAREM algorithm implementation
// - Input: Transition unordered_map, set of initial states, set of final states

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <cstdlib>
#include <ctime>

int main() {
    // Fase 1: Generación del automata determinista
    // Se necesita modificar la clase AFD que utilize la lista de variables
    int initialState = 0;
    std::set<int> finalState;
    std::vector<std::map<char,int>> transitionList;
    // Método de inserción de estados y transicciones
    transitionList.resize(3);
    finalState.insert(2);
    transitionList[0]['0'] = 1;
    transitionList[0]['1'] = 0;
    transitionList[1]['0'] = 1;
    transitionList[1]['1'] = 2;
    transitionList[2]['0'] = 1;
    transitionList[2]['1'] = 0;

    // Fase 2: Generación del input
    const int length = 200;
    std::string binaryString;
    binaryString.reserve(length);
    std::srand(std::time(0));
    for (int i = 0; i < length; ++i) {
        char bit = (rand() % 2) ? '1' : '0';
        binaryString += bit;
    }
    // binaryString = "010101010101010101010101010101";
    std::cout << "Entrada a evaluar: " << binaryString << std::endl;

    // Fase 3: Recorrido del input dentro del automata
    // Aquí se debe paralelizar el proceso con OMP
    int found = 0;
    int currentState = initialState;
    std::string way = std::to_string(currentState);
    if(finalState.find(currentState) != finalState.end())
        found++;
    for (int i = 0; i < binaryString.size(); ++i) {
        char bit = binaryString[i];
        if(transitionList[currentState].find(bit) != transitionList[currentState].end()) {
            currentState = transitionList[currentState][bit];
            way = way + " - " + std::to_string(currentState);
            if(finalState.find(currentState) != finalState.end())
                found++;
        } else {
            way = way + " - X";
            break;
        }
    }
    std::cout << "Número de matches: " << found << std::endl;
    std::cout << "Recorrido en el automata:\n" << way << std::endl; 
    return 0;
};

// - Output: Set of reachable states