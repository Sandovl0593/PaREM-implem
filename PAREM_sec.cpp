// PAREM algorithm implementation
// - Input: Transition map, the initial state and set of final states

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <cstdlib>
#include <ctime>

int main() {
    // Fase 1: Generación del automata determinista
    // Se necesita modificar la clase AFD que utilize la lista de variables
    const int numState = 3;
    const int initialState = 0;
    std::vector<bool> F(numState, false);
    F[2] = true;
    std::vector<std::map<char,int>> Tt;
    Tt.resize(numState);
    Tt[0]['0'] = 1;
    Tt[0]['1'] = 0;
    Tt[1]['0'] = 1;
    Tt[1]['1'] = 2;
    Tt[2]['0'] = 1;
    Tt[2]['1'] = 0;

    // Fase 2: Generación del input
    const int length = 50;
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
    std::vector<int> route = {currentState};
    if(F[currentState])
        found++;
    for (int i = 0; i < binaryString.size(); i++) {
        char bit = binaryString[i];
        if(Tt[currentState].find(bit) != Tt[currentState].end()) {
            currentState = Tt[currentState][bit];
            route.push_back(currentState);
            if(F[currentState])
                found++;
        } else {
            route.push_back(-1);
            break;
        }
    }
    std::cout << "Número de matches: " << found << std::endl;
    std::cout << "Recorrido en el automata:\n";
    for (int val : route)
        std::cout << val << " ";
    std::cout << '\n';
    return 0;
};

// O(n)
// - Output: Set of reachable states