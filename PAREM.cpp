// #include <iostream>
// #include <vector>
// #include <fstream>
// #include <unordered_map>
// #include <omp.h>
// #include "Transition.hh"
// #include "AFD.hh"


// PAREM algorithm implementation
// - Input: Transition unordered_map, set of initial states, set of final states

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <cstdlib>
#include <ctime>

#include <omp.h>

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
    // binaryString = "101010101011110000001111111111111000000";
    // binaryString = "10111010010111000011001000110111000101110100100101001100011010010000011001011101001001001100101111110110000010111111011101110001100111011010010101001101111101011010110111110110111011100010010111000011";
    std::cout << "Entrada a evaluar: " << binaryString << std::endl;

    // Fase 3: Recorrido del input dentro del automata
    // Aquí se debe paralelizar el proceso con OMP
    int total_found = 0;
    const int num_threads = 4;
    int chunk_size = binaryString.size() / num_threads;

    std::vector<std::string> ways(num_threads);

    #pragma omp parallel num_threads(num_threads) reduction(+:total_found)
    {
        int thread_id = omp_get_thread_num();
        int start = thread_id * chunk_size;
        int end = (thread_id == num_threads - 1) ? binaryString.size() : start + chunk_size;

        int currentState = initialState;
        int found = 0;
        std::string way;
    
        for (int i = start; i < end; ++i) {
            int currentState = initialState;
            way = std::to_string(currentState);

            for (int j = i; j < binaryString.size(); ++j) {
                char bit = binaryString[j];
                if (transitionList[currentState].find(bit) != transitionList[currentState].end()) {
                    currentState = transitionList[currentState][bit];
                    way += " - " + std::to_string(currentState);
                    if (finalState.find(currentState) != finalState.end()) {
                        found++;
                    }
                } else {
                    way += " - X";
                    break;
                }
            }
    }   

        total_found += found;
        ways[thread_id] = way;

    }

    std::cout << "Número total de matches: " << total_found << std::endl;

    for (int i = 0; i < num_threads; ++i)
    {
        std::cout << "Hilo " << i << " recorrido en el automata:\n"
                  << ways[i] << std::endl;
    }

    return 0;

    // std::cout << "Número de matches: " << found << std::endl;
    // std::cout << "Recorrido en el automata:\n" << way << std::endl; 
    // return 0;
};

// - Output: Set of reachable states