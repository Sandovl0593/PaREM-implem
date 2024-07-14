// PAREM algorithm implementation
// - Input: Transition unordered_map, set of initial states, set of final states

#include <iostream>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <algorithm>
#include <utility>
#include <cstdlib>
#include <ctime>
#include <omp.h>

int main() {
    // Fase 1: Generación del automata determinista
    // Se necesita modificar la clase AFD que utilize la lista de variables
    const int numState = 5;
    int initialState = 0;
    std::vector<bool> F(numState, false);
    F[3] = true;
    std::vector<std::unordered_map<char,int>> Tt;
    Tt.resize(numState);
    Tt[0]['0'] = 1;
    Tt[0]['1'] = 2;
    Tt[1]['0'] = 1;
    Tt[1]['1'] = 2;
    Tt[2]['0'] = 3;
    Tt[2]['1'] = 2;
    Tt[3]['0'] = 4;
    Tt[3]['1'] = 4;
    Tt[4]['0'] = 1;
    Tt[4]['1'] = 2;

    // Fase 2: Generación del input
    const int length = 120;
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
    std::cout << "Numero de caracteres: " << binaryString.size() << std::endl;
    
    // Fase 3: Recorrido del input dentro del automata
    // Aquí se debe paralelizar el proceso con OMP
    int total_found = 0;
    const int num_threads = 8;
    int start, end, found = 0, currentState;
    std::unordered_map<int,std::vector<int>> route[num_threads];
    route[0][initialState] = std::vector<int>();
    #pragma omp parallel num_threads(num_threads) private(start, end, found, currentState)
    { 
        int id = omp_get_thread_num();
        int chunk_size = binaryString.size() / num_threads;
        int start = chunk_size * id;
        int end = chunk_size*(id+1)-1;
        if(id == num_threads-1) 
            end = binaryString.size()-1;

        if(id != num_threads-1) {
            for(int i=0; i < Tt.size(); i++)
                route[id+1][Tt[i][binaryString[end]]] = std::vector<int>();
        } // O(s); s: Número de estados del automata

        #pragma omp barrier

        for(const auto& pair: route[id]) {
            currentState = pair.first;
            route[id][currentState] = {currentState};
            if(F[currentState])
                found++;
            for (int i = start; i <= end; i++) {
                currentState = Tt[currentState][binaryString[i]];
                route[id][pair.first].push_back(currentState);
                if(F[currentState])
                    found++;
            }
        } // O(s*((n/p)+p))

        #pragma omp critical
        {
            std::cout << '\n' << "Proceso " << id << ":" << '\n';
            for (const auto& pair : route[id]) {
                std::cout << " - Clave: " << pair.first << " -> Valores: ";
                for (int val : pair.second)
                    std::cout << val << " ";
                std::cout << std::endl;
            }
        }
    }
    return 0;
};