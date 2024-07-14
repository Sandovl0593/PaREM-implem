// PAREM algorithm implementation
// - Input: Transition unordered_map, set of initial states, set of final states

#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <omp.h>

int main() {
    // Fase 1: Generación del automata determinista
    // Se necesita modificar la clase AFD que utilize la lista de variables
    int initialState = 0;
    std::set<int> F;
    std::vector<std::map<char,int>> Tt;
    // Método de inserción de estados y transicciones
    Tt.resize(3);
    F.insert(2);
    Tt[0]['0'] = 1;
    Tt[0]['1'] = 0;
    Tt[1]['0'] = 1;
    Tt[1]['1'] = 2;
    Tt[2]['0'] = 1;
    Tt[2]['1'] = 0;

    // Fase 2: Generación del input
    const int length = 40;
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
    const int num_threads = 4;
    std::set<int> S[num_threads], L[num_threads];
    S[0].insert(initialState);
    L[num_threads-1].insert(F.begin(), F.end());
    std::vector<std::string> ways(num_threads);
    

    int start, end, found, currentState;
    std::set<int> R;
    std::map<int, std::vector<int>> route[num_threads];
    #pragma omp parallel num_threads(num_threads) private(start, end, found, currentState, R)
    { 
        int id = omp_get_thread_num();
        int chunk_size = binaryString.size() / num_threads;
        start = chunk_size * id;
        end = chunk_size*(id+1)-1;
        if(id == 3) 
            end = binaryString.size()-1;

        if(id != 0) {
            for(int i=0; i < Tt.size(); i++)
                if(Tt[i].find(binaryString[start]) != Tt[i].end())
                    S[id].insert(i);
        }

        if(id != num_threads-1) {
            for(int i=0; i < Tt.size(); i++)
                if(Tt[i].find(binaryString[end]) != Tt[i].end())
                    L[id].insert(Tt[i][binaryString[end]]);
        }

        #pragma omp barrier
        if(id != 0)
            std::set_intersection(S[id].begin(), S[id].end(), L[id-1].begin(), L[id-1].end(), std::inserter(R, R.begin()));
        else
            R.insert(S[id].begin(),S[id].end());

        for(int init: R) {
            currentState = init;
            route[id][currentState] = {currentState};
            if(F.find(currentState) != F.end())
                found++;
            for (int i = start; i <= end; i++) {
                char bit = binaryString[i];
                if(Tt[currentState].find(bit) != Tt[currentState].end()) {
                    currentState = Tt[currentState][bit];
                    route[id][init].push_back(currentState);
                    if(F.find(currentState) != F.end())
                        found++;
                } else {
                    route[id][init].push_back(-1);
                    break;
                }
            }
        }

        #pragma omp critical
        {
            std::cout << '\n' << "Proceso " << id << ":" << '\n';
            for (const auto& pair : route[id]) {
                std::cout << "Clave: " << pair.first << " -> Valores: ";
                for (int val : pair.second) {
                    std::cout << val << " ";
                }
                std::cout << std::endl;
            }
        }
    }
    return 0;
};

// - Output: Set of reachable states

// Compilar: gcc Ejemp01_master.c -fopenmp
// Ejecutar: .\a.exe