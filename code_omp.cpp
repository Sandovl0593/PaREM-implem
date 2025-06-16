#include "AFD.hpp"
#include <random>
#include <omp.h>
#include <string>
#include <fstream>

pair<vector<int>, int> run_parem(AFD afd, const string &T, int num_threads) {
    int N = T.size();
    int chunk = (N + num_threads - 1) / num_threads;
    vector<int> fullRoute;  // ruta completa del AFD
    int totalMatches = 0;

    struct BlockResult {
        vector<int> route;   // ruta parcial (sin contar el estado inicial global)
        int matches;
        int startState;      // estado inicial de este bloque
        int endState;        // estado final tras procesar el bloque
    };

    vector<BlockResult> results(num_threads);

    // Cada hilo procesa su bloque **para todos** los posibles estados de entrada,
    // pero en nuestra AFD sólo importa el real, así que directamente lo ejecutamos
    #pragma omp parallel for schedule(static)
    for (int t = 0; t < num_threads; t++) {
        int left = t * chunk;
        int right = min(N, left + chunk);
        BlockResult &R = results[t];
        R.matches = 0;

        if (left >= N) {
            R.startState = R.endState = -1;
            continue;
        }

        // Para el primer bloque, partimos del estado inicial; para el resto,
        // se rellenará tras la fase de reducción.
        if (t == 0) R.startState = 0;
        else        R.startState = -1;

        // Procesar el bloque de entrada T[left:right]
        int curr = 0;
        vector<int> localRoute;
        localRoute.reserve(right - left + 1);
        localRoute.push_back(curr);
        if (afd.accept(curr)) R.matches++;

        for (int i = left; i < right; i++) {
            curr = afd.nextState(curr, T[i]);
            localRoute.push_back(curr);
            if (afd.accept(curr)) R.matches++;
        }
        R.route     = std::move(localRoute);
        R.endState  = curr;
    }

    // Reducción secuencial de resultados
    totalMatches = 0;
    for (int t = 0; t < num_threads; t++) {
        if (t == 0) {
            fullRoute.push_back(results[t].startState);
        }
        // saltarnos el estado inicial en bloques >0 para evitar duplicarlo
        for (size_t i = (t == 0); i < results[t].route.size(); i++)
            fullRoute.push_back(results[t].route[i]);
        totalMatches += results[t].matches;

        // Propagar estado de inicio al siguiente bloque
        if (t+1 < num_threads)
            results[t+1].startState = results[t].endState;
    }
    return { fullRoute, totalMatches };
}

int main() {
    const int numStates  = 3;
    std::vector<int> finales = { 2 };
    AFD afd(numStates, finales);

    std::mt19937 rng(static_cast<unsigned int>(time(nullptr)));

    std::vector<char> alphabet = { 'a', 'b', 'c' };
    // alphabet.reserve(26);
    // for (int i = 0; i < 26; i++) {
    //     alphabet.push_back(static_cast<char>('a' + i));
    // }

    afd.addTransition(0, 'a', 1);  // \delta(0, 'a') = 1
    afd.addTransition(0, 'b', 0);  // \delta(0, 'b') = 0
    afd.addTransition(0, 'c', 2);  // \delta(0, 'c') = 2
    afd.addTransition(1, 'a', 1);  // \delta(1, 'a') = 1
    afd.addTransition(1, 'b', 2);  // \delta(1, 'b') = 2
    afd.addTransition(1, 'c', 0);  // \delta(1, 'c') = 0
    afd.addTransition(2, 'a', 1);  // \delta(2, 'a') = 1
    afd.addTransition(2, 'b', 0);  // \delta(2, 'b') = 0
    afd.addTransition(2, 'c', 2);  // \delta(2, 'c') = 2 (transición a sí mismo)

    const int length = 10000;

    // std::uniform_int_distribution<int> distChar(0, alphabet.size() - 1);
    // std::string T;
    // T.reserve(length);
    // for (int i = 0; i < length; i++) {
    //     T.push_back(alphabet[ distChar(rng) ]);
    // }

    // example from "example.txt"
    std::ifstream file("example.txt");
    std::string T; std::getline(file, T);
    file.close();

    double t0 = omp_get_wtime();
    auto [ruta, totMatches] = run_parem(afd, T, 8);
    double t1 = omp_get_wtime();

    std::cout << "Done! en tiempo: " 
              << (t1 - t0) * 1e6 
              << " microsegundos\n\n";

    std::cout << "Número de matches: " << totMatches << "\n";
    // std::cout << "Recorrido en el autómata:\n";
    // for (int estado : ruta) {
    //     std::cout << estado << ' ';
    // }
    std::cout << "\n";
}