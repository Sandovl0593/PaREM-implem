#include "AFD.hpp"
#include <random>
#include <omp.h>
#include <string>
#include <fstream>

struct BlockInfo {
    int endState;
    int matches;
};

// JumpTable signfica que cada bloque precomputará el estado final y el número de matches
vector<vector<BlockInfo>> precompute_jump(AFD afd, const string &T, int num_threads, int chunk){
    int Q = afd.numStates;
    vector<vector<BlockInfo>> jumpTable(num_threads, vector<BlockInfo>(Q));

    #pragma omp parallel for schedule(static)
    for (int t = 0; t < num_threads; ++t) {
        int left  = t * chunk;
        int right = min((int)T.size(), left + chunk);
        for (int q = 0; q < Q; ++q) {
            int curr = q;
            int cnt = afd.accept(curr) ? 1 : 0;
            for (int i = left; i < right; ++i) {
                curr = afd.nextState(curr, T[i]);
                if (afd.accept(curr)) cnt++;
            }
            jumpTable[t][q] = { curr, cnt };
        }
    }
    return jumpTable;
}

// Inicializa el estado inicial de cada bloque
vector<int> compute_start_states(const vector<vector<BlockInfo>> &jumpTable) {
    int p = jumpTable.size();
    vector<int> startState(p);
    startState[0] = 0;
    for (int t = 1; t < p; ++t) {
        int prev = startState[t-1];
        startState[t] = jumpTable[t-1][prev].endState;
    }
    return startState;
}

// Recorre en paralelo y escribe directamente en fullRoute
pair<vector<int>, int> run_parem(AFD afd, const string &T, int num_threads) {
    int N = T.size();
    int chunk = (N + num_threads - 1) / num_threads;

    auto jumpTable = precompute_jump(afd, T, num_threads, chunk);

    auto startState = compute_start_states(jumpTable);

    vector<int> fullRoute(N + 1);
    int totalMatches = 0;

    #pragma omp parallel for reduction(+:totalMatches) schedule(static)
    for (int t = 0; t < num_threads; ++t) {
        int left  = t * chunk;
        int right = min(N, left + chunk);
        if (left >= right) continue;

        int curr = startState[t];
        int localMatches = 0;
        int writePos = left;

        // escribe primer estado
        fullRoute[writePos++] = curr;
        if (afd.accept(curr)) localMatches++;

        // procesa bloque usando salto de bloque a bloque
        // (si prefieres procesar carácter a carácter podrías usar jumpTable[t])
        for (int i = left; i < right; ++i) {
            curr = afd.nextState(curr, T[i]);
            fullRoute[writePos++] = curr;
            if (afd.accept(curr)) localMatches++;
        }

        totalMatches += localMatches;
    }

    return { fullRoute, totalMatches };
}


int main() {
    const int numStates  = 3;
    vector<int> finales = { 2 };
    AFD afd(numStates, finales);

    mt19937 rng(static_cast<unsigned int>(time(nullptr)));

    vector<char> alphabet = { 'a', 'b', 'c' };
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

    // uniform_int_distribution<int> distChar(0, alphabet.size() - 1);
    // string T;
    // T.reserve(length);
    // for (int i = 0; i < length; i++) {
    //     T.push_back(alphabet[ distChar(rng) ]);
    // }

    // example from "example.txt"
    ifstream file("example.txt");
    string T; getline(file, T);
    file.close();

    double t0 = omp_get_wtime();
    auto [ruta, totMatches] = run_parem(afd, T, 8);
    double t1 = omp_get_wtime();

    cout << "Done! en tiempo: " 
              << (t1 - t0) * 1e6 
              << " microsegundos\n\n";

    cout << "Número de matches: " << totMatches << "\n";
    // cout << "Recorrido en el autómata:\n";
    // for (int estado : ruta) {
    //     cout << estado << ' ';
    // }
    cout << "\n";
}