#include "AFD.hpp"
#include <random>
#include <omp.h>
#include <string>

// Procesa la cadena 'T'; devuelve un par con:
//  - vector<int> ruta completa de estados (incluye estado inicial; en caso de fallo en transicion,
//    añade -1 y corta)
//  - numero de veces que se entro en un estado final (incluyendo el inicial si aplica)
pair<vector<int>, int> run(AFD afd, const string &T) {
    vector<int> route;
    route.reserve(T.size() + 1);

    int matches = 0;
    int curr = 0;
    route.push_back(curr);
    if (afd.accept(curr)) {
        matches++;
    }

    for (char c : T) {
        int nxt = afd.nextState(curr, c);
        if (nxt < 0) {
            route.push_back(-1);
            break;
        }
        curr = nxt;
        route.push_back(curr);
        if (afd.accept(curr)) {
            matches++;
        }
    }
    return { route, matches };
}


int main() {
    const int numStates  = 3;
    std::vector<int> finales = { 2 };
    AFD afd(numStates, finales);

    std::mt19937 rng(static_cast<unsigned int>(time(nullptr)));

    std::vector<char> alphabet = { 'a', 'b', 'c' };
    // alphabet.reserve(26);
    // for (int i = 0; i < 26; ++i) {
    //     alphabet.push_back(static_cast<char>('a' + i));
    // }

    std::cout << "Símbolos: ";
    for (char c : alphabet) std::cout << c << ' ';
    std::cout << "\n\n";

    afd.addTransition(0, 'a', 1);  // \delta(0, 'a') = 1
    afd.addTransition(0, 'b', 0);  // \delta(0, 'b') = 0
    afd.addTransition(0, 'c', 2);  // \delta(0, 'c') = 2
    afd.addTransition(1, 'a', 1);  // \delta(1, 'a') = 1
    afd.addTransition(1, 'b', 2);  // \delta(1, 'b') = 2
    afd.addTransition(1, 'c', 0);  // \delta(1, 'c') = 0
    afd.addTransition(2, 'a', 1);  // \delta(2, 'a') = 1
    afd.addTransition(2, 'b', 0);  // \delta(2, 'b') = 0
    afd.addTransition(2, 'c', 2);  // \delta(2, 'c') = 2 (transición a sí mismo)

    const int length = 100;

    std::uniform_int_distribution<int> distChar(0, 25);
    std::string T;
    T.reserve(length);
    for (int i = 0; i < length; ++i) {
        T.push_back(alphabet[ distChar(rng) ]);
    }
    std::cout << "Entrada a evaluar: " << T << "\n\n";

    // Medición con OpenMP (secuencial)
    double t0 = omp_get_wtime();
    auto [ruta, totMatches] = run(afd, T);
    double t1 = omp_get_wtime();

    std::cout << "Done! en tiempo: " 
              << (t1 - t0) * 1e6 
              << " microsegundos\n\n";

    std::cout << "Número de matches: " << totMatches << "\n";
    std::cout << "Recorrido en el autómata:\n";
    for (int estado : ruta) {
        std::cout << estado << ' ';
    }
    std::cout << "\n";
}
