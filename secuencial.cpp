#include <iostream>
#include <vector>
#include <array>
#include <random>
#include <ctime>

class AFD {
private:
    int numStates;
    int initialState;
    std::vector<bool> isFinal;                 // isFinal[s] == true si el estado s es final
    std::vector<std::array<int,2>> trans;      // trans[s][0] = siguiente estado leyendo '0'
                                               // trans[s][1] = siguiente estado leyendo '1'
public:
    // Constructor: recibe número de estados, estado inicial y lista de estados finales
    AFD(int n,int init, const std::vector<int> &finalStates)
        : numStates(n),
          initialState(init),
          isFinal(n, false),
          trans(n) {
        // Marcar estados finales
        for (int f : finalStates) {
            if (f >= 0 && f < n) {
                isFinal[f] = true;
            }
        }
        // Inicializamos todas las transiciones a -1 (sin definir)
        for (int s = 0; s < n; ++s) {
            trans[s][0] = -1;
            trans[s][1] = -1;
        }
    }

    // Agregar/transicionar desde un estado 's' con símbolo '0' o '1'
    void addTransition(int s, char symbol, int nextState) {
        int idx = (symbol == '1') ? 1 : 0;
        if (s >= 0 && s < numStates && nextState >= 0 && nextState < numStates) {
            trans[s][idx] = nextState;
        }
    }

    // Verifica si un estado es final
    bool isAccepting(int s) const {
        return (s >= 0 && s < numStates) ? isFinal[s] : false;
    }

    // Obtiene el siguiente estado desde 's' con el símbolo 'c'; devuelve -1 si no existe
    int nextState(int s, char c) const {
        if (s < 0 || s >= numStates) return -1;
        int idx = (c == '1') ? 1 : 0;
        return trans[s][idx];
    }

    // Procesa la cadena 'input'; devuelve un par con:
    //  - vector<int> ruta completa de estados (incluye estado inicial; en caso de fallo en transición, 
    //    añade -1 y corta)
    //  - número de veces que se entró en un estado final (incluyendo el inicial si aplica)
    std::pair<std::vector<int>, int> run(const std::string &input) const {
        std::vector<int> route;
        route.reserve(input.size() + 1);

        int matches = 0;
        int curr = initialState;
        route.push_back(curr);
        if (isAccepting(curr)) {
            matches++;
        }

        for (char c : input) {
            int nxt = nextState(curr, c);
            if (nxt < 0) {
                route.push_back(-1);
                break;
            }
            curr = nxt;
            route.push_back(curr);
            if (isAccepting(curr)) {
                matches++;
            }
        }
        return { route, matches };
    }
};


int main() {
    // Configuración del autómata determinista ---
    const int numState = 3;
    const int initialState = 0;
    // Estado final = {2}
    std::vector<int> finales = { 2 };
    AFD afd(numState, initialState, finales);

    // Definimos transiciones (igual lógica que en el código original)
    afd.addTransition(0, '0', 1);
    afd.addTransition(0, '1', 0);
    afd.addTransition(1, '0', 1);
    afd.addTransition(1, '1', 2);
    afd.addTransition(2, '0', 1);
    afd.addTransition(2, '1', 0);

    // Generación de la cadena binaria aleatoria ---
    const int length = 50;
    std::string binaryString;
    binaryString.reserve(length);

    // Usamos <random> para mejor calidad de aleatoriedad
    std::mt19937 rng(static_cast<unsigned int>(std::time(nullptr)));
    std::uniform_int_distribution<int> dist(0, 1);

    for (int i = 0; i < length; ++i) {
        binaryString.push_back(dist(rng) == 1 ? '1' : '0');
    }

    std::cout << "Entrada a evaluar: " << binaryString << std::endl;

    // Recorrido de la cadena en el autómata ---
    auto [ruta, totMatches] = afd.run(binaryString);

    std::cout << "Número de matches: " << totMatches << std::endl;
    std::cout << "Recorrido en el automata:\n";
    for (int estado : ruta) {
        std::cout << estado << " ";
    }
    std::cout << std::endl;

    return 0;
}
