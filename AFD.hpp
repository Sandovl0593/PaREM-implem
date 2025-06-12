#include <iostream>
#include <vector>
#include <unordered_map>

using namespace std;

class AFD {
    vector<bool> isFinal;                                // isFinal[s] == true si el estado s es final
    vector<unordered_map<char,int>> trans;               // trans[s][c] = siguiente estado leyendo el símbolo c
public:
    int numStates;
    // Constructor: recibe numero de estados, estado inicial y lista de estados finales
    AFD(int n, const vector<int> &finalStates)
        : numStates(n),
          isFinal(n, false),
          trans(n) {
        // Marcar estados finales
        for (int f : finalStates) {
            if (f >= 0 && f < n) {
                isFinal[f] = true;
            }
        }
    }

    // Agregar/transicionar desde un estado 's' con símbolo cualquiera del alfabeto
    void addTransition(int s, char symbol, int nextState) {
        if (s >= 0 && s < numStates && nextState >= 0 && nextState < numStates) {
            trans[s][symbol] = nextState;
        }
    }

    // Verifica si un estado es final
    bool accept(int s) const {
        return (s >= 0 && s < numStates) ? isFinal[s] : false;
    }

    // Obtiene el siguiente estado desde 's' con el símbolo 'c'; devuelve -1 si no existe
    int nextState(int s, char c) const {
        if (s < 0 || s >= numStates) return -1;
        auto it = trans[s].find(c);
        if (it == trans[s].end()) return -1;
        return it->second;
    }
}; 