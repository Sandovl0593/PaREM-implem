#include "dfa.h"
#include <chrono>


// ✅ FUNCIÓN PARA ANÁLISIS SECUENCIAL
void runSecuencialAnalysis(const DFA& dfa, const string& input, const string& testName) {
    cout << "\n=== ANÁLISIS DETALLADO: " << testName << " ===\n";
    cout << "Entrada: [" << input << "]\n";
    cout << "Tamaño: " << input.size() << " caracteres\n";
    
    auto start = chrono::high_resolution_clock::now();
    auto trackingResult = dfa.run_with_tracking(input);
    bool result = trackingResult.first;
    vector<int> stateSequence = trackingResult.second;
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(end - start);
    
    cout << "Tiempo de ejecución: " << duration.count() << " μs\n";
    cout << "Resultado: " << (result ? "ACEPTADA" : "RECHAZADA") << "\n";
    
    // Mostrar ruta de estados real del tracking
    cout << "Ruta de estados: ";
    for (size_t i = 0; i < stateSequence.size(); i++) {
        cout << stateSequence[i];
        if (i < stateSequence.size() - 1) cout << " -> ";
    }
    cout << "\n";
    
    if (!stateSequence.empty()) {
        int finalState = stateSequence.back();
        cout << "Estado final: " << finalState << "\n";
        cout << "Es estado de aceptación: " << (dfa.isAccepting(finalState) ? "SÍ" : "NO") << "\n";
    }
} 

// Crear DFA que acepta strings que contienen "00"
DFA createTestDFA() {
    // Constructor: DFA(states, initial, accepting, alphabet)
    DFA dfa(3, 0, {2}, "01");
    
    // Agregar transiciones para detectar patrón "00"
    // Estado 0: inicial
    dfa.addTransition(0, '0', 1);  // 0 -> 1 con '0'
    dfa.addTransition(0, '1', 0);  // 0 -> 0 con '1'
    
    // Estado 1: vio un '0'
    dfa.addTransition(1, '0', 2);  // 1 -> 2 con '0' (encontró "00")
    dfa.addTransition(1, '1', 0);  // 1 -> 0 con '1'
    
    // Estado 2: vio "00" (estado final)
    dfa.addTransition(2, '0', 2);  // 2 -> 2 con '0'
    dfa.addTransition(2, '1', 0);  // 2 -> 0 con '1'
    
    return dfa;
}


int main() {
    // Crear DFA de prueba
    DFA dfa = createTestDFA();
    
    // Definir entrada de prueba
    string input = "0010101001100";  // Contiene el patrón "00"
    
    // Ejecutar análisis secuencial
    runSecuencialAnalysis(dfa, input, "Prueba de DFA con patrón '00'");
    
    return 0;
}