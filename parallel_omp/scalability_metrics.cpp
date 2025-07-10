#include "dfa.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <omp.h>
#include <random>

using namespace std;
using namespace std::chrono;

// Función para generar cadenas aleatorias
std::string generateRandomString(const std::string& alphabet, int length) {
    std::string result;
    result.reserve(length);

    std::mt19937 rng(static_cast<unsigned int>(std::chrono::system_clock::now().time_since_epoch().count()));
    std::uniform_int_distribution<int> dist(0, alphabet.length() - 1);

    for (int i = 0; i < length; ++i) {
        result.push_back(alphabet[dist(rng)]);
    }
    return result;
} 

// Crear DFA que acepta strings que contienen "00"
DFA example_DFA() {
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

// Medir tiempo secuencial
double measureSequentialTime(const DFA& dfa, const string& input) {
    auto start = high_resolution_clock::now();
    
    // Usar el método run() público de la clase
    bool result = dfa.run(input);
    
    auto end = high_resolution_clock::now();
    
    // Evitar optimización del compilador
    volatile bool dummy = result;
    (void)dummy;
    
    return duration<double, milli>(end - start).count();
}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cout << "Uso: " << argv[0] << " <num_threads>" << endl;
        cout << "Ejemplo: " << argv[0] << " 4" << endl;
        return 1;
    }
    int num_threads = atoi(argv[1]);
    if (num_threads < 1 || num_threads > 32) {
        cout << "Error: Número de hilos debe estar entre 1 y 32" << endl;
        return 1;
    }
    
    omp_set_num_threads(num_threads);
    
    // Información en secuencial
    cout << "🧮 Generando métricas de escalabilidad para " << num_threads << " hilos" << endl;
    cout << "📊 DFA: |Q| = 3, |Σ| = 2, patrón = \"00\"" << endl;
    cout << "📏 Tamaños: N = {100K, 1M, 10M, 100M} caracteres" << endl;
    
    DFA dfa = example_DFA();
    vector<int> input_sizes = {100000, 1000000, 10000000, 100000000};
    
    // Abrir archivo en modo append
    ofstream outFile;
    outFile.open("scalabilityMetrics.txt", ios::app);
    if (!outFile.is_open()) {
        cout << "Error: No se pudo abrir scalabilityMetrics.txt" << endl;
        return 1;
    }
    
    // Escribir encabezado para este conjunto de pruebas
    outFile << "\n=== MÉTRICAS DE ESCALABILIDAD PARA " << num_threads << " HILOS ===" << endl;
    outFile << "DFA: |Q| = 3, |Σ| = 2, acepta strings con patrón \"00\"" << endl;
    outFile << "\nFormato: N | T_seq(ms) | T_comp(theo) | T_comp(ms) | T_jum_opt(ns) | T_comp_opt(theo) | T_comp_opt(ms) | Speedup | Speedup_opt | "
            << "Efficiency | Efficiency_opt" << endl;
    outFile << string(120, '-') << endl;
    
    for (int N : input_sizes) {
        string input = generateRandomString("01", N);
        double T_seq = 0.0, T_comp = 0.0, T_comp_opt = 0.0, T_jum_opt = 0.0;

        // Medir tiempo secuencial (solo proceso 0)
        T_seq = measureSequentialTime(dfa, input);
        
        // Caso solo para más de un hilo
        if (num_threads != 1) {
            // Medir versión con comunicación bloqueante
            auto result_ = processOMPParallel(dfa, input, num_threads);
            T_comp = result_.computationTime;
          
            auto result_opt = runOptimizedOMPParallel(dfa, input, num_threads);
            T_jum_opt = result_opt.preComputeJumpTime;
            T_comp_opt = result_opt.computationTime;
        }
        
        // Protección contra división por cero
        double speedup = (T_comp > 0.001) ? T_seq / T_comp : 1.0;
        double speedup_opt = (T_comp_opt > 0.001) ? T_seq / T_comp_opt : 1.0;
        double efficiency = (num_threads > 0) ? speedup / num_threads : 0.0;
        double efficiency_opt = (num_threads > 0) ? speedup_opt / num_threads : 0.0;

        int p = num_threads;
        int Q = dfa.getNumStates();
        int sigma = dfa.getAlphabet().size();
        double Qsigma = double(Q) * double(sigma);

        // --- 1) Modelo teórico ---
        double T_theo = (Qsigma + N)/p + N;
        double T_theo_ompt = (Qsigma + N)/p + log2(double(p)) + double(Q) * log2(double(p));

        // Métricas reales ---
        double real_speedup    = T_seq / T_comp;
        double real_efficiency = real_speedup / p * 100.0;

        // Métricas teóricas ---
        double theo_speedup    = T_seq / T_theo;
        double theo_efficiency = theo_speedup / p * 100.0;
        
        // Escribir resultados
        outFile << fixed << setprecision(3);
        outFile << setw(8) << N << " | "
                << setw(9) << T_seq << " | "
                << setw(15) << T_theo << " | "
                << setw(10) << T_comp << " | "
                << setw(15) << T_jum_opt << " | "
                << setw(15) << T_theo_ompt << " | "
                << setw(13) << T_comp_opt << " | "
                << setw(18) << speedup << " | "
                << setw(16) << speedup_opt << " | "
                << setw(12) << efficiency << " | "
                << setw(18) << efficiency_opt << endl;
        cout << "✅ N = " << setw(7) << N 
             << ", Speedup(first) = " << setw(5) << speedup
             << ", Speedup(opt) = " << setw(5) << speedup_opt << endl;
    }
    
    outFile << string(120, '-') << endl;
    outFile.close();
    cout << "📈 Métricas para " << num_threads << " hilos guardadas en scalabilityMetrics.txt" << endl;
    
    return 0;
} 