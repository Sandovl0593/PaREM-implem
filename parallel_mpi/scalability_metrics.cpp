#include "scalability_metrics.h"
#include "optimized_parallel_dfa.h"
#include "utils.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <mpi.h>

using namespace std;
using namespace std::chrono;

// Generar string de entrada con patrón específico
string generateTestInput(int size) {
    string input;
    input.reserve(size);
    
    // Patrón: 010010010... para asegurar transiciones de estado
    for (int i = 0; i < size; i++) {
        if (i % 3 == 0) input += '0';
        else if (i % 3 == 1) input += '1';
        else input += '0';
    }
    
    return input;
}

// Crear DFA que acepta strings que contienen "00"
OptimizedDFA createTestDFA() {
    // Constructor: OptimizedDFA(states, initial, accepting, alphabet)
    OptimizedDFA dfa(3, 0, {2}, "01");
    
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
double measureSequentialTime(const OptimizedDFA& dfa, const string& input) {
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
        cout << "Uso: " << argv[0] << " <numero_procesos>" << endl;
        cout << "Ejemplo: " << argv[0] << " 4" << endl;
        return 1;
    }
    
    int target_processes = atoi(argv[1]);
    if (target_processes < 1 || target_processes > 32) {
        cout << "Error: Número de procesos debe estar entre 1 y 32" << endl;
        return 1;
    }
    
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    if (size != target_processes) {
        if (rank == 0) {
            cout << "Error: Se requieren exactamente " << target_processes 
                 << " procesos MPI, pero se detectaron " << size << endl;
        }
        MPI_Finalize();
        return 1;
    }
    
    // Solo el proceso 0 imprime información
    if (rank == 0) {
        cout << "🧮 Generando métricas de escalabilidad para P = " << target_processes << endl;
        cout << "📊 DFA: |Q| = 3, |Σ| = 2, patrón = \"00\"" << endl;
        cout << "📏 Tamaños: N = {100K, 1M, 10M, 100M, 1B} caracteres" << endl;
    }
    
    OptimizedDFA dfa = createTestDFA();
    vector<int> input_sizes = {100000, 1000000, 10000000, 100000000, 1000000000};
    
    // Abrir archivo en modo append
    ofstream outFile;
    if (rank == 0) {
        outFile.open("MetricsGraphics.txt", ios::app);
        if (!outFile.is_open()) {
            cout << "Error: No se pudo abrir MetricsGraphics.txt" << endl;
            MPI_Finalize();
            return 1;
        }
        
        // Escribir encabezado para este conjunto de pruebas
        outFile << "\n=== MÉTRICAS DE ESCALABILIDAD PARA P = " << target_processes << " ===" << endl;
        outFile << "DFA: |Q| = 3, |Σ| = 2, acepta strings con patrón \"00\"" << endl;
        outFile << "Fecha: " << __DATE__ << " " << __TIME__ << endl;
        outFile << "\nFormato: N | T_seq(ms) | T_comp(ms) | T_comm_block(ms) | T_par_block(ms) | "
                << "T_comm_nonblock(ms) | T_par_nonblock(ms) | Speedup_block | Speedup_nonblock | "
                << "Efficiency_block | Efficiency_nonblock" << endl;
        outFile << string(120, '-') << endl;
    }
    
    for (int N : input_sizes) {
        string input = generateTestInput(N);
        double T_seq = 0.0, T_comp = 0.0, T_comm_block = 0.0, T_comm_nonblock = 0.0;
        double T_par_block = 0.0, T_par_nonblock = 0.0;
        
        // Medir tiempo secuencial (solo proceso 0)
        if (rank == 0) {
            T_seq = measureSequentialTime(dfa, input);
        }
        
        // Caso especial para P = 1 (secuencial)
        if (target_processes == 1) {
            T_comp = T_seq;
            T_comm_block = 0.0;
            T_comm_nonblock = 0.0;
            T_par_block = T_seq;
            T_par_nonblock = T_seq;
        } else {
            // Medir versión con comunicación bloqueante
            MPI_Barrier(MPI_COMM_WORLD);
            auto result_block = processOptimizedDFAParallel(dfa, input, size, rank, false);
            T_comp = max(0.001, result_block.computation_time);  // Mínimo 0.001ms
            T_comm_block = max(0.0, result_block.communication_time);
            T_par_block = T_comp + T_comm_block;
            
            // Medir versión con comunicación no bloqueante  
            MPI_Barrier(MPI_COMM_WORLD);
            auto result_nonblock = processOptimizedDFAParallel(dfa, input, size, rank, true);
            T_comm_nonblock = max(0.0, result_nonblock.communication_time);
            T_par_nonblock = T_comp + T_comm_nonblock;  // Usar el mismo T_comp
        }
        
        // Solo el proceso 0 calcula métricas y escribe resultados
        if (rank == 0) {
            // Protección contra división por cero
            double speedup_block = (T_par_block > 0.001) ? T_seq / T_par_block : 1.0;
            double speedup_nonblock = (T_par_nonblock > 0.001) ? T_seq / T_par_nonblock : 1.0;
            double efficiency_block = (target_processes > 0) ? speedup_block / target_processes : 0.0;
            double efficiency_nonblock = (target_processes > 0) ? speedup_nonblock / target_processes : 0.0;
            
            // Escribir resultados
            outFile << fixed << setprecision(3);
            outFile << setw(8) << N << " | "
                    << setw(9) << T_seq << " | "
                    << setw(10) << T_comp << " | "
                    << setw(15) << T_comm_block << " | "
                    << setw(13) << T_par_block << " | "
                    << setw(18) << T_comm_nonblock << " | "
                    << setw(16) << T_par_nonblock << " | "
                    << setw(12) << speedup_block << " | "
                    << setw(15) << speedup_nonblock << " | "
                    << setw(16) << efficiency_block << " | "
                    << setw(18) << efficiency_nonblock << endl;
            
            cout << "✅ N = " << setw(7) << N 
                 << ", Speedup(block) = " << setw(5) << speedup_block
                 << ", Speedup(nonblock) = " << setw(5) << speedup_nonblock << endl;
        }
    }
    
    if (rank == 0) {
        outFile << string(120, '-') << endl;
        outFile.close();
        cout << "📈 Métricas para P = " << target_processes << " guardadas en MetricsGraphics.txt" << endl;
    }
    
    MPI_Finalize();
    return 0;
} 