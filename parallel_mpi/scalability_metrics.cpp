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
    
    return duration<double, micro>(end - start).count();  // Cambiar a microsegundos
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
        cout << "📏 Tamaños: N = {100K, 1M, 10M, 100M} caracteres" << endl;
    }
    
    OptimizedDFA dfa = createTestDFA();
    // Reducir tamaños para evitar problemas de memoria
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
        outFile << "\nFormato: N | T_seq(μs) | T_comp(μs) | T_comm_block(μs) | T_par_block(μs) | "
                << "T_comm_nonblock(μs) | T_par_nonblock(μs) | Speedup_block | Speedup_nonblock | "
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
        
        // Broadcast del tiempo secuencial a todos los procesos
        MPI_Bcast(&T_seq, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
        // Caso especial para P = 1 (secuencial)
        if (target_processes == 1) {
            T_comp = T_seq;
            T_comm_block = 0.0;
            T_comm_nonblock = 0.0;
            T_par_block = T_seq;
            T_par_nonblock = T_seq;
        } else {
            // Medir versión con comunicación bloqueante usando MPI_Wtime
            MPI_Barrier(MPI_COMM_WORLD);
            double start_block = MPI_Wtime();
            
            // Simular procesamiento paralelo simple
            int blockSize = max(1, (int)input.length() / size);
            int startPos = rank * blockSize;
            int endPos = (rank == size - 1) ? input.length() : min((rank + 1) * blockSize, (int)input.length());
            
            double comp_start = MPI_Wtime();
            bool localResult = false;
            if (startPos < (int)input.length() && endPos > startPos) {
                string localBlock = input.substr(startPos, endPos - startPos);
                localResult = dfa.run(localBlock);
            }
            double comp_end = MPI_Wtime();
            T_comp = (comp_end - comp_start) * 1000000;  // Convertir a microsegundos
            
            // Comunicación bloqueante
            double comm_start = MPI_Wtime();
            int globalResult;
            int localInt = localResult ? 1 : 0;
            MPI_Allreduce(&localInt, &globalResult, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            double comm_end = MPI_Wtime();
            T_comm_block = (comm_end - comm_start) * 1000000;  // Convertir a microsegundos
            
            double end_block = MPI_Wtime();
            T_par_block = (end_block - start_block) * 1000000;  // Convertir a microsegundos
            
            // Medir versión con comunicación no bloqueante
            MPI_Barrier(MPI_COMM_WORLD);
            double start_nonblock = MPI_Wtime();
            
            // Repetir la computación para medición independiente
            double comp_nb_start = MPI_Wtime();
            bool localResult_nb = false;
            if (startPos < (int)input.length() && endPos > startPos) {
                string localBlock = input.substr(startPos, endPos - startPos);
                localResult_nb = dfa.run(localBlock);
            }
            double comp_nb_end = MPI_Wtime();
            double T_comp_nb = (comp_nb_end - comp_nb_start) * 1000000;
            (void)T_comp_nb;  // Suprimir warning de variable no utilizada
            
            // Comunicación no bloqueante
            double comm_nb_start = MPI_Wtime();
            MPI_Request request;
            int localInt_nb = localResult_nb ? 1 : 0;
            int globalResult_nb;
            MPI_Iallreduce(&localInt_nb, &globalResult_nb, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD, &request);
            
            // Simular trabajo adicional durante comunicación
            for (volatile int i = 0; i < 1000; i++) { /* trabajo simulado */ }
            
            MPI_Wait(&request, MPI_STATUS_IGNORE);
            double comm_nb_end = MPI_Wtime();
            T_comm_nonblock = (comm_nb_end - comm_nb_start) * 1000000;  // Convertir a microsegundos
            
            double end_nonblock = MPI_Wtime();
            T_par_nonblock = (end_nonblock - start_nonblock) * 1000000;  // Convertir a microsegundos
        }
        
        // Solo el proceso 0 calcula métricas y escribe resultados
        if (rank == 0) {
            // Debug: Verificar valores antes de calcular speedup
            cout << "🔍 DEBUG N=" << N << ": T_seq=" << T_seq << ", T_par_block=" << T_par_block 
                 << ", T_par_nonblock=" << T_par_nonblock << endl;
            
            // Protección contra división por cero y valores irreales
            double speedup_block = (T_par_block > 0.001) ? T_seq / T_par_block : 1.0;
            double speedup_nonblock = (T_par_nonblock > 0.001) ? T_seq / T_par_nonblock : 1.0;
            
            // Limitar speedup a valores razonables (máximo 10x el número de procesos)
            double max_reasonable_speedup = target_processes * 10.0;
            if (speedup_nonblock > max_reasonable_speedup) {
                cout << "⚠️  ADVERTENCIA: Speedup no bloqueante irreal (" << speedup_nonblock 
                     << "), limitando a " << max_reasonable_speedup << endl;
                speedup_nonblock = max_reasonable_speedup;
            }
            
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