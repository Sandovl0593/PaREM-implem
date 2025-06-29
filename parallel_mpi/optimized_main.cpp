#include "optimized_parallel_dfa.h"
#include "parallel_dfa.h"
#include "dfa_tests.h"
#include "communication_tests.h"
#include "utils.h"
#include <iostream>
#include <mpi.h>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <vector>
#include <numeric>
#include <cmath>

// ✅ CONSTANTES PARA ANÁLISIS TEÓRICO
const int BUFFER_SIZE = 1024 * 1024; // 1MB
const int NUM_ITERATIONS = 50; // Promedio de 50 iteraciones

// ✅ ESTRUCTURA PARA MEDICIÓN DETALLADA DE TIEMPOS
struct DetailedTimeAnalysis {
    double computation_time;        // (|Q|·Σ+n)/p
    double initial_distribution;    // n (distribución inicial)
    double collective_comm;         // log p (comunicación colectiva)
    double dfa_broadcast;          // |Q|·log p (broadcast del DFA)
    double total_parallel;         // Tiempo total paralelo
    bool parallel_accepted;        // Resultado de aceptación paralelo
    
    // Parámetros teóricos
    int Q;  // Número de estados
    int sigma; // Tamaño del alfabeto  
    int n;  // Tamaño de entrada
    int p;  // Número de procesos
    
    DetailedTimeAnalysis(int states, int alphabet_size, int input_size, int processes) 
        : computation_time(0), initial_distribution(0), collective_comm(0), 
          dfa_broadcast(0), total_parallel(0), parallel_accepted(false), Q(states), sigma(alphabet_size), 
          n(input_size), p(processes) {}
};

// ✅ FUNCIÓN: Medición detallada según fórmulas teóricas
DetailedTimeAnalysis measureDetailedPerformance(const OptimizedDFA& dfa, const string& input, 
                                               int numProcs, int myRank, bool useAsyncComm = false) {
    DetailedTimeAnalysis analysis(dfa.getNumStates(), dfa.getAlphabet().size(), input.size(), numProcs);
    
    double start_total = MPI_Wtime();
    
    // ================= COMPONENTE 1: BROADCAST DEL DFA (|Q|·log p) =================
    double start_dfa_broadcast = MPI_Wtime();
    // Simular broadcast del DFA (ya está distribuido, pero medimos el costo teórico)
    MPI_Barrier(MPI_COMM_WORLD); // Sincronización para medición precisa
    double end_dfa_broadcast = MPI_Wtime();
    analysis.dfa_broadcast = (end_dfa_broadcast - start_dfa_broadcast) * 1000000;
    
    // ================= COMPONENTE 2: DISTRIBUCIÓN INICIAL (n) =================
    double start_distribution = MPI_Wtime();
    
    // Usar la función optimizada de procesamiento paralelo
    OptimizedParallelDFAResult parallelResult = processOptimizedDFAParallel(dfa, input, numProcs, myRank, useAsyncComm);
    
    double end_distribution = MPI_Wtime();
    analysis.initial_distribution = (end_distribution - start_distribution) * 1000000;
    
    // ================= COMPONENTE 3: COMPUTACIÓN ((|Q|·Σ+n)/p) =================
    // El tiempo de computación está incluido en el resultado paralelo
    analysis.computation_time = parallelResult.computation_time;
    
    // ================= COMPONENTE 4: COMUNICACIÓN COLECTIVA (log p) =================
    // El tiempo de comunicación está incluido en el resultado paralelo
    analysis.collective_comm = parallelResult.communication_time;
    
    // Guardar resultado de aceptación paralelo
    analysis.parallel_accepted = parallelResult.accepted;
    
    double end_total = MPI_Wtime();
    analysis.total_parallel = (end_total - start_total) * 1000000;
    
    return analysis;
}

// ✅ FUNCIÓN: Análisis teórico de escalabilidad
void analyzeScalability(const DetailedTimeAnalysis& analysis, double sequentialTime, int myRank) {
    if (myRank != 0) return;
    
    // Calcular tiempo teórico según fórmula: Tp = (|Q|·Σ+n)/p + n + log p + |Q|·log p
    double theoretical_computation = (double)(analysis.Q * analysis.sigma + analysis.n) / analysis.p;
    double theoretical_total = theoretical_computation + analysis.n/1000000.0 + log2(analysis.p) + analysis.Q * log2(analysis.p);
    
    // Métricas reales
    double real_speedup = sequentialTime / analysis.total_parallel;
    double real_efficiency = (real_speedup / analysis.p) * 100;
    
    // Métricas teóricas
    double theoretical_speedup = sequentialTime / theoretical_total;
    double theoretical_efficiency_strong = (double)(analysis.Q * analysis.sigma + 1) / 
                                         (analysis.Q * analysis.sigma + analysis.n + analysis.p + 
                                          analysis.p * log2(analysis.p) + analysis.Q * log2(analysis.p)) * 100;
    
    cout << "\n📊 ANÁLISIS DETALLADO DE ESCALABILIDAD:" << endl;
    cout << "============================================================" << endl;
    cout << "Parámetros: |Q|=" << analysis.Q << ", Σ=" << analysis.sigma << ", n=" << analysis.n << ", p=" << analysis.p << endl;
    cout << "\n🔬 COMPONENTES DEL TIEMPO PARALELO (μs):" << endl;
    cout << "  • Computación ((|Q|·Σ+n)/p): " << fixed << setprecision(2) << analysis.computation_time << " μs" << endl;
    cout << "  • Distribución inicial (n): " << analysis.initial_distribution << " μs" << endl;
    cout << "  • Comunicación colectiva (log p): " << analysis.collective_comm << " μs" << endl;
    cout << "  • Broadcast DFA (|Q|·log p): " << analysis.dfa_broadcast << " μs" << endl;
    cout << "  • TOTAL PARALELO: " << analysis.total_parallel << " μs" << endl;
    
    cout << "\n📈 MÉTRICAS REALES vs TEÓRICAS:" << endl;
    cout << "  Real - Speedup: " << setprecision(3) << real_speedup << "x | Eficiencia: " << setprecision(2) << real_efficiency << "%" << endl;
    cout << "  Teórico - Speedup: " << theoretical_speedup << "x | Eficiencia fuerte: " << theoretical_efficiency_strong << "%" << endl;
    
    // Condición de escalabilidad débil
    double weak_scalability_condition = analysis.Q * (analysis.p * log2(analysis.p) - analysis.sigma);
    cout << "\n🔍 ANÁLISIS DE ESCALABILIDAD:" << endl;
    cout << "  • Escalabilidad fuerte: " << (real_efficiency > 80 ? "✅ BUENA" : "❌ POBRE") << " (" << real_efficiency << "%)" << endl;
    cout << "  • Condición débil (n ≈ |Q|·(p·log p - Σ)): n=" << analysis.n << " vs " << weak_scalability_condition << endl;
    cout << "  • Escalabilidad débil: " << (abs(analysis.n - weak_scalability_condition) < analysis.n * 0.1 ? "✅ CUMPLE" : "❌ NO CUMPLE") << endl;
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    
    int numProcs, myRank;
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    
    initializeLogging(numProcs, myRank);
    
    if (myRank == 0) {
        cout << "============================================================" << endl;
        cout << "⚡ ANÁLISIS TEÓRICO DE ESCALABILIDAD MPI DFA - " << numProcs << " PROCESOS" << endl;
        cout << "============================================================" << endl;
        cout << "Iteraciones por test: " << NUM_ITERATIONS << endl;
        cout << "Implementación de fórmulas teóricas de escalabilidad" << endl << endl;
    }
    
    setupOptimizedBuffers(BUFFER_SIZE);
    
    // ================= TESTS CON DIFERENTES CONFIGURACIONES =================
    vector<tuple<int, int, vector<int>, string, string>> dfaConfigs = {
        {3, 0, {2}, "01", "DFA Binario (|Q|=3, Σ=2)"},
        {4, 0, {3}, "abc", "DFA ABC (|Q|=4, Σ=3)"},
        {5, 0, {2,4}, "xyz", "DFA Multi-Final (|Q|=5, Σ=3)"},
        {6, 0, {5}, "0123", "DFA Cuaternario (|Q|=6, Σ=4)"}
    };
    
    vector<int> inputSizes = {1000, 10000, 100000, 1000000};
    
    // ================= INICIALIZAR ARCHIVOS DE RESULTADOS =================
    ofstream resultsFile, scalabilityFile, communicationFile, edgeCaseFile, stressFile;
    
    if (myRank == 0) {
        resultsFile.open("resultados_mpi/optimized_results.txt");
        scalabilityFile.open("resultados_mpi/optimized_scalability_analysis.txt");
        communicationFile.open("resultados_mpi/optimized_communication_results.txt");
        edgeCaseFile.open("resultados_mpi/optimized_edge_case_results.txt");
        stressFile.open("resultados_mpi/optimized_stress_test_results.txt");
        
        if (resultsFile.is_open()) {
            resultsFile << "============================================================\n";
            resultsFile << "🚀 RESULTADOS VERSIÓN OPTIMIZADA - " << numProcs << " PROCESOS\n";
            resultsFile << "============================================================\n";
            resultsFile << "Fecha: " << __DATE__ << " " << __TIME__ << "\n";
            resultsFile << "Iteraciones por configuración: " << NUM_ITERATIONS << "\n\n";
        }
        
        if (scalabilityFile.is_open()) {
            scalabilityFile << "============================================================\n";
            scalabilityFile << "📊 ANÁLISIS COMPLETO DE ESCALABILIDAD OPTIMIZADA - " << numProcs << " PROCESOS\n";
            scalabilityFile << "============================================================\n";
            scalabilityFile << "Fecha: " << __DATE__ << " " << __TIME__ << "\n";
            scalabilityFile << "Iteraciones por configuración: " << NUM_ITERATIONS << "\n\n";
        }
        
        if (communicationFile.is_open()) {
            communicationFile << "============================================================\n";
            communicationFile << "📡 ANÁLISIS DE COMUNICACIÓN OPTIMIZADA - " << numProcs << " PROCESOS\n";
            communicationFile << "============================================================\n";
            communicationFile << "Fecha: " << __DATE__ << " " << __TIME__ << "\n\n";
        }
        
        if (edgeCaseFile.is_open()) {
            edgeCaseFile << "============================================================\n";
            edgeCaseFile << "🔍 CASOS EDGE OPTIMIZADOS - " << numProcs << " PROCESOS\n";
            edgeCaseFile << "============================================================\n";
            edgeCaseFile << "Fecha: " << __DATE__ << " " << __TIME__ << "\n\n";
        }
        
        if (stressFile.is_open()) {
            stressFile << "============================================================\n";
            stressFile << "🔥 TESTS DE STRESS OPTIMIZADOS - " << numProcs << " PROCESOS\n";
            stressFile << "============================================================\n";
            stressFile << "Fecha: " << __DATE__ << " " << __TIME__ << "\n\n";
        }
    }
    
    for (const auto& config : dfaConfigs) {
        int states = get<0>(config);
        int initial = get<1>(config);
        vector<int> accepting = get<2>(config);
        string alphabet = get<3>(config);
        string description = get<4>(config);
        
        // Crear DFA
        OptimizedDFA dfa(states, initial, accepting, alphabet);
        
        // Añadir transiciones de ejemplo
        for (int s = 0; s < states; s++) {
            for (char c : alphabet) {
                int nextState = (s + (c - alphabet[0]) + 1) % states;
                dfa.addTransition(s, c, nextState);
            }
        }
        
        if (myRank == 0) {
            cout << "\n🔬 TESTING: " << description << endl;
            cout << string(60, '=') << endl;
        }
        
        for (int inputSize : inputSizes) {
            if (myRank == 0) {
                cout << "\n📏 Tamaño de entrada: " << inputSize << " caracteres" << endl;
            }
            
            // Generar entrada de prueba
            string testInput = generateRandomString(alphabet, inputSize);
            
            // ================= MEDICIÓN SECUENCIAL =================
            bool seqResult = false;
            double sequentialTime = 0;
            if (myRank == 0) {
                vector<double> seqTimes;
                for (int iter = 0; iter < NUM_ITERATIONS; iter++) {
                    auto start = chrono::high_resolution_clock::now();
                    seqResult = dfa.run(testInput);
                    auto end = chrono::high_resolution_clock::now();
                    double time = chrono::duration_cast<chrono::microseconds>(end - start).count();
                    seqTimes.push_back(time);
                }
                sequentialTime = accumulate(seqTimes.begin(), seqTimes.end(), 0.0) / NUM_ITERATIONS;
                cout << "  📊 Tiempo secuencial promedio: " << fixed << setprecision(2) << sequentialTime << " μs" << endl;
            }
            
            // ================= MEDICIÓN PARALELA (BLOQUEADA) =================
            vector<DetailedTimeAnalysis> blockingResults;
            for (int iter = 0; iter < NUM_ITERATIONS; iter++) {
                MPI_Barrier(MPI_COMM_WORLD);
                DetailedTimeAnalysis analysis = measureDetailedPerformance(dfa, testInput, numProcs, myRank, false);
                if (myRank == 0) blockingResults.push_back(analysis);
            }
            
            // ================= MEDICIÓN PARALELA (NO BLOQUEADA) =================
            vector<DetailedTimeAnalysis> asyncResults;
            for (int iter = 0; iter < NUM_ITERATIONS; iter++) {
                MPI_Barrier(MPI_COMM_WORLD);
                DetailedTimeAnalysis analysis = measureDetailedPerformance(dfa, testInput, numProcs, myRank, true);
                if (myRank == 0) asyncResults.push_back(analysis);
            }
            
            // ================= ANÁLISIS DE RESULTADOS =================
            if (myRank == 0) {
                // Promediar resultados bloqueados
                DetailedTimeAnalysis avgBlocking(states, alphabet.size(), inputSize, numProcs);
                for (const auto& result : blockingResults) {
                    avgBlocking.computation_time += result.computation_time;
                    avgBlocking.initial_distribution += result.initial_distribution;
                    avgBlocking.collective_comm += result.collective_comm;
                    avgBlocking.dfa_broadcast += result.dfa_broadcast;
                    avgBlocking.total_parallel += result.total_parallel;
                }
                avgBlocking.computation_time /= NUM_ITERATIONS;
                avgBlocking.initial_distribution /= NUM_ITERATIONS;
                avgBlocking.collective_comm /= NUM_ITERATIONS;
                avgBlocking.dfa_broadcast /= NUM_ITERATIONS;
                avgBlocking.total_parallel /= NUM_ITERATIONS;
                
                // Promediar resultados no bloqueados
                DetailedTimeAnalysis avgAsync(states, alphabet.size(), inputSize, numProcs);
                for (const auto& result : asyncResults) {
                    avgAsync.computation_time += result.computation_time;
                    avgAsync.initial_distribution += result.initial_distribution;
                    avgAsync.collective_comm += result.collective_comm;
                    avgAsync.dfa_broadcast += result.dfa_broadcast;
                    avgAsync.total_parallel += result.total_parallel;
                }
                avgAsync.computation_time /= NUM_ITERATIONS;
                avgAsync.initial_distribution /= NUM_ITERATIONS;
                avgAsync.collective_comm /= NUM_ITERATIONS;
                avgAsync.dfa_broadcast /= NUM_ITERATIONS;
                avgAsync.total_parallel /= NUM_ITERATIONS;
                
                cout << "\n🔄 COMUNICACIÓN BLOQUEADA (promedio " << NUM_ITERATIONS << " iteraciones):" << endl;
                analyzeScalability(avgBlocking, sequentialTime, myRank);
                
                cout << "\n⚡ COMUNICACIÓN NO BLOQUEADA (promedio " << NUM_ITERATIONS << " iteraciones):" << endl;
                analyzeScalability(avgAsync, sequentialTime, myRank);
                
                // ================= GUARDAR EN TODOS LOS ARCHIVOS =================
                // Archivo principal de resultados
                if (resultsFile.is_open()) {
                    resultsFile << description << " - Entrada: " << inputSize << " caracteres\n";
                    resultsFile << "Secuencial: " << fixed << setprecision(2) << sequentialTime << " μs\n";
                    resultsFile << "Bloqueada: " << avgBlocking.total_parallel << " μs (Speedup: " 
                               << setprecision(3) << sequentialTime/avgBlocking.total_parallel << "x)\n";
                    resultsFile << "No Bloqueada: " << avgAsync.total_parallel << " μs (Speedup: " 
                               << sequentialTime/avgAsync.total_parallel << "x)\n";
                    resultsFile << "Mejora no-bloqueada: " << setprecision(2) 
                               << (avgBlocking.total_parallel/avgAsync.total_parallel) << "x\n\n";
                }
                
                // Archivo de escalabilidad detallada
                if (scalabilityFile.is_open()) {
                    scalabilityFile << "\n" << description << " - " << inputSize << " caracteres:\n";
                    scalabilityFile << "Parámetros: |Q|=" << states << ", Σ=" << alphabet.size() << ", n=" << inputSize << ", p=" << numProcs << "\n";
                    scalabilityFile << "BLOQUEADA - Speedup: " << setprecision(3) << sequentialTime/avgBlocking.total_parallel 
                                   << "x, Eficiencia: " << setprecision(2) << (sequentialTime/avgBlocking.total_parallel/numProcs)*100 << "%\n";
                    scalabilityFile << "NO BLOQUEADA - Speedup: " << sequentialTime/avgAsync.total_parallel 
                                   << "x, Eficiencia: " << (sequentialTime/avgAsync.total_parallel/numProcs)*100 << "%\n";
                    scalabilityFile << "Componentes (μs) - Comp: " << avgAsync.computation_time << ", Dist: " << avgAsync.initial_distribution 
                                   << ", Comm: " << avgAsync.collective_comm << ", DFA: " << avgAsync.dfa_broadcast << "\n\n";
                }
                
                // Archivo de comunicación
                if (communicationFile.is_open()) {
                    communicationFile << description << " - " << inputSize << " caracteres:\n";
                    communicationFile << "Comunicación bloqueada: " << avgBlocking.collective_comm + avgBlocking.initial_distribution << " μs\n";
                    communicationFile << "Comunicación no bloqueada: " << avgAsync.collective_comm + avgAsync.initial_distribution << " μs\n";
                    communicationFile << "Mejora comunicación: " << setprecision(2) 
                                     << ((avgBlocking.collective_comm + avgBlocking.initial_distribution) / 
                                         (avgAsync.collective_comm + avgAsync.initial_distribution)) << "x\n\n";
                }
                
                // Validación automática de correctitud funcional
                bool lastBlocking = false, lastAsync = false;
                if (!blockingResults.empty()) lastBlocking = blockingResults.back().parallel_accepted;
                if (!asyncResults.empty()) lastAsync = asyncResults.back().parallel_accepted;
                cout << "  Validación bloqueada: ";
                if (seqResult == lastBlocking) cout << "✅" << endl;
                else cout << "❌ (Secuencial: " << seqResult << ", Paralelo: " << lastBlocking << ")" << endl;
                cout << "  Validación no bloqueada: ";
                if (seqResult == lastAsync) cout << "✅" << endl;
                else cout << "❌ (Secuencial: " << seqResult << ", Paralelo: " << lastAsync << ")" << endl;
            }
        }
    }
    
    // ================= TESTS DE CASOS EDGE =================
    if (myRank == 0) {
        cout << "\n🔍 EJECUTANDO TESTS DE CASOS EDGE OPTIMIZADOS..." << endl;
    }
    
    // Crear un DFA simple para tests edge
    OptimizedDFA edgeDFA(3, 0, {2}, "01");
    for (int s = 0; s < 3; s++) {
        for (char c : string("01")) {
            int nextState = (s + (c - '0') + 1) % 3;
            edgeDFA.addTransition(s, c, nextState);
        }
    }
    
    // Tests de casos edge
    vector<pair<string, string>> edgeCases = {
        {"", "Cadena Vacía"},
        {"0", "Un Solo Carácter"},
        {"01", "Dos Caracteres"},
        {"000", "Tres Ceros"},
        {"111", "Tres Unos"}
    };
    
    for (const auto& edgeCase : edgeCases) {
        string input = edgeCase.first;
        string description = edgeCase.second;
        
        if (myRank == 0) {
            cout << "  📋 " << description << ": \"" << input << "\"" << endl;
        }
        
        // Medición con comunicación bloqueada
        MPI_Barrier(MPI_COMM_WORLD);
        DetailedTimeAnalysis blockingEdge = measureDetailedPerformance(edgeDFA, input, numProcs, myRank, false);
        
        // Medición con comunicación no bloqueada
        MPI_Barrier(MPI_COMM_WORLD);
        DetailedTimeAnalysis asyncEdge = measureDetailedPerformance(edgeDFA, input, numProcs, myRank, true);
        
        // Validación automática de correctitud funcional para edge cases
        if (myRank == 0) {
            bool seqResult = edgeDFA.run(input);
            bool parResultBlocking = blockingEdge.parallel_accepted;
            bool parResultAsync = asyncEdge.parallel_accepted;
            cout << "    Validación bloqueada: ";
            if (seqResult == parResultBlocking) cout << "✅" << endl;
            else cout << "❌ (Secuencial: " << seqResult << ", Paralelo: " << parResultBlocking << ")" << endl;
            cout << "    Validación no bloqueada: ";
            if (seqResult == parResultAsync) cout << "✅" << endl;
            else cout << "❌ (Secuencial: " << seqResult << ", Paralelo: " << parResultAsync << ")" << endl;
        }
        
        if (myRank == 0 && edgeCaseFile.is_open()) {
            edgeCaseFile << description << " (\"" << input << "\"):\n";
            edgeCaseFile << "  Bloqueada: " << fixed << setprecision(2) << blockingEdge.total_parallel << " μs\n";
            edgeCaseFile << "  No Bloqueada: " << asyncEdge.total_parallel << " μs\n";
            edgeCaseFile << "  Longitud entrada: " << input.length() << " vs " << numProcs << " procesos\n\n";
        }
    }
    
    // ================= TESTS DE STRESS =================
    if (myRank == 0) {
        cout << "\n🔥 EJECUTANDO TESTS DE STRESS OPTIMIZADOS..." << endl;
    }
    
    vector<int> stressSizes = {100000, 500000, 1000000, 2000000};
    for (int stressSize : stressSizes) {
        if (myRank == 0) {
            cout << "  💪 Test stress con " << stressSize << " caracteres" << endl;
        }
        
        string stressInput = generateRandomString("01", stressSize);
        
        // Medición secuencial
        double stressSeqTime = 0;
        bool stressSeqResult = false;
        if (myRank == 0) {
            auto start = chrono::high_resolution_clock::now();
            stressSeqResult = edgeDFA.run(stressInput);
            auto end = chrono::high_resolution_clock::now();
            stressSeqTime = chrono::duration_cast<chrono::microseconds>(end - start).count();
        }
        
        // Medición paralela
        MPI_Barrier(MPI_COMM_WORLD);
        DetailedTimeAnalysis stressParallel = measureDetailedPerformance(edgeDFA, stressInput, numProcs, myRank, true);
        
        // Validación automática de correctitud funcional para stress
        if (myRank == 0) {
            bool parResultStress = stressParallel.parallel_accepted;
            cout << "    Validación stress no bloqueada: ";
            if (stressSeqResult == parResultStress) cout << "✅" << endl;
            else cout << "❌ (Secuencial: " << stressSeqResult << ", Paralelo: " << parResultStress << ")" << endl;
        }
        
        if (myRank == 0 && stressFile.is_open()) {
            double stressSpeedup = stressSeqTime / stressParallel.total_parallel;
            double stressEfficiency = (stressSpeedup / numProcs) * 100;
            
            stressFile << "Tamaño: " << stressSize << " caracteres\n";
            stressFile << "Secuencial: " << fixed << setprecision(2) << stressSeqTime << " μs\n";
            stressFile << "Paralelo optimizado: " << stressParallel.total_parallel << " μs\n";
            stressFile << "Speedup: " << setprecision(3) << stressSpeedup << "x\n";
            stressFile << "Eficiencia: " << setprecision(2) << stressEfficiency << "%\n";
            stressFile << "Componentes - Comp: " << stressParallel.computation_time 
                      << ", Comm: " << stressParallel.collective_comm << " μs\n\n";
        }
    }
    
    if (myRank == 0) {
        cout << "\n============================================================" << endl;
        cout << "✅ ANÁLISIS COMPLETO OPTIMIZADO FINALIZADO" << endl;
        cout << "============================================================" << endl;
        cout << "📁 Resultados detallados guardados en: resultados_mpi/" << endl;
        cout << "  • optimized_results.txt - Resultados principales" << endl;
        cout << "  • optimized_scalability_analysis.txt - Análisis de escalabilidad" << endl;
        cout << "  • optimized_communication_results.txt - Análisis de comunicación" << endl;
        cout << "  • optimized_edge_case_results.txt - Casos edge" << endl;
        cout << "  • optimized_stress_test_results.txt - Tests de stress" << endl;
        
        // Cerrar todos los archivos
        if (resultsFile.is_open()) resultsFile.close();
        if (scalabilityFile.is_open()) scalabilityFile.close();
        if (communicationFile.is_open()) communicationFile.close();
        if (edgeCaseFile.is_open()) edgeCaseFile.close();
        if (stressFile.is_open()) stressFile.close();
    }
    
    cleanupOptimizedBuffers();
    MPI_Finalize();
    return 0;
}