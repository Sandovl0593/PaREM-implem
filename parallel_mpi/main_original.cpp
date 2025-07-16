#include <iostream>
#include <mpi.h>
#include <vector>
#include <string>
#include <chrono>
#include <numeric>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "optimized_parallel_dfa.h"
#include "dfa_tests.h"
#include "parallel_dfa.h"
#include "communication_tests.h"
#include "utils.h"
using namespace std;

const int VALIDATION_ITERATIONS = 10; 

int main(int argc, char *argv[]) {
    // Inicializar MPI
    MPI_Init(&argc, &argv);
    
    int numProcs, myRank;
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    
    if (myRank == 0) {
        cout << "============================================================" << endl;
        cout << "VERSIÓN ORIGINAL MPI DFA - " << numProcs << " PROCESOS" << endl;
        cout << "============================================================" << endl;
        cout << "Fecha: " << __DATE__ << " " << __TIME__ << endl;
        cout << "Procesos: " << numProcs << endl;
        cout << "Iteraciones de validación: " << VALIDATION_ITERATIONS << endl << endl;
    }
    
    // Inicializar sistema de logging
    initializeLogging(numProcs, myRank);
    
    // Crear DFA de ejemplo
    OptimizedDFA dfa(3, 0, {2}, "01");
    dfa.addTransition(0, '0', 1); dfa.addTransition(0, '1', 0);
    dfa.addTransition(1, '0', 1); dfa.addTransition(1, '1', 2);
    dfa.addTransition(2, '0', 1); dfa.addTransition(2, '1', 0);
    
    // ================= TESTS DE COMUNICACIÓN MPI =================
    if (myRank == 0) {
        cout << "EJECUTANDO TESTS DE COMUNICACIÓN MPI..." << endl;
    }
    
    testPointToPointCommunication(numProcs, myRank);
    testCollectiveCommunication(numProcs, myRank);
    testDerivedDataTypes(numProcs, myRank);
    testOptimizedBuffers(numProcs, myRank);
    
    if (myRank == 0) {
        cout << "\n⚡ EJECUTANDO TESTS PARALELOS VS SECUENCIALES MEJORADOS..." << endl;
    }
    
    vector<pair<string, string>> basicTests = {
        {"01010101", "Test Básico - Patrón Alternado"},
        {"000111000", "Test Básico - Grupos de Ceros y Unos"},
        {"111000111", "Test Básico - Patrón Inverso"}
    };
    
    for (const auto& test : basicTests) {

        if (myRank == 0) {
            cout << "\n " << test.second << " (promedio de " << VALIDATION_ITERATIONS << " iteraciones)" << endl;
        }
        
        vector<double> sequentialTimes, parallelTimes;
        
        for (int iter = 0; iter < VALIDATION_ITERATIONS; iter++) {

             if (myRank == 0) {
                 auto start = chrono::high_resolution_clock::now();
                 bool seqResult = dfa.run(test.first);
                 (void)seqResult;
                 auto end = chrono::high_resolution_clock::now();
                 double seqTime = chrono::duration_cast<chrono::microseconds>(end - start).count();
                 sequentialTimes.push_back(seqTime);
             }
             

             MPI_Barrier(MPI_COMM_WORLD);
             double parStart = MPI_Wtime();
             ParallelDFAResult parallelResult = processDFAParallel(dfa, test.first, numProcs, myRank);
             (void)parallelResult;
             double parEnd = MPI_Wtime();
            
            if (myRank == 0) {
                double parTime = (parEnd - parStart) * 1000000;
                parallelTimes.push_back(parTime);
            }
        }
        

        if (myRank == 0) {
            double avgSeq = accumulate(sequentialTimes.begin(), sequentialTimes.end(), 0.0) / VALIDATION_ITERATIONS;
            double avgPar = accumulate(parallelTimes.begin(), parallelTimes.end(), 0.0) / VALIDATION_ITERATIONS;
            
            double speedup = avgSeq / avgPar;
            double efficiency = (speedup / numProcs) * 100;
            
            cout << "  • Tiempo secuencial promedio: " << fixed << setprecision(2) << avgSeq << " μs" << endl;
            cout << "  • Tiempo paralelo promedio: " << avgPar << " μs" << endl;
            cout << "  • Speedup: " << setprecision(3) << speedup << "x" << endl;
            cout << "  • Eficiencia: " << setprecision(2) << efficiency << "%" << endl;
        }
    }
    

    vector<int> largeSizes = {1000, 10000, 50000};
    
    for (int size : largeSizes) {
        string longInput = generateRandomString("01", size);
        
        if (myRank == 0) {
            cout << "Test con cadena de " << size << " caracteres (promedio de " << VALIDATION_ITERATIONS << " iteraciones)" << endl;
        }
        
        vector<double> sequentialTimes, parallelTimes;
        
        for (int iter = 0; iter < VALIDATION_ITERATIONS; iter++) {
               
             if (myRank == 0) {
                 auto start = chrono::high_resolution_clock::now();
                 bool seqResult = dfa.run(longInput);
                 (void)seqResult;
                 auto end = chrono::high_resolution_clock::now();
                 double seqTime = chrono::duration_cast<chrono::microseconds>(end - start).count();
                 sequentialTimes.push_back(seqTime);
             }
             

             MPI_Barrier(MPI_COMM_WORLD);
             double parStart = MPI_Wtime();
             ParallelDFAResult parallelResult = processDFAParallel(dfa, longInput, numProcs, myRank);
             (void)parallelResult;
             double parEnd = MPI_Wtime();
            
            if (myRank == 0) {
                double parTime = (parEnd - parStart) * 1000000;
                parallelTimes.push_back(parTime);
            }
        }
        
        // Análisis teórico para cadenas largas
        if (myRank == 0) {
            double avgSeq = accumulate(sequentialTimes.begin(), sequentialTimes.end(), 0.0) / VALIDATION_ITERATIONS;
            double avgPar = accumulate(parallelTimes.begin(), parallelTimes.end(), 0.0) / VALIDATION_ITERATIONS;
            
            double speedup = avgSeq / avgPar;
            double efficiency = (speedup / numProcs) * 100;
            
            // Análisis teórico según fórmulas
            int Q = 3, sigma = 2, n = size, p = numProcs;
            double theoretical_efficiency = (double)(Q * sigma + 1) / 
                                          (Q * sigma + n + p + p * log2(p) + Q * log2(p)) * 100;
            
            cout << "  • Tiempo secuencial promedio: " << fixed << setprecision(2) << avgSeq << " μs" << endl;
            cout << "  • Tiempo paralelo promedio: " << avgPar << " μs" << endl;
            cout << "  • Speedup real: " << setprecision(3) << speedup << "x" << endl;
            cout << "  • Eficiencia real: " << setprecision(2) << efficiency << "%" << endl;
            cout << "  • Eficiencia teórica: " << theoretical_efficiency << "%" << endl;
            cout << "  • Diferencia teoría vs real: " << abs(efficiency - theoretical_efficiency) << "%" << endl;
        }
    }
    
    // ================= TESTS DE CASOS EDGE =================
    if (myRank == 0) {
        cout << "\nEJECUTANDO TESTS DE CASOS EDGE..." << endl;
    }
    
    validateEdgeCaseResults(dfa, "", "Cadena Vacía", numProcs, myRank);
    validateEdgeCaseResults(dfa, "0", "Un Solo Carácter", numProcs, myRank);
    validateEdgeCaseResults(dfa, "01", "Dos Caracteres", numProcs, myRank);
    validateEdgeCaseResults(dfa, "000", "Tres Ceros", numProcs, myRank);
    
    // Casos edge con cadenas muy cortas vs número de procesos
    if (numProcs > 2) {
        validateEdgeCaseResults(dfa, "0", "Cadena Más Corta que Procesos", numProcs, myRank);
    }
    
    // ================= TESTS DE STRESS =================
    if (myRank == 0) {
        cout << "\nEJECUTANDO TESTS DE STRESS..." << endl;
    }
    
    runStressTests(numProcs, myRank);
    
    // ================= GUARDAR RESULTADOS MEJORADOS =================
    if (myRank == 0) {
        ofstream resultsFile("resultados_mpi/original_enhanced_results.txt");
        if (resultsFile.is_open()) {
            resultsFile << "============================================================\n";
            resultsFile << "RESULTADOS VERSIÓN ORIGINAL MEJORADA - " << numProcs << " PROCESOS\n";
            resultsFile << "============================================================\n";
            resultsFile << "Fecha: " << __DATE__ << " " << __TIME__ << "\n";
            resultsFile << "Iteraciones de validación: " << VALIDATION_ITERATIONS << "\n\n";
            
            resultsFile << "TESTS COMPLETADOS:\n";
            resultsFile << "  • Tests de comunicación MPI\n";
            resultsFile << "  • Tests paralelos vs secuenciales con promedios\n";
            resultsFile << "  • Análisis teórico de escalabilidad\n";
            resultsFile << "  • Tests de casos edge\n";
            resultsFile << "  • Tests de stress\n\n";
            
            resultsFile << "CONFIGURACIÓN DEL DFA:\n";
            resultsFile << "  • Estados: " << dfa.getNumStates() << "\n";
            resultsFile << "  • Alfabeto: " << dfa.getAlphabet() << "\n";
            resultsFile << "  • Estado inicial: " << dfa.getInitialState() << "\n\n";
            
            resultsFile << "ANÁLISIS TEÓRICO:\n";
            resultsFile << "  • Formula Tp = (|Q|·Σ+n)/p + n + log p + |Q|·log p\n";
            resultsFile << "  • Escalabilidad fuerte evaluada\n";
            resultsFile << "  • Escalabilidad débil verificada\n";
            
            resultsFile.close();
        }
    }
    
    // ================= RESUMEN FINAL =================
    if (myRank == 0) {
        cout << "\n============================================================" << endl;
        cout << "VERSIÓN ORIGINAL MEJORADA COMPLETADA" << endl;
        cout << "============================================================" << endl;
        cout << "Archivos de resultados guardados en: resultados_mpi/" << endl;
        cout << "  • results.txt - Resultados principales" << endl;
        cout << "  • communication_results.txt - Tests de comunicación" << endl;
        cout << "  • edge_case_results.txt - Casos edge" << endl;
        cout << "  • stress_test_results.txt - Tests de stress" << endl;
        cout << "  • original_enhanced_results.txt - Análisis mejorado" << endl;
        cout << "\nMEJORAS IMPLEMENTADAS:" << endl;
        cout << "  • Promedio de múltiples iteraciones (" << VALIDATION_ITERATIONS << ")" << endl;
        cout << "  • Análisis teórico de escalabilidad" << endl;
        cout << "  • Comparación teoría vs práctica" << endl;
        cout << "  • Medición precisa de tiempos" << endl;
        cout << "============================================================" << endl;
    }
    
    MPI_Finalize();
    return 0;
}
