#include "dfa.h"
#include <random>
#include <fstream>

// ✅ EJECUTA EL OMP INICIAL 
void testParallelVsSequential(const DFA& dfa, const std::string& input, 
                              const std::string& testName, int numThreads) {
    // Procesamiento secuencial
    int sequentialResult = -1;
    double sequentialTime = 0;

    double startTime = omp_get_wtime();
    sequentialResult = dfa.getInitialState();
    for (char c : input) {
        sequentialResult = dfa.nextState(sequentialResult, c);
    }
    double endTime = omp_get_wtime();
    sequentialTime = (endTime - startTime) * 1000000;
    
    // Procesamiento paralelo
    DFAResult parallelResult = processOMPParallel(dfa, input, numThreads);
    DFAOptimizedResult parallelOptimizedResult = runOptimizedOMPParallel(dfa, input, numThreads);

    logParallelResults(testName, sequentialTime, parallelResult, parallelOptimizedResult, numThreads);
}


// ✅ FUNCIÓN: Log de resultados paralelos con validación
void logParallelResults(const string& testName, double sequentialTime, 
                       const DFAResult& parallelResult, const DFAOptimizedResult& parallelOptResult, 
                       int threads) {    
    double speedup = sequentialTime / parallelResult.computationTime;
    double efficiency = speedup / threads * 100;
    
    string message = "=== TEST CONCURRENCIA: " + testName + " ===\n";
    message += "Procesos: " + to_string(threads) + "\n";
    message += "📊 Secuencial: Estado " + to_string(parallelResult.finalState) + 
               " | Tiempo: " + to_string(sequentialTime) + " μs\n";
    message += "📊 Paralelo OMP: Estado " + to_string(parallelResult.finalState) +
               " | Tiempo total: " + to_string(parallelResult.computationTime) + " μs\n";
    message += "📊 Paralelo OMP Optimizado: Estado " + to_string(parallelOptResult.finalState) +
               " | Tiempo total: " + to_string(parallelOptResult.computationTime) + " μs\n";
    message += "  ⚡ Tiempo computación: " + to_string(parallelResult.computationTime) + " μs\n";
    message += "🚀 Speedup: " + to_string(speedup) + "x | Eficiencia: " + to_string(efficiency) + "%\n\n";
    
    string finalFile = "result_" +  testName + ".txt";
    logToFile(message, finalFile);
}

// ✅ FUNCIÓN: Escribir mensaje a archivo con debugging
void logToFile(const string& message, const string& filename) {
    string resultsDir = "resultados_omp";
    string fullPath = resultsDir + "/" + filename;
    ofstream file(fullPath, ios::app);
    if (file.is_open()) {
        file << message << "\n";
        file.close();
        // Debug: mostrar que se escribió correctamente
        cout << "📝 Log escrito en: " << fullPath << endl;
    } else {
        cout << "❌ ERROR: No se pudo abrir archivo: " << fullPath << endl;
    }
}