#pragma once
#include <iostream>
#include <mpi.h>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>
#include <cmath>
#include "optimized_parallel_dfa.h"
using namespace std;

// ✅ ESTRUCTURA PARA BLOQUE DE ENTRADA
struct InputBlock {
    std::string data;
    int startPos;
    int endPos;
    int processId;
    
    InputBlock(const std::string& d, int start, int end, int pid) 
        : data(d), startPos(start), endPos(end), processId(pid) {}
};

// ✅ ESTRUCTURA PARA RESULTADO PARALELO
struct ParallelDFAResult {
    int finalState;
    bool accepted;
    double computationTime;
    double communicationTime;
    int blockSize;
    
    ParallelDFAResult() : finalState(-1), accepted(false), 
                         computationTime(0.0), communicationTime(0.0), blockSize(0) {}
};

// ✅ ESTRUCTURA OPTIMIZADA PARA COMUNICACIÓN MPI
struct MPICommunicationData {
    MPI_Datatype blockInfoType;    // Tipo derivado para información de bloques
    MPI_Datatype dfaDataType;      // Tipo derivado para datos del DFA
    MPI_Datatype resultType;       // Tipo derivado para resultados
    bool typesInitialized;
    
    MPICommunicationData() : typesInitialized(false) {}
    
    // ✅ Inicializar tipos de datos derivados MPI
    void initializeMPITypes();
    
    // ✅ Liberar tipos de datos derivados MPI
    void cleanupMPITypes();
};

// ✅ ESTRUCTURA PARA INFORMACIÓN DE BLOQUE (optimizada para MPI)
struct BlockInfo {
    int size;
    int startPos;
    int endPos;
    int processId;
};

// ✅ ESTRUCTURA PARA DATOS DEL DFA (optimizada para MPI)
struct DFAData {
    int numStates;
    int initialState;
    int numFinalStates;
    int alphabetSize;
    char alphabet[128];  // Alfabeto fijo para evitar comunicación dinámica
    int finalStates[100]; // Estados finales fijos
    int transitions[100][128]; // Matriz de transiciones fija
};

// ✅ ESTRUCTURA PARA RESULTADO (optimizada para MPI)
struct ResultData {
    int finalState;
    int accepted;
    double computationTime;
    double communicationTime;
    int blockSize;
};

// ✅ FUNCIONES DE DIVISIÓN Y VALIDACIÓN
std::vector<InputBlock> divideInputIntoBlocks(const std::string& input, int numProcesses);
bool validateBlockDivision(const std::string& originalInput, const std::vector<InputBlock>& blocks);

// ✅ FUNCIONES DE PROCESAMIENTO PARALELO
int computeInitialStateForBlock(const OptimizedDFA& dfa, const std::string& input, int blockStartPos);
ParallelDFAResult processBlockParallel(const OptimizedDFA& dfa, const std::string& blockData, 
                                       const std::string& fullInput, int blockStartPos);

// ✅ FUNCIONES DE COMUNICACIÓN OPTIMIZADAS
void broadcastDFAOptimized(OptimizedDFA& dfa, int myRank, MPICommunicationData& commData);
ParallelDFAResult processDFAParallelOptimized(const OptimizedDFA& dfa, const std::string& input, 
                                              int numProcs, int myRank, MPICommunicationData& commData);

// ✅ FUNCIONES DE COMUNICACIÓN ORIGINALES (para compatibilidad)
void broadcastDFA(OptimizedDFA& dfa, int myRank);
ParallelDFAResult processDFAParallel(const OptimizedDFA& dfa, const std::string& input, 
                                     int numProcs, int myRank);

// ✅ FUNCIONES DE TEST Y VALIDACIÓN
void testParallelVsSequential(const OptimizedDFA& dfa, const std::string& input, 
                              const std::string& testName, int numProcs, int myRank);
void runParallelProcessingTests(int numProcs, int myRank);

// ✅ FUNCIONES DE OPTIMIZACIÓN Y ESCALABILIDAD
void runScalabilityTests(int numProcs, int myRank);
void runStressTests(int numProcs, int myRank);
void runEdgeCaseTests(int numProcs, int myRank);

// ✅ FUNCIONES DE MANEJO DE CASOS EDGE EN PARALELO
void handleEmptyStringCase(const OptimizedDFA& dfa, int numProcs, int myRank);
void handleSingleCharacterCase(const OptimizedDFA& dfa, int numProcs, int myRank);
void handleVeryShortStringCase(const OptimizedDFA& dfa, int numProcs, int myRank);
void handleUnevenBlockDistribution(const OptimizedDFA& dfa, int numProcs, int myRank);
void validateEdgeCaseResults(const OptimizedDFA& dfa, const string& input, 
                            const string& testName, int numProcs, int myRank);

// ✅ FUNCIONES DE TESTS DE STRESS CON CADENAS MUY LARGAS
void runStressTests(int numProcs, int myRank);
void testStressWithDifferentSizes(int numProcs, int myRank);
void testStressWithDifferentPatterns(int numProcs, int myRank);
void testStressWithMemoryPressure(int numProcs, int myRank);
void testStressWithConcurrentProcessing(int numProcs, int myRank);
void generateStressTestData(vector<string>& inputs, vector<string>& testNames);
void measureStressPerformance(const OptimizedDFA& dfa, const string& input, 
                             const string& testName, int numProcs, int myRank, int iterations);
void analyzeStressResults(const vector<double>& sequentialTimes, 
                         const vector<double>& parallelTimes, 
                         const vector<double>& speedups, 
                         const string& testName, int numProcs, int myRank);

// ✅ FUNCIONES DE LOGGING A ARCHIVOS
void initializeLogging(int numProcs, int myRank);
void logToFile(const string& message, const string& filename = "results.txt");
void logCommunicationResults(const string& testName, const vector<double>& times, 
                           const string& operation, int numProcs, int myRank);
void logParallelResults(const string& testName, double sequentialTime, 
                       const ParallelDFAResult& parallelResult, int numProcs, int myRank);
void logEdgeCaseResults(const string& testName, double sequentialTime, 
                       const ParallelDFAResult& parallelResult, int numProcs, int myRank);
void logStressResults(const string& testName, const vector<double>& sequentialTimes, 
                     const vector<double>& parallelTimes, const vector<double>& speedups, 
                     int numProcs, int myRank); 