#pragma once
#include <iostream>
#include <mpi.h>
#include <vector>
#include <string>
#include <chrono>
#include <iomanip>
#include <cmath>
#include "parallel_dfa.h"
using namespace std;

// ✅ FUNCIÓN DE TEST DE DIVISIÓN
void testBlockDivision(const string& testName, const string& input, int numProcesses);

// ✅ FUNCIÓN DE TESTS COMPLETOS DE DIVISIÓN
void runBlockDivisionTests();

// ✅ FUNCIÓN MEJORADA: MEDICIÓN CON MÚLTIPLES ITERACIONES Y PROMEDIOS
void measureCommunicationTimeImproved(const string& testName, const string& input, int numProcs, int myRank, int iterations);

// ✅ FUNCIÓN DE TESTS COMPLETOS DE COMUNICACIÓN MPI
void runMPICommunicationTests(int numProcs, int myRank);

// ✅ FUNCIONES DE TESTS DE COMUNICACIÓN
void testPointToPointCommunication(int numProcs, int myRank);
void testCollectiveCommunication(int numProcs, int myRank);
void testDerivedDataTypes(int numProcs, int myRank);
void testOptimizedBuffers(int numProcs, int myRank);
void logCommunicationResults(const string& testName, const vector<double>& times, 
                           const string& operation, int numProcs, int myRank); 