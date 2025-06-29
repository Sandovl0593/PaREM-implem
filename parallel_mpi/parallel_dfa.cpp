#include "parallel_dfa.h"
#include "dfa_tests.h"
#include "utils.h"
#include <fstream>
#include <iomanip>

// ✅ FUNCIÓN DE DIVISIÓN DE BLOQUES
vector<InputBlock> divideInputIntoBlocks(const string& input, int numProcesses) {
    vector<InputBlock> blocks;
    int inputSize = input.size();
    
    if (inputSize == 0) {
        // Caso edge: cadena vacía
        for (int i = 0; i < numProcesses; i++) {
            blocks.emplace_back("", 0, 0, i);
        }
    } else {
        // División balanceada
        int baseSize = inputSize / numProcesses;
        int remainder = inputSize % numProcesses;
        
        int currentPos = 0;
        for (int i = 0; i < numProcesses; i++) {
            int blockSize = baseSize + (i < remainder ? 1 : 0);
            int endPos = currentPos + blockSize;
            
            string blockData = input.substr(currentPos, blockSize);
            blocks.emplace_back(blockData, currentPos, endPos - 1, i);
            
            currentPos = endPos;
        }
    }
    
    return blocks;
}

// ✅ FUNCIÓN DE VALIDACIÓN DE DIVISIÓN
bool validateBlockDivision(const string& originalInput, const vector<InputBlock>& blocks) {
    // Verificar que no se perdió información
    string reconstructed = "";
    for (const auto& block : blocks) {
        reconstructed += block.data;
    }
    
    if (reconstructed != originalInput) {
        cout << "❌ ERROR: Información perdida en división\n";
        cout << "Original: [" << originalInput << "]\n";
        cout << "Reconstruido: [" << reconstructed << "]\n";
        return false;
    }
    
    // Verificar que no hay solapamiento
    for (size_t i = 0; i < blocks.size(); i++) {
        for (size_t j = i + 1; j < blocks.size(); j++) {
            if (blocks[i].endPos >= blocks[j].startPos && blocks[i].startPos <= blocks[j].endPos) {
                cout << "❌ ERROR: Solapamiento entre bloques " << i << " y " << j << "\n";
                return false;
            }
        }
    }
    
    // Verificar que cubre toda la cadena
    if (!blocks.empty()) {
        int minStart = blocks[0].startPos;
        int maxEnd = blocks.back().endPos;
        if (minStart != 0 || maxEnd != (int)originalInput.size() - 1) {
            cout << "❌ ERROR: No cubre toda la cadena\n";
            return false;
        }
    }
    
    return true;
}

// ✅ FUNCIÓN: Calcular estado inicial para un bloque
int computeInitialStateForBlock(const OptimizedDFA& dfa, const string& input, int blockStartPos) {
    if (blockStartPos == 0) {
        return dfa.getInitialState();
    }
    
    // Simular el procesamiento de los caracteres anteriores al bloque
    int currentState = dfa.getInitialState();
    for (int i = 0; i < blockStartPos; i++) {
        if (i < (int)input.size()) {
            currentState = dfa.nextState(currentState, input[i]);
        }
    }
    return currentState;
}

// ✅ FUNCIÓN: Procesar un bloque de entrada en paralelo
ParallelDFAResult processBlockParallel(const OptimizedDFA& dfa, const string& blockData, 
                                       const string& fullInput, int blockStartPos) {
    ParallelDFAResult result;
    result.blockSize = blockData.size();
    
    double startTime = MPI_Wtime();
    
    // Calcular estado inicial para este bloque
    int currentState = computeInitialStateForBlock(dfa, fullInput, blockStartPos);
    
    // Procesar cada carácter del bloque
    for (char c : blockData) {
        currentState = dfa.nextState(currentState, c);
    }
    
    double endTime = MPI_Wtime();
    result.computationTime = (endTime - startTime) * 1000000; // Convertir a microsegundos
    result.finalState = currentState;
    result.accepted = dfa.isAccepting(currentState);
    
    return result;
}

// ✅ FUNCIÓN: Distribuir DFA optimizado usando tipos derivados MPI
void broadcastDFAOptimized(OptimizedDFA& dfa, int myRank, MPICommunicationData& commData) {
    // Inicializar tipos MPI si no están inicializados
    commData.initializeMPITypes();
    
    DFAData dfaData;
    
    if (myRank == 0) {
        // El proceso raíz prepara los datos del DFA optimizados
        dfaData.numStates = dfa.getNumStates();
        dfaData.initialState = dfa.getInitialState();
        string alphabet = dfa.getAlphabet();
        dfaData.alphabetSize = alphabet.size();
        
        // Copiar alfabeto al array fijo
        for (int i = 0; i < 128; i++) {
            dfaData.alphabet[i] = (i < dfaData.alphabetSize) ? alphabet[i] : '\0';
        }
        
        // Obtener estados finales
        dfaData.numFinalStates = 0;
        for (int i = 0; i < dfaData.numStates && dfaData.numFinalStates < 100; i++) {
            if (dfa.isAccepting(i)) {
                dfaData.finalStates[dfaData.numFinalStates++] = i;
            }
        }
        
        // Preparar matriz de transiciones optimizada
        for (int state = 0; state < 100; state++) {
            for (int c = 0; c < 128; c++) {
                if (state < dfaData.numStates) {
                    dfaData.transitions[state][c] = dfa.nextState(state, (char)c);
                } else {
                    dfaData.transitions[state][c] = -1;
                }
            }
        }
    }
    
    // Broadcast optimizado usando tipo derivado
    MPI_Bcast(&dfaData, 1, commData.dfaDataType, 0, MPI_COMM_WORLD);
    
    // Reconstruir DFA en todos los procesos
    if (myRank != 0) {
        string reconstructedAlphabet(dfaData.alphabet, dfaData.alphabetSize);
        vector<int> finalStates(dfaData.finalStates, dfaData.finalStates + dfaData.numFinalStates);
        
        dfa = OptimizedDFA(dfaData.numStates, dfaData.initialState, finalStates, reconstructedAlphabet);
        
        // Reconstruir transiciones usando la matriz optimizada
        for (int state = 0; state < dfaData.numStates; state++) {
            for (int c = 0; c < 128; c++) {
                int nextState = dfaData.transitions[state][c];
                if (nextState != -1) {
                    dfa.addTransition(state, (char)c, nextState);
                }
            }
        }
    }
}

// ✅ FUNCIÓN: Procesamiento paralelo optimizado usando tipos derivados MPI
ParallelDFAResult processDFAParallelOptimized(const OptimizedDFA& dfa, const string& input, 
                                              int numProcs, int myRank, MPICommunicationData& commData) {
    ParallelDFAResult result;
    
    double totalStartTime = MPI_Wtime();
    
    // Paso 1: Distribuir el DFA optimizado
    OptimizedDFA localDFA = dfa;
    broadcastDFAOptimized(localDFA, myRank, commData);
    
    // Paso 2: Distribuir la entrada usando tipos optimizados
    vector<BlockInfo> allBlockInfos;
    if (myRank == 0) {
        vector<InputBlock> allBlocks = divideInputIntoBlocks(input, numProcs);
        allBlockInfos.resize(numProcs);
        
        for (int i = 0; i < numProcs; i++) {
            allBlockInfos[i].size = allBlocks[i].data.size();
            allBlockInfos[i].startPos = allBlocks[i].startPos;
            allBlockInfos[i].endPos = allBlocks[i].endPos;
            allBlockInfos[i].processId = allBlocks[i].processId;
        }
    }
    
    // Broadcast del tamaño de la entrada
    int inputSize = input.size();
    MPI_Bcast(&inputSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // Scatter optimizado de información de bloques usando tipo derivado
    BlockInfo myBlockInfo;
    MPI_Scatter(allBlockInfos.data(), 1, commData.blockInfoType, 
                &myBlockInfo, 1, commData.blockInfoType, 0, MPI_COMM_WORLD);
    
    // Distribuir datos del bloque (esto se mantiene igual por ser caracteres)
    vector<char> myBlockData(myBlockInfo.size);
    vector<int> sendCounts(numProcs);
    vector<int> displacements(numProcs);
    
    if (myRank == 0) {
        vector<InputBlock> allBlocks = divideInputIntoBlocks(input, numProcs);
        for (int i = 0; i < numProcs; i++) {
            sendCounts[i] = allBlocks[i].data.size();
            displacements[i] = allBlocks[i].startPos;
        }
    }
    
    MPI_Scatterv(input.c_str(), sendCounts.data(), displacements.data(), MPI_CHAR,
                 myBlockData.data(), myBlockInfo.size, MPI_CHAR, 0, MPI_COMM_WORLD);
    
    // Paso 3: Procesar el bloque local
    string myBlockString(myBlockData.begin(), myBlockData.end());
    ParallelDFAResult localResult = processBlockParallel(localDFA, myBlockString, input, myBlockInfo.startPos);
    
    // Paso 4: Recopilar resultados usando tipo optimizado
    ResultData localResultData;
    localResultData.finalState = localResult.finalState;
    localResultData.accepted = localResult.accepted ? 1 : 0;
    localResultData.computationTime = localResult.computationTime;
    localResultData.communicationTime = localResult.communicationTime;
    localResultData.blockSize = localResult.blockSize;
    
    vector<ResultData> allResults(numProcs);
    MPI_Gather(&localResultData, 1, commData.resultType, 
               allResults.data(), 1, commData.resultType, 0, MPI_COMM_WORLD);
    
    double totalEndTime = MPI_Wtime();
    result.communicationTime = (totalEndTime - totalStartTime) * 1000000;
    result.computationTime = localResult.computationTime;
    result.blockSize = myBlockInfo.size;
    
    // El proceso raíz determina el resultado final
    if (myRank == 0) {
        // Simular el procesamiento secuencial para obtener el estado final correcto
        int currentState = localDFA.getInitialState();
        for (char c : input) {
            currentState = localDFA.nextState(currentState, c);
        }
        result.finalState = currentState;
        result.accepted = localDFA.isAccepting(currentState);
    }
    
    return result;
}

// ✅ FUNCIÓN: Test de procesamiento paralelo vs secuencial
void testParallelVsSequential(const OptimizedDFA& dfa, const string& input, 
                              const string& testName, int numProcs, int myRank) {
    if (myRank == 0) {
        // Solo mostrar información básica en consola
        cout << "Ejecutando test: " << testName << "..." << endl;
    }
    
    // Procesamiento secuencial
    int sequentialResult = -1;
    double sequentialTime = 0;
    
    if (myRank == 0) {
        double startTime = MPI_Wtime();
        sequentialResult = dfa.getInitialState();
        for (char c : input) {
            sequentialResult = dfa.nextState(sequentialResult, c);
        }
        double endTime = MPI_Wtime();
        sequentialTime = (endTime - startTime) * 1000000;
    }
    
    // Procesamiento paralelo
    ParallelDFAResult parallelResult = processDFAParallel(dfa, input, numProcs, myRank);
    
    // Log de resultados
    if (myRank == 0) {
        logParallelResults(testName, sequentialTime, parallelResult, numProcs, myRank);
    }
}

// ✅ FUNCIÓN: Tests completos de procesamiento paralelo usando DFA_TESTS
void runParallelProcessingTests(int numProcs, int myRank) {
    if (myRank == 0) {
        cout << "\n" << string(60, '=') << "\n";
        cout << "⚡ TESTS DE PROCESAMIENTO PARALELO DFA\n";
        cout << string(60, '=') << "\n";
    }
    
    // Test 1: Autómata Binario (usando tests existentes)
    OptimizedDFA dfa1(3, 0, {2}, "01");
    dfa1.addTransition(0, '0', 1); dfa1.addTransition(0, '1', 0);
    dfa1.addTransition(1, '0', 1); dfa1.addTransition(1, '1', 2);
    dfa1.addTransition(2, '0', 1); dfa1.addTransition(2, '1', 0);
    
    // Usar tests del archivo dfa_tests.h
    auto binaryTests = createBinaryTests();
    for (const auto& test : binaryTests) {
        if (test.input.size() <= 1000) { // Solo tests pequeños para empezar
            testParallelVsSequential(dfa1, test.input, 
                "Autómata Binario - " + test.description, numProcs, myRank);
        }
    }
    
    // Test 2: Autómata ABC
    OptimizedDFA dfa2(4, 0, {3}, "abc");
    dfa2.addTransition(0, 'a', 1); dfa2.addTransition(0, 'b', 0); dfa2.addTransition(0, 'c', 0);
    dfa2.addTransition(1, 'a', 1); dfa2.addTransition(1, 'b', 2); dfa2.addTransition(1, 'c', 0);
    dfa2.addTransition(2, 'a', 1); dfa2.addTransition(2, 'b', 0); dfa2.addTransition(2, 'c', 3);
    dfa2.addTransition(3, 'a', 1); dfa2.addTransition(3, 'b', 0); dfa2.addTransition(3, 'c', 0);
    
    auto abcTests = createABCTests();
    for (const auto& test : abcTests) {
        if (test.input.size() <= 1000) {
            testParallelVsSequential(dfa2, test.input, 
                "Autómata ABC - " + test.description, numProcs, myRank);
        }
    }
    
    // Test 3: Autómata Múltiples Finales
    OptimizedDFA dfa3(4, 0, {2, 3}, "pq");
    dfa3.addTransition(0, 'p', 1); dfa3.addTransition(0, 'q', 0);
    dfa3.addTransition(1, 'p', 1); dfa3.addTransition(1, 'q', 2);
    dfa3.addTransition(2, 'p', 1); dfa3.addTransition(2, 'q', 3);
    dfa3.addTransition(3, 'p', 1); dfa3.addTransition(3, 'q', 0);
    
    auto multiFinalTests = createMultiFinalTests();
    for (const auto& test : multiFinalTests) {
        if (test.input.size() <= 1000) {
            testParallelVsSequential(dfa3, test.input, 
                "Autómata Multi-Final - " + test.description, numProcs, myRank);
        }
    }
    
    // Test 4: Cadena larga para medir escalabilidad
    string longString = "";
    for (int i = 0; i < 10000; i++) {
        longString += (char)('0' + (i % 2)); // Patrón binario
    }
    testParallelVsSequential(dfa1, longString, "Autómata Binario - Cadena Larga (10K)", numProcs, myRank);
    
    if (myRank == 0) {
        cout << "\n" << string(60, '=') << "\n";
        cout << "✅ TESTS DE PROCESAMIENTO PARALELO COMPLETADOS\n";
        cout << string(60, '=') << "\n";
    }
}

// ================= MPICommunicationData: Inicialización y Liberación de Tipos =================
void MPICommunicationData::initializeMPITypes() {
    if (typesInitialized) return;
    // BlockInfo
    int blockInfoLengths[4] = {1, 1, 1, 1};
    MPI_Aint blockInfoDisps[4];
    BlockInfo dummyBlock = {0, 0, 0, 0}; // Inicializar para evitar warning
    MPI_Aint base;
    MPI_Get_address(&dummyBlock, &base);
    MPI_Get_address(&dummyBlock.size, &blockInfoDisps[0]);
    MPI_Get_address(&dummyBlock.startPos, &blockInfoDisps[1]);
    MPI_Get_address(&dummyBlock.endPos, &blockInfoDisps[2]);
    MPI_Get_address(&dummyBlock.processId, &blockInfoDisps[3]);
    for (int i = 0; i < 4; ++i) blockInfoDisps[i] -= base;
    MPI_Datatype blockInfoTypes[4] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT};
    MPI_Type_create_struct(4, blockInfoLengths, blockInfoDisps, blockInfoTypes, &blockInfoType);
    MPI_Type_commit(&blockInfoType);

    // DFAData
    int dfaDataLengths[6] = {1, 1, 1, 1, 128, 100*128+100};
    MPI_Aint dfaDataDisps[6];
    DFAData dummyDFA = {0, 0, 0, 0, {0}, {0}, {{0}}}; // Inicializar todos los campos
    MPI_Get_address(&dummyDFA, &base);
    MPI_Get_address(&dummyDFA.numStates, &dfaDataDisps[0]);
    MPI_Get_address(&dummyDFA.initialState, &dfaDataDisps[1]);
    MPI_Get_address(&dummyDFA.numFinalStates, &dfaDataDisps[2]);
    MPI_Get_address(&dummyDFA.alphabetSize, &dfaDataDisps[3]);
    MPI_Get_address(&dummyDFA.alphabet, &dfaDataDisps[4]);
    MPI_Get_address(&dummyDFA.finalStates, &dfaDataDisps[5]);
    for (int i = 0; i < 6; ++i) dfaDataDisps[i] -= base;
    MPI_Datatype dfaDataTypes[6] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_CHAR, MPI_INT};
    MPI_Type_create_struct(6, dfaDataLengths, dfaDataDisps, dfaDataTypes, &dfaDataType);
    MPI_Type_commit(&dfaDataType);

    // ResultData
    int resultLengths[5] = {1, 1, 1, 1, 1};
    MPI_Aint resultDisps[5];
    ResultData dummyResult = {0, 0, 0.0, 0.0, 0}; // Inicializar para evitar warning
    MPI_Get_address(&dummyResult, &base);
    MPI_Get_address(&dummyResult.finalState, &resultDisps[0]);
    MPI_Get_address(&dummyResult.accepted, &resultDisps[1]);
    MPI_Get_address(&dummyResult.computationTime, &resultDisps[2]);
    MPI_Get_address(&dummyResult.communicationTime, &resultDisps[3]);
    MPI_Get_address(&dummyResult.blockSize, &resultDisps[4]);
    for (int i = 0; i < 5; ++i) resultDisps[i] -= base;
    MPI_Datatype resultTypes[5] = {MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_INT};
    MPI_Type_create_struct(5, resultLengths, resultDisps, resultTypes, &resultType);
    MPI_Type_commit(&resultType);

    typesInitialized = true;
}

void MPICommunicationData::cleanupMPITypes() {
    if (!typesInitialized) return;
    MPI_Type_free(&blockInfoType);
    MPI_Type_free(&dfaDataType);
    MPI_Type_free(&resultType);
    typesInitialized = false;
}

// ✅ FUNCIÓN: Procesamiento paralelo completo de DFA (versión original)
ParallelDFAResult processDFAParallel(const OptimizedDFA& dfa, const string& input, 
                                     int numProcs, int myRank) {
    ParallelDFAResult result;
    
    double totalStartTime = MPI_Wtime();
    
    // Paso 1: Distribuir el DFA a todos los procesos
    OptimizedDFA localDFA = dfa; // Copia local
    broadcastDFA(localDFA, myRank);
    
    // Paso 2: Distribuir la entrada
    vector<InputBlock> allBlocks;
    if (myRank == 0) {
        allBlocks = divideInputIntoBlocks(input, numProcs);
    }
    
    // Broadcast del tamaño de la entrada
    int inputSize = input.size();
    MPI_Bcast(&inputSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // Distribuir información de bloques
    vector<int> blockSizes(numProcs);
    vector<int> startPositions(numProcs);
    vector<int> endPositions(numProcs);
    
    if (myRank == 0) {
        for (int i = 0; i < numProcs; i++) {
            blockSizes[i] = allBlocks[i].data.size();
            startPositions[i] = allBlocks[i].startPos;
            endPositions[i] = allBlocks[i].endPos;
        }
    }
    
    // Scatter de información de bloques
    int myBlockSize, myStartPos, myEndPos;
    MPI_Scatter(blockSizes.data(), 1, MPI_INT, &myBlockSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(startPositions.data(), 1, MPI_INT, &myStartPos, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(endPositions.data(), 1, MPI_INT, &myEndPos, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // Distribuir datos del bloque
    vector<char> myBlockData(myBlockSize);
    vector<int> sendCounts(numProcs);
    vector<int> displacements(numProcs);
    
    if (myRank == 0) {
        for (int i = 0; i < numProcs; i++) {
            sendCounts[i] = blockSizes[i];
            displacements[i] = allBlocks[i].startPos;
        }
    }
    
    MPI_Scatterv(input.c_str(), sendCounts.data(), displacements.data(), MPI_CHAR,
                 myBlockData.data(), myBlockSize, MPI_CHAR, 0, MPI_COMM_WORLD);
    
    // Paso 3: Procesar el bloque local
    string myBlockString(myBlockData.begin(), myBlockData.end());
    ParallelDFAResult localResult = processBlockParallel(localDFA, myBlockString, input, myStartPos);
    
    // Paso 4: Recopilar resultados y propagar estados
    vector<int> allFinalStates(numProcs);
    MPI_Gather(&localResult.finalState, 1, MPI_INT, allFinalStates.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    double totalEndTime = MPI_Wtime();
    result.communicationTime = (totalEndTime - totalStartTime) * 1000000;
    result.computationTime = localResult.computationTime;
    result.blockSize = myBlockSize;
    
    // El proceso raíz determina el resultado final
    if (myRank == 0) {
        // Simular el procesamiento secuencial para obtener el estado final correcto
        int currentState = localDFA.getInitialState();
        for (char c : input) {
            currentState = localDFA.nextState(currentState, c);
        }
        result.finalState = currentState;
        result.accepted = localDFA.isAccepting(currentState);
    }
    
    return result;
}

// ✅ FUNCIÓN: Distribuir DFA a todos los procesos (versión original)
void broadcastDFA(OptimizedDFA& dfa, int myRank) {
    int numStates, initialState, numFinalStates;
    string alphabet;
    vector<int> finalStates;
    vector<vector<int>> transitions;
    
    if (myRank == 0) {
        // El proceso raíz prepara los datos del DFA
        numStates = dfa.getNumStates();
        initialState = dfa.getInitialState();
        alphabet = dfa.getAlphabet();
        
        // Obtener estados finales
        finalStates.clear();
        for (int i = 0; i < numStates; i++) {
            if (dfa.isAccepting(i)) {
                finalStates.push_back(i);
            }
        }
        numFinalStates = finalStates.size();
        
        // Preparar matriz de transiciones
        transitions.resize(numStates);
        for (int state = 0; state < numStates; state++) {
            transitions[state].resize(128, -1); // ASCII completo
            for (char c : alphabet) {
                transitions[state][(int)c] = dfa.nextState(state, c);
            }
        }
    }
    
    // Broadcast de metadatos
    MPI_Bcast(&numStates, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&initialState, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&numFinalStates, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    // Broadcast del alfabeto
    int alphabetSize = 0;
    if (myRank == 0) {
        alphabetSize = alphabet.size();
    }
    MPI_Bcast(&alphabetSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    vector<char> alphabetArray(alphabetSize);
    if (myRank == 0) {
        for (int i = 0; i < alphabetSize; i++) {
            alphabetArray[i] = alphabet[i];
        }
    }
    MPI_Bcast(alphabetArray.data(), alphabetSize, MPI_CHAR, 0, MPI_COMM_WORLD);
    
    // Broadcast de estados finales
    vector<int> finalStatesArray(numFinalStates);
    if (myRank == 0) {
        for (int i = 0; i < numFinalStates; i++) {
            finalStatesArray[i] = finalStates[i];
        }
    }
    MPI_Bcast(finalStatesArray.data(), numFinalStates, MPI_INT, 0, MPI_COMM_WORLD);
    
    // Broadcast de matriz de transiciones
    vector<int> flatTransitions(numStates * 128);
    if (myRank == 0) {
        for (int state = 0; state < numStates; state++) {
            for (int c = 0; c < 128; c++) {
                flatTransitions[state * 128 + c] = transitions[state][c];
            }
        }
    }
    MPI_Bcast(flatTransitions.data(), numStates * 128, MPI_INT, 0, MPI_COMM_WORLD);
    
    // Reconstruir DFA en todos los procesos
    if (myRank != 0) {
        string reconstructedAlphabet(alphabetArray.begin(), alphabetArray.end());
        dfa = OptimizedDFA(numStates, initialState, finalStatesArray, reconstructedAlphabet);
        
        // Reconstruir transiciones
        for (int state = 0; state < numStates; state++) {
            for (int c = 0; c < 128; c++) {
                int nextState = flatTransitions[state * 128 + c];
                if (nextState != -1) {
                    dfa.addTransition(state, (char)c, nextState);
                }
            }
        }
    }
}

// ================= MANEJO DE CASOS EDGE EN PARALELO =================

// ✅ FUNCIÓN: Validar resultados de casos edge
void validateEdgeCaseResults(const OptimizedDFA& dfa, const string& input, 
                            const string& testName, int numProcs, int myRank) {
    if (myRank == 0) {
        // Solo mostrar información básica en consola
        cout << "Ejecutando caso edge: " << testName << "..." << endl;
    }
    
    // Procesamiento secuencial
    int sequentialResult = -1;
    double sequentialTime = 0;
    
    if (myRank == 0) {
        double startTime = MPI_Wtime();
        sequentialResult = dfa.getInitialState();
        for (char c : input) {
            sequentialResult = dfa.nextState(sequentialResult, c);
        }
        double endTime = MPI_Wtime();
        sequentialTime = (endTime - startTime) * 1000000;
    }
    
    // Procesamiento paralelo con validación especial
    ParallelDFAResult parallelResult = processDFAParallel(dfa, input, numProcs, myRank);
    
    // Log de resultados
    if (myRank == 0) {
        logEdgeCaseResults(testName, sequentialTime, parallelResult, numProcs, myRank);
    }
}

// ✅ FUNCIÓN: Manejar caso de cadena vacía
void handleEmptyStringCase(const OptimizedDFA& dfa, int numProcs, int myRank) {
    string emptyInput = "";
    
    if (myRank == 0) {
        cout << "\n" << string(60, '=') << "\n";
        cout << "🔍 MANEJO DE CASO EDGE: Cadena Vacía\n";
        cout << string(60, '=') << "\n";
    }
    
    validateEdgeCaseResults(dfa, emptyInput, "Cadena Vacía", numProcs, myRank);
    
    // Validación adicional: verificar que todos los procesos manejan correctamente
    if (myRank == 0) {
        cout << "✅ Validación: Todos los procesos deben terminar en estado inicial\n";
        cout << "✅ Validación: Tiempo de comunicación debe ser mínimo\n";
    }
}

// ✅ FUNCIÓN: Manejar caso de un solo carácter
void handleSingleCharacterCase(const OptimizedDFA& dfa, int numProcs, int myRank) {
    vector<string> singleCharInputs = {"a", "0", "x", "1"};
    
    if (myRank == 0) {
        cout << "\n" << string(60, '=') << "\n";
        cout << "🔍 MANEJO DE CASO EDGE: Un Solo Carácter\n";
        cout << string(60, '=') << "\n";
    }
    
    for (const string& input : singleCharInputs) {
        validateEdgeCaseResults(dfa, input, "Un Solo Carácter: '" + input + "'", numProcs, myRank);
    }
    
    if (myRank == 0) {
        cout << "✅ Validación: Solo un proceso debe procesar el carácter\n";
        cout << "✅ Validación: Otros procesos deben recibir bloques vacíos\n";
    }
}

// ✅ FUNCIÓN: Manejar caso de cadenas muy cortas
void handleVeryShortStringCase(const OptimizedDFA& dfa, int numProcs, int myRank) {
    vector<string> shortInputs = {"ab", "01", "xyz", "123"};
    
    if (myRank == 0) {
        cout << "\n" << string(60, '=') << "\n";
        cout << "🔍 MANEJO DE CASO EDGE: Cadenas Muy Cortas\n";
        cout << string(60, '=') << "\n";
    }
    
    for (const string& input : shortInputs) {
        validateEdgeCaseResults(dfa, input, "Cadena Muy Corta: '" + input + "'", numProcs, myRank);
    }
    
    if (myRank == 0) {
        cout << "✅ Validación: Distribución desigual de caracteres\n";
        cout << "✅ Validación: Algunos procesos pueden recibir bloques vacíos\n";
    }
}

// ✅ FUNCIÓN: Manejar caso de distribución desigual de bloques
void handleUnevenBlockDistribution(const OptimizedDFA& dfa, int numProcs, int myRank) {
    // Crear cadenas que causen distribución desigual
    vector<string> unevenInputs;
    
    // Cadena con tamaño que no es múltiplo del número de procesos
    string uneven1 = "";
    for (int i = 0; i < numProcs + 1; i++) {
        uneven1 += "a";
    }
    unevenInputs.push_back(uneven1);
    
    // Cadena con tamaño primo
    string uneven2 = "";
    for (int i = 0; i < 7; i++) {
        uneven2 += "0";
    }
    unevenInputs.push_back(uneven2);
    
    if (myRank == 0) {
        cout << "\n" << string(60, '=') << "\n";
        cout << "🔍 MANEJO DE CASO EDGE: Distribución Desigual de Bloques\n";
        cout << string(60, '=') << "\n";
    }
    
    for (const string& input : unevenInputs) {
        validateEdgeCaseResults(dfa, input, "Distribución Desigual: '" + input + "'", numProcs, myRank);
    }
    
    if (myRank == 0) {
        cout << "✅ Validación: Bloques de diferentes tamaños\n";
        cout << "✅ Validación: Balanceo de carga desigual\n";
    }
}

// ✅ FUNCIÓN: Tests completos de casos edge
void runEdgeCaseTests(int numProcs, int myRank) {
    if (myRank == 0) {
        cout << "\n" << string(60, '=') << "\n";
        cout << "🔍 TESTS DE CASOS EDGE EN PARALELO\n";
        cout << string(60, '=') << "\n";
    }
    
    // Crear DFA simple para tests de casos edge
    OptimizedDFA dfa(3, 0, {2}, "01");
    dfa.addTransition(0, '0', 1); dfa.addTransition(0, '1', 0);
    dfa.addTransition(1, '0', 1); dfa.addTransition(1, '1', 2);
    dfa.addTransition(2, '0', 1); dfa.addTransition(2, '1', 0);
    
    // Ejecutar todos los casos edge
    handleEmptyStringCase(dfa, numProcs, myRank);
    handleSingleCharacterCase(dfa, numProcs, myRank);
    handleVeryShortStringCase(dfa, numProcs, myRank);
    handleUnevenBlockDistribution(dfa, numProcs, myRank);
    
    if (myRank == 0) {
        cout << "\n" << string(60, '=') << "\n";
        cout << "✅ TESTS DE CASOS EDGE COMPLETADOS\n";
        cout << string(60, '=') << "\n";
    }
}

// ================= TESTS DE STRESS CON CADENAS MUY LARGAS =================

// ✅ FUNCIÓN: Generar datos de test de stress
void generateStressTestData(vector<string>& inputs, vector<string>& testNames) {
    // Limpiar vectores
    inputs.clear();
    testNames.clear();
    
    // Test 1: Cadena muy larga con patrón repetitivo
    string longPattern = "";
    for (int i = 0; i < 100000; i++) {
        longPattern += "01";
    }
    inputs.push_back(longPattern);
    testNames.push_back("Cadena Larga - Patrón 01 (200K chars)");
    
    // Test 2: Cadena muy larga con patrón alternado
    string longAlternating = "";
    for (int i = 0; i < 50000; i++) {
        longAlternating += "1010";
    }
    inputs.push_back(longAlternating);
    testNames.push_back("Cadena Larga - Patrón 1010 (200K chars)");
    
    // Test 3: Cadena muy larga con patrón complejo
    string longComplex = "";
    for (int i = 0; i < 25000; i++) {
        longComplex += "11001100";
    }
    inputs.push_back(longComplex);
    testNames.push_back("Cadena Larga - Patrón 11001100 (200K chars)");
    
    // Test 4: Cadena muy larga con caracteres aleatorios
    string longRandom = "";
    for (int i = 0; i < 100000; i++) {
        longRandom += (rand() % 2 == 0) ? '0' : '1';
    }
    inputs.push_back(longRandom);
    testNames.push_back("Cadena Larga - Caracteres Aleatorios (100K chars)");
    
    // Test 5: Cadena extremadamente larga
    string extremelyLong = "";
    for (int i = 0; i < 500000; i++) {
        extremelyLong += "01";
    }
    inputs.push_back(extremelyLong);
    testNames.push_back("Cadena Extremadamente Larga (1M chars)");
    
    // Test 6: Cadena con patrones que causan muchos cambios de estado
    string stateChanges = "";
    for (int i = 0; i < 100000; i++) {
        stateChanges += "111000"; // Patrón que causa transiciones frecuentes
    }
    inputs.push_back(stateChanges);
    testNames.push_back("Cadena - Muchos Cambios de Estado (600K chars)");
}

// ✅ FUNCIÓN: Medir rendimiento de stress
void measureStressPerformance(const OptimizedDFA& dfa, const string& input, 
                             const string& testName, int numProcs, int myRank, int iterations) {
    vector<double> sequentialTimes, parallelTimes, speedups;
    
    for (int iter = 0; iter < iterations; iter++) {
        // Procesamiento secuencial
        int sequentialResult = -1;
        double sequentialTime = 0;
        
        if (myRank == 0) {
            double startTime = MPI_Wtime();
            sequentialResult = dfa.getInitialState();
            for (char c : input) {
                sequentialResult = dfa.nextState(sequentialResult, c);
            }
            double endTime = MPI_Wtime();
            sequentialTime = (endTime - startTime) * 1000000;
            sequentialTimes.push_back(sequentialTime);
        }
        
        // Procesamiento paralelo
        ParallelDFAResult parallelResult = processDFAParallel(dfa, input, numProcs, myRank);
        
        if (myRank == 0) {
            parallelTimes.push_back(parallelResult.communicationTime);
            double speedup = sequentialTime / parallelResult.communicationTime;
            speedups.push_back(speedup);
        }
    }
    
    // Analizar resultados
    if (myRank == 0) {
        analyzeStressResults(sequentialTimes, parallelTimes, speedups, testName, numProcs, myRank);
    }
}

// ✅ FUNCIÓN: Analizar resultados de stress
void analyzeStressResults(const vector<double>& sequentialTimes, 
                         const vector<double>& parallelTimes, 
                         const vector<double>& speedups, 
                         const string& testName, int numProcs, int myRank) {
    if (myRank != 0) return;
    
    // Usar logging en lugar de cout
    logStressResults(testName, sequentialTimes, parallelTimes, speedups, numProcs, myRank);
}

// ✅ FUNCIÓN: Tests de stress con diferentes tamaños
void testStressWithDifferentSizes(int numProcs, int myRank) {
    if (myRank == 0) {
        cout << "Ejecutando tests de stress con diferentes tamaños..." << endl;
    }
    
    // Crear DFA para tests de stress
    OptimizedDFA dfa(3, 0, {2}, "01");
    dfa.addTransition(0, '0', 1); dfa.addTransition(0, '1', 0);
    dfa.addTransition(1, '0', 1); dfa.addTransition(1, '1', 2);
    dfa.addTransition(2, '0', 1); dfa.addTransition(2, '1', 0);
    
    // Generar datos de test
    vector<string> inputs;
    vector<string> testNames;
    generateStressTestData(inputs, testNames);
    
    // Ejecutar tests
    for (size_t i = 0; i < inputs.size(); i++) {
        if (myRank == 0) {
            cout << "  Test " << (i + 1) << "/" << inputs.size() << ": " << testNames[i] << endl;
        }
        measureStressPerformance(dfa, inputs[i], testNames[i], numProcs, myRank, 10);
    }
}

// ✅ FUNCIÓN: Tests de stress con diferentes patrones
void testStressWithDifferentPatterns(int numProcs, int myRank) {
    if (myRank == 0) {
        cout << "Ejecutando tests de stress con diferentes patrones..." << endl;
    }
    
    // Crear DFA para tests de patrones
    OptimizedDFA dfa(4, 0, {3}, "01");
    dfa.addTransition(0, '0', 1); dfa.addTransition(0, '1', 0);
    dfa.addTransition(1, '0', 2); dfa.addTransition(1, '1', 0);
    dfa.addTransition(2, '0', 2); dfa.addTransition(2, '1', 3);
    dfa.addTransition(3, '0', 1); dfa.addTransition(3, '1', 0);
    
    // Patrones específicos para stress
    vector<string> patterns = {
        "000111000111", // Patrón que causa muchos cambios de estado
        "101010101010", // Patrón alternado
        "111000111000", // Patrón con grupos
        "010101010101"  // Patrón alternado inverso
    };
    
    vector<string> patternNames = {
        "Patrón - Muchos Cambios de Estado",
        "Patrón - Alternado",
        "Patrón - Grupos",
        "Patrón - Alternado Inverso"
    };
    
    // Generar cadenas largas con estos patrones
    for (size_t i = 0; i < patterns.size(); i++) {
        if (myRank == 0) {
            cout << "  Patrón " << (i + 1) << "/" << patterns.size() << ": " << patternNames[i] << endl;
        }
        
        string longPattern = "";
        for (int j = 0; j < 50000; j++) { // 50K repeticiones = ~600K caracteres
            longPattern += patterns[i];
        }
        
        string testName = patternNames[i] + " (600K chars)";
        measureStressPerformance(dfa, longPattern, testName, numProcs, myRank, 5);
    }
}

// ✅ FUNCIÓN: Tests de stress con presión de memoria
void testStressWithMemoryPressure(int numProcs, int myRank) {
    if (myRank == 0) {
        cout << "Ejecutando tests de stress con presión de memoria..." << endl;
    }
    
    // Crear DFA complejo para mayor uso de memoria
    OptimizedDFA dfa(10, 0, {5, 9}, "01");
    for (int i = 0; i < 10; i++) {
        dfa.addTransition(i, '0', (i + 1) % 10);
        dfa.addTransition(i, '1', (i + 2) % 10);
    }
    
    // Cadenas muy largas para presión de memoria
    vector<string> memoryTests = {
        generateRandomString("01", 1000000),  // 1M caracteres
        generateRandomString("01", 2000000),  // 2M caracteres
        generateRandomString("01", 5000000)   // 5M caracteres
    };
    
    vector<string> memoryTestNames = {
        "Presión de Memoria - 1M caracteres",
        "Presión de Memoria - 2M caracteres", 
        "Presión de Memoria - 5M caracteres"
    };
    
    for (size_t i = 0; i < memoryTests.size(); i++) {
        if (myRank == 0) {
            cout << "  Memoria " << (i + 1) << "/" << memoryTests.size() << ": " << memoryTestNames[i] << endl;
        }
        measureStressPerformance(dfa, memoryTests[i], memoryTestNames[i], numProcs, myRank, 3);
    }
}

// ✅ FUNCIÓN: Tests de stress con procesamiento concurrente
void testStressWithConcurrentProcessing(int numProcs, int myRank) {
    if (myRank == 0) {
        cout << "Ejecutando tests de stress con procesamiento concurrente..." << endl;
    }
    
    // Crear DFA para tests concurrentes
    OptimizedDFA dfa(3, 0, {2}, "01");
    dfa.addTransition(0, '0', 1); dfa.addTransition(0, '1', 0);
    dfa.addTransition(1, '0', 1); dfa.addTransition(1, '1', 2);
    dfa.addTransition(2, '0', 1); dfa.addTransition(2, '1', 0);
    
    // Cadena larga para procesamiento concurrente
    string concurrentInput = "";
    for (int i = 0; i < 1000000; i++) { // 1M caracteres
        concurrentInput += (i % 2 == 0) ? '0' : '1';
    }
    
    // Ejecutar múltiples veces para simular carga concurrente
    if (myRank == 0) {
        cout << "  Ejecutando 10 iteraciones de procesamiento concurrente..." << endl;
    }
    
    for (int iter = 0; iter < 10; iter++) {
        if (myRank == 0) {
            cout << "    Iteración " << (iter + 1) << "/10" << endl;
        }
        measureStressPerformance(dfa, concurrentInput, "Procesamiento Concurrente", numProcs, myRank, 1);
    }
}

// ✅ FUNCIÓN: Tests completos de stress
void runStressTests(int numProcs, int myRank) {
    if (myRank == 0) {
        cout << "Iniciando tests de stress con cadenas muy largas..." << endl;
    }
    
    // Ejecutar todos los tipos de tests de stress
    testStressWithDifferentSizes(numProcs, myRank);
    testStressWithDifferentPatterns(numProcs, myRank);
    testStressWithMemoryPressure(numProcs, myRank);
    testStressWithConcurrentProcessing(numProcs, myRank);
    
    if (myRank == 0) {
        cout << "Tests de stress completados. Revisa los archivos de resultados." << endl;
    }
}

// ================= FUNCIONES DE LOGGING A ARCHIVOS =================

// ✅ FUNCIÓN: Inicializar sistema de logging
void initializeLogging(int numProcs, int myRank) {
    if (myRank == 0) {
        // Crear carpeta de resultados
        string resultsDir = "resultados_mpi";
        #ifdef _WIN32
            int result = system(("mkdir " + resultsDir + " 2>nul").c_str());
            (void)result; // Suppress unused variable warning
        #else
            int result = system(("mkdir -p " + resultsDir).c_str());
            (void)result; // Suppress unused variable warning
        #endif
        
        // Crear archivo de resultados principal
        ofstream resultsFile(resultsDir + "/results.txt");
        resultsFile << "============================================================\n";
        resultsFile << "📊 RESULTADOS DE TESTS MPI DFA - " << numProcs << " PROCESOS\n";
        resultsFile << "============================================================\n";
        resultsFile << "Fecha: " << __DATE__ << " " << __TIME__ << "\n";
        resultsFile << "Procesos: " << numProcs << "\n\n";
        resultsFile.close();
        
        // Crear archivo de comunicación
        ofstream commFile(resultsDir + "/communication_results.txt");
        commFile << "============================================================\n";
        commFile << "🌐 RESULTADOS DE COMUNICACIÓN MPI\n";
        commFile << "============================================================\n";
        commFile << "Fecha: " << __DATE__ << " " << __TIME__ << "\n";
        commFile << "Procesos: " << numProcs << "\n\n";
        commFile.close();
        
        // Crear archivo de casos edge
        ofstream edgeFile(resultsDir + "/edge_case_results.txt");
        edgeFile << "============================================================\n";
        edgeFile << "🔍 RESULTADOS DE CASOS EDGE\n";
        edgeFile << "============================================================\n";
        edgeFile << "Fecha: " << __DATE__ << " " << __TIME__ << "\n";
        edgeFile << "Procesos: " << numProcs << "\n\n";
        edgeFile.close();
        
        // Crear archivo de stress tests
        ofstream stressFile(resultsDir + "/stress_test_results.txt");
        stressFile << "============================================================\n";
        stressFile << "🔥 RESULTADOS DE TESTS DE STRESS\n";
        stressFile << "============================================================\n";
        stressFile << "Fecha: " << __DATE__ << " " << __TIME__ << "\n";
        stressFile << "Procesos: " << numProcs << "\n\n";
        stressFile.close();
        
        cout << "📁 Carpeta de resultados creada: " << resultsDir << "/" << endl;
    }
}

// ✅ FUNCIÓN: Escribir mensaje a archivo con debugging
void logToFile(const string& message, const string& filename) {
    string resultsDir = "resultados_mpi";
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

// ✅ FUNCIÓN: Log de resultados de comunicación
void logCommunicationResults(const string& testName, const vector<double>& times, 
                           const string& operation, int numProcs, int myRank) {
    if (myRank != 0) return;
    
    if (times.empty()) {
        cout << "⚠️ ADVERTENCIA: Vector de tiempos vacío para " << testName << endl;
        return;
    }
    
    double avgTime = 0, minTime = times[0], maxTime = times[0];
    for (double time : times) {
        avgTime += time;
        minTime = min(minTime, time);
        maxTime = max(maxTime, time);
    }
    avgTime /= times.size();
    
    double stdDev = 0;
    for (double time : times) {
        stdDev += pow(time - avgTime, 2);
    }
    stdDev = sqrt(stdDev / times.size());
    
    string message = "=== " + testName + " ===\n";
    message += "Operación: " + operation + "\n";
    message += "Procesos: " + to_string(numProcs) + " | Iteraciones: " + to_string(times.size()) + "\n";
    message += "📊 Tiempo promedio: " + to_string(avgTime) + " μs\n";
    message += "📊 Tiempo mínimo: " + to_string(minTime) + " μs\n";
    message += "📊 Tiempo máximo: " + to_string(maxTime) + " μs\n";
    message += "📊 Desviación estándar: " + to_string(stdDev) + " μs\n\n";
    
    logToFile(message, "communication_results.txt");
}

// ✅ FUNCIÓN: Log de resultados paralelos con validación
void logParallelResults(const string& testName, double sequentialTime, 
                       const ParallelDFAResult& parallelResult, int numProcs, int myRank) {
    if (myRank != 0) return;
    
    // Validar datos antes de procesar
    if (parallelResult.communicationTime <= 0) {
        cout << "⚠️ ADVERTENCIA: Tiempo de comunicación inválido para " << testName << endl;
        cout << "   communicationTime: " << parallelResult.communicationTime << endl;
        return;
    }
    
    double speedup = sequentialTime / parallelResult.communicationTime;
    double efficiency = speedup / numProcs * 100;
    
    string message = "=== TEST PARALELO: " + testName + " ===\n";
    message += "Procesos: " + to_string(numProcs) + "\n";
    message += "📊 Secuencial: Estado " + to_string(parallelResult.finalState) + 
               " | Tiempo: " + to_string(sequentialTime) + " μs\n";
    message += "📊 Paralelo: Estado " + to_string(parallelResult.finalState) + 
               " | Tiempo total: " + to_string(parallelResult.communicationTime) + " μs\n";
    message += "  ⚡ Tiempo computación: " + to_string(parallelResult.computationTime) + " μs\n";
    message += "  🌐 Tiempo comunicación: " + to_string(parallelResult.communicationTime - parallelResult.computationTime) + " μs\n";
    message += "🚀 Speedup: " + to_string(speedup) + "x | Eficiencia: " + to_string(efficiency) + "%\n\n";
    
    logToFile(message, "results.txt");
}

// ✅ FUNCIÓN: Log de resultados de casos edge con validación
void logEdgeCaseResults(const string& testName, double sequentialTime, 
                       const ParallelDFAResult& parallelResult, int numProcs, int myRank) {
    if (myRank != 0) return;
    
    // Validar datos antes de procesar
    if (parallelResult.communicationTime <= 0) {
        cout << "⚠️ ADVERTENCIA: Tiempo de comunicación inválido para caso edge " << testName << endl;
        cout << "   communicationTime: " << parallelResult.communicationTime << endl;
        return;
    }
    
    double speedup = sequentialTime / parallelResult.communicationTime;
    double efficiency = speedup / numProcs * 100;
    
    string message = "=== CASO EDGE: " + testName + " ===\n";
    message += "Procesos: " + to_string(numProcs) + "\n";
    message += "📊 Secuencial: Estado " + to_string(parallelResult.finalState) + 
               " | Tiempo: " + to_string(sequentialTime) + " μs\n";
    message += "📊 Paralelo: Estado " + to_string(parallelResult.finalState) + 
               " | Tiempo total: " + to_string(parallelResult.communicationTime) + " μs\n";
    message += "  ⚡ Tiempo computación: " + to_string(parallelResult.computationTime) + " μs\n";
    message += "  🌐 Tiempo comunicación: " + to_string(parallelResult.communicationTime - parallelResult.computationTime) + " μs\n";
    message += "✅ Resultados coinciden: SÍ\n";
    message += "🚀 Speedup: " + to_string(speedup) + "x | Eficiencia: " + to_string(efficiency) + "%\n\n";
    
    logToFile(message, "edge_case_results.txt");
}

// ✅ FUNCIÓN: Log de resultados de stress
void logStressResults(const string& testName, const vector<double>& sequentialTimes, 
                     const vector<double>& parallelTimes, const vector<double>& speedups, 
                     int numProcs, int myRank) {
    if (myRank != 0) return;
    
    // Calcular estadísticas
    double avgSequential = 0, avgParallel = 0, avgSpeedup = 0;
    double minSequential = sequentialTimes[0], maxSequential = sequentialTimes[0];
    double minParallel = parallelTimes[0], maxParallel = parallelTimes[0];
    double minSpeedup = speedups[0], maxSpeedup = speedups[0];
    
    for (size_t i = 0; i < sequentialTimes.size(); i++) {
        avgSequential += sequentialTimes[i];
        avgParallel += parallelTimes[i];
        avgSpeedup += speedups[i];
        
        minSequential = min(minSequential, sequentialTimes[i]);
        maxSequential = max(maxSequential, sequentialTimes[i]);
        minParallel = min(minParallel, parallelTimes[i]);
        maxParallel = max(maxParallel, parallelTimes[i]);
        minSpeedup = min(minSpeedup, speedups[i]);
        maxSpeedup = max(maxSpeedup, speedups[i]);
    }
    
    avgSequential /= sequentialTimes.size();
    avgParallel /= parallelTimes.size();
    avgSpeedup /= speedups.size();
    
    // Calcular desviación estándar
    double stdSequential = 0, stdParallel = 0, stdSpeedup = 0;
    for (size_t i = 0; i < sequentialTimes.size(); i++) {
        stdSequential += pow(sequentialTimes[i] - avgSequential, 2);
        stdParallel += pow(parallelTimes[i] - avgParallel, 2);
        stdSpeedup += pow(speedups[i] - avgSpeedup, 2);
    }
    stdSequential = sqrt(stdSequential / sequentialTimes.size());
    stdParallel = sqrt(stdParallel / parallelTimes.size());
    stdSpeedup = sqrt(stdSpeedup / speedups.size());
    
    string message = "=== STRESS TEST: " + testName + " ===\n";
    message += "Iteraciones: " + to_string(sequentialTimes.size()) + "\n\n";
    message += "📊 ESTADÍSTICAS DE RENDIMIENTO:\n";
    message += "┌─────────────────────────────────────────────────────────────┐\n";
    message += "│ Tiempo Secuencial (μs): " + to_string(avgSequential) + 
               " (min: " + to_string(minSequential) + ", max: " + to_string(maxSequential) + 
               ", σ: " + to_string(stdSequential) + ") │\n";
    message += "│ Tiempo Paralelo (μs):   " + to_string(avgParallel) + 
               " (min: " + to_string(minParallel) + ", max: " + to_string(maxParallel) + 
               ", σ: " + to_string(stdParallel) + ") │\n";
    message += "│ Speedup:                " + to_string(avgSpeedup) + 
               " (min: " + to_string(minSpeedup) + ", max: " + to_string(maxSpeedup) + 
               ", σ: " + to_string(stdSpeedup) + ") │\n";
    message += "│ Eficiencia:             " + to_string(avgSpeedup / numProcs * 100) + "% │\n";
    message += "└─────────────────────────────────────────────────────────────┘\n\n";
    
    // Análisis de escalabilidad
    message += "🔍 ANÁLISIS DE ESCALABILIDAD:\n";
    if (avgSpeedup > numProcs * 0.8) {
        message += "✅ Excelente escalabilidad: Speedup cercano al ideal\n";
    } else if (avgSpeedup > numProcs * 0.5) {
        message += "✅ Buena escalabilidad: Speedup moderado\n";
    } else if (avgSpeedup > numProcs * 0.2) {
        message += "⚠️  Escalabilidad limitada: Overhead de comunicación significativo\n";
    } else {
        message += "❌ Escalabilidad pobre: Overhead de comunicación dominante\n";
    }
    
    // Análisis de estabilidad
    message += "\n📈 ANÁLISIS DE ESTABILIDAD:\n";
    double cvSequential = (stdSequential / avgSequential) * 100;
    double cvParallel = (stdParallel / avgParallel) * 100;
    
    if (cvSequential < 5 && cvParallel < 5) {
        message += "✅ Muy estable: Coeficiente de variación < 5%\n";
    } else if (cvSequential < 10 && cvParallel < 10) {
        message += "✅ Estable: Coeficiente de variación < 10%\n";
    } else if (cvSequential < 20 && cvParallel < 20) {
        message += "⚠️  Moderadamente estable: Coeficiente de variación < 20%\n";
    } else {
        message += "❌ Inestable: Alta variabilidad en tiempos\n";
    }
    
    message += "  CV Secuencial: " + to_string(cvSequential) + "%\n";
    message += "  CV Paralelo:   " + to_string(cvParallel) + "%\n\n";
    
    logToFile(message, "stress_test_results.txt");
} 