#include "communication_tests.h"
#include "utils.h"
#include <fstream>
#include <iomanip>

// ✅ FUNCIÓN DE TEST DE DIVISIÓN
void testBlockDivision(const string& testName, const string& input, int numProcesses) {
    cout << "\n=== TEST: " << testName << " ===\n";
    cout << "Entrada: [" << input << "] (tamaño: " << input.size() << ")\n";
    cout << "Procesos: " << numProcesses << "\n";
    
    double startTime = MPI_Wtime();
    vector<InputBlock> blocks = divideInputIntoBlocks(input, numProcesses);
    double endTime = MPI_Wtime();
    double divisionTime = (endTime - startTime) * 1000000; // Convertir a microsegundos
    
    cout << "⏱️  Tiempo de división: " << divisionTime << " μs\n";
    cout << "📊 Distribución de bloques:\n";
    
    for (const auto& block : blocks) {
        cout << "  Proceso " << block.processId << ": ";
        cout << "[" << block.data << "] ";
        cout << "(pos " << block.startPos << "-" << block.endPos << ", ";
        cout << "tamaño " << block.data.size() << ")\n";
    }
    
    // Validar la división
    bool isValid = validateBlockDivision(input, blocks);
    cout << (isValid ? "✅ División válida" : "❌ División inválida") << "\n";
    
    // Calcular balanceo
    if (!blocks.empty()) {
        int minSize = blocks[0].data.size();
        int maxSize = blocks[0].data.size();
        for (const auto& block : blocks) {
            minSize = min(minSize, (int)block.data.size());
            maxSize = max(maxSize, (int)block.data.size());
        }
        double balanceRatio = (double)minSize / maxSize;
        cout << "⚖️  Balanceo: " << (balanceRatio * 100) << "% (min: " << minSize << ", max: " << maxSize << ")\n";
    }
}

// ✅ FUNCIÓN DE TESTS COMPLETOS DE DIVISIÓN
void runBlockDivisionTests() {
    cout << "\n" << string(60, '=') << "\n";
    cout << "🧪 TESTS DE DIVISIÓN DE BLOQUES PARA MPI\n";
    cout << string(60, '=') << "\n";
    
    // Test 1: Cadena corta
    testBlockDivision("Cadena Corta", "abc", 2);
    
    // Test 2: Cadena mediana
    testBlockDivision("Cadena Mediana", "abcdefghijklmnop", 4);
    
    // Test 3: Cadena larga
    string longString = "";
    for (int i = 0; i < 1000; i++) {
        longString += (char)('a' + (i % 26));
    }
    testBlockDivision("Cadena Larga (1000 chars)", longString, 8);
    
    // Test 4: Cadena muy larga
    string veryLongString = "";
    for (int i = 0; i < 100000; i++) {
        veryLongString += (char)('0' + (i % 10));
    }
    testBlockDivision("Cadena Muy Larga (100K chars)", veryLongString, 16);
    
    // Test 5: Casos edge
    testBlockDivision("Cadena Vacía", "", 4);
    testBlockDivision("Un Solo Carácter", "x", 8);
    testBlockDivision("Menos Caracteres que Procesos", "ab", 4);
    
    // Test 6: Diferentes números de procesos
    testBlockDivision("Muchos Procesos", "abcdefghijklmnopqrstuvwxyz", 32);
    testBlockDivision("Pocos Procesos", "abcdefghijklmnopqrstuvwxyz", 2);
    
    cout << "\n" << string(60, '=') << "\n";
    cout << "✅ TESTS DE DIVISIÓN COMPLETADOS\n";
    cout << string(60, '=') << "\n";
}

// ✅ FUNCIÓN MEJORADA: MEDICIÓN CON MÚLTIPLES ITERACIONES Y PROMEDIOS
void measureCommunicationTimeImproved(const string& testName, const string& input, int numProcs, int myRank, int iterations = 100) {
    if (myRank == 0) {
        cout << "\n=== MEDICIÓN MEJORADA: " << testName << " ===\n";
        cout << "Entrada: [" << input << "] (tamaño: " << input.size() << ")\n";
        cout << "Procesos: " << numProcs << " | Iteraciones: " << iterations << "\n";
    }
    
    // Arrays para almacenar tiempos de cada iteración
    vector<double> broadcast1Times(iterations);
    vector<double> broadcast2Times(iterations);
    vector<double> scatter1Times(iterations);
    vector<double> scatter2Times(iterations);
    vector<double> scatter3Times(iterations);
    vector<double> scatter4Times(iterations);
    vector<double> totalTimes(iterations);
    
    // Ejecutar múltiples iteraciones
    for (int iter = 0; iter < iterations; iter++) {
        double totalStartTime = MPI_Wtime();
        
        // El proceso raíz divide la entrada en bloques
        vector<InputBlock> allBlocks;
        if (myRank == 0) {
            allBlocks = divideInputIntoBlocks(input, numProcs);
        }
        
        // Broadcast del tamaño de la entrada
        double broadcast1Start = MPI_Wtime();
        int inputSize = input.size();
        MPI_Bcast(&inputSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
        double broadcast1End = MPI_Wtime();
        broadcast1Times[iter] = (broadcast1End - broadcast1Start) * 1000000;
        
        // Broadcast del número de bloques
        double broadcast2Start = MPI_Wtime();
        int numBlocks = numProcs;
        MPI_Bcast(&numBlocks, 1, MPI_INT, 0, MPI_COMM_WORLD);
        double broadcast2End = MPI_Wtime();
        broadcast2Times[iter] = (broadcast2End - broadcast2Start) * 1000000;
        
        // Preparar arrays para Scatter
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
        
        // Scatter de tamaños de bloques
        double scatter1Start = MPI_Wtime();
        int myBlockSize;
        MPI_Scatter(blockSizes.data(), 1, MPI_INT, &myBlockSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
        double scatter1End = MPI_Wtime();
        scatter1Times[iter] = (scatter1End - scatter1Start) * 1000000;
        
        // Scatter de posiciones iniciales
        double scatter2Start = MPI_Wtime();
        int myStartPos;
        MPI_Scatter(startPositions.data(), 1, MPI_INT, &myStartPos, 1, MPI_INT, 0, MPI_COMM_WORLD);
        double scatter2End = MPI_Wtime();
        scatter2Times[iter] = (scatter2End - scatter2Start) * 1000000;
        
        // Scatter de posiciones finales
        double scatter3Start = MPI_Wtime();
        int myEndPos;
        MPI_Scatter(endPositions.data(), 1, MPI_INT, &myEndPos, 1, MPI_INT, 0, MPI_COMM_WORLD);
        double scatter3End = MPI_Wtime();
        scatter3Times[iter] = (scatter3End - scatter3Start) * 1000000;
        
        // Preparar array de caracteres para Scatterv
        vector<char> inputChars;
        if (myRank == 0) {
            inputChars.resize(inputSize);
            for (int i = 0; i < inputSize; i++) {
                inputChars[i] = input[i];
            }
        }
        
        // Scatterv de datos de entrada
        double scatter4Start = MPI_Wtime();
        vector<char> myBlockData(myBlockSize);
        
        vector<int> sendCounts(numProcs);
        vector<int> displacements(numProcs);
        
        if (myRank == 0) {
            for (int i = 0; i < numProcs; i++) {
                sendCounts[i] = allBlocks[i].data.size();
                displacements[i] = allBlocks[i].startPos;
            }
        }
        
        MPI_Scatterv(inputChars.data(), sendCounts.data(), displacements.data(), MPI_CHAR,
                     myBlockData.data(), myBlockSize, MPI_CHAR, 0, MPI_COMM_WORLD);
        
        double scatter4End = MPI_Wtime();
        scatter4Times[iter] = (scatter4End - scatter4Start) * 1000000;
        
        double totalEndTime = MPI_Wtime();
        totalTimes[iter] = (totalEndTime - totalStartTime) * 1000000;
        
        // Validar que los datos se transmitieron correctamente
        if (iter == 0) { // Solo validar en la primera iteración
            string expectedBlock = "";
            if (myRank < (int)input.size()) {
                int start = myStartPos;
                int end = min(myEndPos + 1, (int)input.size());
                expectedBlock = input.substr(start, end - start);
            }
            
            string receivedBlock = string(myBlockData.begin(), myBlockData.end());
            if (expectedBlock != receivedBlock) {
                cout << "❌ ERROR Proceso " << myRank << ": Datos incorrectos!\n";
                cout << "  Esperado: [" << expectedBlock << "]\n";
                cout << "  Recibido: [" << receivedBlock << "]\n";
            }
        }
    }
    
    // Calcular estadísticas
    auto calculateStats = [](const vector<double>& times) -> tuple<double, double, double, double> {
        double sum = 0, min_val = times[0], max_val = times[0];
        for (double t : times) {
            sum += t;
            min_val = min(min_val, t);
            max_val = max(max_val, t);
        }
        double avg = sum / times.size();
        
        // Calcular desviación estándar
        double variance = 0;
        for (double t : times) {
            variance += (t - avg) * (t - avg);
        }
        double std_dev = sqrt(variance / times.size());
        
        return {avg, min_val, max_val, std_dev};
    };
    
    auto [avg_b1, min_b1, max_b1, std_b1] = calculateStats(broadcast1Times);
    auto [avg_b2, min_b2, max_b2, std_b2] = calculateStats(broadcast2Times);
    auto [avg_s1, min_s1, max_s1, std_s1] = calculateStats(scatter1Times);
    auto [avg_s2, min_s2, max_s2, std_s2] = calculateStats(scatter2Times);
    auto [avg_s3, min_s3, max_s3, std_s3] = calculateStats(scatter3Times);
    auto [avg_s4, min_s4, max_s4, std_s4] = calculateStats(scatter4Times);
    auto [avg_total, min_total, max_total, std_total] = calculateStats(totalTimes);
    
    // Mostrar resultados estadísticos
    cout << "Proceso " << myRank << " - Estadísticas (" << iterations << " iteraciones):\n";
    cout << "  📊 Broadcast inputSize: " << fixed << setprecision(2) 
         << avg_b1 << " μs (min: " << min_b1 << ", max: " << max_b1 << ", σ: " << std_b1 << ")\n";
    cout << "  📊 Broadcast numBlocks: " << avg_b2 << " μs (min: " << min_b2 << ", max: " << max_b2 << ", σ: " << std_b2 << ")\n";
    cout << "  📊 Scatter blockSizes: " << avg_s1 << " μs (min: " << min_s1 << ", max: " << max_s1 << ", σ: " << std_s1 << ")\n";
    cout << "  📊 Scatter startPos: " << avg_s2 << " μs (min: " << min_s2 << ", max: " << max_s2 << ", σ: " << std_s2 << ")\n";
    cout << "  📊 Scatter endPos: " << avg_s3 << " μs (min: " << min_s3 << ", max: " << max_s3 << ", σ: " << std_s3 << ")\n";
    cout << "  📊 Scatterv data: " << avg_s4 << " μs (min: " << min_s4 << ", max: " << max_s4 << ", σ: " << std_s4 << ")\n";
    cout << "  ⚡ TOTAL PROMEDIO: " << avg_total << " μs (min: " << min_total << ", max: " << max_total << ", σ: " << std_total << ")\n";
    
    // Sincronización para evitar mezcla de output
    MPI_Barrier(MPI_COMM_WORLD);
}

// ✅ FUNCIÓN DE TESTS COMPLETOS DE COMUNICACIÓN MPI
void runMPICommunicationTests(int numProcs, int myRank) {
    if (myRank == 0) {
        cout << "\n" << string(60, '=') << "\n";
        cout << "🌐 TESTS DE COMUNICACIÓN MPI MEJORADOS\n";
        cout << string(60, '=') << "\n";
    }
    
    // Test 1: Cadena corta con medición mejorada
    measureCommunicationTimeImproved("Cadena Corta", "abc", numProcs, myRank, 50);
    
    // Test 2: Cadena mediana con medición mejorada
    measureCommunicationTimeImproved("Cadena Mediana", "abcdefghijklmnop", numProcs, myRank, 50);
    
    // Test 3: Cadena vacía con medición mejorada
    measureCommunicationTimeImproved("Cadena Vacía", "", numProcs, myRank, 50);
    
    // Test 4: Un solo carácter con medición mejorada
    measureCommunicationTimeImproved("Un Solo Carácter", "x", numProcs, myRank, 50);
    
    if (myRank == 0) {
        cout << "\n" << string(60, '=') << "\n";
        cout << "✅ TESTS DE COMUNICACIÓN MPI COMPLETADOS\n";
        cout << string(60, '=') << "\n";
    }
}

// ✅ FUNCIÓN: Test de comunicación punto a punto
void testPointToPointCommunication(int numProcs, int myRank) {
    if (myRank == 0) {
        cout << "Ejecutando test de comunicación punto a punto..." << endl;
    }
    
    vector<double> times;
    int iterations = 100;
    
    for (int iter = 0; iter < iterations; iter++) {
        double startTime = MPI_Wtime();
        
        if (myRank == 0) {
            // Proceso 0 envía datos a todos los demás procesos
            for (int dest = 1; dest < numProcs; dest++) {
                int data = 42 + iter;
                MPI_Send(&data, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
            }
        } else {
            // Otros procesos reciben datos
            int receivedData;
            MPI_Recv(&receivedData, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        
        double endTime = MPI_Wtime();
        double elapsedTime = (endTime - startTime) * 1000000; // Convertir a microsegundos
        times.push_back(elapsedTime);
    }
    
    // Log de resultados usando la función de parallel_dfa.cpp
    if (myRank == 0) {
        logCommunicationResults("Comunicación Punto a Punto", times, "Send/Recv", numProcs, myRank);
    }
}

// ✅ FUNCIÓN: Test de comunicación colectiva
void testCollectiveCommunication(int numProcs, int myRank) {
    if (myRank == 0) {
        cout << "Ejecutando test de comunicación colectiva..." << endl;
    }
    
    vector<double> times;
    int iterations = 50;
    
    for (int iter = 0; iter < iterations; iter++) {
        double startTime = MPI_Wtime();
        
        // Cada proceso tiene un valor local
        int localValue = myRank + iter;
        int globalSum;
        
        // Operación de reducción (suma)
        MPI_Allreduce(&localValue, &globalSum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        
        double endTime = MPI_Wtime();
        double elapsedTime = (endTime - startTime) * 1000000;
        times.push_back(elapsedTime);
    }
    
    // Log de resultados usando la función de parallel_dfa.cpp
    if (myRank == 0) {
        logCommunicationResults("Comunicación Colectiva", times, "Allreduce", numProcs, myRank);
    }
}

// ✅ FUNCIÓN: Test de comunicación con datos derivados
void testDerivedDataTypes(int numProcs, int myRank) {
    if (myRank == 0) {
        cout << "Ejecutando test de tipos de datos derivados..." << endl;
    }
    
    vector<double> times;
    int iterations = 30;
    
    // Definir estructura de datos
    struct TestData {
        int id;
        double value;
        char flag;
    };
    
    // Crear tipo derivado MPI
    MPI_Datatype MPI_TEST_DATA;
    int blocklengths[3] = {1, 1, 1};
    MPI_Aint displacements[3];
    MPI_Datatype types[3] = {MPI_INT, MPI_DOUBLE, MPI_CHAR};
    
    TestData dummy;
    MPI_Aint base_address;
    MPI_Get_address(&dummy, &base_address);
    MPI_Get_address(&dummy.id, &displacements[0]);
    MPI_Get_address(&dummy.value, &displacements[1]);
    MPI_Get_address(&dummy.flag, &displacements[2]);
    
    displacements[0] = MPI_Aint_diff(displacements[0], base_address);
    displacements[1] = MPI_Aint_diff(displacements[1], base_address);
    displacements[2] = MPI_Aint_diff(displacements[2], base_address);
    
    MPI_Type_create_struct(3, blocklengths, displacements, types, &MPI_TEST_DATA);
    MPI_Type_commit(&MPI_TEST_DATA);
    
    for (int iter = 0; iter < iterations; iter++) {
        double startTime = MPI_Wtime();
        
        TestData localData = {myRank, 3.14159 + iter, static_cast<char>('A' + (iter % 26))};
        TestData receivedData;
        
        if (myRank == 0) {
            // Proceso 0 recibe de todos los demás
            for (int source = 1; source < numProcs; source++) {
                MPI_Recv(&receivedData, 1, MPI_TEST_DATA, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        } else {
            // Otros procesos envían al proceso 0
            MPI_Send(&localData, 1, MPI_TEST_DATA, 0, 0, MPI_COMM_WORLD);
        }
        
        double endTime = MPI_Wtime();
        double elapsedTime = (endTime - startTime) * 1000000;
        times.push_back(elapsedTime);
    }
    
    MPI_Type_free(&MPI_TEST_DATA);
    
    // Log de resultados usando la función de parallel_dfa.cpp
    if (myRank == 0) {
        logCommunicationResults("Tipos de Datos Derivados", times, "Struct Send/Recv", numProcs, myRank);
    }
}

// ✅ FUNCIÓN: Test de comunicación con buffers optimizados
void testOptimizedBuffers(int numProcs, int myRank) {
    if (myRank == 0) {
        cout << "Ejecutando test de buffers optimizados..." << endl;
    }
    
    vector<double> times;
    int iterations = 20;
    int bufferSize = 10000; // 10KB de datos
    
    // Crear buffers
    vector<int> sendBuffer(bufferSize);
    vector<int> recvBuffer(bufferSize);
    
    // Inicializar buffer de envío
    for (int i = 0; i < bufferSize; i++) {
        sendBuffer[i] = myRank * 1000 + i;
    }
    
    for (int iter = 0; iter < iterations; iter++) {
        double startTime = MPI_Wtime();
        
        if (myRank == 0) {
            // Proceso 0 recibe de todos los demás
            for (int source = 1; source < numProcs; source++) {
                MPI_Recv(recvBuffer.data(), bufferSize, MPI_INT, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        } else {
            // Otros procesos envían al proceso 0
            MPI_Send(sendBuffer.data(), bufferSize, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
        
        double endTime = MPI_Wtime();
        double elapsedTime = (endTime - startTime) * 1000000;
        times.push_back(elapsedTime);
    }
    
    // Log de resultados usando la función de parallel_dfa.cpp
    if (myRank == 0) {
        logCommunicationResults("Buffers Optimizados", times, "Large Buffer Send/Recv", numProcs, myRank);
    }
} 