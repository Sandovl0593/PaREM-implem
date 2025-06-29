#include "dfa_tests.h"

// Función para ejecutar todos los tests secuenciales y retornar resultados
vector<TestResult> runSequentialTestsWithResults(const OptimizedDFA& dfa, const vector<DFATest>& tests) {
    vector<TestResult> results;
    cout << "=== EJECUTANDO TESTS SECUENCIALES ===\n\n";

    for (const auto& test : tests) {
        cout << "--- " << test.name << " ---\n";
        cout << "Descripción: " << test.description << "\n";
        cout << "Alfabeto: " << test.alphabet << "\n";
        cout << "Estados: " << test.numStates;
        if (test.finalStates.size() == 1) {
            cout << " (Final: " << test.finalStates[0] << ")\n";
        } else {
            cout << " (Finales: ";
            for (size_t i = 0; i < test.finalStates.size(); i++) {
                cout << test.finalStates[i];
                if (i < test.finalStates.size() - 1) cout << ", ";
            }
            cout << ")\n";
        }
        
        // Descripción compacta de la entrada
        if (test.input.empty()) {
            cout << "Entrada: vacía\n";
        } else if (test.input.length() > 1000) {
            cout << "Entrada: " << test.input.length() << " caracteres\n";
        } else if (test.input.length() > 100) {
            cout << "Entrada: " << test.input.substr(0, 10) << "..." << test.input.substr(test.input.length() - 10) << "\n";
        } else {
            cout << "Entrada: " << test.input << "\n";
        }
        
        cout << "Resultado esperado: " << (test.expectedResult ? "ACEPTADA" : "RECHAZADA") << "\n";

        auto start = chrono::high_resolution_clock::now();
        bool result = dfa.run(test.input);
        auto end = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::microseconds>(end - start);

        bool passed = (result == test.expectedResult);
        
        cout << "Resultado obtenido: " << (result ? "ACEPTADA" : "RECHAZADA") << "\n";
        cout << "Tiempo: " << duration.count() << " μs\n";
        cout << "Estado: " << (passed ? "✅ CORRECTO" : "❌ INCORRECTO") << "\n\n";

        // Guardar resultado
        results.push_back({
            test.name,
            passed,
            test.expectedResult,
            result,
            duration.count(),
            test.description
        });
    }
    
    return results;
}

// Función para imprimir resumen de resultados (VERSIÓN ULTRA-COMPACTA)
void printTestSummary(const vector<TestResult>& results) {
    cout << "=== RESULTADOS ===\n";
    
    int totalTests = results.size();
    int passedTests = 0;
    int failedTests = 0;
    
    for (const auto& result : results) {
        if (result.passed) {
            passedTests++;
        } else {
            failedTests++;
            cout << "❌ " << result.testName << "\n";
        }
    }
    
    cout << "✅ " << passedTests << "/" << totalTests << " tests exitosos\n";
    if (failedTests > 0) {
        cout << "❌ " << failedTests << " tests fallidos\n";
    }
    cout << "\n";
}

// Función para ejecutar todos los tests secuenciales (versión original para compatibilidad)
void runSequentialTests(const OptimizedDFA& dfa, const vector<DFATest>& tests) {
    auto results = runSequentialTestsWithResults(dfa, tests);
    printTestSummary(results);
}

// Función para crear tests específicos
vector<DFATest> createBinaryTests() {
    vector<DFATest> tests;
    
    // Test 1: Cadena que debe ser aceptada (termina en estado 2)
    tests.push_back({
        "Test Binario - Cadena Aceptada",
        3, 0, {2}, "01",
        "1101", // 0->0->0->1->2 (aceptada)
        true,
        "Cadena que termina en estado final 2"
    });

    // Test 2: Cadena que debe ser rechazada (termina en estado 0)
    tests.push_back({
        "Test Binario - Cadena Rechazada",
        3, 0, {2}, "01",
        "00000000", // 0->1->1->1->1->1->1->1->1 (rechazada)
        false,
        "Cadena que termina en estado no final 1"
    });

    // Test 3: Cadena vacía (debe terminar en estado inicial)
    tests.push_back({
        "Test Binario - Cadena Vacía",
        3, 0, {2}, "01",
        "", // Estado inicial 0 (rechazada)
        false,
        "Cadena vacía termina en estado inicial 0"
    });

    // Test 4: Cadena con patrón específico
    tests.push_back({
        "Test Binario - Patrón Específico",
        3, 0, {2}, "01",
        "1010101", // 0->0->1->2->1->2->1->2 (aceptada)
        true,
        "Patrón alternado que termina en estado final"
    });

    return tests;
}

// ✅ NUEVOS TESTS COMPLEJOS Y EDGE CASES
vector<DFATest> createComplexTests() {
    vector<DFATest> tests;
    
    // Test 1: Cadena extremadamente larga (puede causar overflow de memoria)
    tests.push_back({
        "Test Complejo - Cadena Extremadamente Larga",
        3, 0, {2}, "01",
        string(1000000, '1'), // 1 millón de '1's
        false,
        "Cadena de 1 millón de símbolos (debe terminar en estado 0)"
    });

    // Test 2: Cadena con símbolos inválidos (fuera del alfabeto)
    tests.push_back({
        "Test Complejo - Símbolos Inválidos",
        3, 0, {2}, "01",
        "01a2b3c", // Contiene símbolos no válidos
        false,
        "Cadena con símbolos fuera del alfabeto {0,1}"
    });

    // Test 3: Cadena con caracteres especiales y espacios
    tests.push_back({
        "Test Complejo - Caracteres Especiales",
        3, 0, {2}, "01",
        "0 1\t0\n1\r0", // Con espacios, tabs, newlines
        false,
        "Cadena con caracteres de control y espacios"
    });

    // Test 4: Cadena con caracteres Unicode
    tests.push_back({
        "Test Complejo - Caracteres Unicode",
        3, 0, {2}, "01",
        "01αβγδε", // Contiene caracteres Unicode
        false,
        "Cadena con caracteres Unicode fuera del alfabeto"
    });

    return tests;
}

// ✅ TESTS DE STRESS CON CADENAS MUY LARGAS
vector<DFATest> createStressTests() {
    vector<DFATest> tests;
    
    // Test 1: Cadena de 100K caracteres
    tests.push_back({
        "Test Stress - 100K Caracteres",
        3, 0, {2}, "01",
        string(100000, '0'), // 100K ceros
        false,
        "Cadena de 100,000 caracteres (debe terminar en estado 1)"
    });

    // Test 2: Cadena de 500K caracteres
    tests.push_back({
        "Test Stress - 500K Caracteres",
        3, 0, {2}, "01",
        string(500000, '1'), // 500K unos
        false,
        "Cadena de 500,000 caracteres (debe terminar en estado 0)"
    });

    // Test 3: Cadena de 1M caracteres
    tests.push_back({
        "Test Stress - 1M Caracteres",
        3, 0, {2}, "01",
        string(1000000, '0'), // 1M ceros
        false,
        "Cadena de 1,000,000 caracteres (debe terminar en estado 1)"
    });

    // Test 4: Patrón alternado muy largo
    string longPattern;
    for (int i = 0; i < 100000; i++) {
        longPattern += (i % 2 == 0) ? '1' : '0';
    }
    tests.push_back({
        "Test Stress - Patrón Alternado 100K",
        3, 0, {2}, "01",
        longPattern,
        true,
        "Patrón alternado de 100,000 caracteres (debe terminar en estado 2)"
    });

    return tests;
}

// ✅ TESTS PARA AUTÓMATA ABC
vector<DFATest> createABCTests() {
    vector<DFATest> tests;
    
    // Test 1: Cadena que debe ser aceptada (termina en estado 3)
    tests.push_back({
        "Test ABC - Cadena Aceptada",
        4, 0, {3}, "abc",
        "abc", // 0->1->2->3 (aceptada)
        true,
        "Cadena que termina en estado final 3"
    });

    // Test 2: Cadena que debe ser rechazada (termina en estado 0)
    tests.push_back({
        "Test ABC - Cadena Rechazada",
        4, 0, {3}, "abc",
        "ababab", // 0->1->2->1->2->1->2 (rechazada)
        false,
        "Cadena que termina en estado no final 2"
    });

    // Test 3: Cadena vacía
    tests.push_back({
        "Test ABC - Cadena Vacía",
        4, 0, {3}, "abc",
        "", // Estado inicial 0 (rechazada)
        false,
        "Cadena vacía termina en estado inicial 0"
    });

    // Test 4: Cadena con símbolos inválidos
    tests.push_back({
        "Test ABC - Símbolos Inválidos",
        4, 0, {3}, "abc",
        "abc123", // Contiene símbolos no válidos
        false,
        "Cadena con símbolos fuera del alfabeto {a,b,c}"
    });

    return tests;
}

// ✅ TESTS PARA AUTÓMATA CON MÚLTIPLES ESTADOS FINALES
vector<DFATest> createMultiFinalTests() {
    vector<DFATest> tests;
    
    // Test 1: Cadena que termina en estado final 2
    tests.push_back({
        "Test Multi-Final - Estado 2",
        4, 0, {2, 3}, "pq",
        "pq", // 0->1->2 (aceptada)
        true,
        "Cadena que termina en estado final 2"
    });

    // Test 2: Cadena que termina en estado final 3
    tests.push_back({
        "Test Multi-Final - Estado 3",
        4, 0, {2, 3}, "pq",
        "pqq", // 0->1->2->3 (aceptada)
        true,
        "Cadena que termina en estado final 3"
    });

    // Test 3: Cadena que termina en estado no final
    tests.push_back({
        "Test Multi-Final - Estado No Final",
        4, 0, {2, 3}, "pq",
        "ppp", // 0->1->1->1 (rechazada)
        false,
        "Cadena que termina en estado no final 1"
    });

    // Test 4: Cadena vacía
    tests.push_back({
        "Test Multi-Final - Cadena Vacía",
        4, 0, {2, 3}, "pq",
        "", // Estado inicial 0 (rechazada)
        false,
        "Cadena vacía termina en estado inicial 0"
    });

    return tests;
}

// ✅ FUNCIÓN PARA ANÁLISIS DETALLADO
void runDetailedAnalysis(const OptimizedDFA& dfa, const string& input, const string& testName) {
    cout << "\n=== ANÁLISIS DETALLADO: " << testName << " ===\n";
    cout << "Entrada: [" << input << "]\n";
    cout << "Tamaño: " << input.size() << " caracteres\n";
    
    auto start = chrono::high_resolution_clock::now();
    auto trackingResult = dfa.run_with_tracking(input);
    bool result = trackingResult.first;
    vector<int> stateSequence = trackingResult.second;
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(end - start);
    
    cout << "Tiempo de ejecución: " << duration.count() << " μs\n";
    cout << "Resultado: " << (result ? "ACEPTADA" : "RECHAZADA") << "\n";
    
    // Mostrar ruta de estados real del tracking
    cout << "Ruta de estados: ";
    for (size_t i = 0; i < stateSequence.size(); i++) {
        cout << stateSequence[i];
        if (i < stateSequence.size() - 1) cout << " -> ";
    }
    cout << "\n";
    
    if (!stateSequence.empty()) {
        int finalState = stateSequence.back();
        cout << "Estado final: " << finalState << "\n";
        cout << "Es estado de aceptación: " << (dfa.isAccepting(finalState) ? "SÍ" : "NO") << "\n";
    }
} 