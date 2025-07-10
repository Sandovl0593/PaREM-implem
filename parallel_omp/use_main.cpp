#include "dfa.h"
#include <cmath>
#include <omp.h>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include <random>
#include <fstream>

const int NUM_ITERATIONS = 50; // Promedio de 50 iteraciones

// Función para generar cadenas aleatorias
std::string generateRandomString(const std::string& alphabet, int length) {
    std::string result;
    result.reserve(length);

    std::mt19937 rng(static_cast<unsigned int>(std::chrono::system_clock::now().time_since_epoch().count()));
    std::uniform_int_distribution<int> dist(0, alphabet.length() - 1);

    for (int i = 0; i < length; ++i) {
        result.push_back(alphabet[dist(rng)]);
    }
    return result;
} 

// ✅ ESTRUCTURA PARA MEDICIÓN DETALLADA DE TIEMPOS
struct DetailedTimeAnalysis {
    double computation_time;        // (|Q|·Σ+n)/p
    double total_parallel;         // Tiempo total paralelo
    
    // Parámetros teóricos
    int Q;  // Número de estados
    int sigma; // Tamaño del alfabeto  
    int n;  // Tamaño de entrada
    int threads;  // Número de threads
    
    DetailedTimeAnalysis(int states, int alphabet_size, int input_size, int threads) 
        : computation_time(0), total_parallel(0), Q(states), sigma(alphabet_size), 
          n(input_size), threads(threads) {}
};

// ✅ FUNCIÓN: Análisis teórico de escalabilidad
void analyzeScalability(const DetailedTimeAnalysis& A, bool isOptim, double sequentialTime, int num_threads) {
    int p = num_threads;
    double Qsigma = double(A.Q) * double(A.sigma);
    double n  = double(A.n);

    // --- 1) Modelo teórico ---
    double T_theo;
    if (!isOptim) {
        // Modelo simple: T = (QΣ + n)/p   +   n
        T_theo = (Qsigma + n)/p + n;
    }
    else {
        // Modelo Brent optimizado:
        //    T = (QΣ + n)/p + log2(p)   +   Q·log2(p)
        T_theo = (Qsigma + n)/p + log2(double(p)) + double(A.Q) * log2(double(p));
    }

    // Métricas reales ---
    double real_speedup    = sequentialTime / A.total_parallel;
    double real_efficiency = real_speedup / p * 100.0;

    // Métricas teóricas ---
    double theo_speedup    = sequentialTime / T_theo;
    double theo_efficiency = theo_speedup / p * 100.0;

    // string file_threads_name = to_string(p) + ".txt";
    // fstream cout;
    // cout.open("resultados_omp/log_" + file_threads_name, ios::app);

    // Salida por pantalla ---
    cout << fixed << setprecision(5);
    cout << "\n📊 ANÁLISIS DE ESCALABILIDAD "<< (isOptim ? "(OMP optimizado)" : "(OMP estimado)") << ":\n";
    cout << "============================================================\n";
    cout << "Parámetros: |Q|=" << A.Q << ", Σ=" << A.sigma << ", n=" << A.n << ", p=" << p << "\n\n";
    cout << "🔬 Componentes teóricas (unidades de ops):\n";
    if (!isOptim) {
        cout << "  • (|Q|·Σ + n)/p: " << (Qsigma + n)/p << "\n" << "  • + n:          " << n << "\n";
    } else {
        cout << "  • (|Q|·Σ + n)/p: " << (Qsigma + n)/p << "\n" << "  • + log2(p):     " << log2(double(p)) << "\n" << "  • + |Q|·log2(p): " << A.Q * log2(double(p)) << "\n";
    }
    cout << "  → T_theoretical: " << T_theo << "\n\n";
    cout << "⏱️ Tiempos:\n";
    cout << "  • Secuencial: " << sequentialTime << " ms\n";
    cout << "📈 Métricas:\n";
    cout << "  • Real     → Speedup: " << real_speedup << " | Eficiencia: " << real_efficiency << "%\n";
    cout << "  • Teórico  → Speedup: " << theo_speedup << " | Eficiencia: " << theo_efficiency << "%\n\n";
    cout << "🔍 Conclusión fuerte: "<< (real_efficiency >= 80.0 ? "✅ BUENA escalabilidad" 
                                    : "❌ Escalabilidad pobre")<< " (" << real_efficiency << "%)\n";
}


int main(int argc, char* argv[]) {
    if (argc != 2) {
        cout << "Uso: " << argv[0] << " <num_threads>" << endl;
        cout << "Ejemplo: " << argv[0] << " 4" << endl;
        return 1;
    }
    int num_threads = atoi(argv[1]);
    if (num_threads < 1 || num_threads > 8) {
        cout << "Error: Número de hilos debe estar entre 1 y 8" << endl;
        return 1;
    }
    
    omp_set_num_threads(num_threads);

    // ================= TESTS CON DIFERENTES CONFIGURACIONES =================
    vector<tuple<int, int, vector<int>, string, string>> dfaConfigs = {
        {3, 0, {2}, "01", "DFA Binario (|Q|=3, Σ=2)"},
        {4, 0, {3}, "abc", "DFA ABC (|Q|=4, Σ=3)"},
        {5, 0, {2,4}, "xyz", "DFA Multi-Final (|Q|=5, Σ=3)"},
        {6, 0, {5}, "0123", "DFA Cuaternario (|Q|=6, Σ=4)"}
    };
    
    vector<int> inputSizes = {1000, 10000, 100000, 1000000};
    
    for (const auto& config : dfaConfigs) {
        int states = get<0>(config);
        int initial = get<1>(config);
        vector<int> accepting = get<2>(config);
        string alphabet = get<3>(config);
        string description = get<4>(config);
        
        // Crear DFA
        DFA dfa(states, initial, accepting, alphabet);
        
        // Añadir transiciones de ejemplo
        for (int s = 0; s < states; s++) {
            for (char c : alphabet) {
                int nextState = (s + (c - alphabet[0]) + 1) % states;
                dfa.addTransition(s, c, nextState);
            }
        }
        
        cout << "\n🔬 TESTING: " << description << endl;
        cout << string(60, '=') << endl;
        
        for (int inputSize : inputSizes) {
            cout << "\n📏 Tamaño de entrada: " << inputSize << " caracteres" << endl;
            
            // Generar entrada de prueba
            string testInput = generateRandomString(alphabet, inputSize);

            // ================= MEDICIÓN SECUENCIAL =================
            bool seqResult = false;
            double sequentialTime = 0;
            vector<double> seqTimes;
            for (int iter = 0; iter < NUM_ITERATIONS; iter++) {
                auto start = chrono::high_resolution_clock::now();
                seqResult = dfa.run(testInput);
                auto end = chrono::high_resolution_clock::now();
                sequentialTime = chrono::duration_cast<chrono::microseconds>(end - start).count();
                seqTimes.push_back(sequentialTime);
            }
            double avgSeqTime = accumulate(seqTimes.begin(), seqTimes.end(), 0.0) / NUM_ITERATIONS;
            cout << "  📊 Tiempo secuencial promedio: " << fixed << setprecision(5) << avgSeqTime << " μs" << endl;



            // ================= MEDICIÓN PARALELA OMP =================
            DFAResult parallelOMPResult;
            double parallelOMPTime = 0;
            vector<double> parallelTimes;
            for (int iter = 0; iter < NUM_ITERATIONS; iter++) {
                parallelOMPResult = processOMPParallel(dfa, testInput, num_threads);
                parallelTimes.push_back(parallelOMPResult.computationTime);
            }
            double parallelAvgTime = accumulate(parallelTimes.begin(), parallelTimes.end(), 0.0) / NUM_ITERATIONS;
            cout << "  📊 Tiempo OMP promedio: " << fixed << setprecision(5) << parallelAvgTime << " μs" << endl;

            DetailedTimeAnalysis resultOMP(states, alphabet.size(), inputSize, num_threads);
            analyzeScalability(resultOMP, false, sequentialTime, num_threads);



            // ================= MEDICIÓN OMP OPTIMIZADO =================
            DFAOptimizedResult parallelOMPResultOptimized;
            double parallelOMPTimeOptimized = 0;
            vector<double> parallelOptTimes;
            for (int iter = 0; iter < NUM_ITERATIONS; iter++) {
                parallelOMPResultOptimized = runOptimizedOMPParallel(dfa, testInput, num_threads);
                parallelOptTimes.push_back(parallelOMPResultOptimized.computationTime);
            }
            double parallelOptAvgTime = accumulate(parallelOptTimes.begin(), parallelOptTimes.end(), 0.0) / NUM_ITERATIONS;
            cout << "  📊 Tiempo OMP optimizado promedio: " << fixed << setprecision(5) << parallelOptAvgTime << " μs" << endl;

            DetailedTimeAnalysis resultOMPOptimized(states, alphabet.size(), inputSize, num_threads);
            analyzeScalability(resultOMPOptimized, true, sequentialTime, num_threads);

            string file_threads_name = to_string(num_threads) + ".txt";
            fstream resultsFile, scalabilityFile;
            resultsFile.open("resultados_omp/results_" + file_threads_name, ios::app);
            scalabilityFile.open("resultados_omp/scalability_" + file_threads_name, ios::app);

            // ================= GUARDAR EN TODOS LOS ARCHIVOS =================
            if (resultsFile.is_open()) {
                resultsFile << "\n" << description << " - Entrada: " << inputSize << " caracteres\n";
                resultsFile << "Secuencial: " << fixed << setprecision(5) << sequentialTime << " μs\n";
                resultsFile << "OMP Paralelo: " << fixed << setprecision(5) << parallelAvgTime << " μs (Speedup: " 
                           << setprecision(5) << sequentialTime/parallelAvgTime << "x )\n";
                resultsFile << "OMP Optimizado: " << fixed << setprecision(5) 
                           << parallelOptAvgTime << " μs (Speedup: " 
                           << sequentialTime/parallelOptAvgTime << "x)\n";
            }
            
            // Archivo de escalabilidad detallada
            if (scalabilityFile.is_open()) {
                scalabilityFile << "\n" << description << " - " << inputSize << " caracteres:\n";
                scalabilityFile << "Parámetros: |Q|=" << states << ", Σ=" << alphabet.size() << ", n=" << inputSize << ", threads=" << num_threads << "\n";
                scalabilityFile << "OMP - Speedup: " << setprecision(5) << sequentialTime/parallelAvgTime
                               << "x, Eficiencia: " << setprecision(5) << (sequentialTime/parallelAvgTime/num_threads) << "\n";
                scalabilityFile << "OMP Optimizada - Speedup: " << sequentialTime/parallelOptAvgTime
                               << "x, Eficiencia: " << (sequentialTime/parallelOptAvgTime/num_threads) << "\n";
            }
            // Cerrar archivos
            if (resultsFile.is_open()) {
                resultsFile.close();
            }
    
            if (scalabilityFile.is_open()) {
                scalabilityFile.close();
            }
        }
    }
}