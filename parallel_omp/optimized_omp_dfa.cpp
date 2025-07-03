#include "dfa.h"

// ✅ FUNCIÓN PARA PRECOMPUTAR LOS SALTOS DE LOS BLOQUES
std::vector<std::vector<BlockInfo>> precomputeJump(const DFA& afd, const std::string &T, int num_threads, int chunk) {
    int Q = afd.getNumStates();
    std::vector<std::vector<BlockInfo>> jumpTable(num_threads, std::vector<BlockInfo>(Q));

    #pragma omp parallel for schedule(static)
    for (int t = 0; t < num_threads; ++t) {
        int left  = t * chunk;
        int right = std::min((int)T.size(), left + chunk);
        for (int q = 0; q < Q; ++q) {
            int curr = q;
            int cnt  = afd.isAccepting(curr) ? 1 : 0;
            for (int i = left; i < right; ++i) {
                curr = afd.nextState(curr, T[i]);
                if (afd.isAccepting(curr)) cnt++;
            }
            jumpTable[t][q] = { curr, cnt };
        }
    }
    return jumpTable;
}

// ✅ FUNCIÓN PARA CALCULAR ESTADOS INICIALES DE CADA BLOQUE
// Basada en la tabla de saltos precomputada
std::vector<int> computeStartStates(const std::vector<std::vector<BlockInfo>> &jumpTable) {
    int p = jumpTable.size();
    std::vector<int> startState(p);
    startState[0] = 0;
    for (int t = 1; t < p; ++t) {
        int prev = startState[t-1];
        startState[t] = jumpTable[t-1][prev].endState;
    }
    return startState;
}

// ✅ FUNCIÓN PARA EJECUTAR EL DFA EN PARALLEL CON OPENMP
DFAOptimizedResult runOptimizedOMPParallel(const DFA& afd, const std::string &T, int num_threads) {
    int N = T.size();
    int chunk = (N + num_threads - 1) / num_threads;

    DFAOptimizedResult result;
    double t_start = omp_get_wtime();
    auto jumpTable   = precomputeJump(afd, T, num_threads, chunk);
    double t_jump = omp_get_wtime() - t_start;
    auto startStates = computeStartStates(jumpTable);
    result.processStates.resize(N + 1);

    int totalMatches = 0;

    #pragma omp parallel for reduction(+:totalMatches) schedule(static)
    for (int t = 0; t < num_threads; ++t) {
        int left  = t * chunk;
        int right = std::min(N, left + chunk);
        if (left >= right) continue;

        int curr = startStates[t];
        int localMatches = afd.isAccepting(curr) ? 1 : 0;
        int writePos = left;

        result.processStates[writePos++] = curr;
        
        for (int i = left; i < right; ++i) {
            curr = afd.nextState(curr, T[i]);
            result.processStates[writePos++] = curr;
            if (afd.isAccepting(curr)) localMatches++;
        }

        totalMatches += localMatches;
    }
    double t_end = omp_get_wtime();

    // Rellenar campos finales
    result.computationTime = (t_end - t_start) * 1000.0;
    result.preComputeJumpTime = t_jump * 1000.0;
    result.numMatches      = totalMatches;
    if (!result.processStates.empty()) {
        int finalState = result.processStates.back();
        result.finalState = finalState;
        result.accepted = afd.isAccepting(finalState);
    }
    return result;
}
