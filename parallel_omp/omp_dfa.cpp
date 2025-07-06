#include "dfa.h"

// ✅ FUNCIÓN PARA PRECOMPUTAR ESTADOS FINALES Y MATCHES
// CONSIDERANDO EL PRECOMPUTO DEL AFD -> TIENE LA MISMA COMPLEJIDAD DEL ENFOQUE POR BRENT
DFAResult processOMPParallel(const DFA& dfa, const std::string &T, int num_threads) {
    int N = T.size();
    int chunk = (N + num_threads - 1) / num_threads;
    std::vector<int> fullRoute;
    int totalMatches = 0;

    std::vector<BlockResult> results(num_threads);
    
    double t_start = omp_get_wtime();
    #pragma omp parallel for schedule(static)
    for (int t = 0; t < num_threads; t++) {
        int left = t * chunk;
        int right = std::min(N, left + chunk);
        results[t].matches = 0;

        if (left >= N) {
            results[t].startState = results[t].endState = -1;
            continue;
        }
        if (t == 0) results[t].startState = 0;
        else        results[t].startState = -1;

        int curr = 0;
        std::vector<int> localRoute;
        localRoute.reserve(right - left + 1);
        localRoute.push_back(curr);
        if (dfa.isAccepting(curr)) results[t].matches++;

        for (int i = left; i < right; i++) {
            curr = dfa.nextState(curr, T[i]);
            localRoute.push_back(curr);
            if (dfa.isAccepting(curr)) results[t].matches++;
        }
        results[t].route    = std::move(localRoute);
        results[t].endState = curr;
    }
    double t_end = omp_get_wtime();

    // Reducción secuencial
    for (int t = 0; t < num_threads; t++) {
        if (t == 0) {
            fullRoute.push_back(results[t].startState);
        }
        for (size_t i = (t == 0); i < results[t].route.size(); ++i) {
            fullRoute.push_back(results[t].route[i]);
        }
        totalMatches += results[t].matches;

        if (t + 1 < num_threads) {
            results[t+1].startState = results[t].endState;
        }
    }

    // Construir y retornar el resultado
    DFAResult result;
    result.computationTime = (t_end - t_start) * 1000.0;
    result.processStates   = std::move(fullRoute);
    result.numMatches      = totalMatches;
    // El DFA acepta si el estado final es de aceptación
    if (!result.processStates.empty()) {
        int finalState = result.processStates.back();
        result.finalState = finalState;
        result.accepted = dfa.isAccepting(finalState);
    }
    return result;
}
