#include "AFD.hpp"
#include <random>
#include <omp.h>
#include <string>

vector<int> run_parem(AFD afd, const string &T, int num_threads) {
    vector<int> I;
    int size_T = T.size();
    I.reserve(size_T + 1);

    vector<vector<int>> L(num_threads, vector<int>(afd.numStates));
    vector<vector<int>> R(num_threads, vector<int>(afd.numStates, 0));

    omp_set_num_threads(num_threads);
    int sub_length = size_T / num_threads;
    int size = omp_get_num_threads();

    #pragma omp parallel
    {
        int th_id = omp_get_thread_num();
        int pos = th_id * sub_length;
    }
}