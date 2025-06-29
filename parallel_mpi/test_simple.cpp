#include <mpi.h>
#include <iostream>
#include <vector>
#include <string>

using namespace std;

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    if (rank == 0) {
        cout << "Proceso " << rank << " de " << size << " iniciado" << endl;
    }
    
    // Test básico de división
    double a = 1.0;
    double b = 0.0;
    
    cout << "Rank " << rank << ": a=" << a << ", b=" << b << endl;
    
    // Esto debería causar el error si es división por cero
    if (b != 0.0) {
        double result = a / b;
        cout << "Rank " << rank << ": resultado=" << result << endl;
    } else {
        cout << "Rank " << rank << ": evitando división por cero" << endl;
    }
    
    // Test de vector
    vector<int> test_vector(10);
    for (int i = 0; i < 10; i++) {
        test_vector[i] = i;
    }
    
    cout << "Rank " << rank << ": vector creado correctamente" << endl;
    
    // Test de string
    string test_string = "0100100100";
    cout << "Rank " << rank << ": string=" << test_string << ", length=" << test_string.length() << endl;
    
    // Test de cálculos básicos
    int inputSize = test_string.length();
    int numProcs = size;
    
    cout << "Rank " << rank << ": inputSize=" << inputSize << ", numProcs=" << numProcs << endl;
    
    if (numProcs > 0) {
        int blockSize = inputSize / numProcs;
        int remainder = inputSize % numProcs;
        cout << "Rank " << rank << ": blockSize=" << blockSize << ", remainder=" << remainder << endl;
    }
    
    MPI_Finalize();
    return 0;
} 