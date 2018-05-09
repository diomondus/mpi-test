#include <iostream>
#include "mpi.h"

int size, rank;

int main(int argc, char **argv) {

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::cout << "Hello, World!" << std::endl;

    MPI_Finalize();
    return 0;
}