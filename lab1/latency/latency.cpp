//
// Created by Дмитрий Бутилов on 06.05.18.
//
#include <mpi.h>
#include <iostream>

int rank, size;

void Latency(int latencyRepeats) {
    MPI_Status mpiStatus;
    double time = 0;
    for (int i = 0; i < latencyRepeats; i++) {
        if (rank == 0) {
            for (int j = 1; j < size; j++) {
                double startTime = MPI_Wtime();
                MPI_Send(nullptr, 0, MPI_BYTE, j, 1, MPI_COMM_WORLD);
                MPI_Recv(nullptr, 0, MPI_BYTE, j, 2, MPI_COMM_WORLD, &mpiStatus);
                time += MPI_Wtime() - startTime;
            }
        } else {
            MPI_Recv(nullptr, 0, MPI_BYTE, 0, 1, MPI_COMM_WORLD, &mpiStatus);
            MPI_Send(nullptr, 0, MPI_BYTE, 0, 2, MPI_COMM_WORLD);
        }
    }
    double latencyTime = time / (2.0 * latencyRepeats); // T/2N
    if (rank == 0) {
        std::cout << "Latency: " << latencyTime << " s\nRepeats: " << latencyRepeats << "\n";
    }
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int latencyRepeats = 10000;
    if (rank == 0) {
        std::cout << "Start\n";
    }

    double completeTime = MPI_Wtime();
    Latency(latencyRepeats);
    completeTime = MPI_Wtime() - completeTime;

    if (rank == 0) {
        std::cout << "End\nTime: " << completeTime << " seconds\n";
    }
    MPI_Finalize();
    return 0;
}