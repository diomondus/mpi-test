//
// Created by Дмитрий Бутилов on 06.05.18.
//
#include <mpi.h>
#include <cstdio>

int currentRank, size;

void Latency(int repeats) {
    MPI_Status state;
    double time = 0;
    for (int t = 0; t < repeats; t++) {
        if (currentRank == 0) {
            for (int i = 1; i < size; i++) {
                double startTime = MPI_Wtime();
                MPI_Send(nullptr, 0, MPI_BYTE, i, 1, MPI_COMM_WORLD);
                MPI_Recv(nullptr, 0, MPI_BYTE, i, 2, MPI_COMM_WORLD, &state);
                time += MPI_Wtime() - startTime;
            }
        } else {
            MPI_Recv(nullptr, 0, MPI_BYTE, 0, 1, MPI_COMM_WORLD, &state);
            MPI_Send(nullptr, 0, MPI_BYTE, 0, 2, MPI_COMM_WORLD);
        }
    }
    time = time / (2.0 * repeats);
    if (currentRank == 0)
        printf("Latency: %2.8f s; Repeats: %d \n", time, repeats);
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &currentRank);

    double time = 0;

    int latencyRepeats = 100000;
    if (currentRank == 0) {
        printf("Start\n");
    }
    time = MPI_Wtime();
    Latency(latencyRepeats);
    time = MPI_Wtime() - time;
    if (currentRank == 0) {
        printf("Ends in %f seconds\n", time);
    }

    MPI_Finalize();
    return 0;
}