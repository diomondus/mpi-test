//
// Created by Дмитрий Бутилов on 06.05.18.
//
#include <mpi.h>
#include <cstdio>

int currentRank, size;

void Bcast(int repeats) {
    double timeBcast = 0, timeSR = 0;
    MPI_Status status;
    if (currentRank == 0) {
        double startTime = MPI_Wtime();
        int sending;
        int destination;
        for (sending = 0; sending < repeats; ++sending) {
            for (destination = 1; destination < size; ++destination) {
                MPI_Send(nullptr, 0, MPI_BYTE, destination, 1, MPI_COMM_WORLD);
            }
        }
        timeSR += MPI_Wtime() - startTime;
    } else {
        int sending;
        for (sending = 0; sending < repeats; ++sending) {
            MPI_Recv(nullptr, 0, MPI_BYTE, 0, 1, MPI_COMM_WORLD, &status);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    double startTime = 0;
    if (currentRank == 0) {
        startTime = MPI_Wtime();
    }
    int sending;
    for (sending = 0; sending < repeats; ++sending) {
        MPI_Bcast(nullptr, 0, MPI_BYTE, 0, MPI_COMM_WORLD);
    }
    if (currentRank == 0) {
        timeBcast += MPI_Wtime() - startTime;
    }
    printf("Bcast: %f Send_Recv: %f Repeats: %d \n", timeBcast, timeSR, repeats);
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &currentRank);

    double time = 0;

    int BCastReps = 100000;
    if (currentRank == 0) {
        printf("Start\n");
    }
    time = MPI_Wtime();
    Bcast(BCastReps);
    time = MPI_Wtime() - time;
    if (currentRank == 0) {
        printf("End in %f seconds\n", time);
    }

    MPI_Finalize();
    return 0;
}