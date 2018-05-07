//
// Created by Дмитрий Бутилов on 06.05.18.
//
#include <mpi.h>
#include <cstdio>
#include <iostream>

int currentRank, size;

int Reduce(int repeats) {
    MPI_Status status;
    double timeSR = 0, timeReduce = 0;
    if (currentRank == 0) {
        char receive;
        int sending;
        double startTime = MPI_Wtime();
        for (sending = 0; sending < repeats; ++sending) {
            int destination;
            for (destination = 1; destination < size; ++destination) {
                MPI_Recv(&receive, 1, MPI_BYTE, destination, 1, MPI_COMM_WORLD, &status);
            }
        }
        timeSR += MPI_Wtime() - startTime;
    } else {
        int sending;
        for (sending = 0; sending < repeats; ++sending) {
            MPI_Send(&currentRank, 1, MPI_BYTE, 0, 1, MPI_COMM_WORLD);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    double startTime = 0;
    if (currentRank == 0) {
        startTime = MPI_Wtime();
    }
    char receive = 0;
    int sending;
    for (sending = 0; sending < repeats; ++sending) {
        MPI_Reduce(&currentRank, &receive, 1, MPI_BYTE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
    if (currentRank == 0) {
        timeReduce += MPI_Wtime() - startTime;
    }
    std::cout << "Reduce: " << timeReduce << "\nSend_Recv: " << timeSR << "\nRepeats: " << repeats << "\n";
    return 0;
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &currentRank);

    double time = 0;

    int ReduceRepeats = 50;
    if (currentRank == 0) {
        std::cout << "Start\n";
    }
    time = MPI_Wtime();
    Reduce(ReduceRepeats);
    time = MPI_Wtime() - time;
    if (currentRank == 0) {
        std::cout << "End\nTime: " << time << " seconds\n";
    }

    MPI_Finalize();
    return 0;
}