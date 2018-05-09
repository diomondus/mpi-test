//
// Created by Дмитрий Бутилов on 06.05.18.
//
#include <mpi.h>
#include <iostream>

int rank, size;

void Bcast(int repeats) {
    double timeBcast = 0, timeSR = 0;
    MPI_Status status;
    if (rank == 0) {
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
    if (rank == 0) {
        startTime = MPI_Wtime();
    }
    int sending;
    for (sending = 0; sending < repeats; ++sending) {
        MPI_Bcast(nullptr, 0, MPI_BYTE, 0, MPI_COMM_WORLD);
    }
    if (rank == 0) {
        timeBcast += MPI_Wtime() - startTime;
    }
    std::cout << "Bcast: " << timeBcast << " Send_Recv: " << timeSR << " Repeats: " << repeats << "\n";
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double time = 0;

    int BCastReps = 100000;
    if (rank == 0) {
        std::cout << "Start\n";
    }
    time = MPI_Wtime();
    Bcast(BCastReps);
    time = MPI_Wtime() - time;
    if (rank == 0) {
        std::cout << "End\nTime: " << time << " seconds\n";
    }

    MPI_Finalize();
    return 0;
}