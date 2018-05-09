//
// Created by Дмитрий Бутилов on 06.05.18.
//
#include <mpi.h>
#include <iostream>

int rank, size;

void compareSendReceiveAndBcast(int repeats) {
    double sendReceiveTime = 0, bcastTime = 0;
    MPI_Status mpiStatus;
    if (rank == 0) {
        double startTime = MPI_Wtime();
        for (int i = 0; i < repeats; ++i) {
            for (int j = 1; j < size; ++j) {
                MPI_Send(nullptr, 0, MPI_BYTE, j, 1, MPI_COMM_WORLD);
            }
        }
        sendReceiveTime += MPI_Wtime() - startTime;
    } else {
        for (int i = 0; i < repeats; ++i) {
            MPI_Recv(nullptr, 0, MPI_BYTE, 0, 1, MPI_COMM_WORLD, &mpiStatus);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    double startTime = 0;
    if (rank == 0) {
        startTime = MPI_Wtime();
    }
    for (int i = 0; i < repeats; ++i) {
        MPI_Bcast(nullptr, 0, MPI_BYTE, 0, MPI_COMM_WORLD);
    }
    if (rank == 0) {
        bcastTime += MPI_Wtime() - startTime;
    }
    std::cout << "Repeats: " << repeats << "\nSend/Recieve: " << sendReceiveTime << "\nBcast: " << bcastTime << "\n";
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double time = 0;
    int repeats = 100000;
    if (rank == 0) {
        std::cout << "Start\n";
    }

    time = MPI_Wtime();
    compareSendReceiveAndBcast(repeats);
    time = MPI_Wtime() - time;

    if (rank == 0) {
        std::cout << "End\nTime: " << time << " s\n";
    }

    MPI_Finalize();
    return 0;
}