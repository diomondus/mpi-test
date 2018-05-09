//
// Created by Дмитрий Бутилов on 06.05.18.
//
#include <mpi.h>
#include <iostream>

int rank, size;

void compareSendReceiveAndReduce(int repeats) {
    MPI_Status mpiStatus;
    double sendReceiveTime = 0, reduceTime = 0;
    if (rank == 0) {
        char receive;
        double startTime = MPI_Wtime();
        for (int i = 0; i < repeats; ++i) {
            for (int j = 1; j < size; ++j) {
                MPI_Recv(&receive, 1, MPI_BYTE, j, 1, MPI_COMM_WORLD, &mpiStatus);
            }
        }
        sendReceiveTime += MPI_Wtime() - startTime;
    } else {
        for (int i = 0; i < repeats; ++i) {
            MPI_Send(&rank, 1, MPI_BYTE, 0, 1, MPI_COMM_WORLD);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double startTime = 0;
    if (rank == 0) {
        startTime = MPI_Wtime();
    }
    char receive = 0;
    for (int i = 0; i < repeats; ++i) {
        MPI_Reduce(&rank, &receive, 1, MPI_BYTE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
    if (rank == 0) {
        reduceTime += MPI_Wtime() - startTime;
    }
    std::cout << "Repeats: " << repeats << "\nSend/Recieve: " << sendReceiveTime << "\nReduce: " << reduceTime << "\n";
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double time = 0;
    int repeats = 50;
    if (rank == 0) {
        std::cout << "Start\n";
    }

    time = MPI_Wtime();
    compareSendReceiveAndReduce(repeats);
    time = MPI_Wtime() - time;

    if (rank == 0) {
        std::cout << "End\nTime: " << time << " seconds\n";
    }

    MPI_Finalize();
    return 0;
}