#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>

int size, rank;

double *matrix, *vector, *result;
double *matrixPart, *vectorPart, *resultPart;

int solvingStatus, partSizePerProcess; // статус решения задачи, строки для текущего процесса
int *mainRowIndexArray, *mainRowIteration; // ведущие строки для каждой итерации, номера итераций ведущих строк
int *sendcounts; // массив (размер=max_rank), содержащий число элементов, посылаемых каждому процессу
int *displs; // i-ое значение определяет смещение относительно начала sendbuf для данных, посылаемых процессу i

struct {
    double maxValue;
    int currentRank;
} localMax = {}, globalMax = {};

int matrixSize = 5000;

bool usePrint = true;

//----------------------------------------------------------------------------------------------------------------------
double getRandomDouble() {
    return rand() / 10000000.0;
}

void initMatrixAndVector() {
    if (rank == 0) //Заполняем
    {
        srand(static_cast<unsigned int>(time(0))); //псевдослучайные

        matrix = new double[matrixSize * matrixSize];
        vector = new double[matrixSize];
        result = new double[matrixSize];

        for (int i = 0; i < matrixSize * matrixSize; ++i) {
            matrix[i] = getRandomDouble();
            if (usePrint) {
                printf("\nA[%i]=%f", i, matrix[i]);
            }
        };
        for (int i = 0; i < matrixSize; i++) {
            vector[i] = getRandomDouble();
            if (usePrint) {
                printf("\nb[%i]=%f", i, vector[i]);
            }
            result[i] = 0;
        };
    }
}

void initParts() {
    int nonDistributeRowCount = matrixSize;
    for (int i = 0; i < rank; ++i)
        nonDistributeRowCount = nonDistributeRowCount - nonDistributeRowCount / (size - i);
    partSizePerProcess = nonDistributeRowCount / (size - rank);//Определение размера части данных,на конкретном процессе
    matrixPart = new double[partSizePerProcess * matrixSize];
    vectorPart = new double[partSizePerProcess]; // элементы столбца свободных членов
    resultPart = new double[partSizePerProcess];
    mainRowIndexArray = new int[matrixSize]; // массив индексов ведущих строк системы на каждой итерации
    // итерация, на которой соответствующая строка системы, расположенная на процессе, выбрана  ведущей
    mainRowIteration = new int[partSizePerProcess];
    displs = new int[size];
    sendcounts = new int[size];
    for (int i = 0; i < partSizePerProcess; i++) {
        mainRowIteration[i] = -1;
    }
}

void initData() {
    initMatrixAndVector();
    MPI_Bcast(&matrixSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    initParts();
}

//------------------------------------------------------------------------------
void disributeDataBetweenProcesses() { // Распределение исходных данных между процессами
    int *matrElem = new int[size];             // Индекс первого элемента матрицы, передаваемого процессу
    int *matrRang = new int[size];            // Число элементов матрицы, передаваемых процессу
    int previousRowCount = 0, nonDistributeRowCount, previousSize, previousIndex, portion;

    nonDistributeRowCount = matrixSize;
    for (int i = 0; i < size; i++) {  //Определяем, сколько элементов матрицы будет передано каждому процессу
        previousSize = (i == 0) ? 0 : matrRang[i - 1];
        previousIndex = (i == 0) ? 0 : matrElem[i - 1];
        portion = (i == 0) ? 0 : previousRowCount;             //число строк, отданных предыдущему процессу
        nonDistributeRowCount -= portion;
        previousRowCount = nonDistributeRowCount / (size - i);
        matrRang[i] = previousRowCount * matrixSize;
        matrElem[i] = previousIndex + previousSize;
    }

    //Рассылка матрицы
    MPI_Scatterv(matrix, matrRang, matrElem, MPI_DOUBLE, matrixPart, matrRang[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    nonDistributeRowCount = matrixSize;
    for (int i = 0; i < size; i++) {
        previousSize = (i == 0) ? 0 : sendcounts[i - 1];
        previousIndex = (i == 0) ? 0 : displs[i - 1];
        portion = (i == 0) ? 0 : sendcounts[i - 1];
        nonDistributeRowCount -= portion;
        sendcounts[i] = nonDistributeRowCount / (size - i);
        displs[i] = previousIndex + previousSize;
    }

    //Рассылка вектора
    MPI_Scatterv(vector, sendcounts, displs, MPI_DOUBLE, vectorPart, sendcounts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    solvingStatus = 1;
    MPI_Bcast(&solvingStatus, 1, MPI_INT, 0, MPI_COMM_WORLD);
    delete[] matrRang;
    delete[] matrElem;
}

//------------------------------------------------------------------------------

void multiplyByCoeffient(const double *matrix, int numIter) {
    double coeffient;
    for (int i = 0; i < partSizePerProcess; i++) {     // для каждой строки в процессе
        if (mainRowIteration[i] == -1) {
            coeffient = matrixPart[i * matrixSize + numIter] / matrix[numIter];
            for (int j = numIter; j < matrixSize; j++) {
                matrixPart[i * matrixSize + j] -= matrix[j] * coeffient;
            };
            vectorPart[i] -= matrix[matrixSize] * coeffient;

        };
    };
}

//------------------------------------------------------------------------------

void triangulate() {
    int mainRowInCurrentProcess = 0;   // индекс ведущей строки на конкретном процессе

    double *globalMatrix = new double[matrixSize + 1]; // т.е. строка матрицы + значение вектора
    for (int i = 0; i < matrixSize; i++) {
        // Вычисление ведущей строки
        double maxValue = 0;
        int index = -1;
        for (int j = 0; j < partSizePerProcess; j++) {
            index = j;
            if ((mainRowIteration[j] == -1) && (maxValue < fabs(matrixPart[i + matrixSize * j]))) {
                maxValue = fabs(matrixPart[i + matrixSize * j]);
                mainRowInCurrentProcess = j;
            }
        }

        localMax.maxValue = maxValue;
        localMax.currentRank = rank;

        // каждый процесс рассылает свой локально максимальный элемент по всем столцам, все процесы принимают уже глобально максимальный элемент
        MPI_Allreduce(&localMax, &globalMax, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);

        //Вычисление ведущей строки всей системы
        if (rank == globalMax.currentRank) {
            if (globalMax.maxValue == 0) {
                solvingStatus = 2;
                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Send(&solvingStatus, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD);
                mainRowIteration[index] = i;
                mainRowIndexArray[i] = displs[rank] + mainRowInCurrentProcess;
                continue;
            } else {
                // Номер итерации, на которой строка с локальным номером является ведущей для всей системы
                mainRowIteration[mainRowInCurrentProcess] = i;
                //Вычисленный номер ведущей строки системы
                mainRowIndexArray[i] = displs[rank] + mainRowInCurrentProcess;
            }
        }
        MPI_Bcast(&mainRowIndexArray[i], 1, MPI_INT, globalMax.currentRank, MPI_COMM_WORLD);
        if (rank == globalMax.currentRank) {
            for (int j = 0; j < matrixSize; j++) {
                globalMatrix[j] = matrixPart[mainRowInCurrentProcess * matrixSize + j];
            }
            globalMatrix[matrixSize] = vectorPart[mainRowInCurrentProcess];
        }
        //Рассылка ведущей строки всем процессам
        MPI_Bcast(globalMatrix, matrixSize + 1, MPI_DOUBLE, globalMax.currentRank, MPI_COMM_WORLD);
        //Исключение неизвестных в столбце с номером i
        multiplyByCoeffient(globalMatrix, i);
    }
}

/*
stringIndex - номер строки, которая была ведущей на определеной итерации
iterationcurrentRank - процесс, на котором эта строка
IterationItervedindex - локальный номер этой строки (в рамках одного процесса)
*/
void findRankAndRowNumber(int rowIndex, int &iterationRank, int &localRowNumber) {
    //Определяем ранг процесса, содержащего данную строку
    for (int i = 0; i < size - 1; i++) {
        if ((displs[i] <= rowIndex) && (rowIndex < displs[i + 1])) {
            iterationRank = i;
        }
    }
    if (rowIndex >= displs[size - 1]) {
        iterationRank = size - 1;
    }
    localRowNumber = rowIndex - displs[iterationRank];
}


void backStroke() {
    int iterationRank;  // Ранг процесса, хранящего текущую ведущую строку
    int localRowNumber;    // локальный на своем процессе номер текущей ведущ
    double iterationResult;   // значение Xi, найденное на итерации
    // Основной цикл
    for (int i = matrixSize - 1; i >= 0; i--) {
        findRankAndRowNumber(mainRowIndexArray[i], iterationRank, localRowNumber);
        // Определили ранг процесса, содержащего текущую ведущую строку, и номер этой строки на процессе
        // Вычисляем значение неизвестной
        if (rank == iterationRank) {
            if (matrixPart[localRowNumber * matrixSize + i] == 0) {
                if (vectorPart[localRowNumber] == 0) {
                    iterationResult = 0.0;
                } else {
                    solvingStatus = 0;
                    MPI_Barrier(MPI_COMM_WORLD);
                    MPI_Send(&solvingStatus, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD);
                    break;
                }
            } else {
                iterationResult = vectorPart[localRowNumber] / matrixPart[localRowNumber * matrixSize + i];
            }
            //нашли значение переменной
            resultPart[localRowNumber] = iterationResult;
        }
        MPI_Bcast(&iterationResult, 1, MPI_DOUBLE, iterationRank, MPI_COMM_WORLD);
        //подстановка найденной переменной
        double val;
        for (int j = 0; j < partSizePerProcess; j++) {
            if (mainRowIteration[j] < i) {
                val = matrixPart[matrixSize * j + i] * iterationResult;
                vectorPart[j] -= val;
            }
        }
    }
}

//------------------------------------------------------------------------------

void prepareMPI(int argc, char *argv[]) {
    if (rank == 0) {
        printf("Start \n");
    }
    setvbuf(stdout, 0, _IONBF, 0); // режим доступа и размер буфера
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Barrier(MPI_COMM_WORLD);
}


void printExecutionStatus(double time) {
    if (rank == 0) {
        if (solvingStatus == 1) {
            printf("\nDone \n"); // решение найдено
        } else if (solvingStatus == 0) {
            printf("\nNo roots.\n"); // решений нет
        } else if (solvingStatus == 2) {
            printf("\nCan't \n"); // решений бесконечно много
        };
        printf("\nTime: %f\n", time);
    }
}

void calculateWithGaussMethod() {
    MPI_Barrier(MPI_COMM_WORLD);
    double time = MPI_Wtime();
    triangulate();
    backStroke();
    //сбор данных, передача от всех одному (нулевому процессу)
    MPI_Gatherv(resultPart, sendcounts[rank], MPI_DOUBLE, result, sendcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    printExecutionStatus(MPI_Wtime() - time);
}

void printResult() {
    for (int i = 0; i < matrixSize; i++) {
        if (usePrint) {
            printf("\nresult[%i]=%f", i, result[i]);
        }
    }
}

void finalize() {
    if (rank == 0) {
        printf("\nEnd");
    }
    MPI_Finalize();
    delete[] matrix;
    delete[] vector;
    delete[] result;
}

int main(int argc, char *argv[]) {
    prepareMPI(argc, argv);
    initData();
    disributeDataBetweenProcesses();
    calculateWithGaussMethod();
    printResult();
    finalize();
}