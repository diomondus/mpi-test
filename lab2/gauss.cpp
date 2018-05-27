#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>

double *matrix, *vector, *result; //Основные данные. Заполняются через rand()
double *matrixPart, *vectorPart, *resultPart, mess = 0.0; // Нарезка для данного процесса
int matrixSize, status, countStr; //размер матрицы, статус решения задачи, строки для текущего процесса
int *numMainStr, *numMainStrIt; //ведущие строки для каждой итерации, номера итераций ведущих строк
int size, currentRank, *mass1, *range;//Размер, ранг, рассылка, количество на каждый процесс

void initMatrixAndVectors() {
    int balanceStr; //Число строк, ещё не распределённых по процессам
    if (currentRank == 0) //Заполняем
    {
        int show;
        matrixSize = 3;
        show = 1;
        srand(0); //псевдослучайные

        matrix = new double[matrixSize * matrixSize];
        vector = new double[matrixSize];
        result = new double[matrixSize];

        for (int i = 0; i < matrixSize * matrixSize; i++) {
            matrix[i] = rand();
            if (show == 1) printf("\nmatrix[%i]=%f", i, matrix[i]);
        };

        for (int i = 0; i < matrixSize; i++) {
            vector[i] = rand();
            if (show == 1) printf("\nb[%i]=%f", i, vector[i]);
            result[i] = 0;
        };
    }

    MPI_Bcast(&matrixSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

    //Определение размера части данных, расположенных на конкретном процессе

    balanceStr = matrixSize;
    for (int i = 0; i < currentRank; i++)
        balanceStr = balanceStr - balanceStr / (size - i);

    countStr = balanceStr / (size - currentRank);
    matrixPart = new double[countStr * matrixSize];//выделяем память под строки матрицы
    vectorPart = new double[countStr]; //память под элементы столбца свободных членов
    resultPart = new double[countStr]; //память под элементы вектора результата
    numMainStr = new int[matrixSize]; //массив индексов ведущих строк системы на каждой итерации
    numMainStrIt = new int[countStr]; /* итерация, на которой соответствующая строка системы, расположенная на процессе, выбрана  ведущей */
    mass1 = new int[size];
    range = new int[size];

    for (int i = 0; i < countStr; i++)
        numMainStrIt[i] = -1;
}

//------------------------------------------------------------------------------

// Распределение исходных данных между процессами
void disributeDataBetweenProcesses() {
    int *matrElem;             //Индекс первого элемента матрицы, передаваемого процессу
    int *matrRang;            //Число элементов матрицы, передаваемых процессу
    int sizestr;
    int balance;
    int prevSize;
    int prevIndex;
    int portion;
    matrElem = new int[size];
    matrRang = new int[size];

    balance = matrixSize;

    for (int i = 0; i < size; i++)  //Определяем, сколько элементов матрицы будет передано каждому процессу
    {
        prevSize = (i == 0) ? 0 : matrRang[i - 1];
        prevIndex = (i == 0) ? 0 : matrElem[i - 1];
        portion = (i == 0) ? 0 : sizestr;             //число строк, отданных предыдущему процессу
        balance -= portion;
        sizestr = balance / (size - i);
        matrRang[i] = sizestr * matrixSize;
        matrElem[i] = prevIndex + prevSize;
    };

    //Рассылка матрицы
    MPI_Scatterv(matrix, matrRang, matrElem, MPI_DOUBLE, matrixPart, matrRang[currentRank], MPI_DOUBLE, 0,
                 MPI_COMM_WORLD);


    balance = matrixSize;

    for (int i = 0; i < size; i++) {
        int prevSize = (i == 0) ? 0 : range[i - 1];
        int prevIndex = (i == 0) ? 0 : mass1[i - 1];
        int portion = (i == 0) ? 0 : range[i - 1];
        balance -= portion;
        range[i] = balance / (size - i);
        mass1[i] = prevIndex + prevSize;
    };

    //Рассылка вектора
    MPI_Scatterv(vector, range, mass1, MPI_DOUBLE, vectorPart, range[currentRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    status = 1;
    MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&mess, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    delete[] matrRang;
    delete[] matrElem;
}

//------------------------------------------------------------------------------

void Raw(int numIter, double *glStr) {
    double koef;
    //для каждой строки в процессе
    for (int i = 0; i < countStr; i++) {
        if (numMainStrIt[i] == -1) {
            koef = matrixPart[i * matrixSize + numIter] / glStr[numIter];
            for (int j = numIter; j < matrixSize; j++) {
                matrixPart[i * matrixSize + j] -= glStr[j] * koef;
            };
            vectorPart[i] -= glStr[matrixSize] * koef;

        };
    };
}

//------------------------------------------------------------------------------

void gauss() {
    int VedIndex;   // индекс ведущей строки на конкретном процессе
    struct {
        double maxValue;
        int currentRank;
    } localMax, glbMax; //максимальный элемент+номер процесса, у которого он
    double *glbWMatr = new double[matrixSize + 1]; //т.е. строка матрицы+значение вектора
    for (int i = 0; i < matrixSize; i++) {
        // Вычисление ведущей строки
        double maxValue = 0;
        int index = -1;
        for (int j = 0; j < countStr; j++) {
            index = j;
            if ((numMainStrIt[j] == -1) && (maxValue < fabs(matrixPart[i + matrixSize * j]))) {
                maxValue = fabs(matrixPart[i + matrixSize * j]);
                VedIndex = j;
            };
        };

        localMax.maxValue = maxValue;
        localMax.currentRank = currentRank;

        //каждый процесс рассылает свой локально максимальный элемент по всем столцам, все процесы принимают уже глобально максимальный элемент
        MPI_Allreduce(&localMax, &glbMax, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);

        //Вычисление ведущей строки всей системы
        if (currentRank == glbMax.currentRank) {
            if (glbMax.maxValue == 0) {
                status = 2;
                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Send(&status, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD);
                numMainStrIt[index] = i;
                numMainStr[i] = mass1[currentRank] + VedIndex;
                continue;

            } else {
                // Номер итерации, на которой строка с локальным номером является ведущей для всей системы
                numMainStrIt[VedIndex] = i;
                //Вычисленный номер ведущей строки системы
                numMainStr[i] = mass1[currentRank] + VedIndex;
            };
        };
        MPI_Bcast(&numMainStr[i], 1, MPI_INT, glbMax.currentRank, MPI_COMM_WORLD);
        if (currentRank == glbMax.currentRank) {
            for (int j = 0; j < matrixSize; j++)
                glbWMatr[j] = matrixPart[VedIndex * matrixSize + j];
            glbWMatr[matrixSize] = vectorPart[VedIndex];
        };
        //Рассылка ведущей строки всем процессам
        MPI_Bcast(glbWMatr, matrixSize + 1, MPI_DOUBLE, glbMax.currentRank, MPI_COMM_WORLD);
        //Исключение неизвестных в столбце с номером i
        Raw(i, glbWMatr);
    };
}

/*
stringIndex - номер строки, которая была ведущей на определеной итерации
iterationcurrentRank - процесс, на котором эта строка
IterationItervedindex - локальный номер этой строки (в рамках одного процесса)
*/
void Frp(int stringIndex, int &iterationcurrentRank, int &IterationItervedindex) {
    //Определяем ранг процесса, содержащего данную строку
    for (int i = 0; i < size - 1; i++) {
        if ((mass1[i] <= stringIndex) && (stringIndex < mass1[i + 1]))
            iterationcurrentRank = i;
    }
    if (stringIndex >= mass1[size - 1])
        iterationcurrentRank = size - 1;
    IterationItervedindex = stringIndex - mass1[iterationcurrentRank];

}


void gaussRevert() {
    int itCurrentRank;  // Ранг процесса, хранящего текущую ведущую строку
    int indexMain;    // локальный на своем процессе номер текущей ведущ
    double iterRes;   // значение Xi, найденное на итерации
    double val;
    // Основной цикл
    for (int i = matrixSize - 1; i >= 0; i--) {
        Frp(numMainStr[i], itCurrentRank, indexMain);
        // Определили ранг процесса, содержащего текущую ведущую строку, и номер этой строки на процессе
        // Вычисляем значение неизвестной
        if (currentRank == itCurrentRank) {
            if (matrixPart[indexMain * matrixSize + i] == 0) {
                if (vectorPart[indexMain] == 0)
                    iterRes = mess;
                else {
                    status = 0;
                    MPI_Barrier(MPI_COMM_WORLD);
                    MPI_Send(&status, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD);
                    break;
                };
            } else
                iterRes = vectorPart[indexMain] / matrixPart[indexMain * matrixSize + i];
            //нашли значение переменной
            resultPart[indexMain] = iterRes;
        };
        MPI_Bcast(&iterRes, 1, MPI_DOUBLE, itCurrentRank, MPI_COMM_WORLD);
        //подстановка найденной переменной
        for (int j = 0; j < countStr; j++)
            if (numMainStrIt[j] < i) {
                val = matrixPart[matrixSize * j + i] * iterRes;
                vectorPart[j] -= val;
            };
    };
}

//------------------------------------------------------------------------------

int main(int argc, char *argv[]) {
    double startTime, time;

    printf("Start \n");
    setvbuf(stdout, 0, _IONBF, 0);   // режим доступа и размер буфера
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &currentRank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Barrier(MPI_COMM_WORLD);

    initMatrixAndVectors();
    disributeDataBetweenProcesses();

    MPI_Barrier(MPI_COMM_WORLD);

    startTime = MPI_Wtime();

    gauss();
    gaussRevert();

    //сбор данных, передача от всех одному (нулевому процессу)
    MPI_Gatherv(resultPart, range[currentRank], MPI_DOUBLE, result, range, mass1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    time = MPI_Wtime() - startTime;
    if (currentRank == 0) {
        // решение найдено
        if (status == 1) {
            printf("\nDone \n");
        };
        // решений нет
        if (status == 0) {
            printf("\nNo roots.\n");
        }
        // решений бесконечно много
        if (status == 2) {
            printf("\nCan't \n");
        };
        printf("\nTim: %f\n", time);
    };
    if (currentRank == 0) {
        printf("\n End");
    };

    for (int i = 0; i < matrixSize; i++) {
        printf("\nresult[%i]=%f", i, result[i]);
    };

    MPI_Finalize();

    delete[] matrix;
    delete[] vector;
    delete[] result;
}