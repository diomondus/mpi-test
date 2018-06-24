#include <mpi.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <memory.h>

#define LATTICE_DIRECTIONS 9

double weights[] = {4.0 / 9,
                    1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9,
                    1.0 / 36, 1.0 / 36, 1.0 / 36, 1.0 / 36};

const double elementalVectors[LATTICE_DIRECTIONS][2] = {{0,  0},
                                                        {1,  0},
                                                        {0,  1},
                                                        {-1, 0},
                                                        {0,  -1},
                                                        {1,  1},
                                                        {-1, 1},
                                                        {-1, -1},
                                                        {1,  -1}};

typedef struct {
    double density;//макроскопическая плотность
    double velocity[2];    //макроскопическая скорость, 0 - горизонтельно, 1 - вертикально
} MacroData;

typedef struct {
    double particleDistribution[LATTICE_DIRECTIONS];        //распределения частиц по направлениям
    double tmp[LATTICE_DIRECTIONS];
} Cell;

typedef struct {
    int height, width;            //размеры сетки
    double relaxationTime, latticeSpeed;        //время релаксации и скорость сетки
    Cell **nodes;
} Grid;

typedef struct {
    int first;
    int last;
} RowBounds;

void sumVectors(double *first, double *second, double *result);

void multiplyVector(double *vector, double multiplier, double *result);

double modulusOfVector(double *vector);

double scalarMultiplication(double *first, double *second);

double cosBetweenVectors(double *first, double *second);

RowBounds getMyBounds(int gridWidth, int computationalProcessorsCount, int index);

double calculateDensity(double *particleDistribution);

void calculateVelocity(double *particleDistribution, double macroscopicDensity, double latticeSpeed, double *result);

void initializeSimulationParameters(const int argc, char **argv,
                                    const int worldSize,
                                    double *speed, double *relaxationTime, int *totalTime, int *snapshotRate,
                                    int *gridWidth);

void initializeGrid(Grid *pg, int gridSize, RowBounds bounds, double latticeSpeed, double relaxationTime);

void getSnapshot(Grid *pg, MacroData *snapshot);

void saveSnapshots(MacroData *snapshots, int width, int snapshotIndex);

void streaming(Grid *pg, int rank, int worldSize);

void collision(Grid *pg);

int main(int argc, char *argv[]) {

    int rank, worldSize, gridWidth, totalTime, snapshotRate;
    double speed, relaxationTime;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
//    if (worldSize < 2) {
//        printf("world is too small");
//        exit(1);
//    }
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    initializeSimulationParameters(argc, argv, worldSize, &speed, &relaxationTime, &totalTime, &snapshotRate,
                                   &gridWidth);

    Grid grid;
    int isMaster = rank == 0;

    if (isMaster) {
        printf("World size = %d\n", worldSize);
        printf("speed = %lf\n", speed);
        printf("relaxationTime = %lf\n", relaxationTime);
        printf("totalTime = %d\n", totalTime);
        printf("snapshotRate = %d\n", snapshotRate);
        printf("gridWidth = %d\n", gridWidth);
        printf("size in MB is = %ld\n", (long) (sizeof(Cell) * gridWidth * gridWidth / 1024 / 1024));
        printf("size of Cell = %d Grid=%d\n", (int) sizeof(Cell), (int) sizeof(Grid));
    }

    //Служебные данные для передачи снепшотов
    int *snapshotSizes = static_cast<int *>(calloc((size_t) worldSize, sizeof(int)));
    int *snapshotOffsets = static_cast<int *>(calloc((size_t) worldSize, sizeof(int)));
    for (int nonMasterNode = 1; nonMasterNode < worldSize; ++nonMasterNode) {
        RowBounds bounds = getMyBounds(gridWidth, worldSize - 1, nonMasterNode - 1);
        snapshotSizes[nonMasterNode] = (bounds.last - bounds.first + 1) * gridWidth * sizeof(MacroData);
        snapshotOffsets[nonMasterNode] = bounds.first * gridWidth * sizeof(MacroData);
    }

    if (!isMaster) {
        RowBounds rowBounds = getMyBounds(gridWidth, worldSize - 1, rank - 1);
        initializeGrid(&grid, gridWidth, rowBounds, speed, relaxationTime);
    }
    MacroData *snapshot;
    if (isMaster) {
        snapshot = static_cast<MacroData *>(calloc((size_t) gridWidth * gridWidth, sizeof(MacroData)));
    } else {
        snapshot = static_cast<MacroData *>(calloc((size_t) grid.height * grid.width, sizeof(MacroData)));
    }
    MPI_Barrier(MPI_COMM_WORLD);
    double timeBefore = MPI_Wtime();
    for (int i = 0; i < totalTime; i++) {
        if (i % snapshotRate == 0) {
            if (!isMaster) {
                getSnapshot(&grid, snapshot);
            }
            MPI_Gatherv(snapshot, isMaster ? 0 : grid.width * grid.height * sizeof(MacroData), MPI_BYTE, snapshot,
                        snapshotSizes, snapshotOffsets, MPI_BYTE, 0, MPI_COMM_WORLD);
            if (isMaster) {
                saveSnapshots(snapshot, gridWidth, i / snapshotRate);
            }
        }
        if (!isMaster) {
            streaming(&grid, rank, worldSize);
            collision(&grid);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    double timeAfter = MPI_Wtime();
    if (isMaster) {
        printf("Algorithm time is %lf", timeAfter - timeBefore);
    }
    free(snapshot);
//    if (!isMaster) {
//        freeGrid(&grid);
//    }
    MPI_Finalize();
}

void sumVectors(double *first, double *second, double *result) {
    int i;
    for (i = 0; i < 2; ++i) {
        result[i] = first[i] + second[i];
    }
}

void multiplyVector(double *vector, double multiplier, double *result) {
    int i;
    for (i = 0; i < 2; ++i) {
        result[i] = vector[i] * multiplier;
    }
}

double modulusOfVector(double *vector) {
    return sqrt(pow(vector[0], 2) + pow(vector[1], 2));
}

double scalarMultiplication(double *first, double *second) {
    double result = 0;
    int i;
    for (i = 0; i < 2; ++i) {
        result += first[i] * second[i];
    }
    return result;
}

double cosBetweenVectors(double *first, double *second) {
    return scalarMultiplication(first, second) / (modulusOfVector(first) * modulusOfVector(second));
}

/**
* @param gridWidth ширина квадратной сетки
* @param computationalProcessorsCount количество вычислитетей
* @param index номер вычислителя начиная с 0
* @return Индексы первой и последней строк в сетке.
*/
RowBounds getMyBounds(int gridWidth, int computationalProcessorsCount, int index) {
    int remainder = gridWidth % computationalProcessorsCount;
    RowBounds res;
    res.first = gridWidth / computationalProcessorsCount * index + (index < remainder ? index : remainder);
    res.last = gridWidth / computationalProcessorsCount * (index + 1) - 1 + (index < remainder ? index + 1 : remainder);
    return res;
}

int minimumRowCount(int dataTypeSizeInBytes, int numberOfComputationalNodes, int minimumSizeOfSystemPerNode) {
    return (int) ceil((sqrt((minimumSizeOfSystemPerNode * numberOfComputationalNodes)
                            / dataTypeSizeInBytes)));
}

void fillTempFieldForNode(const Grid *grid, const Cell *upperBound, int hasUpperBound, const Cell *lowerBound,
                          int hasLowerBound, int row, int column);

/**
 * @param particleDistribution распределение частиц по направлениям
 * @param macroscopicDensity микроскопическая плотность в точке
 * @param latticeSpeed скорость сетки
 * @param result микроскопическая скорость в точке
 */
void calculateVelocity(double *particleDistribution, double macroscopicDensity, double latticeSpeed, double *result) {
    double temp[2];
    int i;
    for (i = 0; i < 2; ++i) {
        result[i] = 0;
    }
    int direction;
    for (direction = 0; direction < LATTICE_DIRECTIONS; ++direction) {
        multiplyVector((double *) elementalVectors[direction], particleDistribution[direction], temp);
        multiplyVector((double *) temp, latticeSpeed, temp);
        sumVectors(result, temp, result);
    }

    multiplyVector(result, 1. / macroscopicDensity, result);
}

/**
 * @param particleDistribution распределение частиц по направлениям
 * @return микроскопическая плотность в точке
 */
double calculateDensity(double *particleDistribution) {
    double density = 0;
    int direction;
    for (direction = 0; direction < LATTICE_DIRECTIONS; ++direction) {
        density += particleDistribution[direction];
    }
    return density;
}

/**
 * @param direction направление
 * @param c скорость сетки
 * @param velocity микроскопическая скорость
 * @return Коэффициент для вычисления равновесного распределения по направлениям
 */
double s(int direction, double c, double *velocity) {
    double eu = scalarMultiplication((double *) elementalVectors[direction], velocity);
    double u2 = scalarMultiplication(velocity, velocity);
    return 3 * (eu + (3 * pow(eu, 2) - u2) / (c * 2)) / c;
}

/**
 * feq[i]=w[i]*ro(1+3/c*(e[i],u)+3/2c^2*(3(e[i],u)^2-(u,u)))
 * @param latticeSpeed скорость сетки
 * @param density микроскопическая плотность
 * @param velocity микроскопическая скорость
 * @param result равновесное распределение по направлениям (OUT)
 */
void calculateEquilibriumDistribution(double latticeSpeed, double density, double *velocity, double *result) {
    int direction;
    for (direction = 0; direction < LATTICE_DIRECTIONS; ++direction) {
        result[direction] = (1 + s(direction, latticeSpeed, velocity)) * density * weights[direction];
    }
}

void streaming(Grid *pg, int rank, int worldSize) {
    int hasUpperBound = rank != 1;
    int hasLowerBound = rank != (worldSize - 1);
    Cell *upperBound, *lowerBound;
    size_t rowSize = sizeof(Cell) * pg->width;


    //Обмен смежными строками сетки
    if (hasUpperBound) {
        upperBound = static_cast<Cell *>(malloc(rowSize));
        //Копируем то, что нужно передать
        memcpy(upperBound, pg->nodes[0], rowSize);
    }
    if (hasLowerBound) {
        lowerBound = static_cast<Cell *>(malloc(rowSize));
        //Копируем то, что нужно передать
        memcpy(lowerBound, pg->nodes[pg->height - 1], rowSize);
    }
    MPI_Status status;
    for (int i = 0; i < 2; ++i) {
        if (hasLowerBound && (rank % 2 == i)) {
            MPI_Sendrecv_replace(lowerBound, (int) rowSize, MPI_BYTE, rank + 1, 0, rank + 1, 0, MPI_COMM_WORLD,
                                 &status);
        } else if (hasUpperBound & (rank % 2 != i)) {
            MPI_Sendrecv_replace(upperBound, (int) rowSize, MPI_BYTE, rank - 1, 0, rank - 1, 0, MPI_COMM_WORLD,
                                 &status);
        }
    }

    //обработка распространения
    for (int row = 0; row < pg->height; row++) {
        for (int column = 0; column < pg->width; column++) {
            fillTempFieldForNode(pg, upperBound, hasUpperBound, lowerBound, hasLowerBound, row, column);
        }
    }
}

/**
* Заполняет поле tmp в ноде сетки
* @param grid сетка
* @param upperBound верхняя граница
* @param hasUpperBound есть ли верхняя граница у сетки
* @param lowerBound нижняя граница
* @param hasLowerBound есть ли нижняя граница у сетки
* @param row строка ноды
* @param column столбец ноды
*/
void fillTempFieldForNode(const Grid *grid, const Cell *upperBound, int hasUpperBound, const Cell *lowerBound,
                          int hasLowerBound, int row, int column) {
    Cell *currentNode = &grid->nodes[row][column];
    double *tmp = currentNode->tmp;
    tmp[0] = currentNode->particleDistribution[0];

    int firstRow = row == 0;
    int firstColumn = column == 0;
    int lastRow = row == grid->height - 1;
    int lastColumn = column == grid->width - 1;

    if (firstRow) {
        if (hasUpperBound) {
            tmp[4] = upperBound[column].particleDistribution[4];
        } else {
            tmp[4] = currentNode->particleDistribution[2];
        }
    } else {
        tmp[4] = grid->nodes[row - 1][column].particleDistribution[4];
    }
    if (firstColumn) {
        tmp[1] = currentNode->particleDistribution[3];
    } else {
        tmp[1] = grid->nodes[row][column - 1].particleDistribution[1];
    }
    if (lastRow) {
        if (hasLowerBound) {
            tmp[2] = lowerBound[column].particleDistribution[2];
        } else {
            tmp[2] = currentNode->particleDistribution[4];
        }
    } else {
        tmp[2] = grid->nodes[row + 1][column].particleDistribution[2];
    }
    if (lastColumn) {
        tmp[3] = currentNode->particleDistribution[1];
    } else {
        tmp[3] = grid->nodes[row][column + 1].particleDistribution[3];
    }
    if (firstRow || firstColumn) {
        if (!firstColumn && hasUpperBound) {
            tmp[8] = upperBound[column - 1].particleDistribution[8];

        } else {
            tmp[8] = currentNode->particleDistribution[6];
        }
    } else {
        tmp[8] = grid->nodes[row - 1][column - 1].particleDistribution[8];
    }
    if (lastRow || lastColumn) {
        if (!lastColumn && hasLowerBound) {
            tmp[6] = lowerBound[column + 1].particleDistribution[6];
        } else {
            tmp[6] = currentNode->particleDistribution[8];
        }
    } else {
        tmp[6] = grid->nodes[row + 1][column + 1].particleDistribution[6];
    }
    if (firstRow || lastColumn) {
        if (!lastColumn && hasUpperBound) {
            tmp[7] = upperBound[column + 1].particleDistribution[7];
        } else {
            tmp[7] = currentNode->particleDistribution[5];
        }
    } else {
        tmp[7] = grid->nodes[row - 1][column + 1].particleDistribution[7];
    }
    if (lastRow || firstColumn) {
        if (!firstColumn && hasLowerBound) {
            tmp[5] = lowerBound[column - 1].particleDistribution[5];
        } else {
            tmp[5] = currentNode->particleDistribution[7];
        }
    } else {
        tmp[5] = grid->nodes[row + 1][column - 1].particleDistribution[5];
    }
}

/**
* @param tempDistribution значение распределения в точке, полученное во время шага Streaming
* @param equilibriumDistribution равновесное распределение на основе
* @param relaxationTime время релаксации газа
* @param result новое распределение частиц
*/
void updateDistribution(double *tempDistribution,
                        double *equilibriumDistribution,
                        double relaxationTime,
                        double *result) {
    for (int direction = 0; direction < LATTICE_DIRECTIONS; ++direction) {
        result[direction] = tempDistribution[direction] +
                            (equilibriumDistribution[direction] - tempDistribution[direction]) / relaxationTime;
    }
}

void collision(Grid *pg) {
    //обработка столкновений
    for (int row = 0; row < pg->height; ++row) {
        for (int column = 0; column < pg->width; ++column) {
            Cell *currentNode = &pg->nodes[row][column];
            // плотность.
            double density = calculateDensity(currentNode->tmp);
            // скорость в точке
            double velocity[2];
            calculateVelocity(currentNode->tmp, density, pg->latticeSpeed, velocity);
            double equilibriumDistribution[LATTICE_DIRECTIONS];
            calculateEquilibriumDistribution(pg->latticeSpeed, density, velocity, equilibriumDistribution);
            // новое распределение
            updateDistribution(currentNode->tmp, equilibriumDistribution, pg->relaxationTime,
                               currentNode->particleDistribution);
        }
    }
}

void initializeSimulationParameters(const int argc, char **argv,
                                    const int worldSize,
                                    double *speed, double *relaxationTime, int *totalTime, int *snapshotRate,
                                    int *gridWidth) {
    if (argc < 5) {
        printf(
                "usage: boltzman <lattice-speed>  <relaxation-time>  <simulation-time>  <snapshot-rate> (optional <grid-width>)\n");
        exit(1);
    }

    if (sscanf(argv[1], "%lf", speed) != 1) {
        fprintf(stderr, "speed is not double");
        exit(1);
    }
    if (sscanf(argv[2], "%lf", relaxationTime) != 1) {
        fprintf(stderr, "relaxation time is not double");
        exit(1);
    }
    if (sscanf(argv[3], "%i", totalTime) != 1) {
        fprintf(stderr, "simulation time is not integer");
        exit(1);
    }

    if (sscanf(argv[4], "%i", snapshotRate) != 1) {
        fprintf(stderr, "snapshot rate is not integer");
        exit(1);
    }

    if (argc == 6) {
        if (sscanf(argv[5], "%i", gridWidth) != 1) {
            fprintf(stderr, "grid width is not integer");
            exit(1);
        }
    } else {
        // 100мб на каждый вычислитель
        *gridWidth = minimumRowCount(sizeof(Cell), worldSize - 1, 100 * 1024 * 1024);
    }
}

double generateNormalizedRandom() { return rand() / (double) RAND_MAX; }

/**
 * @param from вектор, который проецируется
 * @param to векток, на который нужно спроецировать
 * @return позитивный косинус угла между векторами в кубе, или 0
 */
double tangentProjectionCubed(double *from, double *to) {
    //Так как здесь в качестве вектора to только элементарные вектора,
    //можно просто умножить элементарный вектор на проекцию
    double cos = -cosBetweenVectors(from, to);
    return cos > 0 ? cos : generateNormalizedRandom() / 20;
}

/**
 * Генерирует распределение частиц по направлениям в точке для формирования воронки.
 * @param centerOfGrid центр воронки
 * @param row строка
 * @param column столбец
 * @param result распределение частиц по направлениям в данной точке.
 */
void generateTwisterData(double *centerOfGrid, int row, int column, double *result) {
    double perpendicular[2];
    perpendicular[0] = centerOfGrid[1] - column;
    perpendicular[1] = row - centerOfGrid[0];
    int direction;
    for (direction = 0; direction < LATTICE_DIRECTIONS; ++direction) {
        result[direction] = tangentProjectionCubed(perpendicular, (double *) elementalVectors[direction]);
    }
}

void initializeGrid(Grid *pg, int gridSize, RowBounds bounds, double latticeSpeed, double relaxationTime) {
    double centerOfGrid = (gridSize - 1.) / 2;
    //инициализация решётки
    double center[2] = {centerOfGrid, centerOfGrid};
    pg->width = gridSize;
    pg->height = bounds.last - bounds.first + 1;
    pg->nodes = static_cast<Cell **>(calloc((size_t) pg->height, sizeof(Cell *)));
    pg->latticeSpeed = latticeSpeed;
    pg->relaxationTime = relaxationTime;
    int row;
    for (row = 0; row < pg->height; ++row) {
        pg->nodes[row] = static_cast<Cell *>(calloc((size_t) pg->width, sizeof(Cell)));
        int column;
        for (column = 0; column < pg->width; ++column) {
            Cell *currentNode = &pg->nodes[row][column];
            generateTwisterData(center, bounds.first + row, column, currentNode->particleDistribution);
        }
    }
}

void getSnapshot(Grid *pg, MacroData *snapshot) {
    for (int row = 0; row < pg->height; ++row) {
        for (int column = 0; column < pg->width; ++column) {
            MacroData *currentSnapshot = &snapshot[row * pg->width + column];
            Cell *currentNode = &pg->nodes[row][column];
            currentSnapshot->density = calculateDensity(currentNode->particleDistribution);
            calculateVelocity(currentNode->particleDistribution, currentSnapshot->density, pg->latticeSpeed,
                              currentSnapshot->velocity);
        }
    }
}

void saveSnapshots(MacroData *snapshots, int width, int snapshotIndex) {
    //сохранение снимков
    char fileName[30];
    sprintf(fileName, "snapshot%d.csv", snapshotIndex);
    FILE *file = fopen(fileName, "w");
    fprintf(file, "x,y,Vx,Vy,p\n");
    for (int row = 0; row < width; ++row) {
        for (int column = 0; column < width; ++column) {
            MacroData *macro = &snapshots[row * width + column];
            fprintf(file, "%d,%d,%f,%f,%f\n", row, column, macro->velocity[0], macro->velocity[1], macro->density);
        }
    }
    fclose(file);
}