#include <mpi.h>
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <memory.h>

#define DIRECTIONS 9

typedef struct {
    double density;        // макроскопическая плотность
    double velocity[2];    // макроскопическая скорость, 0 - горизонтельно, 1 - вертикально
} MacroData;

typedef struct {
    double particleDistribution[DIRECTIONS];    // распределения частиц по направлениям
    double tmp[DIRECTIONS];
} Cell;

typedef struct {
    int first;
    int last;
} RowBounds;

typedef struct {
    int height, width;                          // размеры сетки
    double relaxationTime, latticeSpeed;        // время релаксации и скорость сетки
    Cell **nodes;
} Grid;

double weights[] = {4.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 36, 1.0 / 36, 1.0 / 36,
                    1.0 / 36}; // весовые коэфф

const double elementalVectors[DIRECTIONS][2] = {{0,  0},
                                                {1,  0},
                                                {0,  1},
                                                {-1, 0},
                                                {0,  -1},
                                                {1,  1},
                                                {-1, 1},
                                                {-1, -1},
                                                {1,  -1}}; // вектора скоростей

void addVectors(const double *first, const double *second, double *result) {
    int i;
    for (i = 0; i < 2; ++i) {
        result[i] = first[i] + second[i];
    }
}

void mulVector(const double *vector, double multiplier, double *result) {
    int i;
    for (i = 0; i < 2; ++i) {
        result[i] = vector[i] * multiplier;
    }
}

double modulusOfVector(double *vector) {
    return sqrt(pow(vector[0], 2) + pow(vector[1], 2));
}

double scalarMultiplication(const double *first, const double *second) {
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
* @param worldSizeMinesOne количество вычислитетей
* @param index номер вычислителя начиная с 0
* @return Индексы первой и последней строк в сетке.
*/
RowBounds getMyBounds(int gridWidth, int worldSizeMinesOne, int index) {
    int remainder = gridWidth % worldSizeMinesOne;
    RowBounds res;
    res.first = gridWidth / worldSizeMinesOne * index + (index < remainder ? index : remainder);
    res.last = gridWidth / worldSizeMinesOne * (index + 1) - 1 + (index < remainder ? index + 1 : remainder);
    return res;
}

int minimumRowCount(int dataTypeSizeInBytes, int numberOfComputationalNodes, int minimumSizeOfSystemPerNode) {
    return (int) ceil((sqrt((minimumSizeOfSystemPerNode * numberOfComputationalNodes)
                            / dataTypeSizeInBytes)));
}

void fillTempFieldForNode(const Grid *grid, const Cell *upperBound, int hasUpperBound, const Cell *lowerBound,
                          int hasLowerBound, int row, int column);

void initNodes(const Grid *pg, const RowBounds &bounds, const double *center);

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
    for (direction = 0; direction < DIRECTIONS; ++direction) {
        mulVector((double *) elementalVectors[direction], particleDistribution[direction], temp);
        mulVector((double *) temp, latticeSpeed, temp);
        addVectors(result, temp, result);
    }
    mulVector(result, 1. / macroscopicDensity, result);
}

/**
 * @param directionsDistribution распределение частиц по направлениям
 * @return микроскопическая плотность в точке
 */
double calculateDensity(const double *directionsDistribution) {
    double density = 0;
    for (int direction = 0; direction < DIRECTIONS; ++direction) {
        density += directionsDistribution[direction];
    }
    return density;
}

/**
 * @param direction направление
 * @param latticeVelocity скорость сетки
 * @param microVelocity микроскопическая скорость
 * @return Коэффициент для вычисления равновесного распределения по направлениям
 */
double directionCoeffient(int direction, double latticeVelocity, double *microVelocity) {
    double eu = scalarMultiplication((double *) elementalVectors[direction], microVelocity);
    double u2 = scalarMultiplication(microVelocity, microVelocity);
    return 3 * (eu + (3 * pow(eu, 2) - u2) / (latticeVelocity * 2)) / latticeVelocity;
}


// частицы по узлам

/**
 * feq[i]=w[i]*ro(1+3/c*(e[i],u)+3/2c^2*(3(e[i],u)^2-(u,u)))
 * @param latticeSpeed скорость сетки
 * @param density микроскопическая плотность
 * @param velocity микроскопическая скорость
 * @param result равновесное распределение по направлениям (OUT)
 */
void calculateEquilibriumDistribution(double latticeSpeed, double density, double *velocity, double *result) {
    for (int direction = 0; direction < DIRECTIONS; ++direction) {
        result[direction] = (1 + directionCoeffient(direction, latticeSpeed, velocity)) * density * weights[direction];
    }
}

void streaming(Grid *pg, int rank, int worldSize) {
    int hasUpperBound = rank != 1;
    int hasLowerBound = rank != (worldSize - 1);
    Cell *upperBound = nullptr, *lowerBound = nullptr;
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
void updateDistribution(const double *tempDistribution,
                        const double *equilibriumDistribution,
                        double relaxationTime,
                        double *result) {
    for (int direction = 0; direction < DIRECTIONS; ++direction) {
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
            double equilibriumDistribution[DIRECTIONS];
            calculateEquilibriumDistribution(pg->latticeSpeed, density, velocity, equilibriumDistribution);
            // новое распределение
            updateDistribution(currentNode->tmp, equilibriumDistribution, pg->relaxationTime,
                               currentNode->particleDistribution);
        }
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
void generateTwisterData(const double *centerOfGrid, int row, int column, double *result) {
    double perpendicular[2];
    perpendicular[0] = centerOfGrid[1] - column;
    perpendicular[1] = row - centerOfGrid[0];
    int direction;
    for (direction = 0; direction < DIRECTIONS; ++direction) {
        result[direction] = tangentProjectionCubed(perpendicular, (double *) elementalVectors[direction]);
    }
}

void initializeGrid(Grid *pg, int gridSize, RowBounds bounds, double latticeSpeed, double relaxationTime) {
    pg->latticeSpeed = latticeSpeed;
    pg->nodes = static_cast<Cell **>(calloc((size_t) pg->height, sizeof(Cell *)));
    pg->height = bounds.last - bounds.first + 1;
    double center[2] = {(gridSize - 1.0) / 2, (gridSize - 1.0) / 2};
    pg->width = gridSize;
    pg->relaxationTime = relaxationTime;
    initNodes(pg, bounds, center);
}

void initNodes(const Grid *pg, const RowBounds &bounds, const double *center) {
    for (int row = 0; row < pg->height; ++row) {
        pg->nodes[row] = static_cast<Cell *>(calloc((size_t) pg->width, sizeof(Cell)));
        int column;
        for (column = 0; column < pg->width; ++column) {
            Cell *currentNode = &pg->nodes[row][column];
            generateTwisterData(center, bounds.first + row, column, currentNode->particleDistribution);
        }
    }
}

void getState(Grid *pg, MacroData *state) {
    for (int row = 0; row < pg->height; ++row) {
        for (int column = 0; column < pg->width; ++column) {
            MacroData *currentState = &state[row * pg->width + column];
            Cell *currentNode = &pg->nodes[row][column];
            currentState->density = calculateDensity(currentNode->particleDistribution);
            calculateVelocity(currentNode->particleDistribution, currentState->density, pg->latticeSpeed,
                              currentState->velocity);
        }
    }
}

void saveState(MacroData *macrodata, int width, int index) {
    char fileName[30];
    sprintf(fileName, "state%d.csv", index);
    FILE *file = fopen(fileName, "w");
    fprintf(file, "x,y,Vx,Vy,p\n");
    for (int row = 0; row < width; ++row) {
        for (int column = 0; column < width; ++column) {
            MacroData *macro = &macrodata[row * width + column];
            fprintf(file, "%d,%d,%f,%f,%f\n", row, column, macro->velocity[0], macro->velocity[1], macro->density);
        }
    }
    fclose(file);
}

// тестовые данные 1 0.8 100 999 100
void initializeSimulationParameters(char **argv, const int worldSize, double *speed, double *relaxationTime,
                                    int *totalTime, int *stateRate, int *gridWidth) {
    sscanf(argv[1], "%lf", speed);
    sscanf(argv[2], "%lf", relaxationTime);
    sscanf(argv[3], "%i", totalTime);
    sscanf(argv[4], "%i", stateRate);
    sscanf(argv[5], "%i", gridWidth);
    *gridWidth = minimumRowCount(sizeof(Cell), worldSize - 1, 100 * 1024 * 1024); // 100мб на каждый вычислитель
}

int main(int argc, char *argv[]) {

    int rank, worldSize, gridWidth, totalTime, stateRate;
    double speed, relaxationTime;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    initializeSimulationParameters(argv, worldSize, &speed, &relaxationTime, &totalTime, &stateRate, &gridWidth);
    // Служебные данные для передачи состояний
    auto *stateSizes = static_cast<int *>(calloc((size_t) worldSize, sizeof(int)));
    auto *stateOffsets = static_cast<int *>(calloc((size_t) worldSize, sizeof(int)));
    for (int nonMasterNode = 1; nonMasterNode < worldSize; ++nonMasterNode) {
        RowBounds bounds = getMyBounds(gridWidth, worldSize - 1, nonMasterNode - 1);
        stateSizes[nonMasterNode] = (bounds.last - bounds.first + 1) * gridWidth * sizeof(MacroData);
        stateOffsets[nonMasterNode] = bounds.first * gridWidth * sizeof(MacroData);
    }

    Grid grid;
    if (rank != 0) {
        RowBounds rowBounds = getMyBounds(gridWidth, worldSize - 1, rank - 1);
        initializeGrid(&grid, gridWidth, rowBounds, speed, relaxationTime);
    }
    MacroData *state;
    if (rank == 0) {
        state = static_cast<MacroData *>(calloc((size_t) gridWidth * gridWidth, sizeof(MacroData)));
    } else {
        state = static_cast<MacroData *>(calloc((size_t) grid.height * grid.width, sizeof(MacroData)));
    }
    MPI_Barrier(MPI_COMM_WORLD);
    double time = MPI_Wtime();
    for (int i = 0; i < totalTime; i++) {
        if (i % stateRate == 0) {
            if (rank != 0) {
                getState(&grid, state);
            }
            MPI_Gatherv(state, rank == 0 ? 0 : grid.width * grid.height * sizeof(MacroData), MPI_BYTE, state,
                        stateSizes, stateOffsets, MPI_BYTE, 0, MPI_COMM_WORLD);
            if (rank == 0) {
                saveState(state, gridWidth, i / stateRate);
            }
        }
        if (rank != 0) {
            streaming(&grid, rank, worldSize);
            collision(&grid);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        printf("Algorithm time is %lf", MPI_Wtime() - time);
    }
    free(state);
    MPI_Finalize();
}