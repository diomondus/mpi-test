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
} RowLimits;

typedef struct {
    int height, width;                          // размеры сетки
    double relaxationTime, latticeSpeed;        // время релаксации и скорость сетки
    Cell **nodes;
} Grid;

double weights[] = {4.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 36, 1.0 / 36, 1.0 / 36,
                    1.0 / 36}; // весовые коэффициенты

const double unitVectors[DIRECTIONS][2] = {{0,  0},
                                           {1,  0},
                                           {0,  1},
                                           {-1, 0},
                                           {0,  -1},
                                           {1,  1},
                                           {-1, 1},
                                           {-1, -1},
                                           {1,  -1}}; // вектора скоростей

void addVectors(const double *first, const double *second, double *result) {
    for (int i = 0; i < 2; ++i) {
        result[i] = first[i] + second[i];
    }
}

void mulVector(const double *vector, double multiplier, double *result) {
    for (int i = 0; i < 2; ++i) {
        result[i] = vector[i] * multiplier;
    }
}

double modulusOfVector(double *vector) {
    return sqrt(pow(vector[0], 2) + pow(vector[1], 2));
}

double scalarMultiplication(const double *first, const double *second) {
    double result = 0;
    for (int i = 0; i < 2; ++i) {
        result += first[i] * second[i];
    }
    return result;
}

double cosBetweenVectors(double *first, double *second) {
    return scalarMultiplication(first, second) / (modulusOfVector(first) * modulusOfVector(second));
}

RowLimits getMyLimits(int gridWidth, int worldSizeMinesOne, int index) {
    int remainder = gridWidth % worldSizeMinesOne;
    RowLimits res;
    res.first = gridWidth / worldSizeMinesOne * index + (index < remainder ? index : remainder);
    res.last = gridWidth / worldSizeMinesOne * (index + 1) - 1 + (index < remainder ? index + 1 : remainder);
    return res;
}

int minimumRowCount(int dataTypeSizeInBytes, int numberOfComputationalNodes, int minimumSizeOfSystemPerNode) {
    return (int) ceil((sqrt((minimumSizeOfSystemPerNode * numberOfComputationalNodes) / dataTypeSizeInBytes)));
}

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
        mulVector((double *) unitVectors[direction], particleDistribution[direction], temp);
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
    double eu = scalarMultiplication((double *) unitVectors[direction], microVelocity);
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

/**
* Заполняет поле tmp в ноде сетки
* @param grid сетка
* @param upperLimit верхняя граница
* @param hasUpperLimit есть ли верхняя граница у сетки
* @param lowerLimit нижняя граница
* @param hasLowerLimit есть ли нижняя граница у сетки
*/
void fillTempFieldForNode(const Grid *grid, const Cell *upperLimit, int hasUpperLimit, const Cell *lowerLimit,
                          int hasLowerLimit, int nodeRow, int nodeColumn) {
    Cell *currentNode = &grid->nodes[nodeRow][nodeColumn];
    double *tmp = currentNode->tmp;
    tmp[0] = currentNode->particleDistribution[0];

    int firstRow = nodeRow == 0;
    int firstColumn = nodeColumn == 0;
    int lastRow = nodeRow == grid->height - 1;
    int lastColumn = nodeColumn == grid->width - 1;

    if (firstRow) {
        if (hasUpperLimit) {
            tmp[4] = upperLimit[nodeColumn].particleDistribution[4];
        } else {
            tmp[4] = currentNode->particleDistribution[2];
        }
    } else {
        tmp[4] = grid->nodes[nodeRow - 1][nodeColumn].particleDistribution[4];
    }
    if (firstColumn) {
        tmp[1] = currentNode->particleDistribution[3];
    } else {
        tmp[1] = grid->nodes[nodeRow][nodeColumn - 1].particleDistribution[1];
    }
    if (lastRow) {
        if (hasLowerLimit) {
            tmp[2] = lowerLimit[nodeColumn].particleDistribution[2];
        } else {
            tmp[2] = currentNode->particleDistribution[4];
        }
    } else {
        tmp[2] = grid->nodes[nodeRow + 1][nodeColumn].particleDistribution[2];
    }
    if (lastColumn) {
        tmp[3] = currentNode->particleDistribution[1];
    } else {
        tmp[3] = grid->nodes[nodeRow][nodeColumn + 1].particleDistribution[3];
    }
    if (firstRow || firstColumn) {
        if (!firstColumn && hasUpperLimit) {
            tmp[8] = upperLimit[nodeColumn - 1].particleDistribution[8];

        } else {
            tmp[8] = currentNode->particleDistribution[6];
        }
    } else {
        tmp[8] = grid->nodes[nodeRow - 1][nodeColumn - 1].particleDistribution[8];
    }
    if (lastRow || lastColumn) {
        if (!lastColumn && hasLowerLimit) {
            tmp[6] = lowerLimit[nodeColumn + 1].particleDistribution[6];
        } else {
            tmp[6] = currentNode->particleDistribution[8];
        }
    } else {
        tmp[6] = grid->nodes[nodeRow + 1][nodeColumn + 1].particleDistribution[6];
    }
    if (firstRow || lastColumn) {
        if (!lastColumn && hasUpperLimit) {
            tmp[7] = upperLimit[nodeColumn + 1].particleDistribution[7];
        } else {
            tmp[7] = currentNode->particleDistribution[5];
        }
    } else {
        tmp[7] = grid->nodes[nodeRow - 1][nodeColumn + 1].particleDistribution[7];
    }
    if (lastRow || firstColumn) {
        if (!firstColumn && hasLowerLimit) {
            tmp[5] = lowerLimit[nodeColumn - 1].particleDistribution[5];
        } else {
            tmp[5] = currentNode->particleDistribution[7];
        }
    } else {
        tmp[5] = grid->nodes[nodeRow + 1][nodeColumn - 1].particleDistribution[5];
    }
}

void streaming(Grid *pg, int rank, int worldSize) {
    int hasUpperLimit = rank != 1;
    int hasLowerLimit = rank != (worldSize - 1);
    Cell *upperLimit = nullptr, *lowerLimit = nullptr;
    size_t rowSize = sizeof(Cell) * pg->width;

    //Обмен смежными строками сетки
    if (hasUpperLimit) {
        upperLimit = static_cast<Cell *>(malloc(rowSize));
        //Копируем то, что нужно передать
        memcpy(upperLimit, pg->nodes[0], rowSize);
    }
    if (hasLowerLimit) {
        lowerLimit = static_cast<Cell *>(malloc(rowSize));
        //Копируем то, что нужно передать
        memcpy(lowerLimit, pg->nodes[pg->height - 1], rowSize);
    }
    MPI_Status status;
    for (int i = 0; i < 2; ++i) {
        if (hasLowerLimit && (rank % 2 == i)) {
            MPI_Sendrecv_replace(lowerLimit, (int) rowSize, MPI_BYTE, rank + 1, 0, rank + 1, 0, MPI_COMM_WORLD,
                                 &status);
        } else if (hasUpperLimit & (rank % 2 != i)) {
            MPI_Sendrecv_replace(upperLimit, (int) rowSize, MPI_BYTE, rank - 1, 0, rank - 1, 0, MPI_COMM_WORLD,
                                 &status);
        }
    }

    //обработка распространения
    for (int row = 0; row < pg->height; row++) {
        for (int column = 0; column < pg->width; column++) {
            fillTempFieldForNode(pg, upperLimit, hasUpperLimit, lowerLimit, hasLowerLimit, row, column);
        }
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

void processCollision(Grid *pg) {
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

void generateParticleDistribution(const double *center, int row, int column, double *result) {
    double perpendicular[2];
    perpendicular[0] = center[1] - column;
    perpendicular[1] = row - center[0];
    for (int direction = 0; direction < DIRECTIONS; ++direction) {
        result[direction] = tangentProjectionCubed(perpendicular, (double *) unitVectors[direction]);
    }
}

void initNodes(const Grid *pg, const RowLimits &Limits, const double *center) {
    for (int row = 0; row < pg->height; ++row) {
        pg->nodes[row] = static_cast<Cell *>(calloc((size_t) pg->width, sizeof(Cell)));
        int column;
        for (column = 0; column < pg->width; ++column) {
            Cell *currentNode = &pg->nodes[row][column];
            generateParticleDistribution(center, Limits.first + row, column,
                                         currentNode->particleDistribution);
        }
    }
}

void initializeGrid(Grid *pg, int gridSize, RowLimits Limits, double latticeSpeed, double relaxationTime) {
    pg->latticeSpeed = latticeSpeed;
    pg->nodes = static_cast<Cell **>(calloc((size_t) pg->height, sizeof(Cell *)));
    pg->height = Limits.last - Limits.first + 1;
    double center[2] = {(gridSize - 1.0) / 2, (gridSize - 1.0) / 2};
    pg->width = gridSize;
    pg->relaxationTime = relaxationTime;
    initNodes(pg, Limits, center);
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
        RowLimits Limits = getMyLimits(gridWidth, worldSize - 1, nonMasterNode - 1);
        stateSizes[nonMasterNode] = (Limits.last - Limits.first + 1) * gridWidth * sizeof(MacroData);
        stateOffsets[nonMasterNode] = Limits.first * gridWidth * sizeof(MacroData);
    }

    Grid grid;
    if (rank != 0) {
        RowLimits rowLimits = getMyLimits(gridWidth, worldSize - 1, rank - 1);
        initializeGrid(&grid, gridWidth, rowLimits, speed, relaxationTime);
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
            processCollision(&grid);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        printf("time %lf", MPI_Wtime() - time);
    }
    free(state);
    MPI_Finalize();
}