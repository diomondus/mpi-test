#include <mpi.h>
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <memory.h>

#define DIRECTIONS_COUNT 9

typedef struct {
    double density;        // макроскопическая плотность
    double velocity[2];    // макроскопическая скорость, 0 - горизонтельно, 1 - вертикально
} MacroData;

typedef struct {
    double particleDistribution[DIRECTIONS_COUNT];    // распределения частиц по направлениям
    double data[DIRECTIONS_COUNT];
} Cell;

typedef struct {
    int first;
    int last;
} RowLimits;

typedef struct {
    int height, width;                          // размеры сетки
    double relaxationTime, gridSpeed;        // время релаксации и скорость сетки
    Cell **nodes;
} Grid;

double weights[] = {4.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 9, 1.0 / 36, 1.0 / 36, 1.0 / 36,
                    1.0 / 36}; // весовые коэффициенты

void firstRow(const Grid *grid, const Cell *upperLimit, int hasUpperLimit, int nodeRow, int nodeColumn,
              const Cell *currentNode, const double *data);

void firstColumn(const Grid *grid, int nodeRow, int nodeColumn, const Cell *currentNode, const double *data);

void lastRow(const Grid *grid, const Cell *lowerLimit, int hasLowerLimit, int nodeRow, int nodeColumn,
             const Cell *currentNode, const double *data);

void lastColumn(const Grid *grid, int nodeRow, int nodeColumn, const Cell *currentNode, const double *data);

void firstRowOrColomn(const Grid *grid, const Cell *upperLimit, int hasUpperLimit, int nodeRow, int nodeColumn,
                      const Cell *currentNode, const double *data);

void lastRowOrColomn(const Grid *grid, const Cell *lowerLimit, int hasLowerLimit, int nodeRow, int nodeColumn,
                     const Cell *currentNode, const double *data);

void firstRowOrLastColomn(const Grid *grid, const Cell *upperLimit, int hasUpperLimit, int nodeRow, int nodeColumn,
                          const Cell *currentNode, const double *data);

void lastRowOrFirstColumn(const Grid *grid, const Cell *lowerLimit, int hasLowerLimit, int nodeRow, int nodeColumn,
                          const Cell *currentNode, const double *data);

void exchangeAdjacentGridLines(const Grid *localGrid, int rank, int hasUpperLimit, int hasLowerLimit, size_t rowSize,
                               Cell *&upperLimit, Cell *&lowerLimit);

const double unitVectors[DIRECTIONS_COUNT][2] = {{0,  0},
                                                 {1,  0},
                                                 {0,  1},
                                                 {-1, 0},
                                                 {0,  -1},
                                                 {1,  1},
                                                 {-1, 1},
                                                 {-1, -1},
                                                 {1,  -1}}; // вектора скоростей

double vectorModulus(double *vector) {
    return sqrt(pow(vector[0], 2) + pow(vector[1], 2));
}

void addVectors(const double *first, const double *second, double *result) {
    for (int i = 0; i < 2; ++i) {
        result[i] = first[i] + second[i];
    }
}

void mulVectors(const double *vector, double multiplier, double *result) {
    for (int i = 0; i < 2; ++i) {
        result[i] = vector[i] * multiplier;
    }
}

double scalarMul(const double *first, const double *second) {
    double result = 0;
    for (int i = 0; i < 2; ++i) {
        result += first[i] * second[i];
    }
    return result;
}

double cosBetweenVectors(double *first, double *second) {
    return scalarMul(first, second) / (vectorModulus(first) * vectorModulus(second));
}

RowLimits getMyLimits(int gridWidth, int worldSizeMinesOne, int index) {
    int div = gridWidth % worldSizeMinesOne; //остаток от деления сетки на количество вычислиетелей (остаток сетки)
    RowLimits res;
    res.first = gridWidth / worldSizeMinesOne * index + (index < div ? index : div);
    res.last = gridWidth / worldSizeMinesOne * (index + 1) - 1 + (index < div ? index + 1 : div);
    return res;
}

int minimumRowCount(int dataTypeSize, int worldSize, int minimumSizePerNode) {
    return (int) ceil((sqrt((minimumSizePerNode * worldSize) / dataTypeSize)));
}

void
calculateMicroVelocityInPoint(double *particleDistribution, double microDensity, double gridSpeed, double *result) {
    double temp[2];
    for (int i = 0; i < 2; ++i) {
        result[i] = 0;
    }
    for (int direction = 0; direction < DIRECTIONS_COUNT; ++direction) {
        mulVectors((double *) unitVectors[direction], particleDistribution[direction], temp);
        mulVectors((double *) temp, gridSpeed, temp);
        addVectors(result, temp, result);
    }
    mulVectors(result, 1. / microDensity, result);
}

double calculateMicroDensityInPoint(const double *directionsDistribution) {
    double density = 0;
    for (int direction = 0; direction < DIRECTIONS_COUNT; ++direction) {
        density += directionsDistribution[direction];
    }
    return density;
}

// Коэффициент для равновесного распределения по направлениям
double directionCoeffient(int direction, double gridVelocity, double *microVelocity) {
    double eu = scalarMul((double *) unitVectors[direction], microVelocity);
    double u2 = scalarMul(microVelocity, microVelocity);
    return 3 * (eu + (3 * pow(eu, 2) - u2) / (gridVelocity * 2)) / gridVelocity;
}

void calculateEquilibriumDistribution(double gridSpeed, double microdensity, double *microvelocity, double *result) {
    for (int direction = 0; direction < DIRECTIONS_COUNT; ++direction) {
        result[direction] =
                (1 + directionCoeffient(direction, gridSpeed, microvelocity)) * microdensity * weights[direction];
    }
}

void defineCellData(const Grid *grid, const Cell *upperLimit, int hasUpperLimit, const Cell *lowerLimit,
                    int hasLowerLimit, int nodeRow, int nodeColumn) {
    Cell *currentNode = &grid->nodes[nodeRow][nodeColumn];
    double *data = currentNode->data;
    data[0] = currentNode->particleDistribution[0];

    firstRow(grid, upperLimit, hasUpperLimit, nodeRow, nodeColumn, currentNode, data);
    firstColumn(grid, nodeRow, nodeColumn, currentNode, data);
    lastRow(grid, lowerLimit, hasLowerLimit, nodeRow, nodeColumn, currentNode, data);
    lastColumn(grid, nodeRow, nodeColumn, currentNode, data);
    firstRowOrColomn(grid, upperLimit, hasUpperLimit, nodeRow, nodeColumn, currentNode, data);
    lastRowOrColomn(grid, lowerLimit, hasLowerLimit, nodeRow, nodeColumn, currentNode, data);
    firstRowOrLastColomn(grid, upperLimit, hasUpperLimit, nodeRow, nodeColumn, currentNode, data);
    lastRowOrFirstColumn(grid, lowerLimit, hasLowerLimit, nodeRow, nodeColumn, currentNode, data);
}

void lastRowOrFirstColumn(const Grid *grid, const Cell *lowerLimit, int hasLowerLimit, int nodeRow, int nodeColumn,
                          const Cell *currentNode, double *data) {
    if (nodeRow == grid->height - 1 || nodeColumn == 0) {
        if (nodeColumn != 0 && hasLowerLimit) {
            data[5] = lowerLimit[nodeColumn - 1].particleDistribution[5];
        } else {
            data[5] = currentNode->particleDistribution[7];
        }
    } else {
        data[5] = grid->nodes[nodeRow + 1][nodeColumn - 1].particleDistribution[5];
    }
}

void firstRowOrLastColomn(const Grid *grid, const Cell *upperLimit, int hasUpperLimit, int nodeRow, int nodeColumn,
                          const Cell *currentNode, double *data) {
    if (nodeRow == 0 || nodeColumn == grid->width - 1) {
        if (nodeColumn != grid->width - 1 && hasUpperLimit) {
            data[7] = upperLimit[nodeColumn + 1].particleDistribution[7];
        } else {
            data[7] = currentNode->particleDistribution[5];
        }
    } else {
        data[7] = grid->nodes[nodeRow - 1][nodeColumn + 1].particleDistribution[7];
    }
}

void lastRowOrColomn(const Grid *grid, const Cell *lowerLimit, int hasLowerLimit, int nodeRow, int nodeColumn,
                     const Cell *currentNode, double *data) {
    if (nodeRow == grid->height - 1 || nodeColumn == grid->width - 1) {
        if (nodeColumn != grid->width - 1 && hasLowerLimit) {
            data[6] = lowerLimit[nodeColumn + 1].particleDistribution[6];
        } else {
            data[6] = currentNode->particleDistribution[8];
        }
    } else {
        data[6] = grid->nodes[nodeRow + 1][nodeColumn + 1].particleDistribution[6];
    }
}

void firstRowOrColomn(const Grid *grid, const Cell *upperLimit, int hasUpperLimit, int nodeRow, int nodeColumn,
                      const Cell *currentNode, double *data) {
    if (nodeRow == 0 || nodeColumn == 0) {
        if (nodeColumn != 0 && hasUpperLimit) {
            data[8] = upperLimit[nodeColumn - 1].particleDistribution[8];

        } else {
            data[8] = currentNode->particleDistribution[6];
        }
    } else {
        data[8] = grid->nodes[nodeRow - 1][nodeColumn - 1].particleDistribution[8];
    }
}

void lastColumn(const Grid *grid, int nodeRow, int nodeColumn, const Cell *currentNode, double *data) {
    if (nodeColumn == grid->width - 1) {
        data[3] = currentNode->particleDistribution[1];
    } else {
        data[3] = grid->nodes[nodeRow][nodeColumn + 1].particleDistribution[3];
    }
}

void lastRow(const Grid *grid, const Cell *lowerLimit, int hasLowerLimit, int nodeRow, int nodeColumn,
             const Cell *currentNode, double *data) {
    if (nodeRow == grid->height - 1) {
        if (hasLowerLimit) {
            data[2] = lowerLimit[nodeColumn].particleDistribution[2];
        } else {
            data[2] = currentNode->particleDistribution[4];
        }
    } else {
        data[2] = grid->nodes[nodeRow + 1][nodeColumn].particleDistribution[2];
    }
}

void firstColumn(const Grid *grid, int nodeRow, int nodeColumn, const Cell *currentNode, double *data) {
    if (nodeColumn == 0) {
        data[1] = currentNode->particleDistribution[3];
    } else {
        data[1] = grid->nodes[nodeRow][nodeColumn - 1].particleDistribution[1];
    }
}

void firstRow(const Grid *grid, const Cell *upperLimit, int hasUpperLimit, int nodeRow, int nodeColumn,
              const Cell *currentNode, double *data) {
    if (nodeRow == 0) {
        if (hasUpperLimit) {
            data[4] = upperLimit[nodeColumn].particleDistribution[4];
        } else {
            data[4] = currentNode->particleDistribution[2];
        }
    } else {
        data[4] = grid->nodes[nodeRow - 1][nodeColumn].particleDistribution[4];
    }
}

void streaming(Grid *localGrid, int rank, int worldSize) {
    int hasUpperLimit = rank != 1;
    int hasLowerLimit = rank != (worldSize - 1);

    Cell *upperLimit = nullptr, *lowerLimit = nullptr;
    size_t rowSize = sizeof(Cell) * localGrid->width;

    exchangeAdjacentGridLines(localGrid, rank, hasUpperLimit, hasLowerLimit, rowSize, upperLimit, lowerLimit);

    //обработка распространения
    for (int row = 0; row < localGrid->height; row++) {
        for (int column = 0; column < localGrid->width; column++) {
            defineCellData(localGrid, upperLimit, hasUpperLimit, lowerLimit, hasLowerLimit, row, column);
        }
    }
}

void exchangeAdjacentGridLines(const Grid *localGrid, int rank, int hasUpperLimit, int hasLowerLimit, size_t rowSize,
                               Cell *&upperLimit, Cell *&lowerLimit) {
    if (hasUpperLimit) {
        upperLimit = static_cast<Cell *>(malloc(rowSize));
        memcpy(upperLimit, localGrid->nodes[0], rowSize);
    }
    if (hasLowerLimit) {
        lowerLimit = static_cast<Cell *>(malloc(rowSize));
        memcpy(lowerLimit, localGrid->nodes[localGrid->height - 1], rowSize);
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
}

void getNewParticleDistribution(const double *cellDataDistribution,
                                const double *equilibriumDistribution,
                                double relaxationGasTime,
                                double *result) {
    for (int direction = 0; direction < DIRECTIONS_COUNT; ++direction) {
        result[direction] = cellDataDistribution[direction] +
                            (equilibriumDistribution[direction] - cellDataDistribution[direction]) / relaxationGasTime;
    }
}

void processCollision(Grid *pg) {
    for (int row = 0; row < pg->height; ++row) {
        for (int column = 0; column < pg->width; ++column) {
            Cell *currentNode = &pg->nodes[row][column];
            double density = calculateMicroDensityInPoint(currentNode->data);
            double velocityInPoint[2];
            calculateMicroVelocityInPoint(currentNode->data, density, pg->gridSpeed, velocityInPoint);
            double equilibriumDistribution[DIRECTIONS_COUNT];
            calculateEquilibriumDistribution(pg->gridSpeed, density, velocityInPoint, equilibriumDistribution);
            getNewParticleDistribution(currentNode->data, equilibriumDistribution, pg->relaxationTime,
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
    for (int direction = 0; direction < DIRECTIONS_COUNT; ++direction) {
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

void initializeGrid(Grid *pg, int gridSize, RowLimits Limits, double gridSpeed, double relaxationTime) {
    pg->gridSpeed = gridSpeed;
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
            currentState->density = calculateMicroDensityInPoint(currentNode->particleDistribution);
            calculateMicroVelocityInPoint(currentNode->particleDistribution, currentState->density, pg->gridSpeed,
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

    MacroData *state;
    Grid grid;


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

    if (rank != 0) {
        RowLimits rowLimits = getMyLimits(gridWidth, worldSize - 1, rank - 1);
        initializeGrid(&grid, gridWidth, rowLimits, speed, relaxationTime);
    }
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