#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <pthread.h>

int arrayLength;
double **myArray;
double precision;
int numberOfThreads;

int isNotAnEdge(int arrayLength, int row, int column) {
    return (row >= 1 && column != 0 && column != arrayLength - 1 && row != arrayLength - 1);
}

int isPrecisionMet(double **myArray, double **tempArray, int arrayLength) {
    double max = 0;
    for (int i = 0; i < arrayLength; ++i) {
        for (int j = 0; j < arrayLength; ++j) {
            if (isNotAnEdge(arrayLength, i, j)) {
                double currentPrecision = fabs(myArray[i][j] - tempArray[i][j]);
                if (currentPrecision > max) {
                    max = currentPrecision;
                }
            }
        }
    }
    return (max < precision);
}

void print2DArray(int arrayLength, double **myArray) {
    for (int i = 0; i < arrayLength; i++) {
        for (int j = 0; j < arrayLength; j++) {
            printf("%.100f,", myArray[i][j]);
        }
        printf("\n");
    }
}

double generateRandomNumber() {
    double range = (100 - 0);
    double div = RAND_MAX / range;
    return 0 + (rand() / div);
}

void setupArray(int arrayLength) {
    myArray = malloc(arrayLength * sizeof(double *));

    for (int i = 0; i < arrayLength; ++i) {
        myArray[i] = malloc(arrayLength * sizeof(double));
    }

    for (int i = 0; i < arrayLength; ++i) {
        for (int j = 0; j < arrayLength; ++j) {
            myArray[i][j] = generateRandomNumber();
        }
    }
}

void setupNonRandomTestingArray(int arrayLength) {
    myArray = malloc(arrayLength * sizeof(double *));

    for (int i = 0; i < arrayLength; ++i) {
        myArray[i] = malloc(arrayLength * sizeof(double));
    }

    int count = 0;

    for (int i = 0; i < arrayLength; ++i) {
        for (int j = 0; j < arrayLength; ++j) {
            myArray[i][j] = count;
            count++;
        }
    }
}

void average(double **myArray) {
    int ended = 0;
    double **tempArray;

    tempArray = malloc(arrayLength * sizeof(double *));
    for (int i = 0; i < arrayLength; ++i) {
        tempArray[i] = malloc(arrayLength * sizeof(double));
    }

    while (ended != 1) {
        for (int i = 0; i < arrayLength; ++i) {
            for (int j = 0; j < arrayLength; ++j) {
                tempArray[i][j] = myArray[i][j];
            }
        }

        for (int i = 0; i < arrayLength; ++i) {
            for (int j = 0; j < arrayLength; ++j) {
                if (isNotAnEdge(arrayLength, i, j)) {
                    double above = tempArray[i - 1][j];
                    double below = tempArray[i + 1][j];
                    double left = tempArray[i][j - 1];
                    double right = tempArray[i][j + 1];
                    myArray[i][j] = (above + below + left + right) / 4;
                }
            }
        }

        if (isPrecisionMet(myArray, tempArray, arrayLength)) {
            ended = 1;
        }
    }
}

int main() {
    printf("Entered main!");
    clock_t begin = clock();
    srand(time(NULL));


    numberOfThreads = 4;


    arrayLength = 10;
    precision = 0.0000000000000000000000000000000000000000001;

    setupArray(arrayLength);
    print2DArray(arrayLength,myArray);
    average(myArray);


    clock_t end = clock();
    double time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
    printf("%f\n", time_spent);
    print2DArray(arrayLength,myArray);

    return 0;
}
