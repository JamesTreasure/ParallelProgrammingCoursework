#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <pthread.h>

int arrayLength;
double **myArray;
double precision;
int numberOfThreads;
pthread_mutex_t lock;
pthread_t *threads;
pthread_barrier_t barr;
int ended;

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
            printf("%.4f,", myArray[i][j]);
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

void average(int *inc) {

    int thrNum = *inc;
    int start_row, end_row;
    int n = (arrayLength - 2) / numberOfThreads;
    if (thrNum == numberOfThreads - 1) {
        start_row = thrNum * n + 1;
        end_row = arrayLength - 2;
    } else {
        start_row = thrNum * n + 1;
        end_row = start_row + n - 1;
    }

    printf("Thread %d will average row %d to %d\n", thrNum, start_row, end_row);

    double **tempArray;
    ended = 0;

    tempArray = malloc(arrayLength * sizeof(double *));
    for (int i = 0; i < arrayLength; ++i) {
        tempArray[i] = malloc(arrayLength * sizeof(double));
    }

    for (int i = 0; i < arrayLength; ++i) {
        for (int j = 0; j < arrayLength; ++j) {
            tempArray[i][j] = myArray[i][j];
        }
    }

    while (ended != 1) {
        for (int i = start_row; i <= end_row; ++i) {
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

        printf("Ended is %d\n", ended);

        if (isPrecisionMet(myArray, tempArray, arrayLength)) {
            printf("Precision was met\n");
            ended = 1;
        }

        printf("ended is %d\n", ended);

        //printf("Thread %d finished averaging %d to %d\n", thrNum, start_row, end_row);

        int barrierWait = pthread_barrier_wait(&barr);

        if(barrierWait != 0 && barrierWait != PTHREAD_BARRIER_SERIAL_THREAD)
        {
            printf("Could not wait on barrier.\n");
            exit(-1);
        }else{
            printf("Barrier wait has worked");
        }


        for (int i = 0; i < arrayLength; ++i) {
            for (int j = 0; j < arrayLength; ++j) {
                tempArray[i][j] = myArray[i][j];
            }
        }
    }
}

int main() {
    clock_t begin = clock();
    srand(time(NULL));

    numberOfThreads = 2;


    arrayLength = 10;
    precision = 1.0;

    setupArray(arrayLength);
    /* create threads */
    threads = malloc(sizeof(pthread_t)*numberOfThreads);
    for(int i = 0; i < numberOfThreads; ++i)
    {
        int *inc = malloc(sizeof(i));
        *inc = i;
        if(pthread_create(&threads[i], NULL, (void*(*)(void*))average, inc))
        {
            printf("Could not create thread: %d\n", i);
            return -1;
        }

    }

    pthread_barrier_init(&barr, NULL, numberOfThreads);

    for(int i = 0; i< numberOfThreads; ++i)
    {
        if(pthread_join(threads[i], NULL))
        {
            printf("Could not join thread: %d\n", i);
            return -1;
        }
    }

    clock_t end = clock();
    double time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
    printf("%f", time_spent);

    return 0;
}
