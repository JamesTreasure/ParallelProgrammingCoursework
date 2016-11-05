#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <pthread.h>

int arrayLength;
double **myArray;
double precision = 0.000000000001;
int verbose = 0;
int numberOfThreads = 1;
int iterationCount = 0;
int waitCounter = 0;
pthread_mutex_t lock;
pthread_cond_t condition;
pthread_barrier_t barrier;
pthread_t *threads;

int isNotAnEdge(int arrayLength, int row, int column) {
    return (row >= 1 && column != 0 && column != arrayLength - 1 && row != arrayLength - 1);
}

int isPrecisionMet(double **myArray, double **tempArray, int arrayLength) {
    double max = 0;
    for (int i = 0; i < arrayLength; ++i) {
        for (int j = 0; j < arrayLength; ++j) {
            if (isNotAnEdge(arrayLength, i, j)) {
                double currentPrecision = fabs(myArray[i][j] - tempArray[i][j]);//fabs takes the absoulte value of a double
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
            printf("%.10f,", myArray[i][j]);
        }
        printf("\n");
    }
}

double generateRandomNumber() {
    double range = 100;
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

int getNumberOfLinesInTextFile(char *fileName){
    FILE *input;
    input = fopen(fileName, "r"); // reopen file to reset ptr
    int ch = 0;
    int numberOfLines = 0;

    do
    {
        ch = fgetc(input);
        if(ch == '\n')
            numberOfLines++;
    } while (ch != EOF);

    if(ch != '\n' && numberOfLines != 0)
        numberOfLines++;

    fclose(input);

    arrayLength = sqrt(numberOfLines);

    return numberOfLines;
}

int *readFile(char *fileName) {
    FILE *input;
    input = fopen(fileName, "r"); // reopen file to reset ptr
    int numberOfLines = getNumberOfLinesInTextFile(fileName);

    int *arr = malloc(numberOfLines * sizeof(int));

    for (int i = 0; i < numberOfLines; i++)
    {
        fscanf(input, "%d,", &arr[i] );

    }

    return arr;
}

void useFile(char *fileName){
    int *fileRead;

    fileRead = readFile(fileName);

    myArray = malloc(arrayLength * sizeof(double *));
    for (int i = 0; i < arrayLength; ++i) {
        myArray[i] = malloc(arrayLength * sizeof(double));
    }

    for (int i = 0; i < arrayLength; ++i) {
        for (int j = 0; j < arrayLength; ++j) {
            myArray[i][j] = fileRead[i*arrayLength+j];
        }
    }

}

void relax(int *inc) {

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

    printf("Thread %d will relax row %d to %d\n", thrNum, start_row, end_row);

    int ended = 0;
    double **tempArray;

    tempArray = malloc(arrayLength * sizeof(double *));
    for (int i = 0; i < arrayLength; ++i) {
        tempArray[i] = malloc(arrayLength * sizeof(double));
    }

    for (int i = 0; i < arrayLength; ++i) {
        for (int j = 0; j < arrayLength; ++j) {
            tempArray[i][j] = myArray[i][j];
        }
    }

    //potentially just change to a while true and then break once the precision has been met?
    while (ended != 1) {

        pthread_mutex_lock(&lock);
        waitCounter++;

        if(waitCounter == numberOfThreads){
            waitCounter = 0;
            iterationCount++;
            ended = 0;
            pthread_cond_broadcast(&condition); //unblocks all threads currently blocked on &condition
            if(verbose) printf("Thread %d is broadcasting signals\n", thrNum);
            if(verbose) printf("Round %d starts now\n", iterationCount);
        }else{
            if(verbose) printf("Thread number %d is waiting for lock to be lifted\n", thrNum);
            pthread_cond_wait(&condition, &lock);
        }

        pthread_mutex_unlock(&lock);

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

        int barrierWait = pthread_barrier_wait(&barrier);

        if(barrierWait != 0 && barrierWait != PTHREAD_BARRIER_SERIAL_THREAD)
        {
            //printf("Could not wait on barrier.\n");
            exit(-1);
        }else{
            //printf("Barrier wait has worked\n");
        }

        if (isPrecisionMet(myArray, tempArray, arrayLength)) {
            if(verbose) printf("Thread %d says that precision has been met\n", thrNum);
            ended = 1;
        }else{
            if(verbose) printf("Thread %d says that the precision has NOT been met\n", thrNum);
        }

        if(verbose) printf("Thread %d finished averaging %d to %d\n", thrNum, start_row, end_row);

        for (int i = 0; i < arrayLength; ++i) {
            for (int j = 0; j < arrayLength; ++j) {
                tempArray[i][j] = myArray[i][j];
            }
        }
    }
}

int main(int argc, char *argv[]) {
    clock_t begin = clock();
    struct timespec start, finish;
    double elapsed;

    clock_gettime(CLOCK_MONOTONIC, &start);
    srand(time(NULL));


    useFile("/home/james/ClionProjects/ParallelProgrammingCoursework/testArrayLarge200by200.txt");
    //setupArray(arrayLength);
    print2DArray(arrayLength,myArray);

    pthread_mutex_init(&lock, NULL);
    pthread_cond_init(&condition, NULL);
    pthread_barrier_init(&barrier, NULL, numberOfThreads);

    threads = malloc(sizeof(pthread_t)*numberOfThreads);

    for(int i = 0; i < numberOfThreads; ++i)
    {
        int *inc = malloc(sizeof(i));
        *inc = i;
        if(pthread_create(&threads[i], NULL, (void *(*)(void *)) relax, inc))
        {
            printf("Could not create thread: %d\n", i);
            return -1;
        }

    }

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
    printf("CPU time spent is %f\n", time_spent);

    clock_gettime(CLOCK_MONOTONIC, &finish);

    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    printf("Real time is %f\n", elapsed);

    print2DArray(arrayLength,myArray);

    return 0;
}
