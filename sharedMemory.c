#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
#include <unistd.h>
#include <ctype.h>

#define TRUE  (1==1)
#define FALSE (!TRUE)

int arrayLength;
int precisionMet = FALSE;
double **myArray;
double precision;
int verbose = 0;
int numberOfThreads;
pthread_barrier_t barrier;
pthread_t *threads;

/**
 * Checks to see if a given index of a 2d array is on the edge
 * @param arrayLength
 * @param row
 * @param column
 * @return true or false
 */
int isNotAnEdge(int arrayLength, int row, int column) {
    return (row >= 1 && column != 0 && column != arrayLength - 1 &&
            row != arrayLength - 1);
}

void print2DArray(printArray arrayLength, double **myArray) {
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

/**
 * If no file is passed in then a random array will be generated to the given
 * array length
 * @param arrayLength
 */
void setupRandomArray(int arrayLength) {
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

int getNumberOfLinesInTextFile(char *fileName) {
    FILE *input;
    input = fopen(fileName, "r");
    int ch = 0;
    int numberOfLines = 0;

    do {
        ch = fgetc(input);
        if (ch == '\n')
            numberOfLines++;
    } while (ch != EOF);

    if (ch != '\n' && numberOfLines != 0)
        numberOfLines++;

    fclose(input);

    //as text file is 1 number per line, sqrt to get array length
    arrayLength = sqrt(numberOfLines);

    return numberOfLines;
}

int *readFile(char *fileName) {
    FILE *input;
    input = fopen(fileName, "r");
    if (input == NULL) {
        perror("Error");
        exit(0);
    }

    int numberOfLines = getNumberOfLinesInTextFile(fileName);
    int *fileArray = malloc(numberOfLines * sizeof(int));

    for (int i = 0; i < numberOfLines; i++) {
        fscanf(input, "%d,", &fileArray[i]);
    }

    return fileArray;
}

void loadFileInto2dArray(char *fileName) {
    int *file;

    file = readFile(fileName);

    myArray = malloc(arrayLength * sizeof(double *));
    for (int i = 0; i < arrayLength; ++i) {
        myArray[i] = malloc(arrayLength * sizeof(double));
    }

    for (int i = 0; i < arrayLength; ++i) {
        for (int j = 0; j < arrayLength; ++j) {
            myArray[i][j] = file[i * arrayLength + j];
        }
    }

}

void relaxArray(int *thread) {
    //allocates beginning and end row of array to a thread based on thread number
    int threadNumber = *thread;
    int start_row, end_row;
    int n = (arrayLength - 2) / numberOfThreads;
    if (threadNumber == numberOfThreads - 1) {
        start_row = threadNumber * n + 1;
        end_row = arrayLength - 2;
    } else {
        start_row = threadNumber * n + 1;
        end_row = start_row + n - 1;
    }


    if(verbose) printf("Thread %d will work on row %d to %d\n", threadNumber, start_row, end_row);

    //allocates memory for the tempArray and copies in the original array
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

    //keep running until precision is met
    while (!precisionMet) {

        /**
         * Was originally using a pthread_mutex_lock and pthread_cond_wait here to ensure all
         * the threads started together and to keep track of the precision. It works
         * just as well using a barrier and the code is far simpler and there appears to be no
         * difference in performance
         */
        int preAverageBarrier = pthread_barrier_wait(&barrier);
        if(preAverageBarrier != 0 && preAverageBarrier != PTHREAD_BARRIER_SERIAL_THREAD)
        {
            printf("Could not wait on barrier before averaging.\n");
            exit(-1);
        }else if(preAverageBarrier == PTHREAD_BARRIER_SERIAL_THREAD){
            pthread_barrier_destroy(&barrier);
        }

        precisionMet = TRUE;

        /**
         * loop to average each value in given rows and check if any value doesn't meet precision
         * set precision met to be false so loop will re run
         */
        for (int i = start_row; i <= end_row; ++i) {
            for (int j = 0; j < arrayLength; ++j) {
                if (isNotAnEdge(arrayLength, i, j)) {
                    double above = tempArray[i - 1][j];
                    double below = tempArray[i + 1][j];
                    double left = tempArray[i][j - 1];
                    double right = tempArray[i][j + 1];
                    myArray[i][j] = (above + below + left + right) / 4;

                    double currentPrecision = fabs(myArray[i][j] - tempArray[i][j]);
                    if(currentPrecision > precision){
                        precisionMet = FALSE;
                    }
                }
            }
        }
        if (verbose) printf("Thread %d finished averaging %d to %d\n", threadNumber, start_row, end_row);

        int postAverageBarrier = pthread_barrier_wait(&barrier);
        if(postAverageBarrier != 0 && postAverageBarrier != PTHREAD_BARRIER_SERIAL_THREAD)
        {
            printf("Could not wait on barrier after averaging.\n");
            exit(-1);
        }else if(postAverageBarrier == PTHREAD_BARRIER_SERIAL_THREAD){
            pthread_barrier_destroy(&barrier);
        }

        //update tempArray to match the working array
        for (int i = 0; i < arrayLength; ++i) {
            for (int j = 0; j < arrayLength; ++j) {
                tempArray[i][j] = myArray[i][j];
            }
        }
    }
}

int main(int argc, char **argv) {
    clock_t begin = clock();
    struct timespec start, finish;
    double elapsed;

    clock_gettime(CLOCK_MONOTONIC, &start);
    srand(time(NULL));

    //argument handling
    int c;
    int fileFlag = FALSE;
    char *fileLocation;
    while ((c = getopt(argc, argv, "f:")) != -1)
        switch (c) {
            case 'f':
                fileFlag = TRUE;
                fileLocation = optarg;
                break;
            case '?':
                if (optopt == 'f')
                    fprintf(stderr, "Option -%c requires a path to a array text file\n", optopt);
                else if (isprint(optopt))
                    fprintf(stderr, "Unknown option `-%c'.\n", optopt);
                else
                    fprintf(stderr,
                            "Unknown option character `\\x%x'.\n",
                            optopt);
                return 1;
            default:
                abort();
        }

    if (fileFlag & (argc == 5)) {
        loadFileInto2dArray(fileLocation);
        precision = atof(argv[3]);
        numberOfThreads = strtol(argv[4], NULL, 10);
    } else if ((argc < 4) & (fileFlag == FALSE)) {
        printf("Not enough arguments. Running with default 10x10, 0.1 precision, 1 thread.\n");
        arrayLength = 10;
        precision = 0.1;
        numberOfThreads = 1;
        setupRandomArray(arrayLength);
    } else {
        arrayLength = strtol(argv[1], NULL, 10);
        precision = atof(argv[2]);
        numberOfThreads = strtol(argv[3], NULL, 10);
        setupRandomArray(arrayLength);
    }

    printf("----------------------------------------------------------------------------\n");
    printf("Beginning with the following arguments:\n");
    printf("Matrix dimension: %d by %d\n", arrayLength, arrayLength);
    printf("Precision: %.10f\n", precision);
    printf("Number of threads %d\n", numberOfThreads);
    printf("----------------------------------------------------------------------------\n");
    printf("Starting array:\n");
    printArray(arrayLength, myArray);
    printf("----------------------------------------------------------------------------\n");

    pthread_barrier_init(&barrier, NULL, numberOfThreads);
    threads = malloc(sizeof(pthread_t) * numberOfThreads);

    //create threads
    int threadArguments[numberOfThreads];
    for (int i = 0; i < numberOfThreads; ++i) {
        threadArguments[i]=i;
        if (pthread_create(&threads[i], NULL, (void *(*)(void *)) relaxArray, (void *) &threadArguments[i])) {
            printf("Could not create thread: %d\n", i);
            return -1;
        }
    }

    //join threads
    for (int i = 0; i < numberOfThreads; ++i) {
        if (pthread_join(threads[i], NULL)) {
            printf("Could not join thread: %d\n", i);
            return -1;
        }
    }

    clock_t end = clock();
    double time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

    printf("----------------------------------------------------------------------------\n");
    printf("Finished array:\n");
    printArray(arrayLength, myArray);
    printf("----------------------------------------------------------------------------\n");
    printf("Real time elapsed is %f\n", elapsed);
    printf("CPU time elapsed is %f\n", time_spent);
    return 0;
}
