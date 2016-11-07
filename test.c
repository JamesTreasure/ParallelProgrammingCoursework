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
double **myArray;
double precision;
int verbose = 0;
int numberOfThreads;
int waitCounter = 0;
pthread_mutex_t lock;
pthread_cond_t condition;
pthread_barrier_t barrier;
pthread_t *threads;

int isNotAnEdge(int arrayLength, int row, int column) {
    return (row >= 1 && column != 0 && column != arrayLength - 1 && row != arrayLength - 1);
}

void average(int start_row, int end_row, double *const *tempArray) {
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
}

int isPrecisionMet(double **myArray, double **tempArray, int arrayLength) {
    double max = 0;
    for (int i = 0; i < arrayLength; ++i) {
        for (int j = 0; j < arrayLength; ++j) {
            if (isNotAnEdge(arrayLength, i, j)) {
                double currentPrecision = fabs(
                        myArray[i][j] - tempArray[i][j]);//fabs takes the absoulte value of a double
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

int getNumberOfLinesInTextFile(char *fileName) {
    FILE *input;
    input = fopen(fileName, "r"); // reopen file to reset ptr
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

    arrayLength = sqrt(numberOfLines);

    return numberOfLines;
}

int *readFile(char *fileName) {
    FILE *input;
    input = fopen(fileName, "r");
    if(input == NULL){
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

void useFile(char *fileName) {
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

void relax(int *inc) {

    int threadNumber = *inc;
    int start_row, end_row;
    int n = (arrayLength - 2) / numberOfThreads;
    if (threadNumber == numberOfThreads - 1) {
        start_row = threadNumber * n + 1;
        end_row = arrayLength - 2;
    } else {
        start_row = threadNumber * n + 1;
        end_row = start_row + n - 1;
    }


    printf("Thread %d will relax row %d to %d\n", threadNumber, start_row, end_row);

    int ended = FALSE;
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

    while (!ended) {

        pthread_mutex_lock(&lock);
        waitCounter++;

        if (waitCounter == numberOfThreads) {
            waitCounter = 0;
            ended = FALSE;
            pthread_cond_broadcast(&condition); //unblocks all threads currently blocked on &condition
            if (verbose) printf("Thread %d is broadcasting signals\n", threadNumber);
        } else {
            if (verbose) printf("Thread number %d is waiting for lock to be lifted\n", threadNumber);
            pthread_cond_wait(&condition, &lock);
        }

        pthread_mutex_unlock(&lock);

        average(start_row, end_row, tempArray);

        pthread_barrier_wait(&barrier);

        if (isPrecisionMet(myArray, tempArray, arrayLength)) {
            if (verbose) printf("Thread %d says that precision has been met\n", threadNumber);
            ended = TRUE;
        } else {
            if (verbose) printf("Thread %d says that the precision has NOT been met\n", threadNumber);
        }

        if (verbose) printf("Thread %d finished averaging %d to %d\n", threadNumber, start_row, end_row);

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

    int c;
    int fileFlag;
    char *fileLocation;

    while ((c = getopt (argc, argv, "vf:")) != -1)
        switch (c)
        {
            case 'v':
                verbose = TRUE;
                break;
            case 'f':
                fileFlag = TRUE;
                fileLocation = optarg;
                break;
            case '?':
                if (optopt == 'f')
                    fprintf (stderr, "Option -%c requires a path to a array text file\n", optopt);
                else if (isprint (optopt))
                    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
                else
                    fprintf (stderr,
                             "Unknown option character `\\x%x'.\n",
                             optopt);
                return 1;
            default:
                abort ();
        }

    if (fileFlag & (argc == 5)) {
        useFile(fileLocation);
        precision = atof(argv[3]);
        numberOfThreads = strtol(argv[4], NULL, 10);
    } else if ((argc < 4) & (fileFlag == FALSE)) {
        printf("Not enough arguments. Running with default 10x10, 0.1 precision, 1 thread.\n");
        arrayLength = 10;
        precision = 0.1;
        numberOfThreads = 1;
        setupArray(arrayLength);
    } else {
        arrayLength = strtol(argv[1], NULL, 10);
        precision = atof(argv[2]);
        numberOfThreads = strtol(argv[3], NULL, 10);
        setupArray(arrayLength);
    }

    int index;
    for (index = optind; index < argc; index++)
        printf ("Non-option argument number %d is %s\n", index, argv[index]);

    printf("----------------------------------------------------------------------------\n");
    printf("Beginning with the following arguments:\n");
    printf("Matrix dimension: %d by %d\n", arrayLength, arrayLength);
    printf("Precision: %.10f\n", precision);
    printf("Number of threads %d\n", numberOfThreads);
    printf("----------------------------------------------------------------------------\n");
    printf("Starting array:\n");
    print2DArray(arrayLength, myArray);
    printf("----------------------------------------------------------------------------\n");


    pthread_mutex_init(&lock, NULL);
    pthread_cond_init(&condition, NULL);
    pthread_barrier_init(&barrier, NULL, numberOfThreads);

    threads = malloc(sizeof(pthread_t) * numberOfThreads);

    for (int i = 0; i < numberOfThreads; ++i) {
        int *inc = malloc(sizeof(i));
        *inc = i;
        if (pthread_create(&threads[i], NULL, (void *(*)(void *)) relax, inc)) {
            printf("Could not create thread: %d\n", i);
            return -1;
        }
    }

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
    print2DArray(arrayLength, myArray);
    printf("----------------------------------------------------------------------------\n");
    printf("Real time elapsed %f\n", elapsed);
    printf("CPU time elapsed is %f\n", time_spent);
    return 0;
}
