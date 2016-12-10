#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <getopt.h>
#include <ctype.h>

#define TRUE  (1==1)
#define FALSE (!TRUE)

int precisionMet = FALSE;

int isNotAnEdge(int arrayLength, int index) {
    return !((index % arrayLength == 0) || ((index + 1) % arrayLength == 0));
}

void printArray(int arrayLength, double arr[]) {
    for (int i = 0; i < arrayLength; i++) {
        for (int j = 0; j < arrayLength; j++) {
            printf("%.7f,", arr[i * arrayLength + j]);
        }
        printf("\n");
    }
}

/*
 * Only added as an option to use random numbers. No scalability testing was
 * carried out using random numbers.
 */
double generateRandomNumber() {
    double range = 100;
    double div = RAND_MAX / range;
    return 0 + (rand() / div);
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

void loadFileInto2dArray(int arrayLength, char *fileName, double arr[]) {
    int *file;
    file = readFile(fileName);
    for (int i = 0; i < arrayLength * arrayLength; ++i) {
        arr[i] = file[i];
    }
}

void updateArrayValues(int arrayLength, double precision, int start_row, int end_row, double arr[],
                       double tempArray[]) {
    for (int i = start_row * arrayLength;
         i < end_row * arrayLength + arrayLength; ++i) {
        double currentPrecision = fabs(arr[i] - tempArray[i]);
        if (currentPrecision > precision) {
            precisionMet = FALSE;
        }
        arr[i] = tempArray[i];
    }
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    srand(time(NULL));
    double *arr;

     clock_t begin = clock();
     struct timespec start, finish;
     double elapsed;

     clock_gettime(CLOCK_MONOTONIC, &start);

    //argument handling
    int c;
    int fileFlag = FALSE;
    char *fileLocation;
    int arrayLength;
    double precision;
    int master = 0;

    while ((c = getopt(argc, argv, "f:")) != -1)
        switch (c) {
            case 'f':
                fileFlag = TRUE;
                fileLocation = optarg;
                break;
            case '?':
                if (optopt == 'f')
                    fprintf(stderr,
                            "Option -%c requires a path to a array text file\n",
                            optopt);
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
        arrayLength = strtol(argv[3], NULL, 10);
        arr = malloc(arrayLength * arrayLength * sizeof(double));
        loadFileInto2dArray(arrayLength,fileLocation, arr);
        precision = atof(argv[4]);
    } else if ((argc < 3) & (fileFlag == FALSE)) {
        printf("Not enough arguments. Running with default 10x10, 0.1 precision\n");
        arrayLength = 10;
        arr = malloc(arrayLength * arrayLength * sizeof(double));
        precision = 0.1;
        for (int i = 0; i < arrayLength * arrayLength; ++i) {
            arr[i] = generateRandomNumber();
        }
    } else {
        arrayLength = strtol(argv[1], NULL, 10);
        precision = atof(argv[2]);
        arr = malloc(arrayLength * arrayLength * sizeof(double));
        for (int i = 0; i < arrayLength * arrayLength; ++i) {
            arr[i] = generateRandomNumber();
        }
    }

    // Get the number of processes
    int numberOfProcesses;
    MPI_Comm_size(MPI_COMM_WORLD, &numberOfProcesses);

    // Get the rank of the process
    int processRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &processRank);

    if (processRank == master) {
        printf("------------------------------------------------------\n");
        printf("There are %d processes, precision is %f and array size is %d by %d\n",
               numberOfProcesses, precision, arrayLength, arrayLength);
        printf("------------------------------------------------------\n");
        printf("Array before relaxing is:\n");
        for (int i = 0; i < arrayLength; i++) {
            for (int j = 0; j < arrayLength; j++) {
                printf("%f,", arr[i * arrayLength + j]);
            }
            printf("\n");
        }
        printf("------------------------------------------------------\n");
    }

    // use process rank to divide work;
    int start_row, end_row;
    int n = (arrayLength - 2) / numberOfProcesses;
    if (processRank == numberOfProcesses - 1) {
        start_row = processRank * n + 1;
        end_row = arrayLength - 2;
    } else {
        start_row = processRank * n + 1;
        end_row = start_row + n - 1;
    }


    while (!precisionMet) {
        //broadcast array to all slaves
        MPI_Bcast(arr, arrayLength * arrayLength, MPI_DOUBLE, 0,
                  MPI_COMM_WORLD);

        //setup tempArray
        double tempArray[arrayLength * arrayLength];
        for (int i = 0; i < arrayLength * arrayLength; ++i) {
            tempArray[i] = arr[i];
        }

        //average
        for (int i = start_row * arrayLength;
             i < end_row * arrayLength + arrayLength; ++i) {
            if (isNotAnEdge(arrayLength,i)) {
                double above = arr[i - arrayLength];
                double below = arr[i + arrayLength];
                double left = arr[i - 1];
                double right = arr[i + 1];
                tempArray[i] = (above + below + left + right) / 4;
            }
        }

        if (processRank == master) {
            //set precision met to true
            precisionMet = TRUE;

            //master updates the rows it has averaged
            updateArrayValues(arrayLength, precision, start_row, end_row, arr, tempArray);

            //receive from slave loop
            for (int i = 1; i < numberOfProcesses; ++i) {
                double receiveArray[arrayLength * arrayLength];
                int startEndArray[2];

                //receive start and end row from slave
                MPI_Recv(startEndArray, 2, MPI_INT, i, 0,
                         MPI_COMM_WORLD,
                         MPI_STATUS_IGNORE);

                //receive whole array
                MPI_Recv(receiveArray, arrayLength * arrayLength, MPI_DOUBLE, i,
                         0,
                         MPI_COMM_WORLD,
                         MPI_STATUS_IGNORE);

                //replace rows with updated rows from slaves
                updateArrayValues(arrayLength, precision, startEndArray[0], startEndArray[1], arr,
                                  receiveArray);
            }
        } else {
            //slave creates array with start and end row
            int startEndArray[2] = {start_row, end_row};

            //send start and end row to master
            MPI_Send(startEndArray, 2, MPI_INT, 0, 0, MPI_COMM_WORLD);

            //send whole array to master
            MPI_Send(tempArray, arrayLength * arrayLength, MPI_DOUBLE, 0, 0,
                     MPI_COMM_WORLD);
        }

        //barrier to synchronise 
        MPI_Barrier(MPI_COMM_WORLD);

        //master broadcasts if precision has been met
        MPI_Bcast(&precisionMet, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }

    if (processRank == master) {
        printf("Finished array is:\n");
        printArray(arrayLength, arr);
        printf("----------------------------------------------------------------------------\n");
         clock_t end = clock();
         double time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
         clock_gettime(CLOCK_MONOTONIC, &finish);
         elapsed = (finish.tv_sec - start.tv_sec);
         elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
        printf("Real time elapsed is %f\n", elapsed);
        printf("----------------------------------------------------------------------------\n");
    }
    MPI_Finalize();
    return 0;
}
