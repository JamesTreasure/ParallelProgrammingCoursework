#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#define TRUE  (1==1)
#define FALSE (!TRUE)

int arrayLength = 10;
int precisionMet = FALSE;
double precision = 0.01;


int isNotAnEdge(int index) {
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

void loadFileInto2dArray(char *fileName, double arr[]) {
    int *file;
    file = readFile(fileName);
    for (int i = 0; i < arrayLength * arrayLength; ++i) {
        arr[i] = file[i];
    }
}

void updateArrayValues(int start_row, int end_row, double arr[], double tempArray[]) {
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
    srand(time(NULL));

    // Initialize the MPI environment
    MPI_Init(&argc, &argv);
    clock_t begin = clock();
    struct timespec start, finish;
    double elapsed;

    clock_gettime(CLOCK_MONOTONIC, &start);


    // Get the number of processes
    int numberOfProcesses;
    MPI_Comm_size(MPI_COMM_WORLD, &numberOfProcesses);

    // Get the rank of the process
    int processRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &processRank);

    double arr[arrayLength * arrayLength];
    if (processRank == 0) {
        printf("------------------------------------------------------\n");
        printf("There are %d processes\n", numberOfProcesses);
        printf("Array size is %d by %d\n", arrayLength, arrayLength );
        printf("Precision is %f\n", precision);
        //printf("------------------------------------------------------\n");
       //  This is for a random array
       for (int i = 0; i < arrayLength * arrayLength; ++i) {
           arr[i] = generateRandomNumber();
       }
        // loadFileInto2dArray(
        //         "/Users/jamestreasure/GitHub/ParallelProgrammingCoursework/testArray.txt",
        //         arr);

        // for (int i = 0; i < arrayLength; i++) {
        //     for (int j = 0; j < arrayLength; j++) {
        //         printf("%f,", arr[i * arrayLength + j]);
        //     }
        //     printf("\n");
        // }
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
            if (isNotAnEdge(i)) {
                double above = arr[i - arrayLength];
                double below = arr[i + arrayLength];
                double left = arr[i - 1];
                double right = arr[i + 1];
                tempArray[i] = (above + below + left + right) / 4;
            }
        }

        if (processRank == 0) {
            //set precision met to true
            precisionMet = TRUE;

            //replace rows with updated rows from master
            updateArrayValues(start_row, end_row, arr, tempArray);

            //replace rows with updated rows from slave
            for (int i = 1; i < numberOfProcesses; ++i) {
                double receiveArray[arrayLength * arrayLength];
                int startEndArray[2];
                MPI_Recv(startEndArray, 2, MPI_INT, i, 0,
                         MPI_COMM_WORLD,
                         MPI_STATUS_IGNORE);
                MPI_Recv(receiveArray, arrayLength * arrayLength, MPI_DOUBLE, i,
                         0,
                         MPI_COMM_WORLD,
                         MPI_STATUS_IGNORE);
                //replace rows with updated rows from slaves
                updateArrayValues(startEndArray[0], startEndArray[1], arr, receiveArray);
            }
        } else {
            int startEndArray[2] = {start_row, end_row};
            MPI_Send(startEndArray, 2, MPI_INT, 0, 0, MPI_COMM_WORLD);
            MPI_Send(tempArray, arrayLength * arrayLength, MPI_DOUBLE, 0, 0,
                     MPI_COMM_WORLD);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(&precisionMet, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }

    if (processRank == 0) {
        //printf("Finished array:\n");
        //printArray(arrayLength, arr);
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
