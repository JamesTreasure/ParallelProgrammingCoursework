#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <ctype.h>
#include <mpi.h>

#define TRUE  (1==1)
#define FALSE (!TRUE)
#define sendDataTag 2001
#define returnDataTag 2002

int arrayLength = 10;
int precisionMet = FALSE;
double **myArray;
double precision;
int verbose = 0;
double *flattenedArray;
int numberOfThreads;


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
        precisionMet = TRUE;
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

        //update tempArray to match the working array
        for (int i = 0; i < arrayLength; ++i) {
            for (int j = 0; j < arrayLength; ++j) {
                tempArray[i][j] = myArray[i][j];
            }
        }
    }
}

int main(int argc, char **argv) {

    // Initialize the MPI environment
    MPI_Init(NULL, NULL);

    // Get the number of processes
    int numberOfProcesses;
    MPI_Comm_size(MPI_COMM_WORLD, &numberOfProcesses);

    // Get the rank of the process
    int processRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &processRank);

        // use process rank int threadNumber = *thread;
    int start_row, end_row;
    int n = (arrayLength - 2) / numberOfProcesses;
    if (processRank == numberOfProcesses - 1) {
        start_row = processRank * n + 1;
        end_row = arrayLength - 2;
    } else {
        start_row = processRank * n + 1;
        end_row = start_row + n - 1;
    }

    printf("I am process %d, and I am averaging rows %d to %d\n", processRank, start_row, end_row);


    double arr[arrayLength*arrayLength];
    if(processRank == 0){
        printf("Root process now setting up array\n");
        for (int i = 0; i < 100; ++i)
        {
            arr[i] = i;
        }
    }


    if(processRank == 0){
        printf("I am the root process. Total number of processes: %d\n", numberOfProcesses);
        int id;
        for(id = 1; id < numberOfProcesses; id++){
            printf("Sending %d\n", id);
            MPI_Send(arr, arrayLength*arrayLength, MPI_DOUBLE, id, 0, MPI_COMM_WORLD);
        }
    }else{
        MPI_Recv(arr, arrayLength*arrayLength, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);


        printf("Process %d averaging between index %d to %d \n", processRank, start_row*arrayLength, end_row*arrayLength);




    }
    return 0;
}
