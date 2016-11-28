#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <ctype.h>
#include <mpi.h>

#define TRUE  (1==1)
#define FALSE (!TRUE)

int arrayLength = 5;
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
int isNotAnEdge(int index) {
    return !((index % arrayLength == 0) || ((index + 1) % arrayLength == 0));
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
                //if (isNotAnEdge(arrayLength, i, j)) {
                double above = tempArray[i - 1][j];
                double below = tempArray[i + 1][j];
                double left = tempArray[i][j - 1];
                double right = tempArray[i][j + 1];
                myArray[i][j] = (above + below + left + right) / 4;

                double currentPrecision = fabs(myArray[i][j] - tempArray[i][j]);
                if (currentPrecision > precision) {
                    precisionMet = FALSE;
                }
                //}
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

    srand(time(NULL));

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

    double arr[arrayLength * arrayLength];
    if (processRank == 0) {
        printf("Root process now setting up array\n");
        printf("------------------------------------------------------\n");
        for (int i = 0; i < arrayLength * arrayLength; ++i) {
            arr[i] = generateRandomNumber();
        }

        for (int i = 0; i < arrayLength; i++) {
            for (int j = 0; j < arrayLength; j++) {
                printf("%f,", arr[i * arrayLength + j]);
            }
            printf("\n");
        }
        printf("------------------------------------------------------\n");
    }

    // everyone calls bcast and takes data from the root's buffer. No need to receive.
    MPI_Bcast(arr, arrayLength * arrayLength, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // printf("Process %d averaging between index %d to %d \n", processRank, start_row*arrayLength, end_row*arrayLength+arrayLength);
    double tempArray[arrayLength * arrayLength];
    for (int i = 0; i < arrayLength * arrayLength; ++i) {
        tempArray[i] = arr[i];
    }


    for (int i = start_row * arrayLength;
         i < end_row * arrayLength + arrayLength; ++i) {
        if (isNotAnEdge(i)) {
            //printf("%d is not an edge\n", i);
            double above = arr[i - arrayLength];
            double below = arr[i + arrayLength];
            double left = arr[i - 1];
            double right = arr[i + 1];
            //printf("Value %f with above %f, below %f, left %f and right %f\n", arr[i], above, below, left, right );
            tempArray[i] = (above + below + left + right) / 4;
        } else {
            //printf("%d IS AN EDGE and value is %f\n", i, arr[i] );
        }
    }


    if (processRank == 0) {
        //replace rows with updated rows
        for (int i = start_row * arrayLength;
             i < end_row * arrayLength + arrayLength; ++i) {
            arr[i] = tempArray[i];
        }
        printf("Updated rows from root process\n");

        for (int i = 1; i < numberOfProcesses; ++i) {
            double receiveArray[arrayLength * arrayLength];
            int startEndArray[2];
            MPI_Recv(startEndArray, 2, MPI_INT, i, 0,
                     MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
            MPI_Recv(receiveArray, arrayLength * arrayLength, MPI_DOUBLE, i, 0,
                     MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);

            printf("Start is %d, end is %d\n", startEndArray[0],
                   startEndArray[1]);

            for (int i = startEndArray[0] * arrayLength;
                 i < startEndArray[1] * arrayLength + arrayLength; ++i) {
                arr[i] = receiveArray[i];
            }
            printf("Finished\n");
            printf("------------------------------------------------------\n");
            for (int i = 0; i < arrayLength; i++) {
                for (int j = 0; j < arrayLength; j++) {
                    printf("%f,", arr[i * arrayLength + j]);
                }
                printf("\n");
            }
            printf("------------------------------------------------------\n");

        }
    } else {
        int startEndArray[2] = {start_row, end_row};
        MPI_Send(startEndArray, 2, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(tempArray, arrayLength * arrayLength, MPI_DOUBLE, 0, 0,
                 MPI_COMM_WORLD);
    }

    return 0;
}
