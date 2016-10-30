#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <pthread.h>

int arrayLength;
double **myArray;
double precision;
int numberOfThreads = 1;
pthread_mutex_t lock;
pthread_t *threads;

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
            printf("%.10f,", myArray[i][j]);
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

int getNumberOfLinesInFile(char *fileName){
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

    return numberOfLines;
}

double *readFile(char *fileName) {
    FILE *input;
    input = fopen(fileName, "r"); // reopen file to reset ptr
    int numberOfLines = getNumberOfLinesInFile(fileName);

    double *arr = malloc(numberOfLines * sizeof(int));
    int i;
    int temp;

    for (i = 0; i < 100; ++i)
    {
        // Skip whitespace
        do {
            temp = fgetc(input);
        } while (temp == ' ' || temp == '\n');

        // Store number
        arr[i] = temp - '0';
    }
    return arr;
}

void useFile(){
    double *fileRead;

    fileRead = readFile("/Users/jamestreasure/ClionProjects/Parallel/testArray.txt");

    myArray = malloc(arrayLength * sizeof(double *));

    for (int i = 0; i < arrayLength; ++i) {
        myArray[i] = malloc(arrayLength * sizeof(double));
    }

    for (int i = 0; i < arrayLength; ++i) {
        for (int j = 0; j < arrayLength; ++j) {
            printf("Value is %f\n",fileRead[i*arrayLength+j]);
            myArray[i][j] = fileRead[i*arrayLength+j];
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

        if (isPrecisionMet(myArray, tempArray, arrayLength)) {
            ended = 1;
        }

        printf("Thread %d finished averaging %d to %d\n", thrNum, start_row, end_row);


        for (int i = 0; i < arrayLength; ++i) {
            for (int j = 0; j < arrayLength; ++j) {
                tempArray[i][j] = myArray[i][j];
            }
        }
    }
}

int main(int argc, char **argv) {
    clock_t begin = clock();
    srand(time(NULL));

    arrayLength = 100;
    precision = 0.001;

    //useFile();
    setupArray(arrayLength);
    print2DArray(arrayLength,myArray);

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
    printf("%f\n", time_spent);

    print2DArray(arrayLength,myArray);

    return 0;
}
