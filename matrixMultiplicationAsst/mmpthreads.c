#include <stdio.h>

#include <stdlib.h>

#include <time.h>

#include <pthread.h>



#define M 99

#define N 99

#define K 99

#define NUM_THREADS 9



int A[M][K];

int B[K][N];

int C[M][N];



struct ThreadArgs {

    int start;

    int end;

};



void *multiply(void *arg) {

    struct ThreadArgs *args = arg;

    int start = args->start;

    int end = args->end;

    int i, j, k;



    for (i = start; i < end; i++) {

        for (j = 0; j < N; j++) {

            for (k = 0; k < K; k++) {

                C[i][j] += A[i][k] * B[k][j];

            }

        }

    }



    pthread_exit(NULL);

}



int main() {

    double elapsed_time;

    struct timespec start, end;

    int i, j;

    pthread_t threads[NUM_THREADS];

    struct ThreadArgs threadArgs[NUM_THREADS];

    void *status;



    // Initialize matrices A and B

    for (i = 0; i < M; i++) {

        for (j = 0; j < K; j++) {

            A[i][j] = rand() % 10 + 1;

        }

    }



    for (i = 0; i < K; i++) {

        for (j = 0; j < N; j++) {

            B[i][j] = rand() % 10 + 1;

       }

     
    }

    	
    printf("\n Done with initialization \n \n");


    clock_gettime(CLOCK_MONOTONIC, &start); // start timer



    // Create threads

    for (i = 0; i < NUM_THREADS; i++) {

        int start = i * (M / NUM_THREADS);

        int end = (i + 1) * (M / NUM_THREADS);

        threadArgs[i].start = start;

        threadArgs[i].end = end;

        pthread_create(&threads[i], NULL, multiply, (void *)&threadArgs[i]);

    }



    // Wait for threads to finish

    for (i = 0; i < NUM_THREADS; i++) {

        pthread_join(threads[i], &status);

    }

    // Measure execution time

    // Stop the timer

    clock_gettime(CLOCK_MONOTONIC, &end);


    // Calculate the elapsed time

    elapsed_time = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;

    printf("pthread Elapsed time: %f seconds\n\n", elapsed_time);

    return 0;

}

