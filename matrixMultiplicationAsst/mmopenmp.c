
#include <stdio.h>

#include <stdlib.h>

#include <time.h>

#include <omp.h>



#define M 99

#define K 99

#define N 99



int A[M][K];

int B[K][N];

int C[M][N];



int main() {

    double start, end, elapsed;

    int num_procs = omp_get_num_procs();

    

    printf("Number of processors available: %d\n", num_procs);

    // Assign random values to matrices A and B

    for (int i = 0; i < M; i++) {

        for (int j = 0; j < K; j++) {

            A[i][j] = rand() % 10;

            //printf("%d ", A[i][j]);

        }

        //printf("\n");

    }

    printf("\n");

    for (int i = 0; i < K; i++) {

        for (int j = 0; j < N; j++) {

            B[i][j] = rand() % 10;

            //printf("%d ", B[i][j]);

        }

        //printf("\n");

    }

    printf("\n Done with initialization \n \n");



    start = omp_get_wtime(); // start timer

    // Compute matrix product in parallel

    

    #pragma omp parallel for

    for (int i = 0; i < M; i++) {

        for (int j = 0; j < N; j++) {

            int sum = 0;

            for (int k = 0; k < K; k++) {

                sum += A[i][k] * B[k][j];

            }

            C[i][j] = sum;

        }

    }

    // Measure execution time

    end = omp_get_wtime(); // stop timer

    elapsed = end - start; // calculate elapsed time



    printf("OpenMP Execution time: %f seconds\n\n", elapsed);

    return 0;

}



