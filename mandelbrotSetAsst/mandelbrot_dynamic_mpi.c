#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define IMAGE_WIDTH 1920
#define IMAGE_HEIGHT 1080
#define REAL_MAX 2
#define REAL_MIN -2
#define IMAG_MAX 2
#define IMAG_MIN -2
#define MAX_ITERATIONS 256

unsigned char calculate_mandelbrot(int pixel_x, int pixel_y);
void save_ppm_image(const char* filename, unsigned char* image, int width, int height);

int main(int argc, char* argv[]) {
    int num_processes, rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int rows_per_process = IMAGE_HEIGHT / num_processes;
    int start_row = rank * rows_per_process;
    int end_row = (rank + 1) * rows_per_process;

    unsigned char* image = (unsigned char*)malloc(IMAGE_WIDTH * rows_per_process * sizeof(unsigned char));
    if (image == NULL) {
        fprintf(stderr, "Error: Memory allocation failed.\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Start measuring time
    double start_time = MPI_Wtime();

    // Create a buffer to gather results from all processes
    unsigned char* recv_buffer = NULL;
    if (rank == 0) {
        recv_buffer = (unsigned char*)malloc(IMAGE_WIDTH * IMAGE_HEIGHT * sizeof(unsigned char));
        if (recv_buffer == NULL) {
            fprintf(stderr, "Error: Memory allocation failed.\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    // Dynamic task assignment loop
    while (1) {
        // Master process initiates communication with workers
        if (rank == 0) {
            for (int i = 1; i < num_processes; i++) {
                // Send a region request to worker i
                MPI_Send(&start_row, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            }
        } else {
            int requesting_worker_rank = -1;
            // Receive a region request or termination signal from the master
            MPI_Recv(&requesting_worker_rank, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            if (requesting_worker_rank == -1) {
                // Terminate worker processes
                break;
            }
        }

        // Compute the assigned region
        for (int pixel_y = start_row; pixel_y < end_row; pixel_y++) {
            for (int pixel_x = 0; pixel_x < IMAGE_WIDTH; pixel_x++) {
                unsigned char pixel_value = calculate_mandelbrot(pixel_x, pixel_y);
                image[(pixel_y - start_row) * IMAGE_WIDTH + pixel_x] = pixel_value;
            }
        }

        // Gather results from all processes into the recv_buffer
        MPI_Gather(image, IMAGE_WIDTH * rows_per_process, MPI_UNSIGNED_CHAR,
                   recv_buffer, IMAGE_WIDTH * rows_per_process, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

        // Master process processes the gathered results and saves the final image
        if (rank == 0) {
            // Process the gathered results and save the final image as PPM

            // Create a grayscale image array to hold the final result
            unsigned char* final_image = (unsigned char*)malloc(IMAGE_WIDTH * IMAGE_HEIGHT * sizeof(unsigned char));
            if (final_image == NULL) {
                fprintf(stderr, "Error: Memory allocation failed.\n");
                MPI_Abort(MPI_COMM_WORLD, 1);
            }

            // Copy the gathered results into the final image
            for (int i = 0; i < num_processes; i++) {
                int row_offset = i * rows_per_process;
                for (int pixel_y = 0; pixel_y < rows_per_process; pixel_y++) {
                    for (int pixel_x = 0; pixel_x < IMAGE_WIDTH; pixel_x++) {
                        final_image[(row_offset + pixel_y) * IMAGE_WIDTH + pixel_x] =
                            recv_buffer[(i * rows_per_process + pixel_y) * IMAGE_WIDTH + pixel_x];
                    }
                }
            }

            // Save the final grayscale image as a PPM file
            save_ppm_image("mandelbrot3.ppm", final_image, IMAGE_WIDTH, IMAGE_HEIGHT);

            // Cleanup
            free(final_image);

            // End measuring time
            double end_time = MPI_Wtime();
            printf("Time taken to calculate Mandelbrot set: %lf seconds\n", end_time - start_time);

            // Signal worker processes to terminate
            for (int i = 1; i < num_processes; i++) {
                int termination_signal = -1;
                // Send termination signal to worker i
                MPI_Send(&termination_signal, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
            }
            break;
        }
    }

    free(image);
    if (rank == 0) {
        free(recv_buffer);
    }

    MPI_Finalize();
    return 0;
}

unsigned char calculate_mandelbrot(int pixel_x, int pixel_y) {
    double real = 0; // Real part of complex number (c)
    double imag = 0; // Imaginary part of complex number (c)

    double real_coordinate = REAL_MIN + (pixel_x * ((REAL_MAX - REAL_MIN) / (IMAGE_WIDTH * 1.0)));
    double imag_coordinate = IMAG_MIN + (pixel_y * ((IMAG_MAX - IMAG_MIN) / (IMAGE_HEIGHT * 1.0)));

    int iterations = 0;
    double real_squared = 0;
    double imag_squared = 0;

    while (real_squared + imag_squared <= 4.0 && iterations < MAX_ITERATIONS) {
        imag = 2 * real * imag + imag_coordinate;
        real = real_squared - imag_squared + real_coordinate;
        real_squared = real * real;
        imag_squared = imag * imag;
        iterations++;
    }

    // Calculate grayscale value based on the number of iterations
    unsigned char grayscale_value = (unsigned char)(255 * iterations / (double)MAX_ITERATIONS);
    return grayscale_value;
}

void save_ppm_image(const char* filename, unsigned char* image, int width, int height) {
    FILE* file = fopen(filename, "wb");
    if (!file) {
        fprintf(stderr, "Error: Cannot open file for writing.\n");
        return;
    }

    // Write PPM header
    fprintf(file, "P5\n%d %d\n255\n", width, height);

    // Write image data
    fwrite(image, sizeof(unsigned char), width * height, file);

    fclose(file);
}
