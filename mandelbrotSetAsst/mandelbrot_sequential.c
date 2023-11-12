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

int main() {
    unsigned char* image = (unsigned char*)malloc(IMAGE_WIDTH * IMAGE_HEIGHT * sizeof(unsigned char));
    if (image == NULL) {
        fprintf(stderr, "Error: Memory allocation failed.\n");
        return 1;
    }

    // Start measuring time
    clock_t start_time = clock();

    // Compute the Mandelbrot set
    for (int pixel_y = 0; pixel_y < IMAGE_HEIGHT; pixel_y++) {
        for (int pixel_x = 0; pixel_x < IMAGE_WIDTH; pixel_x++) {
            unsigned char pixel_value = calculate_mandelbrot(pixel_x, pixel_y);
            image[pixel_y * IMAGE_WIDTH + pixel_x] = pixel_value;
        }
    }

    // End measuring time
    clock_t end_time = clock();
    double elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("Time taken to calculate Mandelbrot set: %lf seconds\n", elapsed_time);

    // Save the final grayscale image as a PPM file
    save_ppm_image("mandelbrot_sequential.ppm", image, IMAGE_WIDTH, IMAGE_HEIGHT);

    // Cleanup
    free(image);

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
