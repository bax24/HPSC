#include <stdio.h>
#include <stdlib.h>

// Structure used to define a point in 3D-space
typedef struct point {
    float coordinates[3]; // List of 3 floats
} point_struct;

// Structure used to represent a triangle in 3D-space
typedef struct triangle {
    point_struct points[3]; // List of 3 points
} triangle_struct;

// Structure used to represent a scene by its parameters
typedef struct scene_params {
    // Object color
    unsigned char object_r;
    unsigned char object_g;
    unsigned char object_b;

    // Background color
    unsigned char background_r;
    unsigned char background_g;
    unsigned char background_b;

    // Camera coordinates
    float camera_x;
    float camera_y;
    float camera_z;

    // Height and width of the 2D image
    unsigned int height;
    unsigned int width;

    // Vertical field of view
    float vertical_fow;

    // Far clip and near clip point position
    float far_clip;
    float near_clip;

} scene_struct;

triangle_struct *load_triangles(FILE *fp) {
    unsigned int number_triangles;
    fread(&number_triangles, sizeof(unsigned int), 1, fp); // Read the number of triangles
    printf("There are %u triangles\n", number_triangles); // ... and print them for readability

    triangle_struct *triangles = malloc(number_triangles * sizeof(triangle_struct)); // Allocate space for them

    for (int i = 0; i < number_triangles; i++) { // Load those triangles into the allocated list
        if (i % 1000 == 0) {
            printf("Triangles loaded: %d\n", i);
        }

        triangle_struct triangle;

        // Populate a triangle
        for (int j = 0; j < 3; j++) {
            point_struct point;

            // Populate a point
            for (int k = 0; k < 3; k++) {
                float value;
                fread(&value, sizeof(float), 1, fp);
                point.coordinates[k] = value;
            }
            triangle.points[j] = point;
        }
        triangles[i] = triangle;
    }

    return triangles;
}

void load_scene_params(FILE *fp, scene_struct *scene) {

    printf("Le voici: %uc", scene->object_r);
}

int main(int argc, char **argv) {

    // Parameters check
    if (argc != 4) {
        printf("-- Missing parameters --\n"
               "1. Data file (e.g. bun.dat)\n2. Parameter file (e.g. param_bun.txt)\n3. Output file (e.g. bun.ppm)");
        return 1;
    }

    FILE *dataFile = fopen(argv[1], "r");
    FILE *paramsFile = fopen(argv[2], "r");

    if (!dataFile) {
        printf("Data file not readable");
        return 1;
    }

    if (!paramsFile) {
        printf("Parameters file not readable");
        return 1;
    }

//    triangle_struct *triangles = load_triangles(dataFile);
    scene_struct scene;
    load_scene_params(paramsFile, &scene);

    // Terminate program
//    free(triangles);
    fclose(dataFile);
    fclose(paramsFile);

    return 0;
}