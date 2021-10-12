#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

unsigned int number_triangles;

// Structure used to represent a scene by its parameters
typedef struct scene_params {
    // Object color
    unsigned char object_r[4];
    unsigned char object_g[4];
    unsigned char object_b[4];

    // Background color
    unsigned char background_r[4];
    unsigned char background_g[4];
    unsigned char background_b[4];

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

/**
 * Load the data of the triangles into an array of float
 * @param fp The file containing the data
 * @return An array of 9 * T floats
 */
float *load_triangles(FILE *fp) {
    // Get the number of triangles
    fread(&number_triangles, sizeof(unsigned int), 1, fp); // Read the number of triangles
    printf("\nNumber of triangles: %i\n", number_triangles);

    // Allocate space
    float *triangle_vertices = malloc(9 * number_triangles * sizeof(float));

    // Populate the array
    fread(triangle_vertices, sizeof(float), 9 * number_triangles, fp);

    // Close the file
    fclose(fp);

    return triangle_vertices;
}


/**
 * Load the scene file into a structure that represents it
 * @param fp The parameters file
 * @param scene The structure in which the params are saved
 */
void load_scene_params(FILE *fp, scene_struct *scene) {
    fscanf(fp, "%s %s %s", scene->object_r, scene->object_g, scene->object_b);
    fscanf(fp, "%s %s %s", scene->background_r, scene->background_g, scene->background_b);
    fscanf(fp, "%f %f %f", &scene->camera_x, &scene->camera_y, &scene->camera_z);
    fscanf(fp, "%i %i", &scene->height, &scene->width);
    fscanf(fp, "%f", &scene->vertical_fow);
    fscanf(fp, "%f %f", &scene->near_clip, &scene->far_clip);

    // Close the file
    fclose(fp);
}

void print_scene(scene_struct *scene) {
    printf("\nScene parameters\n");
    printf("----------------\n");
    printf("Object color: %s %s %s\n", scene->object_r, scene->object_g, scene->object_b);
    printf("Background color: %s %s %s\n", scene->background_r, scene->background_g, scene->background_b);
    printf("Camera position: %.2f %.2f %.2f\n", scene->camera_x, scene->camera_y, scene->camera_z);
    printf("Height and width: %i %i\n", scene->height, scene->width);
    printf("Vertical fow: %.6f\n", scene->vertical_fow);
    printf("Near and far clip: %.2f %.2f\n", scene->near_clip, scene->far_clip);
}

/**
 * Give the index of the beginning of triangle
 * @param triangle_idx the index of the triangle
 */
int get_triangle(int triangle_idx) {
    assert(triangle_idx >= 0 && triangle_idx < number_triangles - 1);
    return triangle_idx *= 9;
}

/**
 * Give the index of a vertex in a given triangle
 * @param triangle_idx The index of the triangle
 * @param point_idx The vertex (a,b,c)
 */
int get_point(int triangle_idx, int point_idx) {
    assert(point_idx >= 0 && point_idx < 3);
    return get_triangle(triangle_idx) + point_idx * 3;
}

/**
 * Give the index of a coordinate in a given vertex from a given triangle
 * @param triangle_idx The index of the triangle
 * @param point_idx The vertex of choice (a,b,c)
 * @param coord_idx The coordinate of choice (x,y,z)
 */
int get_coord(int triangle_idx, int point_idx, int coord_idx){
    assert(coord_idx >= 0 && coord_idx < 3);
    return get_point(triangle_idx, point_idx) + coord_idx;
}


int main(int argc, char **argv) {

    // Parameters check
    if (argc != 4) {
        printf("-- Missing parameters --\n"
               "1. Data file (e.g. bun.dat)\n2. Parameter file (e.g. param_bun.txt)\n3. Output file (e.g. bun.ppm)");
        return 1;
    }

    FILE *dataFile = fopen(argv[1], "r");
    if (dataFile == NULL) {
        printf("Data file not readable");
        return 1;
    }

    FILE *paramsFile = fopen(argv[2], "r");
    if (paramsFile == NULL) {
        printf("Parameters file not readable");
        fclose(dataFile);
        return 1;
    }

    // Store the triangles
    float *triangles = load_triangles(dataFile);

    // Store the scene parameters
    scene_struct scene;
    load_scene_params(paramsFile, &scene);

    // Shading parameters


    // Terminate program
    free(triangles);

    return 0;
}