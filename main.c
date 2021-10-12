#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

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
float *load_triangles(FILE *fp, unsigned int number_triangles) {
    float *triangles = malloc(9 * number_triangles * sizeof(float)); // Allocate space
    fread(triangles, sizeof(float), 9 * number_triangles, fp); // Populate the array
    fclose(fp); // Populate the array
    return triangles;
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
 * @param number_triangles The number of triangles
 */
int get_triangle(int triangle_idx) {
    return triangle_idx *= 9;
}

/**
 * Give the index of a vertex in a given triangle
 * @param triangle_idx The index of the triangle
 * @param point_idx The vertex (a,b,c)
 * @param number_triangles The number of triangles
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
 * @param number_triangles The number of triangles
 */
int get_coord(int triangle_idx, int point_idx, int coord_idx) {
    assert(coord_idx >= 0 && coord_idx < 3);
    return get_point(triangle_idx, point_idx) + coord_idx;
}

/**
 * Compute the vector AB
 * @param triangles The list of triangles
 * @param components Where to store the components
 * @param point_a The point a
 * @param point_b The point b
 * @return A list of 3 floats representing the components of the vector
 */
void compute_vector(const float *triangles, float *components, const unsigned int point_a, const unsigned int point_b) {
    for (int i = 0; i < 3; i++) {
        components[i] = triangles[point_b + i] - triangles[point_a + i];
    }
}

/**
 * Computes the Euclidian norm of a 3D-vector
 * @param vector
 * @return
 */
double compute_norm(float* vector){
    return sqrt(pow(vector[0], 2) + pow(vector[1], 2) + pow(vector[2], 2));
}

void vector_product(float *result, float *vector_a, float *vector_b){

}

void compute_shading_params(const float *triangles, float *shading_params,
                            unsigned int number_triangles, scene_struct *scene) {
    for (int i = 0; i < number_triangles; i++) {
        // === Compute the normal of the triangle ===

        // Compute AB and AC
        float vector_ab[3] = {0, 0, 0}, vector_ac[3] = {0, 0, 0};
        compute_vector(triangles, vector_ab, get_point(i, 0), get_point(i, 1));
        compute_vector(triangles, vector_ac, get_point(i, 0), get_point(i, 2));

        // Compute ||AB|| * ||AC||
        double norm = compute_norm(vector_ab) * compute_norm(vector_ac);

        // Change of coordinates

        // Compute shading
    }

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

    // Read the number of triangles
    unsigned int number_triangles;
    fread(&number_triangles, sizeof(unsigned int), 1, dataFile);
    printf("\nNumber of triangles: %i\n", number_triangles);

    // Store the triangles
    float *triangles = load_triangles(dataFile, number_triangles);

    // Store the scene parameters
    scene_struct scene;
    load_scene_params(paramsFile, &scene);

    // Shading parameters
    float *shading_params = malloc(number_triangles * sizeof(float));
    compute_shading_params(triangles, shading_params, number_triangles, &scene);

    // Terminate program
    free(triangles);

    return 0;
}