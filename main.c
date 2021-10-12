#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

const double lightbeam_vector[3] = {-0.57735026919, -0.57735026919, -0.57735026919};

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
    float camera_coordinates[3];

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

    float camera_x, camera_y, camera_z;
    fscanf(fp, "%f %f %f", &camera_x, &camera_y, &camera_z);
    scene->camera_coordinates[0] = camera_x;
    scene->camera_coordinates[1] = camera_y;
    scene->camera_coordinates[2] = camera_z;

    fscanf(fp, "%i %i", &scene->height, &scene->width);
    fscanf(fp, "%f", &scene->vertical_fow);
    fscanf(fp, "%f %f", &scene->near_clip, &scene->far_clip);

    // Close the file
    fclose(fp);
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
 * Computes the Euclidian norm of a vector
 * @param vector The vector
 * @param size The size of the vector
 * @return A scalar
 */
double compute_norm(float *vector, size_t size) {
    double result = 0;
    for (size_t i = 0; i < size; i++) {
        result += pow(vector[i], 2);
    }
    return sqrt(result);
}

/**
 * Compute the vector product of 3D vectors
 * @param result The resulting vector
 * @param u The first vector
 * @param v The second vector
 */
void vector_product(float *result, const float *u, const float *v) {
    result[0] = u[1] * v[2] - u[2] * v[1];
    result[1] = u[2] * v[0] - u[0] * v[2];
    result[2] = u[0] * v[1] - u[1] * v[0];
}

/**
 * Compute the scalar product of 2 vectors
 * @param u The first vector
 * @param v The second vector
 * @param size The size of the vectors
 * @return The scalar product
 */
float scalar_product(const float *u, const float *v, size_t size) {
    float result = 0;
    for (size_t i = 0; i < size; i++) {
        result += u[i] + v[i];
    }
    return result;
}

/**
 * Scale a given vector by a scalar value
 * @param v The vector
 * @param scale The scalar
 * @param size The size of the vector
 */
void scale_vector(float *v, float scale, size_t size) {
    for (size_t i = 0; i < size; i++) {
        v[i] *= scale;
    }
}

void compute_shading_params(const float *triangles, float *shading_params,
                            const unsigned int number_triangles) {

    for (int i = 0; i < (int) number_triangles; i++) {
        // ######################################
        // # Compute the normal of the triangle #
        // ######################################

        // Declare variables
        float vector_ab[3] = {0, 0, 0}, vector_ac[3] = {0, 0, 0}, vector_normal[3] = {0, 0, 0};
        float shading_value = 0;

        // Compute both AB and AC
        compute_vector(triangles, vector_ab, get_point(i, 0), get_point(i, 1));
        compute_vector(triangles, vector_ac, get_point(i, 0), get_point(i, 2));

        // Compute the normal vector
        vector_product(vector_normal, vector_ab, vector_ac);
        double norm = compute_norm(vector_ab, 3) * compute_norm(vector_ac, 3);
        scale_vector(vector_normal, 1.0 / norm, 3);

        // #################################
        // # Compute the shading parameter #
        // #################################
        float result = scalar_product(vector_normal, (float *) lightbeam_vector, 3);
        if (result > 0)
            shading_value = result;

        shading_params[i] = shading_value;
    }
}

void compute_world_to_view_coordinate(float *triangles, scene_struct *scene, const unsigned int number_triangles) {

    // Declare variables
    float E_u[3], E_v[3], *E_w;
    float camera_distance;
    float ground_plan_distance;

    // Compute the distances
    camera_distance = compute_norm(scene->camera_coordinates, 3);
    float vector_tmp[2] = {scene->camera_coordinates[0], scene->camera_coordinates[1]};
    ground_plan_distance = compute_norm(vector_tmp, 2);

    // Compute E_u
    E_u[0] = scene->camera_coordinates[2];
    E_u[1] = 0;
    E_u[2] = -scene->camera_coordinates[0];
    scale_vector(E_u, 1 / ground_plan_distance, 3);

    // Compute E_v
    E_v[0] = -scene->camera_coordinates[0] * scene->camera_coordinates[1];
    E_v[1] = pow(scene->camera_coordinates[2], 2) + pow(scene->camera_coordinates[0], 2);
    E_v[2] = -scene->camera_coordinates[1] * scene->camera_coordinates[2];
    scale_vector(E_v, 1 / camera_distance * ground_plan_distance, 3);

    // Compute E_w
    E_w = scene->camera_coordinates;
    scale_vector(E_w, 1 / camera_distance, 3);

    // Save the new system of coordinates for convenience
    float *new_system[3] = {E_u, E_v, E_w};

    // Compute to vector from the origin to the camera position
    float *vector_om = scene->camera_coordinates;

    // Compute the world to view matrix
    float w_to_v[3][4];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            w_to_v[i][j] = new_system[i][j];
        }
        w_to_v[i][3] = - scalar_product(vector_om, new_system[i], 3);
    }

    // Iterate over the triangles
    for (unsigned int i = 0; i < number_triangles; i++) {

        // Iterate over the points of a triangle
        for (int j = 0; j < 3; j++) {
            float vector_column[4] = {get_coord(i, j, 0), get_coord(i, j, 1), get_coord(i, j, 2), 1.0};

            // Iterate over the coordinates of the point
            for(int k = 0; k < 3; k++){
                int index = get_coord(i, j, k);
                triangles[index] = scalar_product(w_to_v[k], vector_column, 4);
            }
        }
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
    compute_shading_params(triangles, shading_params, number_triangles);

    // Change of coordinates
    compute_world_to_view_coordinate(triangles, &scene, number_triangles);

    // Terminate program
    free(triangles);
    free(shading_params);

    return 0;
}