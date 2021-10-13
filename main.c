#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <stdbool.h>

#define MIN(a, b) ( a < b ) ? a : b
#define MAX(a, b) ( a > b ) ? a : b

const double lightbeam_vector[3] = {-0.57735026919, -0.57735026919, -0.57735026919};

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

typedef struct pixel {
    unsigned char r;
    unsigned char g;
    unsigned char b;
    float depth;
} pixel;

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
    unsigned int color_r, color_g, color_b;

    fscanf(fp, "%i %i %i", &color_r, &color_g, &color_b);
    scene->object_r = color_r;
    scene->object_g = color_g;
    scene->object_b = color_b;

    fscanf(fp, "%i %i %i", &color_r, &color_g, &color_b);
    scene->background_r = color_r;
    scene->background_g = color_g;
    scene->background_b = color_b;

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
 */
int get_triangle(int triangle_idx) {
    triangle_idx *= 9;
    return triangle_idx;
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
double euclidian_norm(float *vector, size_t size) {
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
void compute_vector_product(float *result, const float *u, const float *v) {
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

    for (unsigned int i = 0; i < number_triangles; i++) {
        // Declare variables
        float vector_ab[3] = {0, 0, 0}, vector_ac[3] = {0, 0, 0}, vector_normal[3] = {0, 0, 0};
        float shading_value = 0;

        // Compute both AB and AC
        compute_vector(triangles, vector_ab, get_point(i, 0), get_point(i, 1));
        compute_vector(triangles, vector_ac, get_point(i, 0), get_point(i, 2));

        // Compute the normal vector
        compute_vector_product(vector_normal, vector_ab, vector_ac);
        double norm = euclidian_norm(vector_ab, 3) * euclidian_norm(vector_ac, 3);
        scale_vector(vector_normal, 1.0 / norm, 3);

        // Compute the shading parameter
        float result = scalar_product(vector_normal, (float *) lightbeam_vector, 3);
        if (result > 0)
            shading_value = result;

        shading_params[i] = shading_value;
    }
}

void compute_world_to_view_coordinates(float *triangles, scene_struct *scene, const unsigned int number_triangles) {
    // Declare variables
    float E_u[3], E_v[3], *E_w;
    float camera_distance;
    float ground_plan_distance;

    // Compute the distances
    camera_distance = euclidian_norm(scene->camera_coordinates, 3);
    ground_plan_distance = euclidian_norm(scene->camera_coordinates, 2);

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
        w_to_v[i][3] = -scalar_product(vector_om, new_system[i], 3);
    }

    // Iterate over the triangles
    for (unsigned int i = 0; i < number_triangles; i++) {

        // Iterate over the points of a triangle
        for (int j = 0; j < 3; j++) {
            float vector_column[4] = {triangles[get_coord(i, j, 0)], triangles[get_coord(i, j, 1)],
                                      triangles[get_coord(i, j, 2)], 1.0};

            // Iterate over the coordinates of the point
            for (int k = 0; k < 3; k++) {
                int index = get_coord(i, j, k);
                triangles[index] = scalar_product(w_to_v[k], vector_column, 4);
            }
        }
    }
}

void compute_view_to_projection_coordinates(float *triangles, scene_struct *scene, unsigned number_triangles) {
    // Compute S_y and S_x
    float S_y = 1 / tan(scene->vertical_fow / 2);
    float S_x = S_y / (scene->width / scene->height);

    // Compute the view to projection matrix
    float v_to_p[3][4];
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 4; j++)
            v_to_p[i][j] = 0;
    }
    v_to_p[0][0] = S_x;
    v_to_p[1][1] = S_y;
    v_to_p[2][2] = (scene->near_clip + scene->far_clip) / (scene->near_clip - scene->far_clip);
    v_to_p[2][3] = -1;

    // Iterate over the number of triangles
    for (unsigned int i = 0; i < number_triangles; i++) {

        // Iterate over the points of a triangle
        for (unsigned int j = 0; j < 3; j++) {
            float vector_column[4] = {get_coord(i, j, 0), get_coord(i, j, 1), get_coord(i, j, 2), 1.0};
            float scaling_factor = (scene->near_clip - scene->far_clip) / (2 * scene->near_clip * scene->far_clip *
                                                                           get_coord(i, j, 2));

            // Iterate over the coordinates of a point
            for (unsigned int k = 0; k < 3; k++) {
                int index = get_coord(i, j, k);
                triangles[index] = scalar_product(v_to_p[k], vector_column, 4) * scaling_factor;
            }
        }
    }
}

void init_pixels(pixel *pixels, scene_struct *scene) {
    unsigned int number_pixels = scene->width * scene->height;
    for (unsigned int i = 0; i < number_pixels; i++) {
        pixels[i].depth = 1;
        pixels[i].r = scene->background_r;
        pixels[i].g = scene->background_g;
        pixels[i].b = scene->background_b;
    }
}

bool is_inside(float pxl_x, float pxl_y, unsigned int trgl_idx, float *triangles) {
    float vector1[3] = {0, 0, 0}, vector2[3] = {0, 0, 0}, vector3[3] = {0, 0, 0};
    float point_p[2] = {pxl_x, pxl_y};
    float vector_product1[3], vector_product2[3];

    // First inequality
    for (unsigned int i = 0; i < 2; i++) {
        vector1[i] = triangles[get_coord(trgl_idx, 1, i)] - triangles[get_coord(trgl_idx, 0, i)]; // AB
        vector2[i] = triangles[get_coord(trgl_idx, 2, i)] - triangles[get_coord(trgl_idx, 0, i)]; // AC
        vector3[i] = point_p[i] - triangles[get_coord(trgl_idx, 0, i)]; // AP
    }

    compute_vector_product(vector_product1, vector1, vector3);
    compute_vector_product(vector_product2, vector3, vector2);
    if (scalar_product(vector_product1, vector_product2, 3) < 0) {
        return false;
    }

    // Second inequality
    for (unsigned int i = 0; i < 2; i++) {
        vector1[i] = triangles[get_coord(trgl_idx, 2, i)] - triangles[get_coord(trgl_idx, 1, i)]; // BC
        vector2[i] = triangles[get_coord(trgl_idx, 0, i)] - triangles[get_coord(trgl_idx, 1, i)]; // BA
        vector3[i] = point_p[i] - triangles[get_coord(trgl_idx, 1, i)];  // BP
    }

    compute_vector_product(vector_product1, vector1, vector3);
    compute_vector_product(vector_product2, vector3, vector2);
    if (scalar_product(vector_product1, vector_product2, 3) < 0) {
        return false;
    }

    return true;
}

unsigned int get_pxl_x(float value, scene_struct *scene) {
    return 2 * (value + 0.5) / (scene->width - 1);
}

unsigned int get_pxl_y(float value, scene_struct *scene) {
    return -2 * (value + 0.5) / (scene->height + 1);
}

float compute_pixel_depth(const float pxl_x, const float pxl_y, float *triangles, const unsigned int triangle_idx) {
    float triangle_points[3][2], w_a, w_b, w_c;
    for (unsigned int i = 0; i < 3; i++) {
        for (unsigned int j = 0; j < 2; j++) {
            triangle_points[i][j] = triangles[get_coord(triangle_idx, i, j)];
        }
    }

    // Temporary values to improve readability; in the end -> (num1 - num2) / (denum1 - denum2)
    float num1, num2, denum1, denum2;

    num1 = (triangle_points[1][0] - triangle_points[0][0]) * (pxl_y - triangle_points[0][1]);
    num2 = (triangle_points[1][1] - triangle_points[0][1]) * (pxl_x - triangle_points[0][0]);
    denum1 = (triangle_points[1][0] - triangle_points[0][0]) * (triangle_points[2][1] - triangle_points[0][1]);
    denum2 = (triangle_points[1][1] - triangle_points[0][1]) * (triangle_points[2][0] - triangle_points[0][0]);
    w_a = (num1 - num2) / (denum1 - denum2);

    num1 = (triangle_points[0][0] - triangle_points[2][0]) * (pxl_y - triangle_points[2][1]);
    num2 = (triangle_points[0][1] - triangle_points[2][1]) * (pxl_x - triangle_points[2][0]);
    denum1 = (triangle_points[0][0] - triangle_points[2][0]) * (triangle_points[1][1] - triangle_points[2][1]);
    denum2 = (triangle_points[0][1] - triangle_points[2][1]) * (triangle_points[1][0] - triangle_points[2][0]);
    w_b = (num1 - num2) / (denum1 - denum2);

    w_c = 1 - w_a - w_b;

    float depth = w_a * triangles[get_coord(triangle_idx, 0, 2)];
    depth += w_b * triangles[get_coord(triangle_idx, 1, 2)];
    depth += w_c * triangles[get_coord(triangle_idx, 2, 2)];
    return depth;
}

void compute_pixels(float *triangles, pixel *pixels, scene_struct *scene,
                    const unsigned int number_triangles, float *shading_params) {

    // Iterate over the triangles
    for (unsigned int triangle_idx = 0; triangle_idx < number_triangles; triangle_idx++) {

        float min_x = FLT_MAX, max_x = -FLT_MAX;
        float min_y = FLT_MAX, max_y = -FLT_MAX;

        // Iterate over the points of the triangle to define a square of search
        for (unsigned int point_idx = 0; point_idx < 3; point_idx++) {
            min_x = MIN(min_x, triangles[get_coord(triangle_idx, point_idx, 0)]);
            min_y = MIN(min_y, triangles[get_coord(triangle_idx, point_idx, 1)]);
            max_x = MAX(max_x, triangles[get_coord(triangle_idx, point_idx, 0)]);
            max_y = MAX(max_y, triangles[get_coord(triangle_idx, point_idx, 1)]);
        }

        // Notation: pixel[i][j]
        // TODO: Something is wrong with this computation -> The values are too big
        unsigned int min_j = round((min_x * (scene->width - 1) / 2) - 0.5);
        unsigned int max_j = round((max_x * (scene->width - 1) / 2) - 0.5);
        unsigned int min_i = round((-min_y * (scene->height + 1) / 2) - 0.5);
        unsigned int max_i = round((-max_y * (scene->height + 1) / 2) - 0.5);

        // TODO: put a flag to improve computation (out of the triangle -> can skip (=continue) the line)
        for (unsigned int i = min_i; i < max_i; i++) {
            for (unsigned int j = min_j; j < max_j; j++) {
                float pxl_x = get_pxl_x(j, scene);
                float pxl_y = get_pxl_y(i, scene);

                if (!is_inside(pxl_x, pxl_y, triangle_idx, triangles))
                    continue;

                float pixel_depth = compute_pixel_depth(pxl_x, pxl_y, triangles, triangle_idx);
                if(pixel_depth < pixels[i % scene->width + j].depth){
                    pixels[i % scene->width + j].depth = pixel_depth;
                    pixels[i % scene->width + j].r = scene->object_r * shading_params[triangle_idx];
                    pixels[i % scene->width + j].g = scene->object_g * shading_params[triangle_idx];
                    pixels[i % scene->width + j].b = scene->object_b * shading_params[triangle_idx];
                }
            }
        }
    }
}

int main(int argc, char **argv) {

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

    // Store the triangles
    float *triangles = load_triangles(dataFile, number_triangles);

    // Store the scene parameters
    scene_struct scene;
    load_scene_params(paramsFile, &scene);

    // Shading parameters
    float *shading_params = malloc(number_triangles * sizeof(float));
    compute_shading_params(triangles, shading_params, number_triangles);

    printf("[LOG] Shading parameters computed\n");

    // Change of coordinates
    compute_world_to_view_coordinates(triangles, &scene, number_triangles);
    compute_view_to_projection_coordinates(triangles, &scene, number_triangles);

    printf("[LOG] Change of coordinates computed\n");

    // Allocate space for the pixels
    pixel *pixels = malloc(scene.height * scene.width * sizeof(pixel));
    init_pixels(pixels, &scene);

    // Compute pixel color
    compute_pixels(triangles, pixels, &scene, number_triangles, shading_params);

    printf("[LOG] Pixels computed\n");

    // Terminate program
    free(triangles);
    free(shading_params);
    free(pixels);

    return 0;
}
