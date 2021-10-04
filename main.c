#include <stdio.h>
#include <stdlib.h>

typedef struct point{
    float x;
    float y;
    float z;
} point;

typedef struct triangle{
    point a;
    point b;
    point c;
} triangle;

triangle* load_triangles(FILE* fp){
    unsigned int number_triangles = fscanf(fp, "%u");
}

int main(int argc, char **argv) {

    // Parameters check
    if(argc != 4){
        printf("-- Missing parameters --\n"
               "1. Data file (e.g. bun.dat)\n2. Parameter file (e.g. param_bun.txt)\n3. Output file (e.g. bun.ppm)");
        return 1;
    }

    FILE *dataFile = fopen(argv[1], "r");
    FILE *paramsFile = fopen(argv[2], "r");

    if(!dataFile){
        printf("Data file not readable");
        return 1;
    }

    if(!paramsFile){
        printf("Parameters file not readable");
        return 1;
    }

    unsigned int number_triangles = fscanf(dataFile, "%u");
    printf("%s", number_triangles);

    return 0;
}