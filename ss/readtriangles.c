#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>



typedef struct {
    double x, y, z;
} Point3D;

typedef struct {
    Point3D start;
    Point3D end;
} Line3D;

typedef struct {
    Point3D vertices[3];
} Triangle3D;


typedef struct {
    int pointId;
    int intersectionId;
    int rayAngle;
    Point3D intersectionPoint;
} Result;


Triangle3D *readTriangles(const char *filename, int *numTriangles) {
    // 动态数组初始化
    Triangle3D *triangles = NULL;
    int triangleCapacity = 1000;  // 初始容量
    int triangleIndex = 0;
    triangles = (Triangle3D *)malloc(sizeof(Triangle3D) * triangleCapacity);
    if (triangles == NULL) {
        perror("内存分配失败");
        return NULL;
    }

    FILE *trianglefile;
    char line[1024];

    trianglefile = fopen(filename, "rb");
    if (trianglefile == NULL) {
        perror("无法打开CSV文件");
        free(triangles);  // 如果文件打开失败，则释放之前分配的内存
        return NULL;
    }

    while (fgets(line, sizeof(line), trianglefile)) {
        // 扩容检查
        if (triangleIndex >= triangleCapacity) {
            triangleCapacity *= 2;
            triangles = (Triangle3D *)realloc(triangles, sizeof(Triangle3D) * triangleCapacity);
            if (triangles == NULL) {
                perror("内存重新分配失败");
                return NULL;
            }
        }

        // 解析三角面的顶点并存储
        Triangle3D triangle;
        if (sscanf(line, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", 
                   &triangle.vertices[0].x, &triangle.vertices[0].y, &triangle.vertices[0].z,
                   &triangle.vertices[1].x, &triangle.vertices[1].y, &triangle.vertices[1].z,
                   &triangle.vertices[2].x, &triangle.vertices[2].y, &triangle.vertices[2].z) == 9) {
            triangles[triangleIndex] = triangle;
            triangleIndex++;
        } else {
            printf("无法解析三角面行: %s\n", line);
        }
    }

    fclose(trianglefile);
    *numTriangles = triangleIndex;

    // 如果需要，可以缩小数组至实际大小
    triangles = (Triangle3D *)realloc(triangles, sizeof(Triangle3D) * triangleIndex);
    if (triangles == NULL) {
        perror("内存重新分配失败");
        return NULL;
    }

    return triangles;
}



int main() {
    clock_t start_time = clock();

    //从triangles.csv文件中读取三角面集
    int numTriangles;
    Triangle3D *triangles = readTriangles("triangles.csv", &numTriangles);
    if (triangles == NULL) {
        return 1;
    }

    clock_t end_time = clock();
    double cpu_time_used = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;
    printf("运行时间: %f s\n", cpu_time_used);
    return 0;
}

//运行时间: 0.078436 s