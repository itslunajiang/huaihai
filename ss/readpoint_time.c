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

Point3D *readPoints(const char *filename, int *numPoints) {
    FILE *pointfile;
    char line[1024];
    int pointIndex = 0; // 当前的点索引
    int pointCapacity = 10; // 初始容量
    Point3D *points = NULL; // 初始化为空指针

    // 分配初始容量
    points = (Point3D *) malloc(sizeof(Point3D) * pointCapacity);
    if (points == NULL) {
        perror("无法分配内存");
        return NULL;
    }

    pointfile = fopen(filename, "r");
    if (pointfile == NULL) {
        perror("无法打开CSV文件");
        free(points);  
        return NULL;
    }

    while (fgets(line, sizeof(line), pointfile)) {
        // 扩容检查
        if (pointIndex >= pointCapacity) {
            pointCapacity *= 2;  // 两倍扩容
            Point3D *newPoints = (Point3D *) realloc(points, sizeof(Point3D) * pointCapacity);
            if (newPoints == NULL) {
                perror("无法重新分配内存");
                free(points);  
                fclose(pointfile);  
                return NULL;
            }
            points = newPoints;
        }

        double x, y, z;
        if (sscanf(line, "%lf,%lf,%lf", &x, &y, &z) == 3) {
            points[pointIndex].x = x;
            points[pointIndex].y = y;
            points[pointIndex].z = z;
            pointIndex++;
        } else {
            printf("无法解析坐标行: %s\n", line);
        }
    }

    fclose(pointfile); 
    *numPoints = pointIndex;  
    return points;
}


int main() {
    clock_t start_time = clock();

    //从point.csv文件中读取点集    
    int numPoints;
    Point3D *points = readPoints("point.csv", &numPoints);
    if (points == NULL) {
        return 1;
    }


    clock_t end_time = clock();
    double cpu_time_used = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;
    printf("运行时间: %f s\n", cpu_time_used);
    return 0;
}

//运行时间: 0.000134 s
