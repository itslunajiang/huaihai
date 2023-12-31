#include <stdio.h>
#include <math.h>
#include <time.h>

typedef struct {
    double x;
    double y;
    double z;
} Point3D;

typedef struct {
    Point3D vertices[3];
} Triangle3D;

int main() {
    FILE *trianglesfile;
    char line[1024]; // 用于存储每行数据的字符串缓冲区
    Triangle3D triangles[2000]; // 用于存储提取的三角面数据
    int numTriangles = 0;

    // 打开CSV文件以进行读取
    trianglesfile = fopen("triangles.csv", "r");
    if (trianglesfile == NULL) {
        perror("无法打开CSV文件");
        return 1;
    }

    // 逐行读取CSV文件中的坐标点数据并将其存储为三角面
    while (fgets(line, sizeof(line), trianglesfile) && numTriangles < 20) {
        Triangle3D triangle;
        if (sscanf(line, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", 
                   &triangle.vertices[0].x, &triangle.vertices[0].y, &triangle.vertices[0].z,
                   &triangle.vertices[1].x, &triangle.vertices[1].y, &triangle.vertices[1].z,
                   &triangle.vertices[2].x, &triangle.vertices[2].y, &triangle.vertices[2].z) == 9) {
            printf("解析三角面 #%d:\n", numTriangles + 1);
            printf("Vertex 1: (%.2f, %.2f, %.2f)\n", triangle.vertices[0].x, triangle.vertices[0].y, triangle.vertices[0].z);
            printf("Vertex 2: (%.2f, %.2f, %.2f)\n", triangle.vertices[1].x, triangle.vertices[1].y, triangle.vertices[1].z);
            printf("Vertex 3: (%.2f, %.2f, %.2f)\n", triangle.vertices[2].x, triangle.vertices[2].y, triangle.vertices[2].z);
            triangles[numTriangles] = triangle;
            numTriangles++;
        } else {
            printf("无法解析三角面行: %s\n", line);
        }
    }

    // 关闭文件
    fclose(trianglesfile);

    return 0;
}
