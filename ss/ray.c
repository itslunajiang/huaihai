#include <stdio.h>
#include <math.h>
#include <time.h>

typedef struct {
    double x;
    double y;
    double z;
} Point3D;

typedef struct {
    Point3D start;
    Point3D end;
} Ray3D;

int main() {
    FILE *pointfile;
    char line[1024]; // 用于存储每行数据的字符串缓冲区
    Point3D points[2000]; // 用于存储提取的坐标点
    int numPoints = 0;

    // 打开CSV文件以进行读取
    pointfile = fopen("point.csv", "r");
    if (pointfile == NULL) {
        perror("无法打开CSV文件");
        return 1;
    }

    // 逐行读取CSV文件中的坐标点数据
    while (fgets(line, sizeof(line), pointfile) && numPoints < 6) {
        double x, y, z;
        if (sscanf(line, "%lf,%lf,%lf", &x, &y, &z) == 3) {
            printf("point %.2d: x=%.2f, y=%.2f, z=%.2f\n", numPoints, x, y, z);
            points[numPoints].x = x;
            points[numPoints].y = y;
            points[numPoints].z = z;
            numPoints++;
        } else {
            printf("无法解析坐标行: %s\n", line);
        }
    }

    // 关闭文件
    fclose(pointfile);

    clock_t start_time = clock();

    // 生成射线并输出
    for (int i = 0; i < numPoints; i++) {
        Point3D viewpoint = points[i];

        for (int j = 0; j < 36; j++) {
            double angle = (2.0 * M_PI / 36) * j;
            Ray3D ray;
            ray.start = viewpoint;
            double length = 200.0;
            ray.end.x = viewpoint.x + length * cos(angle);
            ray.end.y = viewpoint.y + length * sin(angle);
            ray.end.z = viewpoint.z; // 平视

            printf("Point %d, Ray %d: Start (%.2f, %.2f, %.2f), End (%.2f, %.2f, %.2f)\n",
                   i + 1, j + 1, ray.start.x, ray.start.y, ray.start.z,
                   ray.end.x, ray.end.y, ray.end.z);
        }
    }

    clock_t end_time = clock();

    double cpu_time_used = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;

    printf("运行时间: %f s\n", cpu_time_used);

    return 0;
}
