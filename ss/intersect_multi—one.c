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

// 计算两个向量的叉积
Point3D crossProduct(Point3D v1, Point3D v2) {
    Point3D result;
    result.x = v1.y * v2.z - v1.z * v2.y;
    result.y = v1.z * v2.x - v1.x * v2.z;
    result.z = v1.x * v2.y - v1.y * v2.x;
    return result;
}

// 使用巴里心坐标判断点是否在三角形内部
int isPointInsideTriangle(Point3D P, Triangle3D triangle) {
    Point3D A = triangle.vertices[0];
    Point3D B = triangle.vertices[1];
    Point3D C = triangle.vertices[2];

    Point3D v0 = {B.x - A.x, B.y - A.y, B.z - A.z};
    Point3D v1 = {C.x - A.x, C.y - A.y, C.z - A.z};
    Point3D v2 = {P.x - A.x, P.y - A.y, P.z - A.z};

    double dot00 = v0.x * v0.x + v0.y * v0.y + v0.z * v0.z;
    double dot01 = v0.x * v1.x + v0.y * v1.y + v0.z * v1.z;
    double dot02 = v0.x * v2.x + v0.y * v2.y + v0.z * v2.z;
    double dot11 = v1.x * v1.x + v1.y * v1.y + v1.z * v1.z;
    double dot12 = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;

    double invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01);
    double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    double v = (dot00 * dot12 - dot01 * dot02) * invDenom;

    return (u >= 0) && (v >= 0) && (u + v <= 1);
}


// 计算线段与平面的交点
Point3D linePlaneIntersection(Line3D line, Triangle3D triangle) {
    Point3D intersection;

    // 计算三角面的法向量
    Point3D edge0 = {
        triangle.vertices[1].x - triangle.vertices[0].x,
        triangle.vertices[1].y - triangle.vertices[0].y,
        triangle.vertices[1].z - triangle.vertices[0].z
    };

    Point3D edge1 = {
        triangle.vertices[2].x - triangle.vertices[0].x,
        triangle.vertices[2].y - triangle.vertices[0].y,
        triangle.vertices[2].z - triangle.vertices[0].z
    };

    Point3D normal = crossProduct(edge0, edge1);

    // 计算线段的方向向量
    Point3D direction = {
        line.end.x - line.start.x,
        line.end.y - line.start.y,
        line.end.z - line.start.z
    };

    // 计算点积
    double dotProduct = direction.x * normal.x + direction.y * normal.y + direction.z * normal.z;

    // 定义一个很小的容差值 epsilon
    const double epsilon = 1e-6;

    // 如果点积接近零，线段与平面平行，无交点
    if (fabs(dotProduct) < epsilon) {
        intersection.x = intersection.y = intersection.z = NAN;
        return intersection;
    }

    // 计算线段与平面的交点的参数 t
    double t = ((triangle.vertices[0].x - line.start.x) * normal.x +
                (triangle.vertices[0].y - line.start.y) * normal.y +
                (triangle.vertices[0].z - line.start.z) * normal.z) /
               dotProduct;

    // 使用容差范围来检查 t 值
    if (t >= epsilon && t <= 1.0 + epsilon) {
        // 计算交点坐标
        intersection.x = line.start.x + t * direction.x;
        intersection.y = line.start.y + t * direction.y;
        intersection.z = line.start.z + t * direction.z;

        // 检查交点是否在三角形内部
        int isInside = isPointInsideTriangle(intersection, triangle);

        if (isInside) {
            return intersection;
            }
    }
    intersection.x = intersection.y = intersection.z = NAN;
    return intersection;
}


Point3D *readPoints(const char *filename, int *numPoints) {
    static Point3D points[2000];
    FILE *pointfile;
    char line[1024];
    int pointIndex = 0;

    pointfile = fopen(filename, "r");
    if (pointfile == NULL) {
        perror("无法打开CSV文件");
        return NULL;
    }

    while (fgets(line, sizeof(line), pointfile) && pointIndex < 2000) {
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

void checkLineIntersections(Point3D *points, int numPoints, Triangle3D triangle) {
    int intersectionIndex = 1;  // 为交点分配一个编号
    Line3D lines[10000];
    int lineIndex = 0;

    // 循环生成所有线段
    for (int i = 0; i < numPoints; i++) {
        for (int j = 0; j < 360; j++) {
            double angle = (2.0 * M_PI / 360) * j;
            Line3D line;
            line.start = points[i];
            double length = 200.0;
            line.end.x = points[i].x + length * cos(angle);
            line.end.y = points[i].y + length * sin(angle);
            line.end.z = points[i].z;  // 在XY平面上进行
            Point3D intersection = linePlaneIntersection(line, triangle);
            if (!isnan(intersection.x)) {
                printf("观测点编号: %d, 交点编号: %d, 交点坐标: (%.2f, %.2f, %.2f)\n", i, intersectionIndex++, intersection.x, intersection.y, intersection.z);
                printf("(%.2f, %.2f, %.2f)\n", line.start.x, line.start.y, line.start.z);
                printf("(%.2f, %.2f, %.2f)\n\n", line.end.x, line.end.y, line.end.z);

            }

            lines[lineIndex++] = line;
        }
    }
}


int main() {
    clock_t start_time = clock();

    Triangle3D triangle = {{{1000.90173, 680.127043, 100.501}, {1001.719601, 676.858675, -100.001}, {1001.719601, 676.858675, 9.001}}};
    
    int numPoints;
    Point3D *points = readPoints("point.csv", &numPoints);
    if (points == NULL) {
        return 1;
    }

    checkLineIntersections(points, numPoints, triangle);

    clock_t end_time = clock();
    double cpu_time_used = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;
    printf("运行时间: %f s\n", cpu_time_used);

    return 0;
}