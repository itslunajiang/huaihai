#include "mpi_geometry_processing.h"
#include <math.h>
#include <time.h>
#include <stddef.h>
#include <string.h>
#include <sys/sysctl.h>
#include <unistd.h>
#include <mach/mach.h>

MPI_Datatype MPI_RESULT_TYPE;


// 计算两个向量的叉积
Point3D crossProduct(Point3D v1, Point3D v2) {
    Point3D result;
    result.x = v1.y * v2.z - v1.z * v2.y;
    result.y = v1.z * v2.x - v1.x * v2.z;
    result.z = v1.x * v2.y - v1.y * v2.x;
    return result;
}

// 使用巴里心坐标判断点是否在三角面内部
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
                (triangle.vertices[0].z - line.start.z) * normal.z) /dotProduct;

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

// 计算两个点之间的距离
double distance(Point3D p1, Point3D p2) {
    return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2) + pow(p1.z - p2.z, 2));
}

//读取点数据
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

//读取面数据
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



// 将数据存储进 redis
int storeDataInRedis(Point3D *points, int numPoints, Triangle3D *triangles, int numTriangles) {
    redisContext *conn = redisConnect("127.0.0.1", 6379); // 连接到 Redis 服务器
    if (conn == NULL || conn->err) {
        if (conn) {
            printf("Redis连接错误：%s\n", conn->errstr);
            redisFree(conn);
        } else {
            printf("Redis连接错误：无法分配Redis上下文\n");
        }
        return 0;  // 连接失败
    }

    // 存储点数据到 Redis 列表
    const char *pointListKey = "points_list"; // Redis 列表的键名

    // 清空 "points_list" 列表
    redisReply *delReply = (redisReply *)redisCommand(conn, "DEL points_list");
    if (delReply == NULL) {
        printf("Failed to execute DEL command for points_list\n");
        // 处理错误
    } else {
        freeReplyObject(delReply);
    }

    for (int i = 0; i < numPoints; i++) {
        char pointData[100];
        snprintf(pointData, sizeof(pointData), "%.2lf,%.2lf,%.2lf", points[i].x, points[i].y, points[i].z);
        redisReply *reply = (redisReply *)redisCommand(conn, "RPUSH %s %s", pointListKey, pointData);
        if (reply == NULL) {
            printf("Failed to execute RPUSH command for points_list\n");
        } else {
            freeReplyObject(reply);
        }
    }

    // 存储三角面数据到 Redis 列表
    const char *triangleListKey = "triangles_list"; // Redis 列表的键名


    // 清空 "triangles_list" 列表
    delReply = (redisReply *)redisCommand(conn, "DEL triangles_list"); // 重新使用同一变量名
    if (delReply == NULL) {
        printf("Failed to execute DEL command for triangles_list\n");
        // 处理错误
    } else {
        freeReplyObject(delReply);
    }

    for (int i = 0; i < numTriangles; i++) {
        char triangleData[200];
        snprintf(triangleData, sizeof(triangleData), "%.2lf,%.2lf,%.2lf,%.2lf,%.2lf,%.2lf,%.2lf,%.2lf,%.2lf",
                triangles[i].vertices[0].x, triangles[i].vertices[0].y, triangles[i].vertices[0].z,
                triangles[i].vertices[1].x, triangles[i].vertices[1].y, triangles[i].vertices[1].z,
                triangles[i].vertices[2].x, triangles[i].vertices[2].y, triangles[i].vertices[2].z);
        redisReply *reply = (redisReply *)redisCommand(conn, "RPUSH %s %s", triangleListKey, triangleData);
        if (reply == NULL) {
            printf("Failed to execute RPUSH command for triangles_list\n");
        } else {
            freeReplyObject(reply);
        }
    }

    // 关闭 Redis 连接并释放资源
    redisFree(conn);

    return 1;  // 连接成功并存储完成
}


// 从 Redis 加载点数据(一次性)
int loadPointsFromRedis(Point3D **points, int *numPoints, int rank) {
    // 连接到 Redis 服务器
    redisContext *conn = redisConnect("127.0.0.1", 6379);
    if (conn == NULL || conn->err) {
        if (conn) {
            printf("Redis连接错误：%s\n", conn->errstr);
            redisFree(conn);
        } else {
            printf("Redis连接错误：无法分配Redis上下文\n");
        }
        return 0;  // 连接失败
    }

    const char *pointListKey = "points_list"; // Redis 列表的键名

    // 从 Redis 列表中获取点数据
    redisReply *rangeReply = (redisReply *)redisCommand(conn, "LRANGE %s 0 -1", pointListKey);
    if (rangeReply == NULL) {
        printf("Failed to execute LRANGE command for points_list\n");
        redisFree(conn);
        return 0; // 失败
    }

    *numPoints = rangeReply->elements;
    *points = (Point3D *)malloc(*numPoints * sizeof(Point3D));

    for (size_t i = 0; i < rangeReply->elements; i++) {
        char *pointData = rangeReply->element[i]->str;
        sscanf(pointData, "%lf,%lf,%lf", &(*points)[i].x, &(*points)[i].y, &(*points)[i].z);
    }

    freeReplyObject(rangeReply);
    redisFree(conn);

    return 1; // 加载成功
}

// 从 Redis 加载面数据(批量)
int loadTrianglesFromRedisPartial(Triangle3D **triangles, int startIndex, int count, int rank) {
    // 连接到 Redis 服务器
    redisContext *conn = redisConnect("127.0.0.1", 6379);
    if (conn == NULL || conn->err) {
        if (conn) {
            printf("Redis连接错误：%s\n", conn->errstr);
            redisFree(conn);
        } else {
            printf("Redis连接错误：无法分配Redis上下文\n");
        }
        return 0;  // 连接失败
    }

    const char *triangleListKey = "triangles_list"; // Redis 列表的键名

    // 从 Redis 列表中获取一部分三角面数据
    redisReply *rangeReply = (redisReply *)redisCommand(conn, "LRANGE %s %d %d", triangleListKey, startIndex, startIndex + count - 1);
    if (rangeReply == NULL) {
        printf("Failed to execute LRANGE command for triangles_list\n");
        redisFree(conn);
        return 0; // 失败
    }

    // 释放之前分配的内存
    if (*triangles != NULL) {
        free(*triangles); // 释放之前分配的内存
    }

    // 为三角面数据分配内存
    *triangles = (Triangle3D *)malloc(count * sizeof(Triangle3D));
    if (*triangles == NULL) {
        printf("内存分配失败\n");
        freeReplyObject(rangeReply);
        redisFree(conn);
        return 0; // 内存分配失败
    }

    // 解析并存储三角面数据
    if (rangeReply->elements > count) {
        printf("Error: Received more elements (%lu) than expected (%d)\n", rangeReply->elements, count);
        freeReplyObject(rangeReply);
        redisFree(conn);
        return 0; // 返回错误
    }

    // 解析并存储三角面数据
    for (size_t i = 0; i < count; i++) {
        char *triangleData = rangeReply->element[i]->str;
        sscanf(triangleData, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
               &(*triangles)[i].vertices[0].x, &(*triangles)[i].vertices[0].y, &(*triangles)[i].vertices[0].z,
               &(*triangles)[i].vertices[1].x, &(*triangles)[i].vertices[1].y, &(*triangles)[i].vertices[1].z,
               &(*triangles)[i].vertices[2].x, &(*triangles)[i].vertices[2].y, &(*triangles)[i].vertices[2].z);
    }

    freeReplyObject(rangeReply);
    redisFree(conn);

    return 1; // 加载成功
}


int globalIntersectionIndex = 1; // 全局变量

// 计算交点的函数，返回 Result 数组的指针
Result* checkLineIntersections(Point3D *points, int numPoints, Triangle3D *triangles, int numTriangles, int *totalResultsOut) {
    int rank, size;
    
    // 获取当前进程的 rank 和总进程数 size
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int pointsPerProcess = numPoints / size;
    int remainder = numPoints % size;

    int startPoint = rank * pointsPerProcess + (rank < remainder ? rank : remainder);
    int endPoint = startPoint + pointsPerProcess + (rank < remainder ? 1 : 0);

    int maxResults = numPoints * 5 * 36;
    Result localResults[maxResults];
    int localResultCount = 0;


    // 循环生成所有线段
    for (int i = startPoint; i < endPoint; i++) {
       for (int offset = 1; offset <= 5; offset++){
            for (int j = 0; j < 36; j++) {
                double angle = (2.0 * M_PI / 36) * j;
                Line3D line;
                line.start = points[i];
                double length = 200.0;
                line.end.x = points[i].x + length * cos(angle);
                line.end.y = points[i].y + length * sin(angle);
                line.end.z = points[i].z + offset;  // 在XY平面上进行，并向Z轴每隔一米做偏移

                Point3D nearestIntersection = {NAN, NAN, NAN};
                double minDistance = 1e30;  // 初始设置一个很大的距离值



                // 遍历所有三角面
                for (int k = 0; k < numTriangles; k++) {
                    Point3D intersection = linePlaneIntersection(line, triangles[k]);
                    if (!isnan(intersection.x)) {
                        double dist = distance(line.start, intersection);
                        if (dist < minDistance) {
                            minDistance = dist;
                            nearestIntersection = intersection;
                        }
                    }
                }

                if (localResultCount >= maxResults) {
                    printf("Error: Buffer overflow detected!\n");
                    MPI_Abort(MPI_COMM_WORLD, -1);  // 终止所有MPI进程
                }

                if (!isnan(nearestIntersection.x)) {
                    localResults[localResultCount].pointId = i + 1;
                    localResults[localResultCount].directionId = j * 10 + offset;
                    localResults[localResultCount].intersectionPoint = nearestIntersection;
                    localResults[localResultCount].parameter = minDistance;
                    localResultCount++;
                }
            }
        }
    }

    // 同步所有进程到该点，再运行后续代码
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (rank == 0) {
        int counts[size];  // 每个进程的结果数量
        int displs[size];  // 每个进程的偏移量
        int totalResults = 0;

        MPI_Gather(&localResultCount, 1, MPI_INT, counts, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // 计算偏移量和总结果数量
        for (int i = 0; i < size; i++) {
            displs[i] = totalResults;
            totalResults += counts[i];
        }

        // 动态分配 allResults 数组
        Result *allResults = (Result *)malloc(sizeof(Result) * totalResults);

        MPI_Gatherv(localResults, localResultCount, MPI_RESULT_TYPE,
                    allResults, counts, displs, MPI_RESULT_TYPE, 0, MPI_COMM_WORLD);

        *totalResultsOut = totalResults;  // 将 totalResults 的值传递出去

        return allResults;  // 返回 allResults 数组

    } else {
        MPI_Gather(&localResultCount, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Gatherv(localResults, localResultCount, MPI_RESULT_TYPE,
                    NULL, NULL, NULL, MPI_RESULT_TYPE, 0, MPI_COMM_WORLD);
                    return NULL;
    }
}


// Function to compute normal vector of a triangle
Point3D getTriangleNormal(Triangle3D triangle) {
    Point3D edge1 = {triangle.vertices[1].x - triangle.vertices[0].x,
                     triangle.vertices[1].y - triangle.vertices[0].y,
                     triangle.vertices[1].z - triangle.vertices[0].z};

    Point3D edge2 = {triangle.vertices[2].x - triangle.vertices[0].x,
                     triangle.vertices[2].y - triangle.vertices[0].y,
                     triangle.vertices[2].z - triangle.vertices[0].z};

    // Cross product of edge1 and edge2
    return crossProduct(edge1, edge2);
}

// Function to normalize a vector
Point3D normalizeVector(Point3D v) {
    double length = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    return (Point3D){v.x / length, v.y / length, v.z / length};
}


// 计算交点的函数，返回 Result 数组的指针
Result* checkIntersectionsAngle(Point3D *points, int numPoints, Triangle3D *triangles, int numTriangles, int *totalResultsOut) {
    int rank, size;

    int nearestTriangleIndex = -1;  // Initialize to an invalid value


    // 获取当前进程的 rank 和总进程数 size
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int pointsPerProcess = numPoints / size;
    int remainder = numPoints % size;

    int startPoint = rank * pointsPerProcess + (rank < remainder ? rank : remainder);
    int endPoint = startPoint + pointsPerProcess + (rank < remainder ? 1 : 0);

    int maxResults = numPoints * 5 * 36;
    Result localResults[maxResults];
    int localResultCount = 0;


    // 循环生成所有线段
    for (int i = startPoint; i < endPoint; i++) {
       for (int offset = 1; offset <= 5; offset++){
            for (int j = 0; j < 36; j++) {
                double angle = (2.0 * M_PI / 36) * j;
                Line3D line;
                line.start = points[i];
                double length = 200.0;
                line.end.x = points[i].x + length * cos(angle);
                line.end.y = points[i].y + length * sin(angle);
                line.end.z = points[i].z + offset;  // 在XY平面上进行，并向Z轴每隔一米做偏移

                Point3D nearestIntersection = {NAN, NAN, NAN};
                double minDistance = 1e30;  // 初始设置一个很大的距离值



                // 遍历所有三角面
                for (int k = 0; k < numTriangles; k++) {
                    Point3D intersection = linePlaneIntersection(line, triangles[k]);
                    if (!isnan(intersection.x)) {
                        double dist = distance(line.start, intersection);
                        if (dist < minDistance) {
                            minDistance = dist;
                            nearestIntersection = intersection;
                            nearestTriangleIndex = k; // Update nearestTriangleIndex here

                        }
                    }
                }

                if (localResultCount >= maxResults) {
                    printf("Error: Buffer overflow detected!\n");
                    MPI_Abort(MPI_COMM_WORLD, -1);  // 终止所有MPI进程
                }

                if (!isnan(nearestIntersection.x) && nearestTriangleIndex != -1) {
                    Point3D lineDirection = {line.end.x - line.start.x,
                                            line.end.y - line.start.y,
                                            line.end.z - line.start.z};
                    Point3D normalizedDirection = normalizeVector(lineDirection);
                    Point3D normal = normalizeVector(getTriangleNormal(triangles[nearestTriangleIndex]));  // Use nearestTriangleIndex

                    // Compute angle using dot product
                    double dot = normalizedDirection.x * normal.x +
                                normalizedDirection.y * normal.y +
                                normalizedDirection.z * normal.z;
                    double angle = 90 - acos(dot)* (180.0 / M_PI);

                    localResults[localResultCount].pointId = i + 1;
                    localResults[localResultCount].directionId = j * 10 + offset;
                    localResults[localResultCount].intersectionPoint = nearestIntersection;
                    localResults[localResultCount].parameter = angle;
                    localResultCount++;
                }
            }
        }
    }

    // 同步所有进程到该点，再运行后续代码
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (rank == 0) {
        int counts[size];  // 每个进程的结果数量
        int displs[size];  // 每个进程的偏移量
        int totalResults = 0;

        MPI_Gather(&localResultCount, 1, MPI_INT, counts, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // 计算偏移量和总结果数量
        for (int i = 0; i < size; i++) {
            displs[i] = totalResults;
            totalResults += counts[i];
        }

        // 动态分配 allResults 数组
        Result *allResults = (Result *)malloc(sizeof(Result) * totalResults);

        MPI_Gatherv(localResults, localResultCount, MPI_RESULT_TYPE,
                    allResults, counts, displs, MPI_RESULT_TYPE, 0, MPI_COMM_WORLD);

        *totalResultsOut = totalResults;  // 将 totalResults 的值传递出去

        return allResults;  // 返回 allResults 数组

    } else {
        MPI_Gather(&localResultCount, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Gatherv(localResults, localResultCount, MPI_RESULT_TYPE,
                    NULL, NULL, NULL, MPI_RESULT_TYPE, 0, MPI_COMM_WORLD);
                    return NULL;
    }
}




// 比较函数，用于排序(升序)
int compareResults(const void *a, const void *b) {
    const Result *result1 = (const Result *)a;
    const Result *result2 = (const Result *)b;

    // 首先按 pointId 排序
    if (result1->pointId != result2->pointId) {
        return result1->pointId - result2->pointId;
    }

    // pointId 相同，按 direction.x 排序
    if (result1->directionId != result2->directionId) {
        return (result1->directionId < result2->directionId) ? -1 : 1;
    }

    return 0;  // 完全相同
}

// 用于判断两个结果是否等价（即 pointId 和 direction 完全相同）
int areResultsEquivalent(const Result *a, const Result *b) {
    return (a->pointId == b->pointId) &&
           (a->directionId == b->directionId);
}

// 筛选结果，只保留最近交点
Result* filterResults(Result *allResults, int totalResults, int *filteredCount) {
    *filteredCount = 0;

    // 首先对结果进行排序
    qsort(allResults, totalResults, sizeof(Result), compareResults);

    // 为过滤后的结果分配足够的内存
    Result *filteredResults = (Result *)malloc(sizeof(Result) * totalResults);
    if (filteredResults == NULL) {
        // 处理内存分配失败的情况
        return NULL;
    }

    for (int i = 0; i < totalResults; i++) {
        int found = 0;
        // 检查当前结果是否已经存在于filteredResults中
        for (int j = 0; j < *filteredCount; j++) {
            if (areResultsEquivalent(&allResults[i], &filteredResults[j])) {
                found = 1;
                // 如果当前结果的directionLength更小，则替换
                if (allResults[i].parameter < filteredResults[j].parameter) {
                    filteredResults[j] = allResults[i];
                }
            }
        }
        // 如果当前结果在filteredResults中未找到，则添加
        if (!found) {
            filteredResults[*filteredCount] = allResults[i];
            (*filteredCount)++;
        }
    }

    return filteredResults;
}


long initialMemoryUsage = 0;
long peakMemoryUsage = 0;
long currentMemoryUsage = 0;
// 获取系统可用内存
long long getAvailableMemory() {
    int mib[2];
    int64_t physical_memory;
    size_t length;

    // Total physical memory
    mib[0] = CTL_HW;
    mib[1] = HW_MEMSIZE;
    length = sizeof(physical_memory);
    if (sysctl(mib, 2, &physical_memory, &length, NULL, 0) == -1) {
        perror("获取物理内存失败");
        return -1;
    }

    // Current memory usage
    vm_size_t page_size;
    mach_port_t mach_port = mach_host_self();
    vm_statistics64_data_t vm_stats;
    mach_msg_type_number_t count = sizeof(vm_stats) / sizeof(natural_t);
    if (host_page_size(mach_port, &page_size) != KERN_SUCCESS || host_statistics64(mach_port, HOST_VM_INFO, (host_info64_t)&vm_stats, &count) != KERN_SUCCESS) {
        perror("获取当前内存使用情况失败");
        return -1;
    }

    // Available memory
    int64_t free_memory = (int64_t)vm_stats.free_count * (int64_t)page_size;

    return free_memory / (1024 * 1024); // 可用内存,单位MB
}

// 计算 double 数组中所有元素所需的内存大小
size_t calculateDoubleArrayDataSize(const double *dataArray, size_t numElements) {
    size_t dataSize = numElements * sizeof(double);
    return dataSize;
}


Result* mpiGeometryProcessing(int m, ComputeFunc computeFunc, int *totalResultsOut) {

    MPI_Init(NULL, NULL);
    // 在 MPI_Init 后获取起始时间
    double start_time = MPI_Wtime();
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int numPoints;
    Point3D *points = NULL; // 初始化为 NULL，以便所有进程都能访问
    int numTriangles;
    Triangle3D *triangles = NULL; // 初始化为 NULL，以便所有进程都能访问

    // rank0 加载所有数据
    if (rank == 0) {
        // 根进程读取数据
        points = readPoints("point.csv", &numPoints);
        if (points == NULL) {
            MPI_Finalize(); // 执行失败后需要释放 MPI 资源
            return NULL;
        }
        triangles = readTriangles("triangles.csv", &numTriangles);
        if (triangles == NULL) {
            free(points);
            MPI_Finalize();
            return NULL;
        }
        MPI_Bcast(&numPoints, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&numTriangles, 1, MPI_INT, 0, MPI_COMM_WORLD);
        printf("Rank %d: read data, %f s \n", rank, MPI_Wtime() - start_time);

    } else {
        MPI_Bcast(&numPoints, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&numTriangles, 1, MPI_INT, 0, MPI_COMM_WORLD);
        printf("Rank %d: wait rank 0 reading data, %f s \n", rank, MPI_Wtime() - start_time);
    }


    //   自定义 mpi 数据结构 - Result 结构体
    MPI_Datatype types[4] = { MPI_INT, MPI_INT, MPI_DOUBLE, MPI_DOUBLE};
    int blocklengths[4] = { 1, 1, 3, 1 }; // 3 doubles for intersectionPoint (x, y, z), 1 double for parameter
    MPI_Aint offsets[4];
    offsets[0] = offsetof(Result, pointId);
    offsets[1] = offsetof(Result, directionId);
    offsets[2] = offsetof(Result, intersectionPoint);
    offsets[3] = offsetof(Result, parameter);
    MPI_Type_create_struct(4, blocklengths, offsets, types, &MPI_RESULT_TYPE);
    MPI_Type_commit(&MPI_RESULT_TYPE);

    Result* outputResults = NULL; // 在函数作用域内声明

    // 目前可用内存情况
    int64_t availableMemory = getAvailableMemory();
    // 估算所需内存大小
    size_t pointDataSize = calculateDoubleArrayDataSize((double *)points, numPoints * 3);
    size_t trianglesDataSize = calculateDoubleArrayDataSize((double *)triangles, numTriangles * 9);
    size_t requiredMemorySize = (pointDataSize + trianglesDataSize) / (1024*1024);

    // 根进程打印 requiredMemorySize 和 availableMemory 的值
    if (rank == 0) {
        printf("requiredMemorySize: %ld MB\n", (long)requiredMemorySize);
        printf("availableMemory: %ld MB\n", (long)availableMemory);
    }

    
    //判断内存是否够用
    if (availableMemory >= (long long)requiredMemorySize) {

        // 可用内存足够，可以将数据加载到内存中并处理

        if (rank == 0){
            printf("sufficient memory, local processing\n");
            printf("numpoints: %d, numtriangles: %d\n", numPoints, numTriangles);

            // 广播点数据和三角面数据的数量到所有进程
            MPI_Bcast(&numPoints, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(&numTriangles, 1, MPI_INT, 0, MPI_COMM_WORLD);

            // 广播点数据和三角面数据
            MPI_Bcast(points, numPoints * sizeof(Point3D), MPI_BYTE, 0, MPI_COMM_WORLD);
            MPI_Bcast(triangles, numTriangles * sizeof(Triangle3D), MPI_BYTE, 0, MPI_COMM_WORLD);

            printf("Rank %d: broadcast data, %f s \n", rank, MPI_Wtime() - start_time);
        } else {
            // 接收点数据和三角面数据的数量
            MPI_Bcast(&numPoints, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(&numTriangles, 1, MPI_INT, 0, MPI_COMM_WORLD);

            // 为点数据和三角面数据分配内存
            points = (Point3D *)malloc(numPoints * sizeof(Point3D));
            triangles = (Triangle3D *)malloc(numTriangles * sizeof(Triangle3D));
            if (points == NULL || triangles == NULL) {
                printf("Memory allocation error.\n");
                MPI_Finalize();
                return NULL;
            }

            // 接收广播的点数据和三角面数据
            MPI_Bcast(points, numPoints * sizeof(Point3D), MPI_BYTE, 0, MPI_COMM_WORLD);
            MPI_Bcast(triangles, numTriangles * sizeof(Triangle3D), MPI_BYTE, 0, MPI_COMM_WORLD);
            printf("Rank %d: receive bcast_data, %f s \n", rank, MPI_Wtime() - start_time);

        }

        // 检查线与面的交点
        outputResults = computeFunc(points, numPoints, triangles, numTriangles, totalResultsOut);

        printf("Rank %d: finish calculation, %f s \n", rank, MPI_Wtime() - start_time);

        // 释放点数据和三角面数据内存
        free(points);
        free(triangles);

    } else {

        // 可用内存不足，可以将数据加载到redis并处理

        // rank 0 存储点和三角面数据到Redis
        if (rank == 0){
            printf("insufficient memory, processing with redis\n");
            printf("numpoints: %d, numtriangles: %d\n", numPoints, numTriangles);

            int success = storeDataInRedis(points, numPoints, triangles, numTriangles);
            if (!success) {
            // 存储到 Redis 失败，处理错误或返回
                printf("Failed to store data in Redis.\n");
                free(points);       // 释放点数据的内存
                free(triangles);    // 释放三角面数据的内存
                MPI_Type_free(&MPI_RESULT_TYPE);
                MPI_Finalize();
                return NULL;
            } else {
                printf("Rank 0: store data in Redis, %f s \n\n", MPI_Wtime() - start_time);
            }
            
            // 释放点数据和三角面数据内存
            free(points);
            free(triangles);
        } else {  
            printf("Rank %d: wait data storage, %f s \n", rank, MPI_Wtime() - start_time);
        }

        // 同步所有进程到该点，再运行后续代码
        MPI_Barrier(MPI_COMM_WORLD);

        int numBatches = (numTriangles + m - 1) / m; // 计算需要的循环次数，向上取整
        Result **allBatchesResults = (Result**)malloc(sizeof(Result*) * numBatches);
        int *totalResultsPerBatch = (int*)malloc(sizeof(int) * numBatches);

        
        //分批加载数据
        for (int batch = 0; batch < numBatches; batch++) {

            Point3D *localPoints = NULL;
            int numLocalPoints = 0;
            Triangle3D *localTriangles = NULL;
            int numTrianglesThisBatch = (batch * m + m > numTriangles) ? (numTriangles - batch * m) : m;


            if (rank == 0) {
                // 从 Redis 加载点数据
                int successPoint = loadPointsFromRedis(&localPoints, &numLocalPoints, rank);
                if (!successPoint) {
                    printf("Rank 0: Failed to load points from Redis.\n");
                    MPI_Type_free(&MPI_RESULT_TYPE);
                    MPI_Finalize();
                    return NULL;
                }
                // 广播点数据的数量给其他进程
                MPI_Bcast(&numLocalPoints, 1, MPI_INT, 0, MPI_COMM_WORLD);
                // 广播点数据
                MPI_Bcast(localPoints, numLocalPoints * sizeof(Point3D), MPI_BYTE, 0, MPI_COMM_WORLD);
                

                // 从redis加载三角面数据
                int successTriangle = loadTrianglesFromRedisPartial(&localTriangles, batch * m, numTrianglesThisBatch, rank);
                MPI_Bcast(&successTriangle, 1, MPI_INT, 0, MPI_COMM_WORLD);
                if (successTriangle) {
                    // 广播三角面数据的数量
                    MPI_Bcast(&numTrianglesThisBatch, 1, MPI_INT, 0, MPI_COMM_WORLD);
                    // 广播三角面数据
                    MPI_Bcast(localTriangles, numTrianglesThisBatch * sizeof(Triangle3D), MPI_BYTE, 0, MPI_COMM_WORLD);
                } else {
                    printf("Rank %d: Failed to load triangles.\n", rank);
                }


                printf("Rank %d: broadcast data, %f s \n", rank, MPI_Wtime() - start_time);
                
                allBatchesResults[batch] = computeFunc(localPoints, numLocalPoints, localTriangles, numTrianglesThisBatch, &totalResultsPerBatch[batch]);

                printf("batch: %d, Rank %d: finish calculation, %f s \n", batch, rank, MPI_Wtime() - start_time);
                
                MPI_Barrier(MPI_COMM_WORLD);

            } else {
            
                printf("Rank %d: wait rank0 data, %f s \n", rank, MPI_Wtime() - start_time);

                // 其他进程接收点数据数量和数据
                MPI_Bcast(&numLocalPoints, 1, MPI_INT, 0, MPI_COMM_WORLD);
                Point3D *localPoints = (Point3D *)malloc(numLocalPoints * sizeof(Point3D));
                MPI_Bcast(localPoints, numLocalPoints * sizeof(Point3D), MPI_BYTE, 0, MPI_COMM_WORLD);

                // 其他进程接收面数据数量和数据
                int successTriangle;
                MPI_Bcast(&successTriangle, 1, MPI_INT, 0, MPI_COMM_WORLD);
                if (successTriangle) {
                    // 接收三角面数据的数量
                    MPI_Bcast(&numTrianglesThisBatch, 1, MPI_INT, 0, MPI_COMM_WORLD);
                    // 分配内存以接收数据
                    localTriangles = (Triangle3D *)malloc(numTrianglesThisBatch * sizeof(Triangle3D));
                    // 接收三角面数据
                    MPI_Bcast(localTriangles, numTrianglesThisBatch * sizeof(Triangle3D), MPI_BYTE, 0, MPI_COMM_WORLD);

                    printf("Rank %d: receive bcast_data, %f s \n\n", rank, MPI_Wtime() - start_time);

                } else {
                    printf("Rank %d: Failed to receive triangles.\n", rank);
                }
                
                allBatchesResults[batch] = computeFunc(localPoints, numLocalPoints, localTriangles, numTrianglesThisBatch, &totalResultsPerBatch[batch]);

                printf("batch: %d, Rank %d: finish calculation, %f s \n", batch, rank, MPI_Wtime() - start_time);
                
                MPI_Barrier(MPI_COMM_WORLD);

            }


            // 在循环的每次迭代结束时释放内存
            if (localTriangles != NULL) {
                free(localTriangles);
                localTriangles = NULL;
            }
            if (localPoints != NULL) {
                free(localPoints);
                localPoints = NULL;
            }


        }
        
        // 同步所有进程到该点，再运行后续代码
        MPI_Barrier(MPI_COMM_WORLD);

        //rank 0 排序并筛选有效结果
        if (rank == 0) {

            printf("Rank %d: start filter, %f s \n", rank, MPI_Wtime() - start_time);


            // 计算所有批次结果的总数
            int totalAllBatchesResults = 0;
            for (int i = 0; i < numBatches; i++) {
                totalAllBatchesResults += totalResultsPerBatch[i];
            }

            // 分配一个足够大的数组来存储所有结果
            Result *allCombinedResults = (Result *)malloc(sizeof(Result) * totalAllBatchesResults);

            // 将每个批次的结果复制到这个大数组中
            int currentPosition = 0;
            for (int batch = 0; batch < numBatches; batch++) {
                memcpy(allCombinedResults + currentPosition, allBatchesResults[batch], sizeof(Result) * totalResultsPerBatch[batch]);
                currentPosition += totalResultsPerBatch[batch];
                free(allBatchesResults[batch]); // 释放原始批次数组的内存
            }

            // 排序和筛选 allResults
            int filteredCount = 0;
            outputResults = filterResults(allCombinedResults, totalAllBatchesResults, &filteredCount);
            *totalResultsOut = filteredCount;
    
            // 释放 allCombinedResults
            free(allCombinedResults);

        } else {
            printf("Rank %d: wait rank 0 filter, %f s \n\n", rank, MPI_Wtime() - start_time);
        }

    }


    if (rank == 0){
        printf("运行时间: %f s ", MPI_Wtime() - start_time);
    }

    // 释放 MPI 自定义数据类型
    MPI_Type_free(&MPI_RESULT_TYPE);
    // 结束 MPI，释放 MPI 资源
    MPI_Finalize();

    // 只有在根进程中返回结果指针，其他进程返回 NULL
    if (rank == 0) {
        return outputResults;
    } else {
        return NULL;
    }
    
}
