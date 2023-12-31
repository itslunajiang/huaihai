#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <stddef.h>
#include <sys/resource.h>
#include <string.h>
#include <mach/mach.h>
#include "hiredis.h"

MPI_Datatype MPI_RESULT_TYPE;

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


// 计算两个点之间的距离
double distance(Point3D p1, Point3D p2) {
    return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2) + pow(p1.z - p2.z, 2));
}


void checkLineIntersections(Point3D *points, int numPoints, Triangle3D *triangles, int numTriangles) {
    int rank, size;
    
    redisContext *conn = redisConnect("127.0.0.1", 6379); // 连接到Redis服务器
        if (conn == NULL || conn->err) {
        if (conn) {
            printf("Redis连接错误：%s\n", conn->errstr);
            redisFree(conn);
        } else {
            printf("Redis连接错误：无法分配Redis上下文\n");
        }
        return;
    }

    redisReply *reply; // Redis回复对象

    
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
    int intersectionIndex = 1;


    // 循环生成所有线段
    for (int i = startPoint; i < endPoint; i++) {
       for (int offset = 1; offset <= 5; offset++){
        for (int j = 0; j < 36; j++) {
            double angle = (2.0 * M_PI / 36) * j;
            int ray_angle = -1;
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
                        ray_angle = (int) round (angle * 180 / M_PI);
                    }
                }
            }

            if (localResultCount >= maxResults) {
                printf("Error: Buffer overflow detected!\n");
                MPI_Abort(MPI_COMM_WORLD, -1);  // 终止所有MPI进程
            }

            if (!isnan(nearestIntersection.x)) {
                localResults[localResultCount].pointId = i + 1;
                localResults[localResultCount].intersectionId = intersectionIndex++;
                localResults[localResultCount].rayAngle = ray_angle;
                localResults[localResultCount].intersectionPoint = nearestIntersection;
                localResultCount++;
            }
        }
       }
    }
        

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

        Result allResults[totalResults];

        MPI_Gatherv(localResults, localResultCount, MPI_RESULT_TYPE,
                    allResults, counts, displs, MPI_RESULT_TYPE, 0, MPI_COMM_WORLD);
                
        // 打开文件并写入数据
        FILE *file = fopen("combined_output.txt", "w");
        for (int i = 0; i < totalResults; i++) {
            fprintf(file, "观测点id: %d, 射线角度：%d, 最近交点坐标: (%.2f, %.2f, %.2f)\n",
                    allResults[i].pointId, allResults[i].rayAngle,
                    allResults[i].intersectionPoint.x, allResults[i].intersectionPoint.y, allResults[i].intersectionPoint.z);
        }
        fclose(file);
    } else {
        MPI_Gather(&localResultCount, 1, MPI_INT, NULL, 0, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Gatherv(localResults, localResultCount, MPI_RESULT_TYPE,
                    NULL, NULL, NULL, MPI_RESULT_TYPE, 0, MPI_COMM_WORLD);
    }
}
        



long initialMemoryUsage = 0;
long peakMemoryUsage = 0;
long currentMemoryUsage = 0;



// 函数用于获取当前进程的内存使用情况（单位：千字节）
long getCurrentMemoryUsage() {
    struct mach_task_basic_info info;
    mach_msg_type_number_t size = MACH_TASK_BASIC_INFO_COUNT;
    kern_return_t kerr = task_info(mach_task_self(), MACH_TASK_BASIC_INFO, (task_info_t)&info, &size);
    if (kerr == KERN_SUCCESS) {
        return (long)(info.resident_size / (1024 * 1024)); // 转换为兆字节
    } else {
        fprintf(stderr, "无法获取内存使用情况: %s\n", mach_error_string(kerr));
        return -1;
    }
}


int main() {
    // 在程序开始前记录初始内存使用情况
    initialMemoryUsage = getCurrentMemoryUsage();


    MPI_Init(NULL, NULL);

    MPI_Datatype types[4] = { MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE };
    int blocklengths[4] = { 1, 1, 1, 3 };  // 3 doubles for intersectionPoint (x, y, z)
    MPI_Aint offsets[4];

    offsets[0] = offsetof(Result, pointId);
    offsets[1] = offsetof(Result, intersectionId);
    offsets[2] = offsetof(Result, rayAngle);
    offsets[3] = offsetof(Result, intersectionPoint);

    MPI_Type_create_struct(4, blocklengths, offsets, types, &MPI_RESULT_TYPE);
    MPI_Type_commit(&MPI_RESULT_TYPE);



    clock_t start_time = clock();

    // 连接到 Redis 服务器
    redisContext *conn = redisConnect("127.0.0.1", 6379);
    if (conn == NULL || conn->err) {
        if (conn) {
            printf("Connection error: %s\n", conn->errstr);
            redisFree(conn);
        } else {
            printf("Connection error: can't allocate redis context\n");
        }
        exit(1);
    }

    // 读取点数据并将其写入 Redis 列表
    int numPoints;
    Point3D *points = readPoints("point.csv", &numPoints);
    if (points == NULL) {
        redisFree(conn); // 关闭连接
        exit(1);
    }

    const char *pointListKey = "points_list"; // Redis 列表的键名

    // 清空 "points_list" 列表
    redisReply *delReply = redisCommand(conn, "DEL points_list");
    if (delReply == NULL) {
        printf("Failed to execute DEL command for points_list\n");
        // 处理错误
    } else {
        freeReplyObject(delReply);
    }

    for (int i = 0; i < numPoints; i++) {
        char pointData[100];
        snprintf(pointData, sizeof(pointData), "%.2lf,%.2lf,%.2lf", points[i].x, points[i].y, points[i].z);
        redisReply *reply = redisCommand(conn, "RPUSH %s %s", pointListKey, pointData);
        if (reply == NULL) {
            printf("Failed to execute RPUSH command\n");
        } else {
            freeReplyObject(reply);
        }
    }

    // 从 triangles.csv 文件中读取三角面数据
    int numTriangles;
    Triangle3D *triangles = readTriangles("triangles.csv", &numTriangles);
    if (triangles == NULL) {
        return 1;
    }

    const char *triangleListKey = "triangles_list"; // Redis 列表的键名
    
    
    // 清空 "triangles_list" 列表
    delReply = redisCommand(conn, "DEL triangles_list"); // 重新使用同一变量名
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
            redisReply *reply = redisCommand(conn, "RPUSH %s %s", triangleListKey, triangleData);
            if (reply == NULL) {
                printf("Failed to execute RPUSH command for triangles\n");
            } else {
                freeReplyObject(reply);
            }
        }

    // 对于每个读取到的三角面，检查线与面的交点
    checkLineIntersections(points, numPoints, triangles, numTriangles);

    clock_t end_time = clock();
    
    double cpu_time_used = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;
    printf("运行时间: %f s, numpoints: %d, numtriangles: %d\n", cpu_time_used, numPoints, numTriangles);
    
    // 在程序结束前记录最终内存使用情况
    peakMemoryUsage = getCurrentMemoryUsage();

    // 断开连接并释放资源
    redisFree(conn);
    free(points);
    free(triangles);


    MPI_Type_free(&MPI_RESULT_TYPE);
    MPI_Finalize();

    // 在程序结束时记录最终内存使用情况
    currentMemoryUsage = getCurrentMemoryUsage();

    // 计算最高峰内存使用情况和差值
    long peakMemoryDelta = peakMemoryUsage - initialMemoryUsage;
    long differenceMemory = currentMemoryUsage - initialMemoryUsage;
    printf("程序运行期间的最高峰内存使用情况: %ld KB\n", peakMemoryDelta);
    printf("程序运行前后的内存差值: %ld KB\n", differenceMemory);


    return 0;
}
