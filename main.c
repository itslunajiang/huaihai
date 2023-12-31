#include "mpi_geometry_processing.h"

int main() {

    int m = 50000; // 每次加载的面数量
    int totalResults;
    Result* results = mpiGeometryProcessing(m, checkLineIntersections, &totalResults);

    // 打开文件用于写入
    if (results != NULL) {
        FILE *file = fopen("intersections.csv", "w");
        if (file == NULL) {
            printf("无法创建文件\n");
            // 释放内存
            free(results);
            MPI_Finalize();
            return -1;
        }
        // 写入标题行
        fprintf(file, "PointID,DirectionID,IntersectionX,IntersectionY,IntersectionZ,Parameter\n");
        // 写入数据
        for (int i = 0; i < totalResults; i++) {
            fprintf(file, "%d,%d,%.2f,%.2f,%.2f,%.2f\n",
                    results[i].pointId,
                    results[i].directionId,
                    results[i].intersectionPoint.x,
                    results[i].intersectionPoint.y,
                    results[i].intersectionPoint.z,
                    results[i].parameter);
        }
        // 关闭文件
        fclose(file);
        // 释放内存
        free(results);
    }

    // 建立 Redis 连接，写入redis
    redisContext *conn = redisConnect("127.0.0.1", 6379);
    if (conn == NULL || conn->err) {
        if (conn) {
            printf("Redis连接错误: %s\n", conn->errstr);
            redisFree(conn);
        } else {
            printf("无法分配 Redis 上下文\n");
        }
        return -1; // 退出程序
    }

    if (results != NULL) {
        // 在进程0上，将所有结果写入 Redis
        for (int i = 0; i < totalResults; i++) {
            char redisKey[128];
            snprintf(redisKey, sizeof(redisKey), "pointId:%d,directionId:%d", results[i].pointId, results[i].directionId);
            redisReply *reply = (redisReply *)redisCommand(conn, "HMSET %s IntersectionX %.2f IntersectionY %.2f IntersectionZ %.2f Parameter %.2f",
                                redisKey, results[i].intersectionPoint.x, results[i].intersectionPoint.y, results[i].intersectionPoint.z, results[i].parameter);

            if (reply != NULL) {
                freeReplyObject(reply);
            } else {
                printf("Redis错误：无法在Redis中设置数据\n");
            }
        }

        // 关闭 Redis 连接
        redisFree(conn);
    }
    
    return 0;
}