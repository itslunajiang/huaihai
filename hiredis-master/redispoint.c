#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hiredis.h"

// 读取点数据并将其写入 Redis 列表
void readAndWritePointsToRedis(redisContext *conn) {
    FILE *pointfile = fopen("point.csv", "r");
    if (pointfile == NULL) {
        perror("无法打开CSV文件");
        return;
    }

    const char *pointListKey = "points_list"; // Redis 列表的键名
    char line[1024];

    while (fgets(line, sizeof(line), pointfile)) {
        double x, y, z;
        if (sscanf(line, "%lf,%lf,%lf", &x, &y, &z) == 3) {
            // 构建点数据字符串
            char pointData[100];
            snprintf(pointData, sizeof(pointData), "%.2lf,%.2lf,%.2lf", x, y, z);

            // 将点数据写入 Redis 列表
            redisReply *reply = redisCommand(conn, "RPUSH %s %s", pointListKey, pointData);
            if (reply == NULL) {
                printf("Failed to execute RPUSH command\n");
            } else {
                freeReplyObject(reply);
            }
        } else {
            printf("无法解析坐标行: %s\n", line);
        }
    }

    fclose(pointfile);
}

int main() {
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
    readAndWritePointsToRedis(conn);

    // 断开连接并释放资源
    redisFree(conn);

    return 0;
}
