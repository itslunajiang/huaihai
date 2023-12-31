#include <stdio.h>
#include <stdlib.h>
#include <string.h> 

// 函数用于计算文件中数据的大小
size_t calculateDataSize(const char* filename) {
    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        perror("无法打开文件");
        return 0;
    }

    size_t dataSize = 0;
    char line[1024];

    while (fgets(line, sizeof(line), file) != NULL) {
        dataSize += strlen(line); // 计算每行的大小
    }

    fclose(file);
    return dataSize;
}

int main() {
    const char* pointFilename = "point.csv";
    const char* trianglesFilename = "triangles.csv";

    size_t pointDataSize = calculateDataSize(pointFilename);
    size_t trianglesDataSize = calculateDataSize(trianglesFilename);

    printf("文件 %s 中数据的大小：%zu 字节\n", pointFilename, pointDataSize);
    printf("文件 %s 中数据的大小：%zu 字节\n", trianglesFilename, trianglesDataSize);

    return 0;
}
