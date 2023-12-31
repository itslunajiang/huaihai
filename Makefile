# 定义编译器
CC = mpicc

# 定义编译选项
CFLAGS = -Wall

# 定义链接库
LIBS = -lhiredis

# 定义目标可执行文件
TARGET = my_program

# 定义源文件
SOURCES = main.c mpi_geometry_processing.c

# 默认目标
all: $(TARGET)

# 定义如何编译目标文件
$(TARGET): $(SOURCES)
	$(CC) $(CFLAGS) -o $(TARGET) $(SOURCES) $(LIBS)

# 定义清理规则
clean:
	rm -f $(TARGET)

# 定义运行规则
run: $(TARGET)
	mpirun -np 4 ./$(TARGET)
