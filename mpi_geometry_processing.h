// mpi_geometry_processing.h

#ifndef MPI_GEOMETRY_PROCESSING_H
#define MPI_GEOMETRY_PROCESSING_H

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "hiredis-master/hiredis.h"

// 定义结构体，因为它们被公共函数用到

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
    int directionId;
    Point3D intersectionPoint;
    double parameter;
} Result;

// 声明全局变量
extern int globalIntersectionIndex;

// 声明公共函数
typedef Result* (*ComputeFunc)(Point3D *, int, Triangle3D *, int, int *);

Result* mpiGeometryProcessing(int m, ComputeFunc computeFunc, int *totalResultsOut);

Result* checkLineIntersections(Point3D *points, int numPoints, Triangle3D *triangles, int numTriangles, int *totalResultsOut);

Result* checkIntersectionsAngle(Point3D *points, int numPoints, Triangle3D *triangles, int numTriangles, int *totalResultsOut);

// 如果有其他公共函数，也在此处声明

#endif // MPI_GEOMETRY_PROCESSING_H
