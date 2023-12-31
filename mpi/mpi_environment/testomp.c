#include <stdio.h>
#include <omp.h>
int main()
{
#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        printf("Hello, OpenMP! Thread ID: %d\n", tid);
    }
    return 0;
}