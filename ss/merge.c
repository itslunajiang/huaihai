#include <stdio.h>
#include <stdlib.h>
#include <time.h>


void merge_files(const char *output, const char *input1, const char *input2, const char *input3, const char *input4) {
    FILE *out = fopen(output, "w");
    if (!out) {
        perror("Error opening output file");
        exit(1);
    }

    const char *files[] = {input1, input2, input3, input4};
    for (int i = 0; i < 4; i++) {
        FILE *in = fopen(files[i], "r");
        if (!in) {
            perror("Error opening input file");
            exit(1);
        }

        char ch;
        while ((ch = fgetc(in)) != EOF) {
            fputc(ch, out);
        }

        fclose(in);
    }

    fclose(out);
}

int main() {
    clock_t start_time = clock();

    merge_files("combined.txt", "output_z_0.txt", "output_z_1.txt", "output_z_2.txt", "output_z_3.txt");

    clock_t end_time = clock();

    double cpu_time_used = ((double) (end_time - start_time)) / CLOCKS_PER_SEC;        
    
    printf("运行时间: %f s", cpu_time_used);

    return 0;
}
