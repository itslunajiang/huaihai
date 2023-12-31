#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_LINE_LENGTH 1024

int compare_files(const char *file1, const char *file2) {
    char line1[MAX_LINE_LENGTH], line2[MAX_LINE_LENGTH];
    FILE *f1 = fopen(file1, "r");
    FILE *f2 = fopen(file2, "r");
    if (!f1 || !f2) {
        printf("Error opening files.\n");
        return -1;
    }

    while (fgets(line1, sizeof(line1), f1) && fgets(line2, sizeof(line2), f2)) {
        if (strcmp(line1, line2) != 0) {
            fclose(f1);
            fclose(f2);
            return 1;  // files are different
        }
    }

    // Check if either file has more lines
    if (fgets(line1, sizeof(line1), f1) || fgets(line2, sizeof(line2), f2)) {
        fclose(f1);
        fclose(f2);
        return 1;  // files are different
    }

    fclose(f1);
    fclose(f2);
    return 0;  // files are identical
}

int main() {
    const char *tempFile = "temp_combined.txt";
    FILE *temp = fopen(tempFile, "w");
    if (!temp) {
        printf("Error creating temporary file.\n");
        return 1;
    }

    for (int i = 2; i <= 5; i++) {
        char filename[10];
        sprintf(filename, "file%d.txt", i);

        FILE *f = fopen(filename, "r");
        if (!f) {
            printf("Error opening %s.\n", filename);
            return 1;
        }

        char line[MAX_LINE_LENGTH];
        while (fgets(line, sizeof(line), f)) {
            fputs(line, temp);
        }
        fclose(f);
    }

    fclose(temp);

    // Assuming the combined file and file1 are sorted
    int result = compare_files("file1.txt", tempFile);
    if (result == 0) {
        printf("The files are identical.\n");
    } else {
        printf("The files are different.\n");
    }

    // Optionally delete temp file
    remove(tempFile);

    return 0;
}
