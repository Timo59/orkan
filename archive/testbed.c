//
// Created by Timo Ziegler on 26.08.24.
//
#include <stdio.h>
#include <stdlib.h>
#include "qlib.h"

char* itobin(uint64_t integer, uint64_t length) {
    char* out = malloc((length + 1) * sizeof(char));
    out[length] = '\0';
    for (uint64_t i = 0; i < length; ++i) {
        if (integer & (1ULL << i)) {
            out[length - i - 1] = '1';
        }
        else {
            out[length - i - 1] = '0';
        }
    }
    return out;
}

void itobinPrint(uint64_t integer, uint64_t length) {
    for (uint64_t i = 0; i < length; ++i) {
        if (integer & (1ULL << i)) {
            printf("1");
        }
        else {
            printf("0");
        }
    }
}

bool binsrch(double arr[], unsigned long long arrc, double val, unsigned long long* index) {
    unsigned long denom = 2;
    *index = arrc / denom;                                                      // j initialized to the floor of arrc/2
    while (*index > 0 && *index < arrc - 1) {
        denom *= 2;
        if (val < arr[*index]) {
            if (arr[*index - 1] <= val) {                                   // If arr[i - 1] <= val < arr[i], index of
                return (arr[*index - 1] != val);                            // val is i
            }
            *index -= (arrc + denom - 1) / denom;
        }
        else if (val > arr[*index]) {
            if (arr[*index + 1] >= val) {                                   // If arr[i] < val <= arr[i + 1], index of
                return (arr[++(*index)] != val);                            // val is i + 1
            }
            *index += (arrc + denom - 1) / denom;
        }
        else {
            return 0;
        }
    }
    if (val == arr[*index]) {
        return 0;
    }
    else if (val > arr[*index]) {
        *index = arrc;
    }
    return 1;
}

unsigned long long valord(double val[], unsigned long long valc, double** arr) {
    unsigned long long arrc = 1;
    if ((*arr = malloc(sizeof(double))) == NULL) {
        fprintf(stderr, "Arr allocation in valord\n");
        return NAN;
    }
    **arr = val[0];

    for (unsigned long long i = 1; i < valc; ++i) {
        unsigned long long index;

        if (binsrch(*arr, arrc, val[i], &index)) {
            if((*arr = realloc(*arr, (++arrc) * sizeof(double))) == NULL) {
                fprintf(stderr, "Arr reallocation failed in valord");
                free(*arr);
                return NAN;
            }

            for (unsigned long long j = arrc - 1; j > index; --j) {
                *(*arr + j) = *(*arr + (j - 1));
            }

            *(*arr + index) = val[i];
        }
        printf("Current sorted list: ");
        for (unsigned long long j = 0; j < arrc - 1 ; ++j) {
            printf("%f, ", *(*arr + j));
        }
        printf("%f\n", *(*arr + arrc - 1));
    }
    printf("arrc = %lld\n", arrc);
    return arrc;
}



int main() {
    double values[16] = {4, 8,5,1, 4, 5, 1, 2, 6, 6, 4, 5, 8};
    double* arr;

    unsigned long long length = valord(values, 13, &arr);
    for (unsigned long long i = 0; i < length; ++i) {
        printf("%f, ", arr[i]);
    }
    printf("\n");
    free(arr);
}