#include "quicksort.h"
#include "partition.h"
#include <omp.h>

void quicksort(float * A, int lo, int hi) {
    unsigned int p;

    if (lo < hi) {
        p = partition(A, lo, hi);

        #pragma omp task
        quicksort(A, lo, p-1);
        #pragma omp task
        quicksort(A, p+1, hi);

    };
}